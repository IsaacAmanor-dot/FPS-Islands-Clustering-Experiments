#!/usr/bin/env python3
"""Generate grids, select residues and validate HB component mappings."""

from pathlib import Path
import array
import json
import csv
import importlib.util
import math
import os
import re
import shlex
import struct
import subprocess
import sys
import time

spec = importlib.util.spec_from_file_location(
    "experiment",
    Path(__file__).with_name("001.prepare_experiments.py"),
)
experiment = importlib.util.module_from_spec(spec)
spec.loader.exec_module(experiment)
ROOT = experiment.ROOT


def grid_input(directory, receptor="temp.mol2", prefix="temp.rec"):
    relative = lambda name: os.path.relpath(
        ROOT / "inputs" / name, directory
    )

    text = f"""allow_non_integral_charges yes
compute_grids yes
grid_spacing 0.3
output_molecule yes
contact_score no
chemical_score no
energy_score yes
energy_cutoff_distance 9999
atom_model a
attractive_exponent 6
repulsive_exponent 9
distance_dielectric yes
dielectric_factor 4
bump_filter yes
bump_overlap 0.75
receptor_file {receptor}
box_file {relative('box.pdb')}
vdw_definition_file {relative('vdw_AMBER_parm99.defn')}
chemical_definition_file {relative('chem.defn')}
score_grid_prefix {prefix}
receptor_out_file {prefix}.grid.mol2
"""

    for line in text.splitlines():
        key, value = line.split()
        if (
            ("file" in key or "prefix" in key)
            and len(value.encode()) > 80
        ):
            raise ValueError(f"Legacy GRID filename limit: {line}")

    (directory / "grid.in").write_text(text)


def validate_grids(directory, prefixes):
    geometry = None
    checksums = []

    for prefix in prefixes:
        energy = directory / f"{prefix}.nrg"
        bump = directory / f"{prefix}.bmp"

        with energy.open("rb") as stream:
            header = stream.read(16)

        with bump.open("rb") as stream:
            box = stream.read(32)

        size, model, attractive, repulsive = struct.unpack("=4i", header)
        bsize, spacing, ox, oy, oz, nx, ny, nz = struct.unpack(
            "=if3f3i", box
        )

        if (
            size <= 0
            or size != bsize
            or size != nx * ny * nz
            or min(nx, ny, nz) < 4
            or model != ord("a")
            or (attractive, repulsive) != (6, 9)
            or not math.isclose(spacing, 0.3, abs_tol=1e-6)
            or not all(
                math.isfinite(value)
                for value in (spacing, ox, oy, oz)
            )
            or energy.stat().st_size != 16 + 12 * size
            or bump.stat().st_size < 32 + size
        ):
            raise ValueError(f"Invalid grid: {energy}")

        current = (header, box)

        if geometry is not None and geometry != current:
            raise ValueError("Component grid geometries differ.")

        geometry = current

        with energy.open("rb") as stream:
            stream.seek(16)

            for chunk in iter(
                lambda: stream.read(4 * 262144), b""
            ):
                values = array.array("f")
                values.frombytes(chunk)

                if not all(math.isfinite(value) for value in values):
                    raise ValueError(f"Nonfinite grid values: {energy}")

        checksums.extend(
            (path.name, path.stat().st_size, experiment.sha(path))
            for path in (energy, bump)
        )

    experiment.table(
        directory / "grid_checksums.tsv",
        ("file", "bytes", "sha256"),
        checksums,
    )

    return geometry


def checked_generator():
    original = (
        ROOT / "inputs/software/multigrid_fp_gen.py"
    ).read_text()

    pattern = r"(?m)^([ \t]*)gridpath[ \t]*=.*$"

    if (
        len(re.findall(pattern, original)) != 1
        or "os.system(" not in original
    ):
        raise ValueError("Unexpected native generator; inspect before adapting.")

    text = re.sub(
        pattern,
        lambda match: (
            match[1]
            + "gridpath = "
            + repr(
                shlex.quote(str(ROOT / "inputs/software/grid"))
            )
        ),
        original,
    )

    text = text.replace("os.system(", "checked_system(")
    anchor = "from __future__ import print_function"

    if text.count(anchor) != 1:
        raise ValueError("Unexpected native future import.")

    helper = """
import subprocess

def checked_system(command):
    print("Executing:", command, flush=True)
    subprocess.run(command, shell=True, check=True)
"""

    text = text.replace(anchor, anchor + "\n" + helper, 1)
    path = ROOT / "inputs/software/multigrid_fp_gen_checked.py"

    compile(text, str(path), "exec")
    path.write_text(text)

    return path


def export(method, name, ligand, population):
    directory = ROOT / method / name
    directory.mkdir(exist_ok=False)
    prefix = experiment.output_prefix(method, "1e12")
    (directory / f"{prefix}.in").write_text(
        experiment.parameters(method, "1e12", True, ligand, data_mode="combined")
    )
    experiment.execute(directory, population, method)
    return directory / f"{prefix}.dat"


def footprint_rows(path):
    header_found = False

    with path.open() as stream:
        for line in stream:
            if line.startswith("# Island Rank Molecule Role "):
                header_found = True

            if not line.strip() or line.startswith("#"):
                continue

            fields = line.split()

            if not header_found or len(fields) != 16:
                raise ValueError(f"Unexpected footprint format: {path}")

            # Member VDW, ES and HB columns from the combined export.
            values = [float(fields[index]) for index in (7, 10, 13)]

            if not all(math.isfinite(value) for value in values):
                raise ValueError("Nonfinite footprint.")

            yield (int(fields[0]), int(fields[1])), int(fields[4]), values

    if not header_found:
        raise ValueError(f"Missing footprint header: {path}")


def select_residues(path, residues, population):
    means = [[0.0] * 3 for _ in residues]
    counts = [0] * len(residues)
    seen = set()

    for name, index, values in footprint_rows(path):
        key = (name, index)

        if key in seen or not 1 <= index <= len(residues):
            raise ValueError("Duplicate or invalid residue footprint.")

        seen.add(key)
        counts[index - 1] += 1

        for component in range(3):
            means[index - 1][component] += (
                abs(values[component]) / population
            )

    if counts != [population] * len(residues):
        raise ValueError("Incomplete Cartesian selection footprints.")

    denominators = [
        sum(values[index] for values in means)
        for index in range(3)
    ]

    ranked = []

    for residue, values in zip(residues, means):
        normalized = [
            value / denominator if denominator else 0.0
            for value, denominator in zip(values, denominators)
        ]

        ranked.append((
            sum(normalized),
            int(residue["mol2_residue_id"]),
            residue["residue_name"],
            values,
            normalized,
        ))

    ranked.sort(key=lambda row: (-row[0], row[1]))

    experiment.table(
        ROOT / "selection/residue_importance.tsv",
        (
            "rank", "mol2_residue_id", "residue_name",
            "mean_abs_vdw", "mean_abs_es", "mean_abs_hb",
            "normalized_vdw", "normalized_es", "normalized_hb",
            "importance",
        ),
        [
            (
                rank, residue_id, name,
                *values, *normalized, score,
            )
            for rank, (
                score, residue_id, name, values, normalized
            ) in enumerate(ranked, 1)
        ],
    )

    return sorted(row[1] for row in ranked[:19])


def component_setup(method, selected, residue_ids):
    directory = ROOT / method / "grids"
    directory.mkdir(exist_ok=False)

    remaining = [
        residue_id for residue_id in residue_ids
        if residue_id not in selected
    ]

    if not remaining or len(set(selected)) != len(selected):
        raise ValueError("Selection requires a nonempty remainder.")

    groups = [[residue_id] for residue_id in selected] + [remaining]
    prefixes = [
        f"resid_{residue_id:03d}" for residue_id in selected
    ] + ["resid_remaining"]

    grid_input(directory)

    experiment.table(
        ROOT / "selection" / f"{method}_hb.tsv",
        ("grid_index", "mol2_residue_id"),
        [
            (index, residue_id)
            for index, group in enumerate(groups)
            for residue_id in group
        ],
    )

    experiment.table(
        ROOT / "selection" / f"{method}_components.tsv",
        ("grid_index", "grid_prefix", "mol2_residue_id"),
        [
            (index, prefix, residue_id)
            for index, (prefix, group) in enumerate(zip(prefixes, groups))
            for residue_id in group
        ],
    )

    experiment.write_json(
        ROOT / method / "components.json",
        {"prefixes": prefixes, "groups": groups},
    )

    for label, values in (
        ("selected", selected), ("remaining", remaining)
    ):
        (
            ROOT / "selection" / f"{method}_{label}.txt"
        ).write_text("".join(f"{value}\n" for value in values))

    return directory, prefixes, groups


def validate_partition(directory, prefixes, groups):
    original = experiment.atoms(
        (ROOT / "inputs/receptor.mol2").read_text()
    )

    # The native splitter renumbers fragment residue IDs.
    signature = lambda atom: (
        atom[1],
        *map(float, atom[2:5]),
        atom[5],
        atom[7],
        float(atom[8]),
    )

    for prefix, group in zip(prefixes, groups):
        expected = sorted(
            signature(atom)
            for atom in original if int(atom[6]) in group
        )

        actual = sorted(
            signature(atom)
            for atom in experiment.atoms(
                (directory / f"{prefix}.mol2").read_text()
            )
        )

        if expected != actual:
            raise ValueError(
                f"Partition altered atoms, coordinates or charges: {prefix}"
            )


def compare_footprints(cartesian_path, grid_path, groups, residue_ids, count):
    cartesian = {
        (name, residue_ids[index - 1]): values
        for name, index, values in footprint_rows(cartesian_path)
    }

    actual = {
        (name, index): values
        for name, index, values in footprint_rows(grid_path)
    }

    names = sorted({name for name, _ in cartesian})

    if (
        len(names) != count
        or len(cartesian) != count * len(residue_ids)
        or len(actual) != count * len(groups)
    ):
        raise ValueError("Incomplete footprint validation data.")

    squared = [0.0, 0.0]
    maximum = [0.0, 0.0]

    for name in names:
        for index, group in enumerate(groups, 1):
            target = [
                sum(cartesian[name, residue_id][component]
                    for residue_id in group)
                for component in range(3)
            ]

            value = actual[name, index]

            if value[2] != target[2]:
                raise ValueError(
                    f"HB mapping differs: {name}, component {index}"
                )

            for component in range(2):
                difference = value[component] - target[component]
                squared[component] += difference * difference
                maximum[component] = max(
                    maximum[component], abs(difference)
                )

    return {
        "hb_exact_match": True,
        "molecules": count,
        "nonzero_cartesian_hb_components": sum(
            values[2] != 0 for values in cartesian.values()
        ),
        "vdw_rmse": math.sqrt(squared[0] / len(actual)),
        "es_rmse": math.sqrt(squared[1] / len(actual)),
        "vdw_max_abs_difference": maximum[0],
        "es_max_abs_difference": maximum[1],
    }


def main():
    manifest = experiment.frozen()
    settings = manifest["config"]

    if any(character.isspace() for character in str(ROOT)):
        raise ValueError("The native generator requires paths without whitespace.")

    if (ROOT / "comparison/setup_complete.json").exists():
        raise ValueError("Setup already complete; refusing regeneration.")

    with (
        ROOT / "selection/receptor_residue_map.tsv"
    ).open() as stream:
        residues = list(csv.DictReader(stream, delimiter="\t"))

    residue_ids = [
        int(row["mol2_residue_id"]) for row in residues
    ]

    if len(residue_ids) != 323:
        raise ValueError("This experiment requires 323 receptor residues.")

    started = time.perf_counter()
    times = {}
    whole = ROOT / "inputs/whole_grid"

    if list(whole.glob("*.nrg")):
        raise ValueError("Existing grids found; use a fresh experiment.")

    grid_input(whole, "../receptor.mol2", "whole")
    step = time.perf_counter()

    with (whole / "execution.log").open("w") as stream:
        subprocess.run(
            [
                str(ROOT / "inputs/software/grid"),
                "-i", "grid.in", "-o", "grid.out",
            ],
            cwd=whole,
            stdout=stream,
            stderr=subprocess.STDOUT,
            check=True,
        )

    geometry = validate_grids(whole, ["whole"])
    times["whole_grid_seconds"] = time.perf_counter() - step

    step = time.perf_counter()

    selection_data = export(
        "cartesian", "selection_export",
        "inputs/ligands.mol2",
        settings["EXPECTED_MOLECULES"],
    )

    selected = select_residues(
        selection_data, residues, settings["EXPECTED_MOLECULES"]
    )

    times["cartesian_selection_seconds"] = time.perf_counter() - step
    generator = checked_generator()

    cartesian_validation = export(
        "cartesian", "validation",
        "inputs/validation_ligands.mol2",
        settings["validation_n"],
    )

    reports = {}

    for method, chosen in (
        ("multigrid_323", residue_ids[:-1]),
        ("multigrid_20", selected),
    ):
        directory, prefixes, groups = component_setup(
            method, chosen, residue_ids
        )

        command = [
            sys.executable, "-u", str(generator),
            str(ROOT / "inputs/receptor.mol2"),
            "resid", "grid.in", *map(str, chosen),
        ]

        step = time.perf_counter()

        with (directory / "generation.log").open("w") as stream:
            subprocess.run(
                command,
                cwd=directory,
                stdout=stream,
                stderr=subprocess.STDOUT,
                check=True,
            )

        times[method + "_generation_seconds"] = (
            time.perf_counter() - step
        )

        validate_partition(directory, prefixes, groups)

        if validate_grids(directory, prefixes) != geometry:
            raise ValueError("Component geometry differs from the whole grid.")

        output = export(
            method, "validation",
            "inputs/validation_ligands.mol2",
            settings["validation_n"],
        )

        reports[method] = compare_footprints(
            cartesian_validation,
            output,
            groups,
            residue_ids,
            settings["validation_n"],
        )

        print(
            f"{method}: generation, partition and HB checks passed.",
            flush=True,
        )

    times["setup_total_seconds"] = time.perf_counter() - started

    experiment.write_json(
        ROOT / "comparison/setup_times.json", times
    )
    experiment.write_json(
        ROOT / "comparison/footprint_validation.json", reports
    )

    tracked = list((ROOT / "selection").glob("*"))

    tracked += [
        ROOT / method / "components.json"
        for method in experiment.METHODS[1:]
    ]

    tracked += [
        ROOT / method / "grids/grid_checksums.tsv"
        for method in experiment.METHODS[1:]
    ]

    tracked.append(ROOT / "inputs/whole_grid/grid_checksums.tsv")

    tracked += [
        ROOT / name
        for name in (
            "001.prepare_experiments.py",
            "002.run_experiments.slurm",
            "003.compare_results.py",
            "004.generate_multigrids.py",
        )
    ]

    experiment.write_json(
        ROOT / "comparison/setup_complete.json",
        {
            "sha256": {
                str(path.relative_to(ROOT)): experiment.sha(path)
                for path in tracked
            }
        },
    )

    experiment.write_runs()
    print("Setup complete. Submit 002.run_experiments.slurm.")


def validate_existing():
    settings = experiment.frozen()["config"]
    if (ROOT / "comparison/setup_complete.json").exists():
        raise ValueError("Setup already complete; preserve it before revalidation.")
    with (ROOT / "selection/receptor_residue_map.tsv").open() as stream:
        residues = list(csv.DictReader(stream, delimiter="\t"))
    ids = [int(row["mol2_residue_id"]) for row in residues]
    if len(ids) != 323:
        raise ValueError("Expected 323 receptor residues.")
    started = time.perf_counter()
    geometry = validate_grids(ROOT / "inputs/whole_grid", ["whole"])
    specifications = {}
    for method, count in (("multigrid_323", 323), ("multigrid_20", 20)):
        specification = json.loads((ROOT / method / "components.json").read_text())
        prefixes, groups = specification["prefixes"], specification["groups"]
        if len(prefixes) != count or len(groups) != count:
            raise ValueError(f"Expected {count} components in {method}.")
        if sorted(rid for group in groups for rid in group) != sorted(ids):
            raise ValueError(f"Components do not partition the receptor: {method}")
        mapped = [(i, rid) for i, group in enumerate(groups) for rid in group]
        with (ROOT / "selection" / f"{method}_hb.tsv").open() as stream:
            saved = [(int(row["grid_index"]), int(row["mol2_residue_id"]))
                     for row in csv.DictReader(stream, delimiter="\t")]
        if saved != mapped:
            raise ValueError(f"HB mapping differs from grid components: {method}")
        folder = ROOT / method / "grids"
        validate_partition(folder, prefixes, groups)
        if validate_grids(folder, prefixes) != geometry:
            raise ValueError(f"Grid geometry differs: {method}")
        specifications[method] = specification
        print(f"{method}: existing grid files and partitions validated.", flush=True)

    step = time.perf_counter()
    selection_data = export("cartesian", "selection_export", "inputs/ligands.mol2",
                            settings["EXPECTED_MOLECULES"])
    ranking_path = ROOT / "selection/residue_importance.tsv"
    previous_ranking = ranking_path.read_bytes()
    selected = select_residues(selection_data, residues, settings["EXPECTED_MOLECULES"])
    saved_selected = [group[0] for group in specifications["multigrid_20"]["groups"][:-1]]
    if selected != saved_selected:
        ranking_path.write_bytes(previous_ranking)
        raise ValueError("Scored footprint selection changed; the 20-grid set must be reviewed before benchmarking.")
    selection_seconds = time.perf_counter() - step
    cartesian = export("cartesian", "validation", "inputs/validation_ligands.mol2", settings["validation_n"])
    reports = {}
    for method in experiment.METHODS[1:]:
        output = export(method, "validation", "inputs/validation_ligands.mol2", settings["validation_n"])
        reports[method] = compare_footprints(cartesian, output,
                                           specifications[method]["groups"], ids, settings["validation_n"])
        print(f"{method}: scored footprint validation passed.", flush=True)
    experiment.write_json(ROOT / "comparison/footprint_validation.json", reports)
    times_file = ROOT / "comparison/setup_times.json"
    times = json.loads(times_file.read_text()) if times_file.exists() else {}
    times.update(rescored_selection_seconds=selection_seconds,
                 reuse_validation_seconds=time.perf_counter() - started)
    experiment.write_json(times_file, times)
    tracked = list((ROOT / "selection").glob("*"))
    tracked += [ROOT / method / "components.json" for method in experiment.METHODS[1:]]
    tracked += [ROOT / method / "grids/grid_checksums.tsv" for method in experiment.METHODS[1:]]
    tracked += [ROOT / "inputs/whole_grid/grid_checksums.tsv"]
    tracked += [ROOT / name for name in ("001.prepare_experiments.py", "001a.generate_multigrids.slurm",
                                         "002.run_experiments.slurm", "003.compare_results.py", "004.generate_multigrids.py")]
    experiment.write_json(ROOT / "comparison/setup_complete.json",
                          {"sha256": {str(path.relative_to(ROOT)): experiment.sha(path) for path in tracked}})
    experiment.write_runs()
    print("Existing grids reused. Scored setup and all benchmark inputs are ready.", flush=True)


if __name__ == "__main__":
    if sys.argv[1:] == ["validate"]:
        validate_existing()
    elif not sys.argv[1:]:
        main()
    else:
        raise SystemExit("Usage: 004.generate_multigrids.py [validate]")
