#!/usr/bin/env python3
"""Prepare fixed inputs and run one method's cutoff experiments."""

from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import Counter
import csv
import hashlib
import json
import math
import os
import shutil
import socket
import subprocess
import sys
import time
import re

ROOT = Path(__file__).resolve().parent
METHODS = ("cartesian", "multigrid_323", "multigrid_20")


def sha(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2) + "\n")


def table(path, columns, rows):
    with Path(path).open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(columns)
        writer.writerows(rows)


def config():
    # Let Bash interpret its own arrays and variable expansions.
    exports = (
        "DOCK_ROOT DOCK_BIN GRID_BIN RECEPTOR_FILE BOX_FILE LIGAND_FILE "
        "MAX_PARALLEL EXPECTED_MOLECULES EXPECTED_RESIDUES "
        "VALIDATION_MOLECULES"
    )

    result = subprocess.run(
        [
            "bash", "-c",
            'source "$1"; printf "%s\\n" "${FPS_CUTOFFS[@]}"',
            "_", str(ROOT / "000.config.sh"),
        ],
        check=True, capture_output=True, text=True,
    )

    raw = subprocess.check_output([
        "bash", "-c",
        f'source "$1"; export {exports}; env -0',
        "_", str(ROOT / "000.config.sh"),
    ])

    values = dict(
        entry.decode().split("=", 1)
        for entry in raw.split(b"\0") if b"=" in entry
    )

    paths = (
        "DOCK_ROOT", "DOCK_BIN", "GRID_BIN",
        "RECEPTOR_FILE", "BOX_FILE", "LIGAND_FILE",
    )
    settings = {key: values[key] for key in paths}

    for key in (
        "MAX_PARALLEL", "EXPECTED_MOLECULES",
        "EXPECTED_RESIDUES", "VALIDATION_MOLECULES",
    ):
        settings[key] = int(values[key])
        if settings[key] < 1:
            raise ValueError(f"{key} must be positive.")

    settings["cutoffs"] = result.stdout.splitlines()
    numbers = [float(value) for value in settings["cutoffs"]]

    if (
        not numbers
        or len(set(numbers)) != len(numbers)
        or any(not math.isfinite(value) or value < 0 for value in numbers)
    ):
        raise ValueError("Cutoffs must be distinct, finite and nonnegative.")

    return settings


def mols(path):
    lines = Path(path).read_text().splitlines(keepends=True)
    starts = [
        index for index, line in enumerate(lines)
        if line.strip() == "@<TRIPOS>MOLECULE"
    ]
    starts.append(len(lines))

    return [
        "".join(lines[start:end])
        for start, end in zip(starts, starts[1:])
    ]


def atoms(text):
    section = None
    records = []

    for line in text.splitlines():
        if line.startswith("@<TRIPOS>"):
            section = line.strip()
        elif section == "@<TRIPOS>ATOM" and line.strip():
            fields = line.split()
            if len(fields) < 9:
                raise ValueError("MOL2 atoms require residue IDs and charges.")
            records.append(fields)

    return records


def frozen():
    manifest = json.loads((ROOT / "inputs/manifest.json").read_text())

    for name, expected in manifest["sha256"].items():
        if sha(ROOT / name) != expected:
            raise ValueError(f"Frozen input changed: {name}")

    return manifest


def output_prefix(method, cutoff):
    label = str(cutoff).replace("-", "minus_").replace(".", "p")
    return f"FPS_Island_{method}_{label}"


def parameters(method, cutoff, data=True, ligand="inputs/ligands.mol2",
               data_mode="individual_islands"):
    engine = "cartesian" if method == "cartesian" else "multigrid"
    prefix = output_prefix(method, cutoff)
    values = [
        ("conformer_search_type", "analysis"),
        ("ligand_atom_file", ROOT / ligand),
        ("limit_max_ligands", "no"),
        ("skip_molecule", "no"),
        ("read_mol_solvation", "no"),
        ("calculate_rmsd", "no"),
        ("cluster_by_similarity_island", "yes"),
        ("similarity_island_clustering_type", "single_score"),
        ("similarity_island_scoring_type", "fps"),
        ("similarity_island_fps_type", engine),
        ("similarity_island_fps_compare_type", "euclidean"),
        ("similarity_island_fps_score_type", "total"),
        ("similarity_island_fps_cutoff", cutoff),
        ("similarity_island_fps_vdw_rep_rad_scale", "1.0"),
        ("similarity_island_fps_receptor_filename", ROOT / "inputs/receptor.mol2"),
    ]
    if engine == "cartesian":
        values += [
            ("similarity_island_fps_grid_prefix", ROOT / "inputs/whole_grid/whole"),
            ("similarity_island_fps_use_distance_dependent_dielectric", "yes"),
            ("similarity_island_fps_dielectric", "4.0"),
        ]
    else:
        components = json.loads((ROOT / method / "components.json").read_text())
        values += [
            ("similarity_island_fps_number_of_grids", len(components["prefixes"])),
            ("similarity_island_fps_hbond_grid_map", ROOT / "selection" / f"{method}_hb.tsv"),
        ]
        values += [
            (f"similarity_island_fps_grid_prefix{i}", ROOT / method / "grids" / name)
            for i, name in enumerate(components["prefixes"])
        ]
    values.append(("similarity_island_fps_write_data", "yes" if data else "no"))
    if data:
        values += [
            ("similarity_island_fps_data_mode", data_mode),
            ("similarity_island_fps_data_file", prefix + ".dat" if data_mode == "combined" else prefix),
        ]
    # Match descriptor-based grid scoring in the original FPS workflow.
    values += [
        ("similarity_island_retention_mode", "all"),
        ("score_molecules", "yes"),
        ("contact_score_primary", "no"),
        ("grid_score_primary", "no"),
        ("gist_score_primary", "no"),
        ("multigrid_score_primary", "no"),
        ("dock3.5_score_primary", "no"),
        ("continuous_score_primary", "no"),
        ("footprint_similarity_score_primary", "no"),
        ("pharmacophore_score_primary", "no"),
        ("hbond_score_primary", "no"),
        ("internal_energy_score_primary", "no"),
        ("descriptor_score_primary", "yes"),
        ("descriptor_use_grid_score", "yes"),
        ("descriptor_use_grid_lig_efficiency", "no"),
        ("descriptor_use_pharmacophore_score", "no"),
        ("descriptor_use_tanimoto", "no"),
        ("descriptor_use_hungarian", "no"),
        ("descriptor_use_volume_overlap", "no"),
        ("descriptor_use_gist", "no"),
        ("descriptor_use_dock3.5", "no"),
        ("descriptor_grid_score_rep_rad_scale", "1"),
        ("descriptor_grid_score_vdw_scale", "1"),
        ("descriptor_grid_score_es_scale", "1"),
        ("descriptor_grid_score_grid_prefix", ROOT / "inputs/whole_grid/whole"),
        ("descriptor_weight_grid_score", "1"),
        ("atom_model", "all"),
        ("vdw_defn_file", ROOT / "inputs/vdw_AMBER_parm99.defn"),
        ("flex_defn_file", ROOT / "inputs/flex.defn"),
        ("ligand_outfile_prefix", prefix),
    ]
    return "".join(f"{key:<65} {value}\n" for key, value in values)


def read_report(path, expected_n, method=None):
    text = Path(path).read_text(errors="replace")

    def number(pattern, cast=int):
        matches = re.findall(pattern, text)
        if len(matches) != 1:
            raise ValueError(f"Missing or ambiguous report field: {pattern}")
        return cast(matches[0])

    if "Similarity Island clustering complete." not in text:
        raise ValueError("Clustering did not complete.")

    if (
        method
        and method != "cartesian"
        and "HB calculated in Cartesian space" not in text
    ):
        raise ValueError("The expected MultiGrid HB implementation is inactive.")

    population = number(r"Molecules Clustered:\s*(\d+)")
    island_count = number(r"Total Islands:\s*(\d+)")
    sizes = [
        int(value)
        for value in re.findall(r"^Island Size:\s*(\d+)", text, re.M)
    ]

    assignments = {}
    current = None

    for line in text.splitlines():
        match = re.match(r"FPS ISLAND\s+(\d+)\s*$", line)
        if match:
            current = int(match[1])

        fields = line.split()

        if (
            len(fields) == 9
            and fields[0].isdigit()
            and fields[-1] in ("HEAD", "MEMBER")
        ):
            if current is None or fields[1] in assignments:
                raise ValueError("Duplicate molecule or missing island label.")

            scores = [float(value) for value in fields[2:8]]

            if not all(math.isfinite(value) for value in scores):
                raise ValueError("Nonfinite FPS score.")

            if (
                abs(scores[5] - sum(scores[1:4])) > 0.0003
                or abs(scores[0] - scores[5]) > 0.0002
            ):
                raise ValueError("FPS total differs from VDW + ES + HB.")

            assignments[fields[1]] = current

    if (
        population != expected_n
        or sum(sizes) != population
        or len(sizes) != island_count
        or len(assignments) != population
    ):
        raise ValueError("Molecule, island or membership counts differ.")

    if sorted(Counter(assignments.values()).values()) != sorted(sizes):
        raise ValueError("Membership counts differ from reported island sizes.")

    return {
        "islands": island_count,
        "singletons": sizes.count(1),
        "singleton_fraction": sizes.count(1) / island_count,
        "largest": max(sizes),
        "sizes": sizes,
        "assignments": assignments,
        "comparisons": number(r"Total FPS Calculations:\s*(\d+)"),
        "analysis_seconds": number(
            r"Total elapsed time:\s*([\d.]+)", float
        ),
    }


def execute(directory, population, method):
    directory = Path(directory)
    if (directory / "execution.json").exists():
        raise ValueError(f"Existing execution record: {directory}")
    inputs = list(directory.glob("FPS_Island_*.in"))
    if len(inputs) != 1:
        raise ValueError(f"Expected exactly one FPS input: {directory}")
    input_file = inputs[0]
    snapshot = input_file.with_name(input_file.name + ".requested")
    if not snapshot.exists():
        if directory.parent.name == "runs":
            raise ValueError(f"Missing submitted-input snapshot: {snapshot}")
        shutil.copyfile(input_file, snapshot)
    # DOCK rewrites its input; preserve submitted parameters separately.
    shutil.copyfile(snapshot, input_file)
    output_file = input_file.with_suffix(".out")
    values = dict(line.split(None, 1) for line in input_file.read_text().splitlines() if line.strip())
    if (values.get("score_molecules") != "yes"
        or values.get("descriptor_score_primary") != "yes"
        or values.get("descriptor_use_grid_score") != "yes"
        or values.get("grid_score_primary") != "no"):
        raise ValueError("Expected descriptor-based grid rescoring in every input.")
    expected_prefix = input_file.stem
    if values["ligand_outfile_prefix"] != expected_prefix:
        raise ValueError("Output prefix differs from experiment and cutoff name.")
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    command = ["/usr/bin/time", "-f", "%e\t%M", "-o", "resources.tsv",
               str(ROOT / "inputs/software/dock6"),
               "-i", input_file.name, "-o", output_file.name]
    cpu = next((line.split(":", 1)[1].strip()
                for line in Path("/proc/cpuinfo").read_text().splitlines()
                if line.startswith("model name")), "unknown")
    record = dict(host=socket.gethostname(), cpu_model=cpu, command=command,
                  input_file=input_file.name, output_file=output_file.name,
                  status="failed", slurm_job_id=os.environ.get("SLURM_JOB_ID"),
                  slurm_array_task_id=os.environ.get("SLURM_ARRAY_TASK_ID"))
    started = time.perf_counter()
    try:
        with (directory / "execution.log").open("w") as stream:
            subprocess.run(command, cwd=directory, env=environment,
                           stdin=subprocess.DEVNULL, stdout=stream,
                           stderr=subprocess.STDOUT, check=True)
        report = read_report(output_file, population, method)
        rescored = re.findall(r"Molecules rescored:\s*(\d+)", output_file.read_text())
        if rescored != [str(population)]:
            raise ValueError("Rescoring did not cover every input record.")
        if values.get("similarity_island_fps_write_data") == "yes":
            prefix = directory / values["similarity_island_fps_data_file"]
            if values["similarity_island_fps_data_mode"] == "individual_islands":
                files = [prefix.parent / f"{prefix.name}_{i:04d}.dat"
                         for i in range(1, report["islands"] + 1)]
            else:
                files = [prefix]
            if not files or any(not file.is_file() or file.stat().st_size == 0 for file in files):
                raise ValueError("Expected FPS footprint output is missing or empty.")
        wall, memory = (directory / "resources.tsv").read_text().strip().split("\t")
        record.update(status="passed", process_wall_seconds=float(wall), peak_rss_kib=int(memory))
    except Exception as error:
        record["error"] = str(error)
        raise
    finally:
        record["driver_wall_seconds"] = time.perf_counter() - started
        write_json(directory / "execution.json", record)


def prepare():
    settings = config()

    if not Path("/usr/bin/time").is_file():
        raise ValueError("GNU time is required at /usr/bin/time.")

    if (ROOT / "inputs/manifest.json").exists():
        raise ValueError("Already prepared; use a fresh experiment directory.")

    for name in (
        "inputs/software", "inputs/whole_grid",
        "selection", "comparison", *METHODS,
    ):
        (ROOT / name).mkdir(parents=True, exist_ok=True)

    sources = {
        "receptor.mol2": settings["RECEPTOR_FILE"],
        "box.pdb": settings["BOX_FILE"],
        "ligands.mol2": settings["LIGAND_FILE"],
        "software/dock6": settings["DOCK_BIN"],
        "software/grid": settings["GRID_BIN"],
    }

    for name in ("vdw_AMBER_parm99.defn", "flex.defn", "chem.defn"):
        sources[name] = str(
            Path(settings["DOCK_ROOT"]) / "parameters" / name
        )

    for name in ("multigrid_fp_gen.py", "mol2.py"):
        sources["software/" + name] = str(
            Path(settings["DOCK_ROOT"]) / "bin" / name
        )

    for destination, source in sources.items():
        source = Path(source)
        if not source.is_file():
            raise FileNotFoundError(source)

        shutil.copy2(source, ROOT / "inputs" / destination)

    executable = ROOT / "inputs/software/dock6"

    if b"HB calculated in Cartesian space" not in executable.read_bytes():
        raise ValueError("The executable lacks the confirmed MultiGrid HB code.")

    for name in ("dock6", "grid"):
        if not os.access(ROOT / "inputs/software" / name, os.X_OK):
            raise ValueError(f"Not executable: {name}")

    ligand_file = ROOT / "inputs/ligands.mol2"
    ligands = mols(ligand_file)
    if len(ligands) != settings["EXPECTED_MOLECULES"]:
        raise ValueError("Unexpected molecule count.")
    with ligand_file.open(newline="") as stream:
        library = stream.read()
    records = []
    def identify(match):
        index = len(records) + 1
        original = match.group(2).strip()
        label = re.sub(r"\s+", "_", original) or "molecule"
        identifier = f"{label}__SI{index:06d}"
        records.append((index, identifier, original))
        return match.group(1) + identifier + match.group(3)
    library, count = re.subn(
        r"(?m)^(@<TRIPOS>MOLECULE[ \t]*\r?\n)([^\r\n]*)(\r?\n|$)", identify, library,
    )
    if count != len(ligands):
        raise ValueError("Could not identify every MOL2 record.")
    with ligand_file.open("w", newline="") as stream:
        stream.write(library)
    record_map = ROOT / "inputs/ligand_records.tsv"
    table(record_map, ("record_index", "experiment_id", "original_title"), records)
    ligands = mols(ligand_file)

    receptors = mols(ROOT / "inputs/receptor.mol2")
    if len(receptors) != 1:
        raise ValueError("Expected one receptor molecule.")

    receptor_atoms = atoms(receptors[0])
    residues = {}
    order = []

    for atom in receptor_atoms:
        residue_id = int(atom[6])

        if residue_id not in residues:
            residues[residue_id] = [atom[7], 0]
            order.append(residue_id)

        if residues[residue_id][0] != atom[7]:
            raise ValueError("One residue ID has multiple residue names.")

        residues[residue_id][1] += 1

    atom_ids = [int(atom[6]) for atom in receptor_atoms]

    if (
        atom_ids != sorted(atom_ids)
        or len(residues) != settings["EXPECTED_RESIDUES"]
    ):
        raise ValueError("Unexpected residue count or unsorted residue IDs.")

    mapping = ROOT / "selection/receptor_residue_map.tsv"

    table(
        mapping,
        ("footprint_position", "mol2_residue_id", "residue_name", "atom_count"),
        [
            (position, residue_id, *residues[residue_id])
            for position, residue_id in enumerate(order, 1)
        ],
    )

    validation_n = min(settings["VALIDATION_MOLECULES"], len(ligands))
    validation_file = ROOT / "inputs/validation_ligands.mol2"
    validation_file.write_text("".join(ligands[:validation_n]))
    settings["validation_n"] = validation_n

    tracked = [ROOT / "000.config.sh", mapping, validation_file, record_map]
    tracked += [ROOT / "inputs" / name for name in sources]

    write_json(
        ROOT / "inputs/manifest.json",
        {
            "config": settings,
            "sources": sources,
            "source_ligand_sha256": sha(Path(settings["LIGAND_FILE"])),
            "sha256": {
                str(path.relative_to(ROOT)): sha(path) for path in tracked
            },
        },
    )

    print(
        f"Prepared {len(ligands)} molecules and {len(residues)} residues. "
        "No grids generated."
    )


def write_runs():
    settings = frozen()["config"]
    setup_path = ROOT / "comparison/setup_complete.json"
    if not setup_path.exists():
        raise ValueError("Complete grid generation and validation first.")
    setup = json.loads(setup_path.read_text())
    planned = {}
    for method in METHODS:
        for cutoff in settings["cutoffs"]:
            directory = ROOT / method / "runs" / f"cutoff_{cutoff}"
            if directory.exists():
                raise ValueError(f"Existing run directory: {directory}")
            prefix = output_prefix(method, cutoff)
            path = directory / f"{prefix}.in"
            planned[path] = parameters(method, cutoff)
    for path, text in planned.items():
        path.parent.mkdir(parents=True, exist_ok=False)
        path.write_text(text)
        snapshot = path.with_name(path.name + ".requested")
        snapshot.write_text(text)
        setup["sha256"][str(snapshot.relative_to(ROOT))] = sha(snapshot)
    setup["stage"] = "ready"
    setup["benchmark_runs"] = len(planned)
    write_json(setup_path, setup)
    print(f"Prepared and checksummed {len(planned)} runs: one run per method and cutoff.")


def run_method(method):
    settings = frozen()["config"]
    setup = json.loads((ROOT / "comparison/setup_complete.json").read_text())
    if setup.get("stage") != "ready":
        raise ValueError("Grid setup and input generation are not complete.")
    for name, expected in setup["sha256"].items():
        if sha(ROOT / name) != expected:
            raise ValueError(f"Frozen setup file changed: {name}")
    folder = ROOT / ("inputs/whole_grid" if method == "cartesian" else method + "/grids")
    with (folder / "grid_checksums.tsv").open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if sha(folder / row["file"]) != row["sha256"]:
                raise ValueError(f"Grid checksum differs: {row['file']}")
    allocated = int(os.environ.get("SLURM_CPUS_PER_TASK", settings["MAX_PARALLEL"]))
    if allocated < settings["MAX_PARALLEL"]:
        raise ValueError("Allocated CPUs must be at least MAX_PARALLEL.")
    with ThreadPoolExecutor(max_workers=settings["MAX_PARALLEL"]) as pool:
        jobs = {pool.submit(execute, ROOT / method / "runs" / f"cutoff_{cutoff}",
                            settings["EXPECTED_MOLECULES"], method): cutoff
                for cutoff in settings["cutoffs"]}
        for future in as_completed(jobs):
            future.result()
            print(f"{method} cutoff {jobs[future]}: PASS", flush=True)


if __name__ == "__main__":
    if sys.argv[1:] == ["prepare"]:
        prepare()
    elif sys.argv[1:] == ["write-runs"]:
        write_runs()
    elif (
        len(sys.argv) == 3
        and sys.argv[1] == "run"
        and sys.argv[2] in METHODS
    ):
        run_method(sys.argv[2])
    else:
        raise SystemExit(
            "Usage: 001.prepare_experiments.py prepare|write-runs|run METHOD"
        )

