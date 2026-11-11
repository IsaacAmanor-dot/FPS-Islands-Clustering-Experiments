#!/usr/bin/env python3
from pathlib import Path
import csv
import math
import re
import sys

ROOT = Path(__file__).resolve().parent
METHODS = ("cartesian", "multigrid_20", "multigrid_323")
CUTOFFS = range(10)
EXPECTED_MOLECULES = 5282

rows = []
errors = []

def integer_field(text, label):
    matches = re.findall(
        rf"^{re.escape(label)}\s*(\d+)\s*$", text, re.M
    )
    if len(matches) != 1:
        raise ValueError(f"Missing or ambiguous field: {label}")
    return int(matches[0])

for cutoff in CUTOFFS:
    for method in METHODS:
        prefix = f"FPS_Island_{method}_{cutoff}"
        path = ROOT / method / "runs" / f"cutoff_{cutoff}" / f"{prefix}.out"

        try:
            text = path.read_text(errors="replace")
            if "Similarity Island clustering complete." not in text:
                raise ValueError("Clustering did not complete.")

            sizes = [
                int(value) for value in re.findall(
                    r"^Island Size:\s*(\d+)\s*$", text, re.M
                )
            ]
            population = integer_field(text, "Molecules Clustered:")
            islands = integer_field(text, "Total Islands:")

            reported_cutoff = re.findall(
                r"^Similarity Cutoff:\s*([+\-\d.eE]+)\s*$", text, re.M
            )
            if len(reported_cutoff) != 1 or not math.isclose(
                float(reported_cutoff[0]), cutoff, abs_tol=1e-9
            ):
                raise ValueError("Reported cutoff differs from directory cutoff.")

            if population != EXPECTED_MOLECULES:
                raise ValueError(f"Unexpected population: {population}")
            if not sizes or min(sizes) < 1:
                raise ValueError("Missing or invalid island sizes.")
            if len(sizes) != islands:
                raise ValueError("Island-size entries differ from total islands.")
            if sum(sizes) != population:
                raise ValueError("Island sizes do not sum to the population.")

            rows.append({
                "method": method,
                "cutoff": cutoff,
                "molecules": population,
                "islands": islands,
                "singletons": sizes.count(1),
                "nonsingleton_islands": sum(size > 1 for size in sizes),
                "molecules_in_nonsingletons": sum(
                    size for size in sizes if size > 1
                ),
                "singleton_fraction": sizes.count(1) / islands,
            })
        except Exception as error:
            errors.append(f"{method}, cutoff {cutoff}: {error}")

print("\nSingleton islands counted directly from DOCK reports")
print(f"{'Cutoff':>6} {'Cartesian':>12} {'MultiGrid 20':>14} {'MultiGrid 323':>15}")

lookup = {(row["method"], row["cutoff"]): row for row in rows}
for cutoff in CUTOFFS:
    counts = [
        str(lookup[method, cutoff]["singletons"])
        if (method, cutoff) in lookup else "ERROR"
        for method in METHODS
    ]
    print(f"{cutoff:6d} {counts[0]:>12} {counts[1]:>14} {counts[2]:>15}")

print("\nCutoff-zero validation")
for method in METHODS:
    row = lookup.get((method, 0))
    if row:
        print(
            f"{method}: {row['islands']} islands; "
            f"{row['singletons']} singletons; "
            f"{row['molecules_in_nonsingletons']} molecules in nonsingleton islands"
        )

# Compare independently counted values with the plotting data.
csv_paths = [
    ROOT / "comparison/all_runs.csv",
    ROOT / "comparison/summary.csv",
    *(ROOT / method / "statistics.csv" for method in METHODS),
]
for path in csv_paths:
    if not path.is_file():
        errors.append(f"Missing statistics file: {path.relative_to(ROOT)}")
        continue

    with path.open(newline="") as stream:
        records = list(csv.DictReader(stream))

    expected_methods = METHODS if path.parent.name == "comparison" else (
        path.parent.name,
    )
    seen = set()
    mismatches = 0

    for record in records:
        try:
            key = (record["method"], float(record["cutoff"]))
            if key in seen:
                raise ValueError(f"Duplicate statistics row: {key}")
            seen.add(key)
            direct = lookup[key]

            for field in ("molecules", "islands", "singletons"):
                if float(record[field]) != direct[field]:
                    raise ValueError(
                        f"{key}: {field} is {record[field]}, "
                        f"direct count is {direct[field]}"
                    )
            if not math.isclose(
                float(record["singleton_fraction"]),
                direct["singleton_fraction"],
                rel_tol=1e-9, abs_tol=1e-12,
            ):
                raise ValueError(f"{key}: singleton fraction differs.")
        except Exception as error:
            errors.append(f"{path.relative_to(ROOT)}: {error}")
            mismatches += 1

    expected_keys = {
        (method, cutoff)
        for method in expected_methods for cutoff in CUTOFFS
    }
    if seen != expected_keys:
        errors.append(f"{path.relative_to(ROOT)}: statistics row set differs.")
        mismatches += 1

    print(f"{path.relative_to(ROOT)}: {'FAIL' if mismatches else 'PASS'}")

destination = ROOT / "comparison/singleton_validation.csv"
if rows:
    with destination.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

if errors:
    print("\nValidation errors:")
    for error in errors:
        print(f"  {error}")
    sys.exit(1)

print(f"\nPASS: all {len(rows)} reports and statistics rows agree.")
print(f"Independent counts: {destination}")
