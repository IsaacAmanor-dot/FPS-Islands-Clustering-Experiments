#!/usr/bin/env python3
"""Compare runtime, singleton statistics and clustering membership."""

from pathlib import Path
from collections import Counter
from itertools import combinations
import csv
import importlib.util
import json
import statistics

spec = importlib.util.spec_from_file_location(
    "experiment",
    Path(__file__).with_name("001.prepare_experiments.py"),
)
experiment = importlib.util.module_from_spec(spec)
spec.loader.exec_module(experiment)
ROOT = experiment.ROOT


def ari(first, second):
    if first.keys() != second.keys():
        raise ValueError("ARI requires the same molecules.")

    choose2 = lambda value: value * (value - 1) // 2
    cells = Counter((first[name], second[name]) for name in first)
    rows = Counter(first.values())
    columns = Counter(second.values())

    pairs = choose2(len(first))
    if not pairs:
        return 1.0

    observed = sum(map(choose2, cells.values()))
    first_sum = sum(map(choose2, rows.values()))
    second_sum = sum(map(choose2, columns.values()))

    expected = first_sum * second_sum / pairs
    maximum = (first_sum + second_sum) / 2

    if maximum == expected:
        return 1.0

    return (observed - expected) / (maximum - expected)


def save_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    settings = experiment.frozen()["config"]
    rows = []
    errors = []
    memberships = {}

    for method in experiment.METHODS:
        method_rows = []

        for cutoff in settings["cutoffs"]:
            previous = None

            for repeat in range(1, settings["REPEATS"] + 1):
                directory = (
                    ROOT / method / "runs"
                    / f"repeat_{repeat:02d}" / f"cutoff_{cutoff}"
                )

                try:
                    record = json.loads(
                        (directory / "execution.json").read_text()
                    )

                    if record["status"] != "passed":
                        raise ValueError("Execution did not pass.")

                    report = experiment.read_report(
                        directory / "run.out",
                        settings["EXPECTED_MOLECULES"],
                        method,
                    )

                    mapping = report.pop("assignments")
                    sizes = report.pop("sizes")

                    if previous is not None and mapping != previous:
                        raise ValueError("Repeated cluster membership changed.")

                    previous = mapping
                    memberships[method, cutoff] = mapping

                    row = {
                        "method": method,
                        "cutoff": float(cutoff),
                        "repeat": repeat,
                        "molecules": settings["EXPECTED_MOLECULES"],
                        **report,
                        "mean_island_size": statistics.mean(sizes),
                        "median_island_size": statistics.median(sizes),
                        "wall_seconds": record["process_wall_seconds"],
                        "peak_rss_kib": record["peak_rss_kib"],
                        "host": record["host"],
                        "cpu_model": record["cpu_model"],
                        "max_parallel": settings["MAX_PARALLEL"],
                    }

                    row["analysis_seconds_per_comparison"] = (
                        row["analysis_seconds"] / row["comparisons"]
                    )
                    method_rows.append(row)

                except Exception as error:
                    errors.append(
                        f"{method} cutoff {cutoff} repeat {repeat}: {error}"
                    )

        if method_rows:
            save_csv(ROOT / method / "statistics.csv", method_rows)
            rows.extend(method_rows)

    (ROOT / "comparison/errors.txt").write_text(
        "\n".join(errors) + ("\n" if errors else "")
    )

    if errors:
        raise SystemExit(
            f"{len(errors)} missing or failed runs. "
            "Inspect comparison/errors.txt."
        )

    save_csv(ROOT / "comparison/all_runs.csv", rows)
    summary = []

    for method in experiment.METHODS:
        for cutoff in settings["cutoffs"]:
            group = [
                row for row in rows
                if row["method"] == method
                and row["cutoff"] == float(cutoff)
            ]

            summary.append({
                "method": method,
                "cutoff": float(cutoff),
                "repeats": len(group),
                "median_analysis_seconds": statistics.median(
                    row["analysis_seconds"] for row in group
                ),
                "min_analysis_seconds": min(
                    row["analysis_seconds"] for row in group
                ),
                "max_analysis_seconds": max(
                    row["analysis_seconds"] for row in group
                ),
                "median_wall_seconds": statistics.median(
                    row["wall_seconds"] for row in group
                ),
                "median_peak_rss_kib": statistics.median(
                    row["peak_rss_kib"] for row in group
                ),
                "islands": group[0]["islands"],
                "singletons": group[0]["singletons"],
                "singleton_fraction": group[0]["singleton_fraction"],
                "comparisons": group[0]["comparisons"],
                "largest": group[0]["largest"],
                "mean_island_size": group[0]["mean_island_size"],
                "median_island_size": group[0]["median_island_size"],
            })

    save_csv(ROOT / "comparison/summary.csv", summary)
    comparisons = []

    for cutoff in settings["cutoffs"]:
        for first, second in combinations(experiment.METHODS, 2):
            comparisons.append({
                "cutoff": float(cutoff),
                "method_a": first,
                "method_b": second,
                "adjusted_rand_index": ari(
                    memberships[first, cutoff],
                    memberships[second, cutoff],
                ),
            })

    save_csv(ROOT / "comparison/membership_ari.csv", comparisons)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        labels = {
            "cartesian": "Cartesian (323 residues)",
            "multigrid_323": "MultiGrid 323",
            "multigrid_20": "MultiGrid 20",
        }

        plots = [
            (
                "median_analysis_seconds",
                "Total Analysis runtime (s)",
                "runtime_vs_cutoff",
            ),
            ("singletons", "Singleton islands", "singletons_vs_cutoff"),
            (
                "singleton_fraction",
                "Singleton islands / all islands",
                "singleton_fraction_vs_cutoff",
            ),
            (
                "comparisons",
                "FPS pair calculations",
                "comparisons_vs_cutoff",
            ),
        ]

        for field, ylabel, filename in plots:
            figure, axis = plt.subplots(
                figsize=(6.5, 4.2), constrained_layout=True
            )

            for method in experiment.METHODS:
                data = [
                    row for row in summary if row["method"] == method
                ]

                axis.plot(
                    [row["cutoff"] for row in data],
                    [row[field] for row in data],
                    marker="o",
                    markersize=3,
                    label=labels[method],
                )

                if field == "median_analysis_seconds":
                    axis.fill_between(
                        [row["cutoff"] for row in data],
                        [row["min_analysis_seconds"] for row in data],
                        [row["max_analysis_seconds"] for row in data],
                        alpha=0.12,
                    )

            axis.set(
                xlabel="Euclidean total FPS cutoff",
                ylabel=ylabel,
            )
            axis.legend(frameon=False)
            axis.grid(alpha=0.2)

            for extension in ("pdf", "png"):
                figure.savefig(
                    ROOT / "comparison" / f"{filename}.{extension}",
                    dpi=300,
                )

            plt.close(figure)

    except ImportError:
        print(
            "CSV comparisons complete. "
            "Install matplotlib and rerun for figures."
        )

    print(
        f"Compared {len(rows)} successful runs. "
        f"Results: {ROOT / 'comparison'}"
    )


if __name__ == "__main__":
    main()
