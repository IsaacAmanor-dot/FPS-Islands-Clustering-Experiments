#!/usr/bin/env python3
"""Compare one cutoff sweep for each of the three FPS methods."""
from pathlib import Path
from collections import Counter
from itertools import combinations
import csv
import importlib.util
import json
import statistics

spec = importlib.util.spec_from_file_location("experiment", Path(__file__).with_name("001.prepare_experiments.py"))
e = importlib.util.module_from_spec(spec)
spec.loader.exec_module(e)
R = e.ROOT


def ari(a, b):
    if a.keys() != b.keys():
        raise ValueError("ARI requires matching molecule records.")
    choose2 = lambda x: x * (x - 1) // 2
    pairs = choose2(len(a))
    if not pairs:
        return 1.0
    cells = Counter((a[name], b[name]) for name in a)
    first = sum(map(choose2, Counter(a.values()).values()))
    second = sum(map(choose2, Counter(b.values()).values()))
    expected = first * second / pairs
    maximum = (first + second) / 2
    return (sum(map(choose2, cells.values())) - expected) / (maximum - expected) if maximum != expected else 1.0


def save_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    c = e.frozen()["config"]
    setup = json.loads((R / "comparison/setup_complete.json").read_text())
    if setup.get("stage") != "ready":
        raise ValueError("Setup is not ready.")
    for name, expected in setup["sha256"].items():
        if e.sha(R / name) != expected:
            raise ValueError(f"Frozen setup file changed: {name}")
    rows, errors, memberships = [], [], {}
    for method in e.METHODS:
        method_rows = []
        for cutoff in c["cutoffs"]:
            directory = R / method / "runs" / f"cutoff_{cutoff}"
            prefix = e.output_prefix(method, cutoff)
            try:
                record = json.loads((directory / "execution.json").read_text())
                if record["status"] != "passed":
                    raise ValueError("Execution did not pass.")
                report = e.read_report(directory / f"{prefix}.out", c["EXPECTED_MOLECULES"], method)
                sizes = report.pop("sizes")
                memberships[method, cutoff] = report.pop("assignments")
                row = dict(method=method, cutoff=float(cutoff), molecules=c["EXPECTED_MOLECULES"],
                           **report, mean_island_size=statistics.mean(sizes),
                           median_island_size=statistics.median(sizes),
                           wall_seconds=record["process_wall_seconds"], peak_rss_kib=record["peak_rss_kib"],
                           host=record["host"], cpu_model=record["cpu_model"],
                           max_parallel=record.get("concurrent_cutoffs", c["MAX_PARALLEL"]),
                           launcher=record.get("launcher", "direct"))
                row["analysis_seconds_per_comparison"] = row["analysis_seconds"] / row["comparisons"]
                method_rows.append(row)
            except Exception as error:
                errors.append(f"{method} cutoff {cutoff}: {error}")
        if method_rows:
            save_csv(R / method / "statistics.csv", method_rows)
            rows.extend(method_rows)
    (R / "comparison/errors.txt").write_text("\n".join(errors) + ("\n" if errors else ""))
    if errors:
        raise SystemExit(f"{len(errors)} missing or failed runs; inspect comparison/errors.txt.")
    save_csv(R / "comparison/all_runs.csv", rows)
    save_csv(R / "comparison/summary.csv", rows)
    comparisons = []
    for cutoff in c["cutoffs"]:
        for a, b in combinations(e.METHODS, 2):
            comparisons.append(dict(cutoff=float(cutoff), method_a=a, method_b=b,
                                    adjusted_rand_index=ari(memberships[a, cutoff], memberships[b, cutoff])))
    save_csv(R / "comparison/membership_ari.csv", comparisons)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        labels = dict(cartesian="Cartesian (323 residues)", multigrid_323="MultiGrid 323", multigrid_20="MultiGrid 20")
        for field, ylabel, name in [
            ("analysis_seconds", "Total Analysis runtime (s)", "runtime_vs_cutoff"),
            ("singletons", "Singleton islands", "singletons_vs_cutoff"),
            ("singleton_fraction", "Singleton islands / all islands", "singleton_fraction_vs_cutoff"),
            ("comparisons", "FPS pair calculations", "comparisons_vs_cutoff")]:
            fig, ax = plt.subplots(figsize=(6.5, 4.2), constrained_layout=True)
            for method in e.METHODS:
                data = [row for row in rows if row["method"] == method]
                ax.plot([row["cutoff"] for row in data], [row[field] for row in data],
                        marker="o", markersize=3, label=labels[method])
            ax.set(xlabel="Euclidean total FPS cutoff", ylabel=ylabel)
            ax.legend(frameon=False)
            ax.grid(alpha=0.2)
            for extension in ("pdf", "png"):
                fig.savefig(R / "comparison" / f"{name}.{extension}", dpi=300)
            plt.close(fig)
    except ImportError:
        print("CSV files complete. Install matplotlib and rerun for figures.")
    print(f"Compared {len(rows)} runs: one run per method and cutoff.")


if __name__ == "__main__":
    main()
