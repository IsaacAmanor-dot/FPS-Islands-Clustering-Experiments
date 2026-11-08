# FPS Islands clustering experiments

Compare Cartesian FPS, 323-component MultiGrid FPS and 20-component MultiGrid
FPS for the same ordered 2ITY ligand library.

## Experiments

- `cartesian/`: direct VDW, ES and HB footprints for all 323 residues.
- `multigrid_323/`: 322 individual residues plus a one-residue remainder.
- `multigrid_20/`: 19 selected residues plus one 304-residue remainder.

Each method runs cutoffs 0–9 once: 30 runs total, with no repeated sweeps.
A three-task SLURM array assigns one exclusive node to each method.
All 10 cutoffs run concurrently within each method.
The `%a` in experiment logs identifies the array task; no job-ID suffix is used.

## Scoring and output

Every generated DOCK input enables `score_molecules yes`,
`descriptor_score_primary yes`, and `descriptor_use_grid_score yes`.
Direct `grid_score_primary` stays `no`. Descriptor grid scales and weight are
one, using the same `inputs/whole_grid/whole` prefix for all three methods.
FPS comparisons use Euclidean `total = dVDW + dES + dHB`, with retention `all`.
MultiGrid interpolates VDW/ES and calculates geometric HB in Cartesian space.

A Cartesian cutoff-10 run is written under `cartesian/runs/cutoff_5/`:

- `FPS_Island_cartesian_5.in`
- `FPS_Island_cartesian_5.out`
- `FPS_Island_cartesian_5_clustered.mol2`
- `FPS_Island_cartesian_5_0001.dat`, and one `.dat` per additional island.

The other prefixes are `FPS_Island_multigrid_323_10` and
`FPS_Island_multigrid_20_10`. All benchmark inputs use `individual_islands`.
Internal selection/validation exports use combined files for parsing, while
retaining descriptor rescoring and experiment-specific output names.
Execution checks completion, rescored record counts and footprint file presence.

## Run

Edit paths in `000.config.sh`. Use Python 3.9+, Bash, GNU time, SLURM and the
HB-enabled DOCK installation. Matplotlib is optional for figures.
For a fresh experiment, run from this directory:

```bash
python3 001.prepare_experiments.py prepare
sbatch 001a.generate_multigrids.slurm
```

The setup script generates new grids or validates an existing complete grid set.
It freezes the residue selection, validates HB mapping and writes benchmark inputs.
After the setup job finishes successfully and reports 30 ready runs:

```bash
sbatch 002.run_experiments.slurm
```

When every array task finishes:

```bash
python3 003.compare_results.py
```

Do not rerun `prepare` on prepared inputs. After an archived-run reset, submit
`001a.generate_multigrids.slurm` directly; existing grids are reused.
If changing concurrency, request the same number of CPUs in `002`.
Use the same concurrency for all methods. A 323-component process needs roughly
10 GiB for grid energy data, plus other program data; account for concurrent processes.

## Reproducibility

All 5,282 molecule records are retained in original order, including repeated
original titles. Original molecule titles are preserved. Record positions are stored
in `inputs/ligand_records.tsv`. The original library is unchanged.

Residue selection uses Cartesian footprints across the full fixed library.
For each residue, calculate mean absolute VDW, ES and HB contributions.
Normalize each component by its sum across residues, treating a zero denominator
as zero. Sum the three normalized contributions with equal weights, select the
19 highest scores, and break ties by increasing MOL2 residue ID. Save the complete
ranking and remaining-residue list under `selection/`. Selection uses DOCK's
six-decimal export and is specific to this library, not an independent test set.

Grids share the box, 0.3 Å spacing, all-atom 6–9 VDW model, distance-dependent
dielectric 4 and cutoff distance 9999 Å. Short relative GRID filenames avoid its
legacy buffer limit. Grid generation follows the native workflow at
https://github.com/IsaacAmanor-dot/Mult-Grid-Generation-in-Dock6 .

Revalidation checks grid headers, geometry, finite values, partitions and HB maps.
It recomputes scored Cartesian selection and requires the same selected residues
before reusing the 20-grid set. Both MultiGrid mappings must reproduce Cartesian
HB counts on the validation subset; VDW/ES differences are reported separately.
All prepared data, executable snapshots and benchmark inputs are checksummed.

## Results

Each method writes `statistics.csv`. `comparison/` contains all runs, summary,
membership ARI, validation reports, setup costs and runtime/singleton figures.
Singleton fraction means singleton islands divided by all islands.

Runtime is total Analysis time, including initialization, descriptor rescoring
and file output. External wall time and peak RSS are also recorded.
There is one observation per method/cutoff; no repeat variability is estimated.
Hosts and CPU models are recorded. Concurrent workload and hardware differences
can affect timings. Aggregation/interpolation can change membership and pair
counts, so equal cutoffs do not guarantee equal clustering work or resolution.
Grid generation and selection costs are reported separately.

Track scripts, inputs, selection records, manifests and final statistics/figures.
Generated grids, copied executables, raw runs and archives remain ignored.
Preserve raw data locally for auditing; Git LFS is unnecessary.

Submitted inputs are frozen as *.in.requested. DOCK runs on the matching
*.in working copy, which it may rewrite with resolved parameters.
Checksums verify submitted snapshots, not mutable runtime input files.

## Concurrent cutoff allocation

Each method runs cutoffs 0 through 9 once: 30 benchmark runs total.
Each method receives one exclusive rn-long-40core node with 10 tasks,
one CPU per task, and the node's full memory allocation (--mem=0).
All 10 cutoffs run concurrently through separate exclusive srun steps.
There are no repeated sweeps.

Memory checks remain enabled. Each process loads its own grid arrays.
Previous results are archived; grids and residue selection are preserved.
Submit 002.run_experiments.slurm after applying the updater.
Preparation and grid generation do not need to be repeated.


## Original molecule names

Original ZINC names and all duplicate molecule records are retained.
No __SI suffixes are added to molecule titles.
Runtime, island counts, singleton counts and footprint statistics remain available.
Name-based membership ARI is left blank when duplicate titles prevent
unambiguous matching between experiments.

Removing earlier name suffixes is a text-only cleanup; recorded timings,
coordinates, scores and grid files are unchanged.
