# 2ITY Similarity Island FPS comparison

Compare Cartesian FPS, 323-component MultiGrid FPS and 20-component MultiGrid
FPS using the same ordered library, receptor, executable and cutoff schedule.

## Methods

- `cartesian/`: direct VDW, ES and HB footprints for 323 residues.
- `multigrid_323/`: 322 individual residue grids plus a one-residue remainder.
- `multigrid_20/`: 19 selected residue grids plus a 304-residue remainder.

All methods use Euclidean `total = dVDW + dES + dHB`, unit weights,
retention `all`, and `score_molecules no`.

MultiGrid interpolates VDW/ES and calculates geometric HB interactions in
Cartesian space, grouped into the corresponding grid components.
The total is a sum of component distances, not one concatenated-vector norm.

## Files

- `000.config.sh`: paths, cutoffs, repeats and concurrency.
- `001.prepare_experiments.py`: freeze inputs, prepare runs and execute a method.
- `001a.generate_multigrids.slurm`: submit grid generation and validation.
- `002.run_experiments.slurm`: three array tasks, one exclusive node per method.
- `003.compare_results.py`: statistics, membership comparisons and figures.
- `004.generate_multigrids.py`: grid generation, selection and HB validation.
- `inputs/`: copied molecular inputs, parameters, checksums and local tools.
- `selection/`: residue ranking, selected/remainder lists and component mappings.
- `comparison/`: setup times, validation results, aggregate CSVs and figures.

## Requirements

Python 3.9 or newer, Bash, GNU time at `/usr/bin/time`, SLURM, and the
HB-enabled DOCK executable are required. Matplotlib is optional for figures.

The native `grid`, `multigrid_fp_gen.py` and `mol2.py` utilities must be
installed under the configured DOCK `bin/` directory.

Use a working directory without whitespace and allow at least 15–20 GiB
of free disk space for generated grids and run records.

## Run

Edit `000.config.sh` first. Run all commands from this directory.

Prepare the fixed inputs:

```bash
python3 001.prepare_experiments.py prepare
```

Submit grid generation, residue selection and validation:

```bash
sbatch 001a.generate_multigrids.slurm
```

Monitor `comparison/grid_generation.out` and
`comparison/grid_generation.err`.

Wait for successful completion and `comparison/setup_complete.json`.
Then submit the three-method benchmark:

```bash
sbatch 002.run_experiments.slurm
```

When all three array tasks finish:

```bash
python3 003.compare_results.py
```

Defaults: cutoffs 0–25 in unit steps, three repeats and two concurrent
cutoffs per method, giving 234 benchmark runs.

If changing `MAX_PARALLEL`, also change `--cpus-per-task` before preparation.
Use the same concurrency for every method.

The array requests one exclusive node and all node memory per method.
Each 323-component process loads roughly 10 GiB of energy grids, plus
ligands and HB data. Adjust partition and time to local scheduling policy.

## Residue selection

A separate Cartesian export calculates footprints for the fixed library.

For each residue:

1. Calculate mean absolute VDW, ES and HB contributions across all molecules.
2. Normalize each component by its sum across all receptor residues.
3. Add the three normalized values with equal weights.
4. Select the top 19 residues; break ties by increasing MOL2 residue ID.

A component with a zero normalization denominator contributes zero.
Selection uses DOCK's six-decimal footprint export. Selected IDs are then
ordered numerically for the native generator. All remaining residues form
component 20.

The complete ranking and selected/remainder lists are saved under
`selection/` and remain fixed across cutoffs and repeats.

Selection uses the same library being clustered. This is a library-specific
representation experiment, not an independent test-set evaluation.

## Grid generation and validation

Settings: common box, 0.3 Å spacing, all-atom model, 6–9 VDW potential,
distance-dependent dielectric 4, and energy cutoff distance 9999 Å.

Short relative GRID filenames avoid its legacy 81-byte filename buffers.
A local checked copy of the native generator stops when a command fails.

Generation follows the validated workflow at:
https://github.com/IsaacAmanor-dot/Mult-Grid-Generation-in-Dock6

Utilities are copied from the configured installation and checksummed.

Setup verifies receptor partitions, atomic coordinates/charges, grid
headers, common geometry, finite values and file checksums.

The first ten ligands must reproduce Cartesian HB counts exactly for both
MultiGrid mappings. VDW/ES interpolation differences are reported separately.
Exact Cartesian identity of interpolated energies is not required.

Preparation and execution refuse to overwrite completed inputs or run records.
Use a fresh experiment directory for a changed configuration or interrupted
setup, preserving failed logs for diagnosis.

## Results

Each method writes `statistics.csv`.

`comparison/` contains:

- `all_runs.csv`: individual repeat results.
- `summary.csv`: median runtime, repeat range and clustering statistics.
- `membership_ari.csv`: pairwise adjusted Rand indices.
- `setup_times.json`: grid generation and selection costs.
- `footprint_validation.json`: HB checks and VDW/ES differences.
- Runtime, singleton, singleton-fraction and comparison-count figures.

Singleton fraction means singleton islands divided by all islands.
Repeated molecule titles are supported. Every input record is retained.
Only the experiment copy receives unique record identifiers;
inputs/ligand_records.tsv maps them to the original titles.
Repeated membership must be identical within each method and cutoff.

DOCK's reported time is total Analysis runtime, including initialization,
input and output. It is not an isolated clustering timer.
External process wall time and peak RSS are also recorded.

Measurements use the configured concurrent workload. Hosts and CPU models
are recorded. Repeats on one node measure within-node variation and do not
remove hardware differences between methods. A stricter follow-up should
rotate methods across identical nodes or run them sequentially on the same
reserved node with concurrency one.

Residue aggregation and grid interpolation change footprint distances.
Equal cutoffs therefore need not produce equal memberships or calculation
counts. Interpret runtime alongside cluster statistics and comparison counts.
Generation and selection costs are reported separately.

## Git

Track scripts, molecular/parameter inputs, selection records, manifests and
final statistics/figures.

Generated grids, raw runs, temporary footprint exports and copied software
are ignored. Grid checksums are retained without uploading grid binaries.
No Git LFS is required. Preserve raw runs locally for auditing.

## Footprint output

Benchmark runs use `individual_islands`, writing one file per island,
such as `footprints_0001.dat`, inside each cutoff and repeat directory.
Execution checks that every island has a nonempty footprint file.

Internal residue selection and numerical validation use combined files
because their parsers require that format.

Reported runtime includes footprint export and file writing.
Raw run directories and their footprint files remain excluded from Git.
