#!/bin/bash

# Edit paths before preparing the experiment.
DOCK_ROOT="/gpfs/projects/rizzo/iamanor/DOCK6_Development/Similarity_Island/dock6_dev"
DOCK_BIN="${DOCK_ROOT}/bin/dock6"
GRID_BIN="${DOCK_ROOT}/bin/grid"

RECEPTOR_FILE="/gpfs/projects/AMS536/2026/group3/tutorial_isaac/Individual_Project/001_structure/2ITY_protein_wH_wcH.mol2"
BOX_FILE="/gpfs/projects/AMS536/2026/group3/tutorial_isaac/Individual_Project/003_gridbox/2ITY.box.pdb"

LIGAND_FILE="/gpfs/projects/rizzo/iamanor/DOCK6_Development/Similarity_Island/TesTING_SMI_2/000_VS_For_DOCKING_SCORE/Dock_Scored_Molecules_for_SMI_Clustering_scored.mol2"

FPS_CUTOFFS=(0 1 2 3 4 5 6 7 8 9)

MAX_PARALLEL=${#FPS_CUTOFFS[@]}

EXPECTED_MOLECULES=5282
EXPECTED_RESIDUES=323
VALIDATION_MOLECULES=10
