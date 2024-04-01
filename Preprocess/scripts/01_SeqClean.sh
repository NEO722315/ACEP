#!/bin/bash

python ./01_SeqClean.py \
	  --dataset_location ../02_output \
	  --rawAAdir ../01_data/aa_genes \
	  --rawCDSdir ../01_data/nt_genes \
	  --case_sps1 Physeter_catodon Lipotes_vexillifer Delphinapterus_leucas Orcinus_orca Tursiops_truncatus \
	  --case_sps2 Miniopterus_natalensis Myotis_davidii Myotis_brandtii Myotis_lucifugus Eptesicus_fuscus Hipposideros_armiger Rhinolophus_sinicus Rousettus_aegyptiacus

# --case_sps1 Physeter_catodon Lipotes_vexillifer Delphinapterus_leucas Orcinus_orca Tursiops_truncatus \
# --case_sps2 Miniopterus_natalensis Myotis_davidii Myotis_brandtii Myotis_lucifugus Eptesicus_fuscus Hipposideros_armiger Rhinolophus_sinicus Rousettus_aegyptiacus \