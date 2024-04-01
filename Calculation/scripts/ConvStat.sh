#!/bin/bash

python ./ConvStatistic.py \
        --dataset_location ../02_output/01_AA_data \
        --case_sps1 Physeter_catodon Lipotes_vexillifer Delphinapterus_leucas Orcinus_orca Tursiops_truncatus \
        --case_sps2 Miniopterus_natalensis Myotis_davidii Myotis_brandtii Myotis_lucifugus Eptesicus_fuscus Hipposideros_armiger Rhinolophus_sinicus Rousettus_aegyptiacus \
        --include bottleneck \
        --simtime 100 \
        --embed_single True \
        --repr_layers 12