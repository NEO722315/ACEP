#!/bin/bash

python ./Embedding.py \
      --dataset_location ../02_output/01_AA_data \
      --combine_dir ../02_output/Combinedir \
      --model_location ../03_model/esm_msa1b_t12_100M_UR50S.pt \
      --repr_layers 12 \
      --bottleneck_weight ../04_model/linear1_weight.pt \
      --include bottleneck \
      --case_sps1 Physeter_catodon Lipotes_vexillifer Delphinapterus_leucas Orcinus_orca Tursiops_truncatus \
      --case_sps2 Miniopterus_natalensis Myotis_davidii Myotis_brandtii Myotis_lucifugus Eptesicus_fuscus Hipposideros_armiger Rhinolophus_sinicus Rousettus_aegyptiacus \
      --embed_single True \
      --device cuda:0
