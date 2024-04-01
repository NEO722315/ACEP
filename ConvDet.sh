#!/bin/bash

python ./ConvDet_main.py \
      --abspath /home/czq/Projects/ACEP \
      --dataset_location ./Test_data/Test_data_output/clean \
      --rawAAdir ./Test_data/Test_data_input/raw/aa_genes \
	    --rawCDSdir ./Test_data/Test_data_input/raw/nt_genes \
	    --codeml_str ./PhyloInfer/01_data/codemlstr.txt \
      --tree_file ./PhyloInfer/01_data/abs_mammal.nwk \
      --dat_file ./PhyloInfer/01_data/wag.dat \
      --combine_dir ./Test_data/Test_data_output/CombineFiles \
      --freq_mode gene \
      --model_location ./Embeddings/03_model/esm_msa1b_t12_100M_UR50S.pt \
      --repr_layers 12 \
      --bottleneck_weight ./Embeddings/03_model/linear1_weight.pt \
      --include bottleneck \
      --device cuda:0 \
      --case_sps1 Physeter_catodon Lipotes_vexillifer Delphinapterus_leucas Orcinus_orca Tursiops_truncatus \
      --case_sps2 Miniopterus_natalensis Myotis_davidii Myotis_brandtii Myotis_lucifugus Eptesicus_fuscus Hipposideros_armiger Rhinolophus_sinicus Rousettus_aegyptiacus \
      --embed_single True \
      --simtime 100