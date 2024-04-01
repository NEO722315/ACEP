#!/bin/bash

if [ ! -d ../02_output ];then
  mkdir ../02_output
fi

cp -rf ../01_data/01_AA_data ../02_output/01_AA_data

python ./01_pamlparas.py \
	  --dataset_location ../02_output/01_AA_data \
    --codeml_str ../01_data/codemlstr.txt \
    --dat_file ../01_data/wag.dat \
    --tree_file ../01_data/abs_mammal.nwk

bash ./02_pamlinfer.sh 1 ../02_output/01_AA_data