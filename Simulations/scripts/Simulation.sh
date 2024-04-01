#!/bin/bash

tmp_fifofile="/tmp/$$.fifo"
mkfifo $tmp_fifofile
exec 6<>$tmp_fifofile
rm $tmp_fifofile

start_time=`date +%s`

thread_num=$1 # 定义最大线程数

for ((i=0;i<${thread_num};i++));do
    echo
done >&6

datapath=$2

for seqid in `ls $datapath`
do
	echo $seqid
	  # cpu simulate and combine process
  read -u6
       {
        python ./Simulation.py \
            --dataset_location ../02_output/01_AA_data \
            --codeml_str ../01_data/codemlstr.txt \
            --freq_mode gene \
            --dat_file ../01_data/wag.dat \
            --tree_file ../01_data/abs_mammal.nwk \
            --simtime 100 \         # default value is 100
            --seqid $seqid
        echo >&6
       } &

done

wait
stop_time=`date +%s` # 定义脚本运行的结束时间
echo "TIME:`expr $stop_time - $start_time`" # 输出脚本运行时间

exec 6>&- # 关闭FD6
echo "over" # 表示脚本运行结束
