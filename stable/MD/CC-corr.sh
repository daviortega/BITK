#! /bin/bash

mkdir results

direct=`pwd`

content=(`cat .input.dat`)

JOB_NAME=${content[0]}
NUM_TRY=${content[1]}
JOB_NUM=${content[2]}


for i in $(seq 1 1 $NUM_TRY)
	do
	cd ./$i
	echo $JOB_NAME-$i-sim-all-nw-aligned.dcd
	ca_ca_corr_matrix 0 1 $JOB_NAME-$i-sim-all-nw-aligned.pdb
	wait
	cp caca_corr_DT_0_0.dat ../results/$JOB_NAME-$i-calcal-corr.dat
	cd ..
	done
exit 0
