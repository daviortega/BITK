#! /env/bin/bash

# Sintax: caca_movie.sh input_file.pdb num_residues number_of_processors DT_MIN DT_MAX DT_STEP (optional) 

infile=$1
N_res=$2
N=$3
DT_MIN=$4
DT_MAX=$5
if [ $6 ]; then
	DT_STEP=$6
else
	DT_STEP=$((($DT_MAX - $DT_MIN) / 30 ))
	echo "Autoselected step: $DT_STEP"
fi

i=0
n=0
DT=$DT_MIN

while [ $DT -le $DT_MAX ]; do
	while [ $n -lt $N -a $DT -le $DT_MAX ]; do
		echo "Running job $i/$((($DT_MAX - $DT_MIN) / $DT_STEP)) with $DT delay"
		python ./ca_ca_corr_matrix.py $DT 1 $infile -gnuplot > output_$i.log &
		DT=$(($DT + $DT_STEP))
		i=$(( $i + 1 ))
		n=$(( $n + 1 ))
	done
	wait
	n=0
done


sh caca_make_movie.sh $2 $4 $5 $6


