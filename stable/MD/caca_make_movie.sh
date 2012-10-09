#! /env/bin/bash

# Sintax: caca_make_movie.sh num_residues DT_MIN DT_MAX DT_STEP (optional) 

N_res=$1
DT_MIN=$2
DT_MAX=$3
if [ $4 ]; then
	DT_STEP=$4
else
	DT_STEP=$((($DT_MAX - $DT_MIN) / 30 ))
	echo "Autoselected step: $DT_STEP"
fi

cd ./gnuplot_image


DT=$DT_MIN
i=0

while [ $DT -le $DT_MAX ]; do
echo "
set pm3d map;
set size square;
set cbrange [-1:1];
set zrange [-1:1];
set title 't=$DT';
set palette define (-1 'blue', 0 'white', 1 'red');
unset key;
set xrange [1:$N_res];
set yrange [1:$N_res];
set cbrange [-1:1];
set terminal png small;
set output 'gnuplot_$DT.png';
splot 'gnuplot_"$DT"_0.dat' with pm3d" | gnuplot
new=$(printf "%06d.png" ${i})
mv gnuplot_$DT.png $new
DT=$(($DT + $DT_STEP))
i=$(( $i + 1 ))
done

convert -delay 20 *.png movie.mp4
echo "Converting to flv"
convert movie.mp4 movie.flv

		
	

