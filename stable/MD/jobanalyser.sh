#! /bin/bash

mkdir results

direct=`pwd`

content=(`cat .input.dat`)

JOB_NAME=${content[0]}
NUM_TRY=${content[1]}
JOB_NUM=${content[2]}


for i in $(seq 1 1 $NUM_TRY)
	do
	#Make the index_nw.dat file
	awk '/ P1 /{print $2-1}' $JOB_NAME.pdb > index_nw.dat
	#Make the nw pdb file
	awk '/ P1 /{print}' $JOB_NAME.pdb > $JOB_NAME-nw.pdb
	#Making the commands to join the dcds, clean water and build pdb and run order_par3
	cd ./$i
	all_command="catdcd -o $JOB_NAME-$i-sim-all.dcd -otype dcd -dcd"
	nw_command="catdcd -o $JOB_NAME-$i-sim-all-nw.dcd -i ../index_nw.dat -dcd $JOB_NAME-$i-sim-all.dcd"
	#Make the vmd alignment program
	echo "mol new $direct/$JOB_NAME-nw.pdb type pdb" > ./vmd_script_tmp.tcl
	echo "mol addfile $direct/$i/$JOB_NAME-$i-sim-all-nw.dcd type dcd first 0 last -1 step 1 waitfor -1 molid 0" >> ./vmd_script_tmp.tcl
	echo "set n [molinfo top get frame]" >> ./vmd_script_tmp.tcl
	echo "set ref [atomselect top \"backbone and resid 15 to 43 46 47 53 to 61 64 to 81 85 to 119 126 to 136 140 to 156\" frame 0]" >> ./vmd_script_tmp.tcl
	echo "set cur [atomselect top \"backbone and resid 15 to 43 46 47 53 to 61 64 to 81 85 to 119 126 to 136 140 to 156\" ]" >> ./vmd_script_tmp.tcl
	echo "set prot [atomselect top all]" >> ./vmd_script_tmp.tcl
	echo "for { set i 1 } { \$i < \$n } { incr i } {" >> ./vmd_script_tmp.tcl
	echo "	\$cur frame \$i" >> ./vmd_script_tmp.tcl
	echo "  \$prot frame \$i" >> ./vmd_script_tmp.tcl
	echo "	\$prot move [measure fit \$cur \$ref]" >> ./vmd_script_tmp.tcl
	echo "	}" >> ./vmd_script_tmp.tcl
	echo "animate write dcd {$direct/$i/$JOB_NAME-$i-sim-all-nw-aligned.dcd} beg 0 end -1" >> ./vmd_script_tmp.tcl
	echo "exit" >> ./vmd_script_tmp.tcl
	vmd_command="vmd -dispdev none -e ./vmd_script_tmp.tcl"
	pdb_command="catdcd -o $JOB_NAME-$i-sim-all-nw-aligned.pdb -otype pdb -s ../$JOB_NAME-nw.pdb -dcd $JOB_NAME-$i-sim-all-nw-aligned.dcd"
	order_command="/home/ortega/PhD/programing/C/order_par/order_par_3 $JOB_NAME-$i-sim-all-nw-aligned.pdb > $JOB_NAME-$i-OP.dat"
	for j in $(seq 1 1 $JOB_NUM)
		do
		all_command="$all_command $JOB_NAME-$i-sim-$j.dcd"
		done
	echo "Writing file $i"
	echo "#! /bin/bash" > ./temp$i.sh
	echo "$all_command" >> ./temp$i.sh
	echo "$nw_command" >> ./temp$i.sh
	echo "$vmd_command" >> ./temp$i.sh
	echo "$pdb_command" >> ./temp$i.sh
	echo "$order_command" >> ./temp$i.sh
	echo "mv ./$JOB_NAME-$i-OP.dat ../results" >> ./temp$i.sh
	echo "rm ./temp$i.sh" >> ./temp$i.sh
	chmod +x ./temp$i.sh
	echo "Submitting try number $i on background"
	nohup ./temp$i.sh > ./temp$i.out &
	cd ..
	done
exit 0
