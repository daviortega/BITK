#!/bin/bash
#########################
#  by Davi Ortega 2010  #
# davi.ortega@gmail.com #
#########################


if [ ! -f ./$2.pdb ] || [ ! -f ./$2.psf ]  || [ ! -f ./par_all27_prot_lipid_na.inp ]
	then
	echo "The directory must have the following files: ./$2.pdb ./$2.psf par_all27_prot_lipid_na.inp"
	exit 1
fi

echo $2 $1 $3 > ./.input.dat

NUM_TRY=$1
JOB_NAME=$2
NUM_JOBS=$3

NP=256 #Number of processors

echo $NUM_TRY
for i in $(seq 1 1 $NUM_TRY)
do
	mkdir ./$i
	cp ./$2.pdb ./$i/$2-$i.pdb
	cp ./$2.psf ./$i/$2-$i.psf
	cp ./*.inp ./$i/

	for j in $(seq 1 1 $NUM_JOBS)
	do
		cat > ./$i/$2-$i-$j.conf << DELIM
#############################################################
## JOB DESCRIPTION                                         ##
#############################################################

# Production of the mutant $2-$i
# 298K, 30mM NaCl and water box.


#############################################################
## ADJUSTABLE PARAMETERS                                   ##
#############################################################

structure          ./$2-$i.psf
coordinates        ./$2-$i.pdb

set temperature    298
set outputname     $2-$i-sim-$j

DELIM
		
		if [ $j -eq 1 ]; then
			cat >> ./$i/$2-$i-$j.conf << DELIM
firsttimestep      0
#############################################################
## SIMULATION PARAMETERS                                   ##
#############################################################

# Input
paraTypeCharmm      on
parameters          ./par_all27_prot_lipid_na.inp
temperature         \$temperature

DELIM
		else
			cat >> ./$i/$2-$i-$j.conf << DELIM
if {1} {
set inputname      $2-$i-sim-$(($j-1))
binCoordinates     \$inputname.restart.coor
binVelocities      \$inputname.restart.vel  ;# remove the "temperature" entry if you use this!
extendedSystem     \$inputname.restart.xsc

proc get_first_ts { xscfile } {
set fd [open \$xscfile r]
gets \$fd
gets \$fd
gets \$fd line
set ts [lindex \$line 0]
close \$fd
return \$ts
}

set firsttime [get_first_ts \$inputname.restart.xsc]
}

firsttimestep      \$firsttime


#############################################################
## SIMULATION PARAMETERS                                   ##
#############################################################

# Input
paraTypeCharmm      on
parameters          ./par_all27_prot_lipid_na.inp

DELIM
		fi
		cat >> ./$i/$2-$i-$j.conf << DELIM
# Force-Field Parameters
exclude             scaled1-4
1-4scaling          1.0
cutoff              12.
switching           on
switchdist          10.
pairlistdist        13.5


# Integrator Parameters
timestep            2.0  ;# 2fs/step
rigidBonds          all  ;# needed for 2fs steps
nonbondedFreq       1
fullElectFrequency  2
stepspercycle       10


# Constant Temperature Control
langevin            on    ;# do langevin dynamics
langevinDamping     5     ;# damping coefficient (gamma) of 5/ps
langevinTemp        \$temperature
langevinHydrogen    off    ;# don't couple langevin bath to hydrogens


#Periodic Boundary Conditions
cellBasisVector1    64.    0.   0.
cellBasisVector2     0.   91.   0.
cellBasisVector3     0.    0   70.
cellOrigin          0.78  6.29  -1.31

wrapAll             on


# PME (for full-system periodic electrostatics)
PME                 yes

PMEGridSizeX 128
PMEGridSizeY 128
PMEGridSizeZ 128


# Constant Pressure Control (variable volume)
useGroupPressure      yes ;# needed for rigidBonds
useFlexibleCell       no
useConstantArea       no

langevinPiston        on
langevinPistonTarget  1.01325 ;#  in bar -> 1 atm
langevinPistonPeriod  100.
langevinPistonDecay   50.
langevinPistonTemp    \$temperature


# Output
outputName          \$outputname

restartfreq         1000    ;# 1000steps = every 2ps
dcdfreq             1000
xstFreq             1000
outputEnergies      100
outputPressure      100



#############################################################
## EXTRA PARAMETERS                                        ##
#############################################################


#############################################################
## EXECUTION SCRIPT                                        ##
#############################################################

DELIM
		
		if [ $j -eq 1 ]; then
			cat >> ./$i/$2-$i-$j.conf << DELIM
# Minimization
reinitvels          \$temperature
run  45000000 ;# 90ns
DELIM
		else
			cat >> ./$i/$2-$i-$j.conf << DELIM
# Minimization
run  [expr 45000000 - \$firsttime] ;# 90ns
DELIM
		fi

	cat > ./$i/$2-$i-$j.sge << DELIM
#$ -N ${2:0:6}-$i-$j
#$ -cwd
#$ -j y
#$ -q medium_phi


mpirun -np \$NSLOTS -mca btl_openib_ib_timeout 30 /data/apps/NAMD/2.6-intel-openmpi/namd2 ./$i/$2-$i-$j.conf
DELIM

	if [ $j -eq 1 ]; then
		qsub -pe openmpi* $NP ././$i/$2-$i-$j.sge
	else
		qsub -pe openmpi* $NP -hold_jid ${2:1:5}-$i-$(($j-1)) ././$i/$2-$i-$j.sge
	fi

	done
done





