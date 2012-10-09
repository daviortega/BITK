#! /usr/bin/env python 
###################################
#    Davi Ortega 6/18/2012 
###################################
import sys
import numpy
import MDAnalysis as MDA
import MDAnalysis.analysis.align as MDAaa

if '-h' in sys.argv:
	print 'Input pdb file and dcd file and change it will spit out rmsd.\n Options: \n\t-a\tAlign to pdb as reference frame'
	sys.exit()


name_file_dcd = sys.argv[2][:-3]
while "/" in name_file_dcd:
        name_file_dcd = name_file_dcd[name_file_dcd.index("/")+1:]

pdb_file = sys.argv[1]
dcd_file = sys.argv[2]

print "Loading trajectory"
ref = MDA.Universe(pdb_file)
mol = MDA.Universe(pdb_file, dcd_file)

mol_all = mol.selectAtoms('all')
ref_all = ref.selectAtoms('all')

rmsd_list = []

if "-a" in sys.argv:
	dataout = open(name_file_dcd + 'rmsd.A.txt', 'w')
	for frame in mol.trajectory:
		MDAaa.alignto(mol, ref, select = 'backbone', mass_weighted = True)
		rmsd_list.append(MDAaa.rmsd(mol_all.coordinates(),ref_all.coordinates()))
else:
	dataout = open(name_file_dcd + 'rmsd.txt', 'w')
	for frame in mol.trajectory:
                rmsd_list.append(MDAaa.rmsd(mol_all.coordinates(),ref_all.coordinates()))
	
output = ''

for i in range(len(rmsd_list)):
	output += str(i) + '\t' + '%.5f' % rmsd_list[i] + '\n'

dataout.write(output)
dataout.close()

