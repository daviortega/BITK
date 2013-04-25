#! /usr/bin/env python 
###################################
#    Davi Ortega 5/9/2012 
###################################
import os
import sys
import copy
import MDAnalysis as MD
import MDAnalysis.analysis.align as MDA
import MDAnalysis.selections.base as MDS
import multiprocessing as multip
import numpy
import scipy.optimize
import matplotlib.pylab as mp
import time
import glob
from numpy.core.umath_tests import inner1d

if len(sys.argv) != 6:
	print "This is a raw program, it will spit out the RMSD matrix for each layer (+- 4) for MCP trajectory file.\nScript takes exactly 4 arguments: pdb_file dcd_file Np skip res_tip"

	sys.exit()

pdb_file = sys.argv[1]
dcd_file = sys.argv[2]
Np = int(sys.argv[3])
skip = int(sys.argv[4])
res_tip = int(sys.argv[5])


#print "Loading trajectory"
#mol = MD.Universe(pdb_file, dcd_file)


name_file_dcd = sys.argv[2][:-3]
while "/" in name_file_dcd:
	name_file_dcd = name_file_dcd[name_file_dcd.index("/")+1:]

def _rmsd(layer):
	#mol is an md universe
	print("Loading Trajectory 2x")
	beg = time.time()
	mol_ref = MD.Universe(pdb_file, dcd_file)
	mol_trg = MD.Universe(pdb_file, dcd_file)
	#res,name_file_dcd,mol=args
	print(' This took ' + str(time.time()-beg) + ' seconds')

	#### building selection ####
	res_list = range(res_tip - layer - 4, res_tip - layer + 4)
	res_list += range(res_tip + layer - 4, res_tip + layer + 4)
	sel_ref = mol_ref.selectAtoms('backbone and resid ' + str(res_list[0]))
	sel_trg = mol_trg.selectAtoms('backbone and resid ' + str(res_list[0]))
	for res in res_list[1:]:
		sel_ref += mol_ref.selectAtoms('backbone and resid ' + str(res))
		sel_trg += mol_trg.selectAtoms('backbone and resid ' + str(res))
	
	############################


	rmsd_matrix = []

	for f1 in range(0,mol_ref.trajectory.numframes,skip):
		print('working on frame ' + str(f1))
		mol_ref.trajectory[f1]
#		ref = sel_ref.coordinates()
#		ref = ref.atoms.CA.coordiantes() - ref.atoms.CA.CenterOfMass()
		rmsd = []
		for f2 in range(f1+1,mol_trg.trajectory.numframes, skip):
			mol_trg.trajectory[f2]
			rmsd.append(MDA.alignto(sel_trg, sel_ref)[1])
		rmsd_matrix.append(rmsd)

	for rmsd in rmsd_matrix:
		print rmsd
			
			

#	mol = MD.Universe(pdb_file, dcd_file)
#	print "["+str(res)+"] Aligning Trajectory"
#	ref_sel = ref.selectAtoms('backbone and resid ' + str(res),'backbone and around ' + cutoff + ' (backbone and resid ' + str(res) + ')')
#	print ref_sel.atoms
#	ref_sel_str = ''
#	for i in list(ref_sel):
#		useful =str(i).split(":")[1]
#		ref_sel_str += '(' + useful.replace('>',' ').replace("'",'').replace(',',' and').replace('of','and') + ') or '
#		#ref_sel_str += "(" + useful.replace('>',' ').replace("'",'') + ') or '
#	ref_sel_str = ref_sel_str[:-4] 
#	#print ref_sel_str
#	filenamedcd = "LA" + cutoff + "rmsd.bb.R" + str(res) + '.dcd'
#	MDA.rms_fit_trj(mol, ref, select=ref_sel_str, filename = filenamedcd ) #'name N and type N and resname MET resid 1 and segid P1')
#	output = ''

#	print "["+str(res)+"] Calculating RMSD for residue " + str(res) + "..." + str( 1 + range(res_beg,res_end).index(res) ) + " of " + str(res_num)
#	mol = MD.Universe(pdb_file, filenamedcd)
#	cf_timer = time.time()
#	atoms_mol = mol.selectAtoms("resid " + str(res) + " and backbone ")
#	mol.trajectory[ref_frame]
#	atoms_ref = copy.deepcopy(mol.selectAtoms("resid " + str(res) + " and backbone ").coordinates())
#	rmsd = []
#	for ts1 in mol.trajectory[0:-1:skip]:
#		rmsd.append(MDA.rmsd(atoms_mol.coordinates(), atoms_ref))
#	
#	rmsd = numpy.array(rmsd)
#	numpy.savetxt(name_file_dcd + "rmsd.bb.LA" + cutoff +".%04d.txt" % res , rmsd)
#	os.remove(filenamedcd)
	
arg_list = []
#for resid in range(3,res_num):
#	arg_list.append([resid,name_file_dcd,mol])
_rmsd(80)
sys.exit()


res_list = range(res_beg,res_end)
pool = multip.Pool(processes=Np)
pool.map(_rmsd, res_list)

print "Grouping and cleaning"

rmsd_all = numpy.loadtxt(name_file_dcd + "rmsd.bb.LA" + cutoff +".%04d.txt" % res_list[0] )
os.system('rm ' + name_file_dcd + "rmsd.bb.LA" + cutoff +".%04d.txt" % res_list[0])
for res in res_list[1:]:
	rmsd_all = numpy.vstack((rmsd_all, numpy.loadtxt(name_file_dcd + "rmsd.bb.LA" + cutoff +".%04d.txt" % res)))
	os.system('rm ' + name_file_dcd + "rmsd.bb.LA" + cutoff +".%04d.txt" % res)

numpy.savetxt(name_file_dcd + 'all.rmsd.bb.LA' + cutoff + '.txt', rmsd_all.T)

print "Done... with capital D"

