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

if len(sys.argv) != 10:
	print "This is a raw program, it will spit out the RMSD of backbone for each residue (column) for each frame in terms of the reference frame (row).\nScript takes exactly 9 arguments: pdb_file dcd_file Np res_beg res_end time_step skip cutoff reference_frame"

	sys.exit()

pdb_file = sys.argv[1]
dcd_file = sys.argv[2]
Np = int(sys.argv[3])
res_beg = int(sys.argv[4])
res_end = int(sys.argv[5])
step = float(sys.argv[6])
skip = int(sys.argv[7])
cutoff = sys.argv[8]
ref_frame = int(sys.argv[9])


res_num = res_end - res_beg


print "Loading trajectory"
ref = MD.Universe(pdb_file)
mol = MD.Universe(pdb_file, dcd_file)



name_file_dcd = sys.argv[2][:-3]
while "/" in name_file_dcd:
	name_file_dcd = name_file_dcd[name_file_dcd.index("/")+1:]

list_Ct = {}

output = 'residue\trmsf\terr_rmsf'

def _rmsd(res):
	#mol is an md universe
	#res,name_file_dcd,mol=args
	mol = MD.Universe(pdb_file, dcd_file)
	print "["+str(res)+"] Aligning Trajectory"
	ref_sel = ref.selectAtoms('backbone and resid ' + str(res),'backbone and around ' + cutoff + ' (backbone and resid ' + str(res) + ')')
	print ref_sel.atoms
	ref_sel_str = ''
	for i in list(ref_sel):
		useful =str(i).split(":")[1]
		ref_sel_str += '(' + useful.replace('>',' ').replace("'",'').replace(',',' and').replace('of','and') + ') or '
		#ref_sel_str += "(" + useful.replace('>',' ').replace("'",'') + ') or '
	ref_sel_str = ref_sel_str[:-4] 
	#print ref_sel_str
	filenamedcd = "LA" + cutoff + "rmsd.noh.R" + str(res) + '.dcd'
	MDA.rms_fit_trj(mol, ref, select=ref_sel_str, filename = filenamedcd ) #'name N and type N and resname MET resid 1 and segid P1')
	output = ''

	print "["+str(res)+"] Calculating RMSD for residue " + str(res) + "..." + str( 1 + range(res_beg,res_end).index(res) ) + " of " + str(res_num)
	mol = MD.Universe(pdb_file, filenamedcd)
	cf_timer = time.time()
	atoms_mol = mol.selectAtoms("resid " + str(res) + " and not name H")
	mol.trajectory[ref_frame]
	atoms_ref = copy.deepcopy(mol.selectAtoms("resid " + str(res) + " and not name H").coordinates())
	rmsd = []
	for ts1 in mol.trajectory[0:-1:skip]:
		rmsd.append(MDA.rmsd(atoms_mol.coordinates(), atoms_ref))
	
	rmsd = numpy.array(rmsd)
	numpy.savetxt(name_file_dcd + "rmsd.noh.LA" + cutoff +".%04d.txt" % res , rmsd)
	os.remove(filenamedcd)
	
arg_list = []
#for resid in range(3,res_num):
#	arg_list.append([resid,name_file_dcd,mol])
#_rmsd(304)
#sys.exit()


res_list = range(res_beg,res_end)
pool = multip.Pool(processes=Np)
pool.map(_rmsd, res_list)

print "Grouping and cleaning"

rmsd_all = numpy.loadtxt(name_file_dcd + "rmsd.noh.LA" + cutoff +".%04d.txt" % res_list[0] )
os.system('rm ' + name_file_dcd + "rmsd.noh.LA" + cutoff +".%04d.txt" % res_list[0])
for res in res_list[1:]:
	rmsd_all = numpy.vstack((rmsd_all, numpy.loadtxt(name_file_dcd + "rmsd.noh.LA" + cutoff +".%04d.txt" % res)))
	os.system('rm ' + name_file_dcd + "rmsd.noh.LA" + cutoff +".%04d.txt" % res)

numpy.savetxt(name_file_dcd + 'all.rmsd.noh.LA' + cutoff + '.txt', rmsd_all.T)

print "Done... with capital D"

