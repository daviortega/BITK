#! /usr/bin/env python 
###################################
#    Davi Ortega 5/9/2012 
###################################
import os
import sys
import MDAnalysis as MD
import MDAnalysis.analysis.align as MDA
import MDAnalysis.selections.base as MDS
import multiprocessing as multip
import numpy
import rpy2
import scipy.optimize
import matplotlib.pylab as mp
import time
from numpy.core.umath_tests import inner1d

if len(sys.argv) != 10:
	print "Spits out the B-factor of Calpha for each residue.\nScript takes exactly 9 arguments: pdb_file dcd_file Np res_beg res_end chain time_step skip cutoff"
	sys.exit()

pdb_file = sys.argv[1]
dcd_file = sys.argv[2]
Np = int(sys.argv[3])
res_beg = int(sys.argv[4])
res_end = int(sys.argv[5])
chain = sys.argv[6]
step = float(sys.argv[7])
skip = int(sys.argv[8])
cutoff = sys.argv[9]
res_num = res_end - res_beg

print "Loading trajectory"
mol = MD.Universe(pdb_file, dcd_file)



name_file_dcd = sys.argv[2][:-3]
while "/" in name_file_dcd:
	name_file_dcd = name_file_dcd[name_file_dcd.index("/")+1:]


list_Ct = {}

output = 'residue\tb-factor\terr_b-factor'

def _bfactor(res):
	#mol is an md universe
#	res,name_file_dcd,mol=args
	ref = MD.Universe(pdb_file)
	mol = MD.Universe(pdb_file, dcd_file)
	if 1 == 1:
		print "["+str(res)+"] Aligning Trajectory"
		ref_sel = ref.selectAtoms('backbone and resid ' + str(res),'backbone and around ' + cutoff + ' (backbone and resid ' + str(res) + ')')
		print ref_sel.atoms
		for i in list(ref_sel):
			print i
	#	print ref_sel.indices
		ref_sel_str = ''
		for i in list(ref_sel):
			useful =str(i).split(":")[1]
			ref_sel_str += '(' + useful.replace('>',' ').replace("'",'').replace(',',' and').replace('of','and') + ') or '
			#ref_sel_str += "(" + useful.replace('>',' ').replace("'",'') + ') or '
		ref_sel_str = ref_sel_str[:-4] 
	#	print ref_sel_str
		filenamedcd = "LA" + cutoff + "bf.R" + str(res) + '.' + chain + '.dcd'
		MDA.rms_fit_trj(mol, ref, select=ref_sel_str, filename = filenamedcd ) #'name N and type N and resname MET resid 1 and segid P1')
		output = ''
		print "["+str(res)+"] Calculating C(t) for residue " + str(res) + " of " + str(res_num)
		mol = MD.Universe(pdb_file, filenamedcd)
	cf_timer = time.time()
	atoms_mol = mol.selectAtoms("resid " + str(res) + " and name CA and segid " +str(chain))
	atoms_ref = ref.selectAtoms("resid " + str(res) + " and name CA and segid " +str(chain))


	if 1 == 1:
		var = []
		
		for ts1 in mol.trajectory[0:-1:skip]:
			
			
			var.append((numpy.linalg.norm(numpy.subtract(atoms_mol.coordinates(), atoms_ref.coordinates())))**2)
		
		var = numpy.array(var)
		rmsf = numpy.average(var)
		err_rmsf = numpy.std(var)
		bf = (8*numpy.pi*numpy.pi/3)*rmsf
		err_bf = (8*numpy.pi**2/3)*err_rmsf

	filename = name_file_dcd + "b-factor.LA" + cutoff +".%03d." % res + chain + ".txt"
	dataout = open(filename, 'w')
	dataout.write(str(res) + "\t%.5f" % bf + "\t%.5f" % err_bf + '\n')
	dataout.close()
	os.remove(filenamedcd)

arg_list = []
#for resid in range(3,res_num):
#	arg_list.append([resid,name_file_dcd,mol])

res_list = range(res_beg,res_end)
pool = multip.Pool(processes=Np)
pool.map(_bfactor, res_list)

print "Grouping and cleaning"
os.system('cat ' + name_file_dcd + 'b-factor.LA' + cutoff + '.*.' + chain + '.txt > ' + name_file_dcd + 'all.bf.LA' + cutoff + '.' + chain + '.txt')
os.system('rm ' + name_file_dcd + 'b-factor.LA' + cutoff + '.*.' + chain + '.txt')
print "Done... with capital D"


print "done"

#print len(dotprod)
#print len(c_func)
#print mol.trajectory
