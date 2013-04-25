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
	print "script takes exactly 9 arguments: pdb_file dcd_file Np res_beg res_end chain time_step skip cutoff"
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

output = ''

def _ct_part(res):
	#mol is an md universe
#	res,name_file_dcd,mol=args
	ref = MD.Universe(pdb_file)
	mol = MD.Universe(pdb_file, dcd_file)
	print "["+str(res)+"] Aligning Trajectory"
	ref_sel = ref.selectAtoms('backbone and resid ' + str(res),'backbone and around ' + cutoff + ' (backbone and resid ' + str(res) + ')')
#	print ref_sel.atoms
#	for i in list(ref_sel):
#		print i
#	print ref_sel.indices
	ref_sel_str = ''
	for i in list(ref_sel):
		useful =str(i).split(":")[1]
		ref_sel_str += '(' + useful.replace('>',' ').replace("'",'').replace(',',' and').replace('of','and') + ') or '
		#ref_sel_str += "(" + useful.replace('>',' ').replace("'",'') + ') or '
	ref_sel_str = ref_sel_str[:-4] 
#	print ref_sel_str
	filenamedcd = "LA" + cutoff + ".R" + str(res) + '.dcd'
	MDA.rms_fit_trj(mol, ref, select=ref_sel_str, filename = filenamedcd ) #'name N and type N and resname MET resid 1 and segid P1')
	output = ''
	print "["+str(res)+"] Calculating C(t) for residue " + str(res) + " of " + str(res_num)
	mol = MD.Universe(pdb_file, filenamedcd)
	cf_timer = time.time()
	N_atoms = mol.selectAtoms("resid " + str(res) + " and (name N or name HN) and segid " +str(chain))

	if len(N_atoms.coordinates()) == 2:
		bond_vec = []
		ct_numpy = []
		ct_std_numpy = []
		ct_python = []
		
		for ts1 in mol.trajectory[0:-1:skip]:
        		coords = N_atoms.coordinates()
	        	coords = numpy.subtract(coords[1], coords[0])
	       		bond_vec.append(coords/numpy.linalg.norm(coords))
		print "["+str(res)+"] Calculating Correlation function as Chen 2004"
		
		for t in range(len(bond_vec)/2):
			T = t*step
		#        	x_chen.append(T)
			if t == 0:
                                a = numpy.array(bond_vec)
                                b = numpy.array(bond_vec)
                        else:
                                a = numpy.array(bond_vec[:-t])
                                b = numpy.array(bond_vec[t:])
			np_dot = inner1d(a, b)
			np_dot = (3*np_dot**2 - 1)/2
                        ct = numpy.mean(np_dot)
			ct_std = numpy.std(np_dot)
#                       ct = (3*(X)/a.shape[0]- 1)/2
                        ct_numpy.append(ct)
			ct_std_numpy.append(ct_std)
#	                NH_dot = 0
#        		for tal in range(len(bond_vec)-t):
 #             			NH_dot += (3*(numpy.dot(bond_vec[tal], bond_vec[tal+t])**2)-1)/2
#			NH_dot = NH_dot/(len(bond_vec)-t)
#		        ct_python.append(NH_dot)
#			output += str(res) + "\t" + str(T) + "\t" + str(NH_dot) + "\n"
			output += str(res) + "\t" + str(T) + "\t%.5f" % ct + "\t%.5f" % ct_std + "\n" 
#			output += str(res) + "\t" + str(T) + "\t" + str(NH_dot) + "\t" + str(ct) + "\n"
#			print str(res) + "\t" + str(T) + "\t" + str(NH_dot) + "\t" + str(ct) + "\n"
#			time.sleep(5)
		ct_numpy = numpy.around(ct_numpy, 5)
	print "["+str(res)+"] Done in " + str(time.time()-cf_timer) + " seconds.\n"
	filename = name_file_dcd + "ct_numpy.LA" + cutoff +".%03d." % res + chain + ".txt" 
	outfile = open(filename, 'w')
	outfile.write(output)
	outfile.close()
	print "["+str(res)+"] Collecting Trash"
	os.remove(filenamedcd)
	print "["+str(res)+"] Done for real"

arg_list = []
#for resid in range(3,res_num):
#	arg_list.append([resid,name_file_dcd,mol])



res_list = range(res_beg,res_end+1)
pool = multip.Pool(processes=Np)
pool.map(_ct_part, res_list)

print "Grouping and cleaning"
os.system('cat ' + name_file_dcd + 'ct_numpy.LA' + cutoff +'.*.' + chain + '.txt > ' + name_file_dcd + 'all.ct.' + chain +'.LA' + cutoff + '.txt')
os.system('rm ' + name_file_dcd + 'ct_numpy.LA' + cutoff + '.*.' + chain + '.txt')
print "Done... with capital D"


print "done"

#print len(dotprod)
#print len(c_func)
#print mol.trajectory
