#! /usr/bin/env python 
###################################
#    Davi Ortega 7/26/2012 
###################################
import sys
import MDAnalysis as MD
import MDAnalysis.analysis.distances as MDdist
import multiprocessing as multip
import numpy
import scipy.optimize
import time
import os
from numpy.core.umath_tests import inner1d


pdb_file = sys.argv[1]
dcd_file = sys.argv[2]
Np = int(sys.argv[3])
res_beg = int(sys.argv[4])
res_end = int(sys.argv[5])
step = float(sys.argv[6])
skip = int(sys.argv[7])
res_num = res_end - res_beg

print "Loading trajectory"
mol = MD.Universe(pdb_file, dcd_file)



name_file_dcd = sys.argv[2][:-3]
while "/" in name_file_dcd:
	name_file_dcd = name_file_dcd[name_file_dcd.index("/")+1:]


list_Ct = {}

output = ''

def cbcb(res):
	#mol is an md universe
#	res,name_file_dcd,mol=args
	mol = MD.Universe(pdb_file, dcd_file)
	output = ''
	print "["+str(res)+"] Calculating CB-CB for residue " + str(res) + ". " + str(res - res_beg) + " of " + str(res_num)
	cf_timer = time.time()
	try:
		CB_A = mol.selectAtoms("resid " + str(res) + " and name CB and segid A")
		CB_B = mol.selectAtoms("resid " + str(res) + " and name CB and segid B")
	except:
		CB_A = mol.selectAtoms("resid " + str(res) + " and name CA and segid A")
	        CB_B = mol.selectAtoms("resid " + str(res) + " and name CA and segid B")
		
	if 1 == 1:
		d_list = []
		frame = 0
		output = ''
		for ts1 in mol.trajectory[0:-1:skip]:
			d = MDdist.dist(CB_A, CB_B)
			d_list.append(d[2])
			output += str(res) + "\t" + str(frame) + "\t%.5f" % d[2] + "\n"
			
			#print str(frame) + '\t' + str(d[2])
			frame += 1
		dlist = numpy.array(d_list)
		print "["+str(res)+"]" + '\t' + str(numpy.average(d_list)) + '\t' + str(numpy.std(d_list)) + '\tin ' + str(time.time()-cf_timer) + " seconds.\n"
		filename = name_file_dcd + "cbcb_numpy.%03d.txt" % res
		outfile = open(filename, 'w')
		outfile.write(output)
		outfile.close()

arg_list = []
#for resid in range(3,res_num):
#	arg_list.append([resid,name_file_dcd,mol])
#_ct_part(79)
#sys.exit()

#cbcb(280)
#sys.exit()

res_list = range(res_beg,res_end)
pool = multip.Pool(processes=Np)
pool.map(cbcb, res_list)

print "Grouping and cleaning"
os.system('cat ' + name_file_dcd + 'cbcb_numpy.*txt > ' + name_file_dcd + 'all.cbcb.txt')
os.system('rm ' + name_file_dcd + 'cbcb_numpy.*txt')
print "Done... with capital D"

#print len(dotprod)
#print len(c_func)
#print mol.trajectory
