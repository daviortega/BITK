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
import scipy.optimize
import matplotlib.pylab as mp
import time
import glob
from numpy.core.umath_tests import inner1d

if len(sys.argv) != 12:
	print "Spits out the RMSF of Calpha for each residue.\nScript takes exactly 10 arguments: pdb_file dcd_file Np res_beg res_end chain time_step skip cutoff ave full_output\n \
	With ave = 0 the script will use x(0) to calculate < x(t) - x(0) > ^ 2. With ave = 1 it will compute <x(t) - <x(t)>>^2.\n \
	Unnecessary to say that it will take double the time with ave = 1\n \
	With full_output = 1 produces extra file with the time series of the value for each residue. Filename is the same but with .fulldata extension"

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
ave_flag = int(sys.argv[10])
fop_flag = int(sys.argv[11])
res_num = res_end - res_beg


if ave_flag not in [0,1]:
	print "ave flag out of context. Pick 0 or 1"
	sys.exit()

print "Loading trajectory"
mol = MD.Universe(pdb_file, dcd_file)



name_file_dcd = sys.argv[2][:-3]
while "/" in name_file_dcd:
	name_file_dcd = name_file_dcd[name_file_dcd.index("/")+1:]

list_Ct = {}

output = 'residue\trmsf\terr_rmsf'

def _rmsf(res):
	global ave_flag
	global fop_flag
	
	#mol is an md universe
	#res,name_file_dcd,mol=args
	ref = MD.Universe(pdb_file)
	mol = MD.Universe(pdb_file, dcd_file)
	print "["+str(res)+"] Aligning Trajectory"
	ref_sel = ref.selectAtoms('backbone and resid ' + str(res),'backbone and around ' + cutoff + ' (backbone and resid ' + str(res) + ')')
	print ref_sel.atoms
	#for i in list(ref_sel):
	#	print i
	#	print ref_sel.indices
	ref_sel_str = ''
	for i in list(ref_sel):
		useful =str(i).split(":")[1]
		ref_sel_str += '(' + useful.replace('>',' ').replace("'",'').replace(',',' and').replace('of','and') + ') or '
		#ref_sel_str += "(" + useful.replace('>',' ').replace("'",'') + ') or '
	ref_sel_str = ref_sel_str[:-4] 
	#print ref_sel_str
	filenamedcd = "LA" + cutoff + "rmsf.R" + str(res) + '.' + chain + '.dcd'
	MDA.rms_fit_trj(mol, ref, select=ref_sel_str, filename = filenamedcd ) #'name N and type N and resname MET resid 1 and segid P1')
	output = ''

	print "["+str(res)+"] Calculating RMSF for residue " + str(res) + "..." + str( 1 + range(res_beg,res_end).index(res) ) + " of " + str(res_num)
	mol = MD.Universe(pdb_file, filenamedcd)
	cf_timer = time.time()
	atoms_mol = mol.selectAtoms("resid " + str(res) + " and name CA and segid " +str(chain))
	atoms_ref = ref.selectAtoms("resid " + str(res) + " and name CA and segid " +str(chain))
	var = []
	if ave_flag == 0:
		for ts1 in mol.trajectory[0:-1:skip]:
			var.append((numpy.linalg.norm(numpy.subtract(atoms_mol.coordinates(), atoms_ref.coordinates())))**2)
	elif ave_flag:
		ref_coord = []
		for ts1 in mol.trajectory[0:-1:skip]:
			ref_coord.append(atoms_mol.coordinates()[0])
		ref_coord = numpy.array(ref_coord)
		#print ref_coord.shape
		#ref_coord = ref_coord[1:5]
		#print numpy.average( ref_coord[:,0])
		#print numpy.average( ref_coord[:,1])
		#print numpy.average( ref_coord[:,2])
		ave = [numpy.average( ref_coord[:,0]) , numpy.average( ref_coord[:,1]), numpy.average( ref_coord[:,2])]
		for ts1 in mol.trajectory[0:-1:skip]:
                        var.append((numpy.linalg.norm(numpy.subtract(atoms_mol.coordinates(), ave)))**2)

		
	var = numpy.array(var)
	rmsf = numpy.average(var)
	err_rmsf = numpy.std(var)

	filename = name_file_dcd + "rmsf.LA" + cutoff +".%04d." % res + chain + ".txt"
	dataout = open(filename, 'w')
	dataout.write(str(res) + "\t%.5f" % rmsf + "\t%.5f" % err_rmsf + '\n')
	dataout.close()
	os.remove(filenamedcd)
	
	if fop_flag == 1:
		output = str(res) +'\n'
		for i in range(len(var)):
			output += '%.5f' % var[i] + '\n'
		fopfile = open(name_file_dcd + "rmsf.LA" + cutoff +".%04d." % res + chain + ".fulldata",'w')
		fopfile.write(output)
		fopfile.close()

arg_list = []
#for resid in range(3,res_num):
#	arg_list.append([resid,name_file_dcd,mol])
#_rmsf(263)
#sys.exit()


res_list = range(res_beg,res_end)
pool = multip.Pool(processes=Np)
pool.map(_rmsf, res_list)

print "Grouping and cleaning"
#if ave_flag == 0:
#	os.system('cat ' + name_file_dcd + 'rmsf.LA' + cutoff + '.*.' + chain + '.txt > ' + name_file_dcd + 'all.rmsf.ref_ini.LA' + cutoff + '.' + chain + '.txt')
#	if fop_flag == 1:
#		os.system('touch '  + name_file_dcd + 'all.rmsf.ref_ini.LA' + cutoff + '.' + chain + '.fulldata')
#		files = glob.glob( name_file_dcd + "rmsf.LA" + cutoff + ".*." + ".fulldata" )
#		files.sort()
#		for f in files:
#			print f
#			os.system('mv ' + name_file_dcd + 'all.rmsf.ref_ini.LA' + cutoff + '.' + chain + '.fulldata temp.dat')
#			os.system('paste temp.dat ' + f + ' > ' + name_file_dcd + 'all.rmsf.ref_ini.LA' + cutoff + '.' + chain + '.fulldata' )
#else:

if ave_flag == 1:
	ave_str = 'ave'
else:
	ave_str = 'ini'

os.system('cat ' + name_file_dcd + 'rmsf.LA' + cutoff + '.*.' + chain + '.txt > ' + name_file_dcd + 'all.rmsf.ref_' + ave_str + '.LA' + cutoff + '.' + chain + '.txt')
if fop_flag == 1:
#	os.system('touch '  + name_file_dcd + 'all.rmsf.ref_ave.LA' + cutoff + '.' + chain + '.fulldata')
	files = glob.glob(name_file_dcd + "rmsf.LA" + cutoff + ".*." + ".fulldata")
	files.sort()
	for f in range(len(files)-1):
		if f == 0:
			print 'Pasting ' + files[f] + ' and ' + files[f+1] 
			os.system('paste ' + files[f] + ' ' + files[f+1] + ' > ' + name_file_dcd + 'all.rmsf.ref_' + ave_str + '.LA' + cutoff + '.' + chain + '.fulldata' )
		else:
			print 'Pasting pasted and ' + files[f+1]
			os.system('mv ' + name_file_dcd + 'all.rmsf.ref_' + ave_str + '.LA' + cutoff + '.' + chain + '.fulldata temp.dat')
			os.system('paste temp.dat ' + files[f+1] + ' > ' + name_file_dcd + 'all.rmsf.ref_' + ave_str + '.LA' + cutoff + '.' + chain + '.fulldata' )

os.system('rm ' + name_file_dcd + 'rmsf.LA' + cutoff + '.*.' + chain + '.txt')
os.system('rm ' + name_file_dcd + 'rmsf.LA' + cutoff + '.*.' + chain + '.fulldata')


print "Done... with capital D"

