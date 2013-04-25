#! /usr/bin/env python 
###################################
#    Davi Ortega 3/2/2013
###################################
import sys
import MDAnalysis as MD
import numpy

if len(sys.argv) != 8:
	print "This is a raw program, it will spit out the X1 dihedral angle (N - Ca - Cb - Cg) of backbone for each residue of the specified chain (column) for each frame.\nScript takes exactly 7 arguments: pdb_file dcd_file chain res_beg res_end time_step skip"

	sys.exit()

pdb_file = sys.argv[1]
dcd_file = sys.argv[2]
chain = sys.argv[3]
res_beg = int(sys.argv[4])
res_end = int(sys.argv[5])
step = float(sys.argv[6])
skip_ts = int(sys.argv[7])


res_num = res_end - res_beg


print "Loading trajectory"
ref = MD.Universe(pdb_file)
mol = MD.Universe(pdb_file, dcd_file)
dihe_col = MD.collection


name_file_dcd = sys.argv[2][:-3]
while "/" in name_file_dcd:
	name_file_dcd = name_file_dcd[name_file_dcd.index("/")+1:]

list_Ct = {}

output = 'time'

for res in range(res_beg, res_end+1):
	try:
		dihe_col.addTimeseries(MD.Timeseries.Dihedral(mol.selectAtoms("segid " + chain +" and resid " + str(res) + " and name N", "segid " + chain +" and resid " + str(res) + " and name CA", "segid " + chain +" and resid " + str(res) + " and name CB", "segid " + chain +" and resid " + str(res) + " and name CG")))
		output += '\tRES' + str(res) 
	except MD.NoDataError:
		try:
			dihe_col.addTimeseries(MD.Timeseries.Dihedral(mol.selectAtoms("segid " + chain +" and resid " + str(res) + " and name N", "segid " + chain +" and resid " + str(res) + " and name CA", "segid " + chain +" and resid " + str(res) + " and name CB", "segid " + chain +" and resid " + str(res) + " and name CG1")))
			output += '\tRES' + str(res)
		except MD.NoDataError:
			try:
				dihe_col.addTimeseries(MD.Timeseries.Dihedral(mol.selectAtoms("segid " + chain +" and resid " + str(res) + " and name N", "segid " + chain +" and resid " + str(res) + " and name CA", "segid " + chain +" and resid " + str(res) + " and name CB", "segid " + chain +" and resid " + str(res) + " and name CG2")))
                        	output += '\tRES' + str(res)
			except MD.NoDataError:
				print "No deal residue " + str(res) + " is a " + mol.selectAtoms('segid ' + chain + ' and resid ' + str(res)).resnames()[0]
				pass
output += '\n'
print dihe_col

dihe_data = mol.trajectory.correl(dihe_col, skip = skip_ts)*180./numpy.pi


for t in range(len(dihe_data.T)):
	output += str(t*step) 
	for value in dihe_data.T[t]:
		output +=  "\t%.4f" % value
	output += '\n'

dataout = open(name_file_dcd + "X1dihedral." + chain +".dat",'w')
dataout.write(output)
dataout.close() 

print "Done... with capital D"

