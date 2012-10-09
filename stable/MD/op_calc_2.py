#! /usr/bin/env python 
###################################
#    Davi Ortega 6/13/2012 
###################################
import sys
import numpy
if '-h' in sys.argv:
	print 'put your explanation here'
	sys.exit()

cutoff = 0.005
filename = sys.argv[1]

data = open(filename, 'r')

cts = {}
sim_list = []
res_list = []

print "Loading the data"

for line in data:
	if "simulation" in line:
		header = line
	else:
		field = line.split('\t')
		if field[0] not in cts.keys():
			cts[field[0]] = {}
			sim_list.append(field[0])
			print "Working on simulation " + field[0]
		if field[1] not in cts[field[0]].keys():
			cts[field[0]][field[1]] = [[],[],[]]
			if field[1] not in res_list:
				res_list.append(field[1])
		cts[field[0]][field[1]][0].append(numpy.double(field[2])) #t
		cts[field[0]][field[1]][1].append(numpy.double(field[3])) #ct
		cts[field[0]][field[1]][2].append(numpy.double(field[4])) #sd(ct)

print res_list
print sim_list

data.close()

#output = 'residue	Coo_1	Ctail_1	|Coo_1 - Ctail_1|	sd(Ctail_1)	Coo_2	Ctail_2	|Coo_2 - Ctail_2|	sd(Ctail_2)	Coo_3	Ctail_3	|Coo_3 - Ctail_3|	sd(Ctail_3)	Coo_4	Ctail_4	|Coo_4 - Ctail_4|	sd(Ctail_4)	Coo_5	Ctail_5	|Coo_5 - Ctail_5|	sd(Ctail_5)	Coo_6	Ctail_6	|Coo_6 - Ctail_6|	sd(Ctail_6)	Coo_7	Ctail_7	|Coo_7 - Ctail_7|	sd(Ctail_7)	Coo_8	Ctail_8	|Coo_8 - Ctail_8|	sd(Ctail_8)	Coo_9	Ctail_9	|Coo_9 - Ctail_9|	sd(Ctail_9)	Coo_10	Ctail_10	|Coo_10 - Ctail_10|	sd(Ctail_10)	Coo_MEAN	Ctail_MEAN	|Coo_MEAN - Ctail_MEAN|	sd(Ctail_MEAN)	1	2	3	4	5	6	7	8	9	10	Converged	S^2	Err.S^2\n'

#generating header

output = ''
for sim in sim_list:
	output += 'residue\tCoo_'+str(sim)+'\tCtail_'+str(sim)+'\t|Coo_'+str(sim)+' - Ctail_'+str(sim)+'|'
for sim in sim_list:
	output += '\t'+str(sim)
output += '\tConverged\tS^2\tErr.S^2\n'


#output = 'residue	Coo_1	Ctail_1	|Coo_1 - Ctail_1|	Coo_2	Ctail_2	|Coo_2 - Ctail_2|	Coo_3	Ctail_3	|Coo_3 - Ctail_3|	Coo_4	Ctail_4	|Coo_4 - Ctail_4|	Coo_5	Ctail_5	|Coo_5 - Ctail_5|	Coo_6	Ctail_6	|Coo_6 - Ctail_6|	Coo_7	Ctail_7	|Coo_7 - Ctail_7|	Coo_8	Ctail_8	|Coo_8 - Ctail_8|	Coo_9	Ctail_9	|Coo_9 - Ctail_9|	Coo_10	Ctail_10	|Coo_10 - Ctail_10|	Coo_Average	Ctail_Average	|Coo_MEAN - Ctail_Average|	1	2	3	4	5	6	7	8	9	10	Converged	S^2	Err.S^2\n'

for res in res_list:
	print "Calculations for residue " + str(res)
	results = []
	output += str(res)
	for sim in sim_list:
		ct = cts[sim][res][1]
		coo = numpy.average(ct)
		ctail = numpy.average(ct[-500:])
#		std_tail = numpy.std(ct[-500:])
		test = abs(coo - ctail)
		output += '\t%.5f' % coo + '\t%.5f' % ctail + '\t%.5f' % test
		if sim != 'AVERAGE':
			if test < cutoff:
				results.append(coo)
			else:
				results.append('')
	for S in results:
		if S == '':
			output += '\tNA'
		else:
			output += '\t%.5f' % S
	while '' in results:
		results.remove('')
	if len(results) != 0:
		results = numpy.array(results)
		ave = numpy.average(results)
		std = numpy.std(results)
		output += '\t' + str(len(results)) + '\t%.5f' % ave + '\t%.5f' % std
	else:
		output += '\t0\tNA\tNA'
	output += '\n'

outfile = open(filename[:-3] + 'op_calc.txt', 'w')
outfile.write(output)
outfile.close()
			
	
	
	
