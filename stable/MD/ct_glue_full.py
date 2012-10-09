#! /usr/bin/env python 
###################################
#    Davi Ortega 5/9/2012 
###################################
import sys
import MDAnalysis as MD
import numpy
import rpy2
import scipy.optimize
import matplotlib.pylab as mp
import time

if '-h' in sys.argv:
	print 'Glue all *.all.ct.txt together for easy use in R and with op_calc_2.* (R or py)... you are welcome... simulations should be out of jobmaker.sh. If not, then please add -c (common) for single trajectory'
	sys.exit()

#sim_num = len(sys.argv[1:])

ct_data = []
sim_list = []
res_list = []
time_list = []

output = 'simulation\tresidue\tt\tct\tErr(ct)\n'
ct_values = []
first = True

if '-c' in sys.argv:
	file_list = [sys.argv[1]]
else:
	file_list = sys.argv[1:]

print file_list

for filename in file_list:
	print "reading " + filename
	
	if len(file_list) == 1:
		sim_num = str(1)
	else:
		sim_num = filename.split('-')[2]
#	if sim_num.isdigit():
#		sim_num = filename0]
#	else:
#		sim_num = filename[8]
		
		
	sim_list.append(sim_num)
	data_file = open(filename,'r')
	old_time = time.time()
	i = 0
	for line in data_file:
		if line != '':
			output += sim_num + '\t' + line
			fields = line.split('\t')
			if float(fields[1])*1000 % 45000 == 0:
                                print "Dealing with simulation " + sim_num + " residue " + fields[0] + " and t = " + fields[1] + ' in ' + str(time.time() - old_time)
				old_time = time.time()

			if float(fields[1])*1000 % 1 == 0:
				if first:
					ct_data.append(numpy.array([float(fields[2])]))
					if fields[0] not in res_list:
	                                        res_list.append(fields[0])
        	                        if fields[1] not in time_list:
                	                        time_list.append(fields[1])
				else:
					try:
						ct_data[i] = numpy.append(ct_data[i],[float(fields[2])])
					except ValueError:
						print fields
						sys.exit()

#					print "not first"
#					try:
#						numpy.append(ct_data[i],[float(fields[2])])
#					except:
#						print "ops... something is wrong"
#						print fields
#						print len(ct_data)
#						print i
#						sys.exit()
				i += 1

#			elif fields[1] not in ct_data[fields[0]].keys():
#				ct_data[fields[0]][fields[1]] = [fields[2]] #numpy.array(float(fields[2]))
#				if fields[1] not in time_list:
#					time_list.append(fields[1])
#			else:
#				numpy.append(ct_data[fields[0]][fields[1]], float(fields[2]))
	
	first = False
	print len(ct_data)
	data_file.close()

#print ct_data[100]
print len(ct_data)
print len(res_list)
print len(time_list)

if '-c' not in sys.argv:
	print "Including Average C(t)"
	for res in range(len(res_list)):
		for t in range(len(time_list)):
			output += 'AVERAGE\t' + str(res_list[res]) + '\t' + str(time_list[t]) + '\t%.5f' % ct_data[res*len(time_list) + t].mean() + '\t%.5f' % ct_data[res*len(time_list) + t].std() + '\n'

dataout = open('glued_data.txt','w')
dataout.write(output)
dataout.close()

sys.exit()







#			fields = line.split('\t')
#			if fields[0] not in ct_data[sim_num].keys():
#				ct_data[sim_num][fields[0]] = [[float(fields[1])], [float(fields[2])], [float(fields[3])]]
#				res_list[sim_num].append(fields[0])
#			else:
#				ct_data[sim_num][0].append(float(fields[0]))
#				ct_data[sim_num][fields[0]][0].append(float(fields[1]))
#				ct_data[sim_num][fields[0]][1].append(float(fields[2]))
#				ct_data[sim_num][fields[0]][2].append(float(fields[3]))
#	data_file.close()

#organize simulations
sim_list.sort()

output = 't'

#print ct_data[sim_list[0]]

for t in ct_data[sim_list[0]][ct_data[sim_num].keys()[0]][0]:
	print t

sys.exit() 


for res in res_list:
	output += '\tR' + res + '\tD' + res
	output += '\n'
	for i in range(len(ct_data[res_list[0]][0])):
		output += str(ct_data[res_list[0]][0][i])
		for res in res_list:
			output += '\t%.5f' %ct_data[res][1][i] + '\t%.5f' % ct_data[res][2][i]
		output += '\n'
	export_file = open(sys.argv[1][:-3] + 'table.txt','w')
	export_file.write(output)
	export_file.close()

sys.exit()

pdb_file = sys.argv[1]
dcd_file = sys.argv[2]
res = sys.argv[3]
step = float(sys.argv[4])
skip = int(sys.argv[5])

print "Loading trajectory"

mol = MD.Universe(pdb_file, dcd_file)
N_atoms = mol.selectAtoms("resid " + res + " and (name N or name HN)")

print "Done... Calculating N-H bond"

bond_vec = []

for ts1 in mol.trajectory[0:-1:skip]:
	coords = N_atoms.coordinates()
        coords = numpy.subtract(coords[1], coords[0])
	bond_vec.append(coords/numpy.linalg.norm(coords))

x = []
c_func = []

x_chen = []
c_func_chen = [] #Chen et all 2004 with Brooks 1992 Eq (2)

#print "Done... Calculating Correlation Function"
#cf_timer = time.time()
#for t in range(int((len(bond_vec)*0.85))/2):
#	x.append(t*step)
#	NH_dot = 0
#	for i in range(len(bond_vec)-t):
#		cos = numpy.dot(bond_vec[i], bond_vec[i+t])
#		NH_dot += (3*cos*cos - 1)/2
#	NH_dot = NH_dot/(len(bond_vec)-t)
#	c_func.append(NH_dot)
#print "Done in " + str(time.time()-cf_timer) + "s... Calculating Correlation function as Chen 2004"
print "Done... Calculating Correlation function as Chen 2004"

cf_timer = time.time()
for t in range(len(bond_vec)):
	x_chen.append(t*step)
	NH_dot = 0
	for tal in range(len(bond_vec)-t):
		NH_dot += (3*(numpy.dot(bond_vec[tal], bond_vec[tal+t])**2)-1)/2
	c_func_chen.append(NH_dot/(len(bond_vec)-t))


print "Done in " + str(time.time()-cf_timer) + "s... Calculating S2 LP approximation"

#print "Done... Calculating S2 LP approximation"

S2 = 0 

for t1 in range(len(bond_vec)/2):
	S2_par = 0
	for t2 in range(len(bond_vec)/2):
		cos = numpy.dot(bond_vec[t1], bond_vec[t1+t2])
		S2_par += cos*cos
	S2 += S2_par/((len(bond_vec)/2)**2)
S2 = 3*S2/2 - 0.5

print S2	

print "Done... Results"


#x = numpy.array(x,dtype='float')
#c_func = numpy.array(c_func, dtype='float')

x_chen = numpy.array(x_chen[0:int(len(bond_vec)*0.9)],dtype='float')
c_func_chen = numpy.array(c_func_chen[0:int(len(bond_vec)*0.9)], dtype='float')

	
#print c_func[:-1]
#print x[:-1]
print numpy.average(c_func_chen[-100:])
print numpy.average(c_func_chen[-10:])

#print len(c_func)
#print len(x)

def _CF_gl_md4(p, x): #Glore et all 1990 with S == Sf as in model 4. or Eq 9 on Chen 2004
        te,S=p
        y = S**2 + (1 - S**2) * numpy.exp(-x/te)
        return y

def _CF_gl(p, x): #Glore et all 1990
	tm,S,Sf,te=p
	y = S**2 + (1 - Sf**2) * numpy.exp(-x/tm) + (Sf**2 - S**2)*numpy.exp(-x/te)
	return y

def _CF_ex(p, x): #Shaw et all 2008
        tm,S,Sf,te=p
        y = numpy.exp(-x/tm) * ( S**2 + (Sf**2 - S**2)*numpy.exp(-x/te))
        return y

def _CF_or(p, x): #Shaw et all 2008
	tm,S,te=p
	y = numpy.exp(-x/tm) * ( S**2 + (1 - S**2)*numpy.exp(-x/te))
        return y

def _CF_davi(p, x): #Shaw et all 2008
        tm,S,Sf,tf,ts, =p
        y = numpy.exp(-x/tm) * (S**2 + (1 - Sf**2) * numpy.exp(-x/tf) + (Sf**2 - S**2)*numpy.exp(-x/ts))
        return y


def _CF_gl_md4_residuals(p,x,y):
        return _CF_gl_md4(p, x)-y

def _CF_gl_residuals(p,x,y):
        return _CF_gl(p, x)-y

def _CF_ex_residuals(p,x,y):
        return _CF_ex(p, x)-y

def _CF_or_residuals(p,x,y):
        return _CF_or(p, x)-y

def _CF_davi_residuals(p,x,y):
        return _CF_davi(p, x)-y


def _get_OP_(x, y):
	print "CF_gl_md4"
 	p_guess = (0.01, 1)
        p_gl_md4, s_gl_md4 = scipy.optimize.leastsq(_CF_gl_md4_residuals, p_guess, args=(x, y))
        ts_gl_md4,S_gl_md4=p_gl_md4
	print "CF_gl"
	p_guess = (1, 0.5, 1, 0.1)
        p_gl, s_gl = scipy.optimize.leastsq(_CF_gl_residuals, p_guess, args=(x, y))
	tm_gl,S_gl,Sf_gl,te_gl=p_gl
	print "CF_ex"
	p_guess = (4, 0.5, 1, 1)
	p_ex, s_ex = scipy.optimize.leastsq(_CF_ex_residuals, p_guess, args=(x, y))
	tm_ex,S_ex,Sf_ex,te_ex=p_ex
	print "CF_or"
	p_guess = (4, 0.5, 1)
        p_or, s_or = scipy.optimize.leastsq(_CF_or_residuals, p_guess, args=(x, y))
        tm_or,S_or, te_or=p_or
	print "CF_davi"
	p_guess = (4,1,1,0.001,1)
	p_davi, s_davi = scipy.optimize.leastsq(_CF_davi_residuals, p_guess, args=(x, y), maxfev = 12000)
        tm_davi,S_davi,Sf_davi,tf_davi,ts_davi=p_davi
	print "Ci(t) for simple cases (Sf^2 == 1 ) by Eq (9) on Chen - black - S2 = " + str(p_gl_md4[1]**2)
	print p_gl_md4
	r = _CF_gl_md4_residuals(p_gl_md4,x,y)
	r = r**2
#       print sum(r)
	print numpy.corrcoef(y,_CF_gl_md4(p_gl_md4,x))[0][1]
	print "Ci(t) for complex cases or Eq 11 on Chen - magenta - S2 = " + str(p_gl[1]**2)
	print p_gl
#	print sum(_CF_gl_residuals(p_gl,x,y)**2)
	print numpy.corrcoef(y,_CF_gl(p_gl,x))[0][1]
	print "C(t) for complex residues as Shaw 2005 - green - S2 = " + str(p_ex[1]**2)
	print p_ex
#	print sum(_CF_ex_residuals(p_ex,x,y)**2)
	print numpy.corrcoef(y,_CF_ex(p_ex,x))[0][1]
	print "C(t) for simple residues as Shaw 2005 - red - S2 = " + str(p_or[1]**2)
	print p_or
#	print sum(_CF_or_residuals(p_or,x,y)**2)
	print numpy.corrcoef(y,_CF_or(p_or,x))[0][1]
	print "C(t) for complex residues Davi 2012 - blue - S2 = " + str(p_davi[1]**2)
	print p_davi
#	print sum(_CF_davi_residuals(p_davi,x,y)**2)
	print numpy.corrcoef(y,_CF_davi(p_davi,x))[0][1]
#	print "S = " + str(S)

#	xp = numpy.linspace(x.min(), x.max(), numpy.floor((x.max()-x.min())*100)
	p_cf_gl_md4 = _CF_gl_md4(p_gl_md4, x)
	p_cf_gl = _CF_gl(p_gl, x)
	p_cf_ex = _CF_ex(p_ex, x)
	p_cf_davi = _CF_davi(p_davi, x)
	p_cf_or = _CF_or(p_or, x)
	
	mp.plot(x, y, 'o', x, p_cf_ex, 'g-', x, p_cf_or, 'r-', x, p_cf_davi, 'b-', x, p_cf_gl, '-m', x, p_cf_gl_md4, '-k')
	mp.grid(True)
	mp.xscale('log')
	mp.show()

	

print "\nfirst 10 ns function"
_get_OP_(x_chen[1:int(10/step)], c_func_chen[1:int(10/step)])
print "\nfirst 20 ns function"
_get_OP_(x_chen[1:int(20/step)], c_func_chen[1:int(20/step)])
print "\nfirst 30 ns function"
_get_OP_(x_chen[1:int(30/step)], c_func_chen[1:int(30/step)])
print "\nfirst 40 ns function"
_get_OP_(x_chen[1:int(40/step)], c_func_chen[1:int(40/step)])
print "\nFull function"
_get_OP_(x_chen, c_func_chen)


print "done"

#print len(dotprod)
#print len(c_func)
#print mol.trajectory
