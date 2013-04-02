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
from numpy.core.umath_tests import inner1d

if '-h' in sys.argv:
	print 'put your explanation here'
	sys.exit()


pdb_file = sys.argv[1]
dcd_file = sys.argv[2]
res = sys.argv[3]
chain = sys.argv[4]
step = float(sys.argv[5])
skip = int(sys.argv[6])

print "Loading trajectory"

mol = MD.Universe(pdb_file, dcd_file)
#N_atoms = mol.selectAtoms("resid " + res + " and (name N or name HN)")
N_atoms = mol.selectAtoms("resid " + str(res) + " and (name N or name HN) and segid " + str(chain))

print "Done... Calculating N-H bond"

if len(N_atoms.coordinates()) == 2:
	bond_vec = []
        ct_numpy = []
        ct_std_numpy = []
        ct_python = []
	x_chen = []

        for ts1 in mol.trajectory[0:-1:skip]:
        	coords = N_atoms.coordinates()
                coords = numpy.subtract(coords[1], coords[0])
                bond_vec.append(coords/numpy.linalg.norm(coords))
        print "["+str(res)+"] Calculating Correlation function as Chen 2004"
	cf_timer = time.time()
        for t in range(len(bond_vec)/2):
                T = t*step
                #               x_chen.append(T)
                if t == 0:
        	        a = numpy.array(bond_vec)
                        b = numpy.array(bond_vec)
                else:
                        a = numpy.array(bond_vec[:-t])
                        b = numpy.array(bond_vec[t:])
                np_dot = inner1d(a, b)
                np_dot = (3*np_dot**2 - 1)/2
                ct = numpy.mean(np_dot)
#               ct_std = numpy.std(np_dot)
                ct_numpy.append(ct)
#               ct_std_numpy.append(ct_std)
        ct_numpy = numpy.around(ct_numpy, 5)
	x_chen.append(t*step)
#x = []
#c_func = []

c_func_chen = ct_numpy #Chen et all 2004 with Brooks 1992 Eq (2)

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
#	print numpy.corrcoef(y,_CF_gl_md4(p_gl_md4,x))[0][1]
	print "Ci(t) for complex cases or Eq 11 on Chen - magenta - S2 = " + str(p_gl[1]**2)
	print p_gl
#	print sum(_CF_gl_residuals(p_gl,x,y)**2)
#	print numpy.corrcoef(y,_CF_gl(p_gl,x))[0][1]
	print "C(t) for complex residues as Shaw 2005 - green - S2 = " + str(p_ex[1]**2)
	print p_ex
#	print sum(_CF_ex_residuals(p_ex,x,y)**2)
#	print numpy.corrcoef(y,_CF_ex(p_ex,x))[0][1]
	print "C(t) for simple residues as Shaw 2005 - red - S2 = " + str(p_or[1]**2)
	print p_or
#	print sum(_CF_or_residuals(p_or,x,y)**2)
#	print numpy.corrcoef(y,_CF_or(p_or,x))[0][1]
	print "C(t) for complex residues Davi 2012 - blue - S2 = " + str(p_davi[1]**2)
	print p_davi
#	print sum(_CF_davi_residuals(p_davi,x,y)**2)
#	print numpy.corrcoef(y,_CF_davi(p_davi,x))[0][1]
#	print "S = " + str(S)

#	xp = numpy.linspace(x.min(), x.max(), numpy.floor((x.max()-x.min())*100)
	p_cf_gl_md4 = _CF_gl_md4(p_gl_md4, x)
	p_cf_gl = _CF_gl(p_gl, x)
	p_cf_ex = _CF_ex(p_ex, x)
	p_cf_davi = _CF_davi(p_davi, x)
	p_cf_or = _CF_or(p_or, x)
	
#	mp.plot(x, y, 'o', x, p_cf_ex, 'g-', x, p_cf_or, 'r-', x, p_cf_davi, 'b-', x, p_cf_gl, '-m', x, p_cf_gl_md4, '-k')
#	mp.grid(True)
#	mp.xscale('log')
#	mp.show()

	

#print "\nfirst 10 ns function"
#_get_OP_(x_chen[1:int(10/step)], c_func_chen[1:int(10/step)])
#print "\nfirst 20 ns function"
#_get_OP_(x_chen[1:int(20/step)], c_func_chen[1:int(20/step)])
#print "\nfirst 30 ns function"
#_get_OP_(x_chen[1:int(30/step)], c_func_chen[1:int(30/step)])
#print "\nfirst 40 ns function"
#_get_OP_(x_chen[1:int(40/step)], c_func_chen[1:int(40/step)])
print "\nFull function"
_get_OP_(x_chen, c_func_chen)


print "done"

#print len(dotprod)
#print len(c_func)
#print mol.trajectory
