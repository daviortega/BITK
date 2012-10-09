import sys
import random
import time
from math import *
from numpy import *
#from numpy import numarray
#from numarray import *
import numpy.numarray.linear_algebra as la
import numpy.linalg as la_new
from os import listdir
from os import getcwd
import os

#FROM python cookbock Memory monitor!
import os

_proc_status = '/proc/%d/status' % os.getpid()

_scale = {'kB': 1024.0, 'mB': 1024.0*1024.0,
          'KB': 1024.0, 'MB': 1024.0*1024.0}

def _VmB(VmKey):
    '''Private.
    '''
    global _proc_status, _scale
     # get pseudo file  /proc/<pid>/status
    try:
        t = open(_proc_status)
        v = t.read()
        t.close()
    except:
        return 0.0  # non-Linux?
     # get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
    i = v.index(VmKey)
    v = v[i:].split(None, 3)  # whitespace
    if len(v) < 3:
        return 0.0  # invalid format?
     # convert Vm value to bytes
    return float(v[1]) * _scale[v[2]]


def memory(since=0.0):
    '''Return memory usage in bytes.
    '''
    return _VmB('VmSize:') - since


def resident(since=0.0):
    '''Return resident memory usage in bytes.
    '''
    return _VmB('VmRSS:') - since


def stacksize(since=0.0):
    '''Return stack size in bytes.
    '''
    return _VmB('VmStk:') - since

#my own codes


def datagen(data_size=10, num_real_class=2):
	data = []
	counts = 0
	data_size = int(data_size)
	print str(num_real_class)
	for i in range(num_real_class):
		sigma = random.uniform(0,1)
        	mu = random.uniform(0,1)
		k = 0
		while len(data) < data_size * num_real_class:
			sample = random.gauss(mu,sigma)
			data.append(sample)
			k = len(data)
			print str(k) + ' --- ' + str(data_size)
	return data
	
def datagenG(data_size=10, num_real_class=2, num_dim=1):
        data = {}
	num_dim = int(num_dim)
	data_size = int(data_size)
	#for i in range(data_size*num_real_class):
	#	data[i] = []
	counts = 0
        print str(num_real_class)
	for i in range(num_real_class):
		sigma = random.uniform(0,2)
         	mu = random.uniform(-10,10)
		print 'New class'
        	while len(data.keys()) < data_size * (i + 1):
			sample = []
			for t in range(num_dim):
				samplei = random.gauss(mu,sigma)
				sample.append(samplei)
			print str(len(data.keys()))
			data[len(data.keys())] = sample
	return data




def dist_matrix1D(data=[]):
	if data == []:
		print 'Null dataset\n Quiting...'
		sys.exit()
	matrix = []
	for i in range(len(data)):
		matrix_line = []
		for j in range(len(data)-i-1):
			k = j + i + 1
			dist = abs(data[i] - data[k])
			matrix_line.append(dist)
		matrix.append(matrix_line)
	return matrix.append

def Umatrix(data,C,m):
	U = []
	m = float(m)
	for i in range(len(data)):
	        Ui = []
	        for Cj in C:
			Uij = 0
	                for Ck in C:
	                        Uij = Uij + (float(abs(data[i] - Cj))/float(abs(data[i] - Ck)))**(float(2)/float(m-1))
		        Uij = float(1)/float(Uij)
			Ui.append(Uij)
		U.append(Ui)
	return U

def UmatrixG(names,data,Cnames,C,m):
        U = []
	m = float(m)
        for xi in names:
		Ui = []
		for Cj in Cnames:
			dij = 0
			for t in range(len(data[xi])):
				dij = dij + (data[xi][t] - C[Cj][t])**2
			dij = dij ** (float(1)/float(m-1))
			Dik = 0
			for Ck in Cnames:
				dik = 0
				for t in range(len(data[xi])):
					dik = dik + (data[xi][t] - C[Ck][t])**2
				dik = dik ** (float(1)/float(m-1))
				Dik = float(Dik) + float(1)/float(dik)
			Uij = float(dij)*float(Dik)
			Uij = float(1)/float(Uij)
			Ui.append(Uij)
                U.append(Ui)
        return U


def centroides(data,U,m):
	C = []
	m = float(m)
	for j in range(len(U[0])):
		Ci_sup = 0
		Ci_inf = 0
		for i in range(len(data)):
			Ci_sup = Ci_sup + (data[i]*(U[i][j] ** m))
			Ci_inf = Ci_inf + (U[i][j] ** m)
		Ci = float(Ci_sup)/float(Ci_inf)
		C.append(Ci)
	return C

def centroidesG(names,data,U,m):
        C = {}
	Cnames = []
	for i in range(len(U[0])):
		C[i+1] = []
		Cnames.append(i+1)
	
	m = float(m)
	j = 0
	for c in Cnames:
		Cj = []
		for t in range(len(data[names[0]])):	
			Cj_sup = 0
			Cj_inf = 0
			for i in range(len(names)):
				Cj_sup = Cj_sup + (data[names[i]][t]*(U[i][j] ** m))
				Cj_inf = Cj_inf + (U[i][j] ** m)
			Cj_total = float(Cj_sup)/float(Cj_inf)
			Cj.append(Cj_total)
		C[c] = Cj
		j = j+1
	
	return C

def reorderDm(All):
	All_new = {}
	master = All[All.keys()[0]]
	All_new[All.keys()[0]] = master
	for file in All.keys()[1:]:
		slave = All[file]
		for i in range(len(master[0])):
			where = slave[0].index(master[0][i])
			tmp = slave[0][where]
			slave[0].pop(where)
			slave[0].insert(i,tmp)
			sample = []
			for name in slave[0]:
				tmp = array([slave[1][name][where]])
				slave[1][name] = concatenate((slave[1][name][:where],slave[1][name][where+1:]))
				slave[1][name] = concatenate((slave[1][name][:i],tmp,slave[1][name][i:]))
		All_new[file] = slave
	print All
	print All_new
	return All_new





def writeDmphy(names, Dmatrix, filename):
	output = '  ' + str(len(Dmatrix.keys())) + '\n'
	for name in names:
		output = output + name + ' '*(11-len(name))
		k = 0
		for val in Dmatrix[name]:
			val = str(val)[:7]
			if len(val) != 7:
				val = val + '0'*(7-len(val))
			output = output + ' ' + val
			k = k+1
			if k == 7:
				k = 0
				output = output + '\n' + ' '*11
		if k != 0:
			output = output + '\n'
		else:
			output = output[:-11]
	datafile = open(filename, 'w')
	datafile.write(output)
	datafile.close()
	
def readDm(datafile):
	Dmat_str = open(datafile,'r')
	Dmatrix = {}
	name = ''
	names = []
	for line in Dmat_str:
		if line[0] != ' ':
			if name != '':
				Dmatrix[name] = val
			name = line[0:11]
			name = name.replace(' ','')
			names.append(name)
			Dmatrix[name] = []
			val = []
		if name != '':
			list = line.split(' ') 
			for i in list:
				try:
					num = float(i)
					val.append(num)
				except ValueError:
					num = ''
	Dmatrix[name] = val
	Dmat_str.close()
	return names, Dmatrix

def readDmGEN(datafile, lim = '\t', trig_sup = 'n'):
        Dmat_str = open(datafile,'r')
        Dmatrix = {}
        name = ''
        names = []
        for line in Dmat_str:
                if line[0] != ' ':
                        if name != '':
                                Dmatrix[name] = val
                        name = line.split(lim)[0]
                        name = name.replace(' ','')
                        names.append(name)
                        Dmatrix[name] = []
                        val = []
                if name != '':
                        list = line.split(lim)
                        for i in list:
                                try:
                                        num = float(i)
                                        val.append(num)
                                except ValueError:
                                        num = ''
        Dmatrix[name] = val
        Dmat_str.close()
        return names, Dmatrix

def readDmGENlist(datafile, lim = '\t', trig_sup = 'n'):
        Dmat_str = open(datafile,'r')
        Dmatrix = []
        names = []
        for line in Dmat_str:
		line = line.replace('\n','')
		fields = line.split(lim)
#		fields.remove('')
		names.append(fields[0])
		vec = []
		for val in fields[1:]:
			try:
				vec.append(float(val))
			except ValueError:
				print fields[1:]
		Dmatrix.append(vec)
        Dmat_str.close()
        return names, Dmatrix

def read_coords(filename, lim = '\t'):
	datafile = open(filename,'r')
	names = []
	coords = []
	mac_format = 0
	for line in datafile:
		if '\r' in line:
			new_lines = line.split('\r')
			mac_format = 1
			break			
		fields = line.replace('\n','').split(lim)
		names.append(fields[0])
#		print names
		try:
			coords.append([float(s) for s in fields[1:]])
		except ValueError:
			print fields
			sys.exit()
	if mac_format == 1:
		for line in new_lines:
			fields = line.replace('\n','').split(lim)
	                names.append(fields[0])
#       	        print names
                	try:
                        	coords.append([float(s) for s in fields[1:]])
	                except ValueError:
        	                print fields
                	        sys.exit()
	
	coords = array(coords)
	return names, coords
	

def coord2(names, Dmatrix, datafile, export = 'n'):
        start = time.time()
	print '\n\n-> Computing coordinates <-\n\n'
#       A = array(zeros((1,len(names))))
	A = []
        Atime = time.time()
        print "Computing A"
        mA = memory()
        for x in range(len(names)):
#		Ai = array([])
		Ai = []
                for i in range(len(Dmatrix[x])):
                	#Ai = append(Ai,float(-0.5) * Dmatrix[x][i] * Dmatrix[x][i])
			Ai.append(float(-0.5) * Dmatrix[x][i] * Dmatrix[x][i])
		#A = vstack((A,Ai))
		A.append(Ai)
#	A = A[1:,]
	try:
		A = array(A)
	except ValueError:
		print len(A)
		for i in A:
			print len(i)
		sys.exit()
	mD = memory()
	print "Released by Dmatrix: " + str(memory(mD))
#	for i in range(len(names)):
#		print A[i][i]
#	sys.exit()
        Aio = A.sum(axis=1)/float(len(A))
        Aoj = A.sum(axis=0)/float(len(A))
        Aoo = A.sum()/float(len(A)*len(A))

	print "A took: " + str(memory(mA)) + " in bytes"
#	B = array(zeros((1,len(names))))
	B = []
#	print "-->Done\t" + str(time.time() - Atime) + "\nComputing B"
        Btime = time.time()
        for i in range(len(A)):
#		Bi = array([])
		Bi = []
                for j in range(len(A[i])):
			#Bij = round(float(A[i][j] - Aio[i] - Aoj[j] + Aoo),5)
                        #Bi = append(Bi, round(float(A[i][j] - Aio[i] - Aoj[j] + Aoo),5))
			Bi.append(float(A[i][j] - Aio[i] - Aoj[j] + Aoo))
		B.append(Bi) 
		#B = vstack((B,Bi))
	#B=B[1:,]
	B = array(B)
	if export == 'y':
		savetxt('Bmatrix.dat', B, delimiter = '\t',fmt = '%.4f')
#		fileout = open('Bmatrix.dat','w')
#		fileout.write(str(B))
#		fileout.close()
	
        #print str(B)
	print "--> Done Calculating B\t" + str(time.time() - Btime) + "\nI don't need A anymore"
        mA = memory()
        del A
        print "A released " + str(memory(mA))
        eigentime = time.time()
        print "-->Done\n\nSolving Eigenvector and eigenvalues"
        #L , G = la.eigenvectors(B)
	L, G = la_new.eigh(B)
        print "-->Done\t" + str(time.time() - eigentime) + "\n\nEigenvalues:"
	#mB = memory()
	#print "B released " + str(memory(mB))
	print L
	idx = L.argsort()
	L = L[idx]
	G = G[:,idx]
	L = L[::-1]
        G = G[:,::-1]
	print L[:10]
	#print sort(L)[::-1]
	print "Calculating percentage of variance"
	Tr = L[L>0].sum()
	Lp = L[L>0]/Tr
	print Lp[:10]*100
        for i in range(L.shape[0]):
	        if L[i] < 10**(-20):
                	L[i] = 0
        print "Transposing"
        #F = transpose(G)*(L**0.5)
	H = G*(L**0.5)
        print "-->Done"
	print "Done Calculating"
        Final = {}
        for i in range(len(names)):
                 Final[names[i]] = H[i]
        Final_3D = {}
        for i in range(len(names)):
                Final_3D[names[i]] = H[i][:3]
#       print Final
        return Final, Final_3D




def coord(names, Dmatrix, datafile):
	start = time.time()
	coordfile = '.' + datafile[:-4] + 'coord'

#	filenames = listdir(getcwd())
#	compute = 'N'
#	if coordfile in filenames:
#		names_file = []
#		G = []
#		input = open(coordfile,'r')
#		for line in input:
#			if line[0] == '*':
#				dev_file = float(line[1:])
#			else:
#				Gi = []
#				names_file.append(line[:line.find(',')])
#				rest = line[line.find(',')+1:]
#				while rest.find(',') != -1:
#					Gi.append(float(rest[:rest.find(',')]))
#					rest = rest[rest.find(',')+1:]
#				G.append(Gi)
#		
#		G = array(G)
#		print G
#		if names != names_file:
#			compute = 'Y'
#			print 'Wrong names'
#			print names
#			print names_file
#		else:
#			gg = DistMtx(G)
#			dev = 0
###			print gg
#			for i in range(len(names)):
#				for j in range(len(Dmatrix[names[i]])):
#					dev = dev + abs(gg[i][j] - Dmatrix[names[i]][j])
#			if dev +dev*0.1 < dev_file:
#				compute = 'Y'
#				print 'Too much error'
#				print str(dev)
#				print str(dev_file)
#	else:
#		compute = 'Y'
#	
	vai = 1
	if vai == 0:
		print '\n\n-> No need to recompute coordinates!!! <-\n\n'
	else:
		print str(time.time() - start)
		print '\n\n-> Computing coordinates <-\n\n'
		A = []
		Atime = time.time()
		print "Computing A"
		print 
		mA = memory()
		for name in names:
			Ai = []
			for i in range(len(Dmatrix[name])):
				Ai.append(float(-0.5) * Dmatrix[name][i] * Dmatrix[name][i])
			A.append(Ai)
		#print A
#		for i in range(len(names)):
#	                print A[i][i]
#		sys.exit()
		Anp = array(A)
		Aio = Anp.sum(axis=1)/float(len(A))
		Aoj = Anp.sum(axis=0)/float(len(A))
		Aoo = Anp.sum()/float(len(A)*len(A))
		del Anp
		print "A took: " + str(memory(mA))
		B = []
		#Aoo = 0
		#for i in range(len(A)):
		#	for j in range(len(A[i])):
		#		Aoo = Aoo + A[i][j]
		#Aoo = float(Aoo)/(len(A)**2)
		print "-->Done\t" + str(time.time() - Atime) + "\nComputing B"
		Btime = time.time()
		for i in range(len(A)):
#			print "\t" + str(i) + "/" + str(len(A)) + "\t" + str(time.time() - Btime)
			Bi = []
			for j in range(len(A[i])):
				#Aio = 0
				#Aoj = 0
				#for k in range(len(A[i])):
				#	Aio = Aio + A[i][k]
				#Aio = float(Aio)/float(len(A[i]))
				#print 'Aio for i =' + str(i) + ' : ' + str(Aio)
				#for k in range(len(A)):
				#	Aoj = Aoj + A[k][j]
				#Aio = float(Aio)/float(len(A[i]))
				#Aoj = float(Aoj)/float(len(A))
				#print 'Aoj for j =' + str(j) + ' : ' + str(Aoj)
				Bij = round(float(A[i][j] - Aio[i] - Aoj[j] + Aoo),5)
				Bi.append(Bij)
			B.append(Bi)
		print "--> Done Calculating B\t" + str(time.time() - Btime) + "\nI don't need A anymore"
		mA = memory()
		del A
		print "A released " + str(memory(mA))
		print "Let's make B an array"
		Bnp = array(B)
#		Dmat = zeros ((1,len(Dmatrix.keys())))
#		for Ba in B:
#			Dmat =  vstack((Dmat,Ba))
#		B = Dmat[1:,]
		print Bnp.size
		print "-->Done\n\nReleasing memory"
		del B
		eigentime = time.time()
		print "-->Done\n\nSolving Eigenvector and eigenvalues"	
		L , G = la.eigenvectors(Bnp)
		print "-->Done\t" + str(time.time() - eigentime) + "\n\nEigenvalues:"
		del Bnp
		L.sort()
		L = L[::-1]
		print L[:10]
		for i in range(L.shape[0]):
			if L[i] < 10**(-20):
				L[i] = 0

		
		print "Transposing"
		G = transpose(G)*(L**0.5)
		print "-->Done"	

#		Calculates the total deviation from the original matrix		
#		gg = DistMtx(G)
#		dev = 0
#		for i in range(len(names)):
#			for j in range(len(Dmatrix[names[i]])):
#				dev = dev + abs(gg[i][j] - Dmatrix[names[i]][j])
#		checking for matrix dimensionality possible reduction


		if '-red_dim' in sys.argv:
	
			old_dim = len(G[0])

			print 'Initiating dimensionality reduction'

			num_red = 0
			for j in range(len(G[0])):
				try:
					G[0][j]
				except IndexError:
					break
				k_num = 0
				for i in range(len(G)):
					if G[i][j] == 0.0:
						k_num = k_num + 1
				if k_num == len(G):
					print 'Reducing ' + str(j) + 'th coordinate'
					new_G = transpose(concatenate((transpose(G)[:j],transpose(G)[j+1:]),0))
					G = new_G
					num_red = num_red + 1
			print 'Total reduction of ' + str(num_red) +' from ' + str(old_dim) + ' to ' + str(len(G[0]))
			#Writing coordfile

		


#		data = '*'+ str(dev) + '\n'
#		for i in range(len(names)):
#			data = data + names[i]
#			for j in range(len(G[i])):
#				if j == 0:
#					data = data + str(G[i][j])
#				else:
#				data = data + ',' + str(G[i][j].real)
#			data = data + '\n'
#		output = open(coordfile,'w')
#		output.write(data)
#		output.close()
	
#	print '\nCoordinates'
#	print G
#	print '\nRecovering Distance Matrix'	
#	print '\n' + str(gg)
	print "Done Calculating"
	Final = {}
	for i in range(len(names)):
		Final[names[i]] = G[i]
	Final_2D = {}
	for i in range(len(names)):
		Final_2D[names[i]] = G[i][:3]
#	print Final
	return Final, Final_2D

def DistMtx(G):
	gg = zeros((G.shape[0],G.shape[0]))
	for i in range(G.shape[0]):
		for j in range(G.shape[0]):
			ggij = 0
			for k in range(G.shape[1]):
				ggij  = ggij + (G[i][k] - G[j][k]) ** 2
			gg[i][j] = ggij ** 0.5
	return gg

def readtable(datafile):
	table = {}
	data = open(datafile,'r')
	for line in data:
		line = line.split('\t')
		tax_name = line[0]
		class_name = line[1][:-1]
		if class_name not in table.keys():
			table[class_name] = [[tax_name]]
		else:
			table[class_name].append([tax_name])
	if 'Class' in table.keys():
		del table['Class']
	data.close()
	return table

def readresults(datafile, comp = 'N', file = 'Y'):
	results_comp = {}
	results = {}

	if file == 'Y':
		data = open(datafile,'r')
	else:
		data = datafile.split('\n')

	for line in data:
		if line == '':
			print 'No line'
		else:
			line = line.split(' ')
			k = 0
			while k == 0:
				try:
					line.remove('')
				except ValueError:
					k = 1
			k = 0
			while k == 0:
			        try:
					line.remove('---')
				except ValueError:
					k = 1
			tax_name = line[0]
			max_score = max(line[1:len(line)-1])
			max_index = line.index(max(line[1:len(line)-1]))
			if max_index not in results.keys():
				results_comp[max_index] = [{tax_name:max_score}]
				results[max_index] = [tax_name]
			else:
				results_comp[max_index].append({tax_name:max_score})
				results[max_index].append(tax_name)

	if file == 'Y':
		data.close()

	if comp == 'N':
		return results
	else:
		return results_comp

