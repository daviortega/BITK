#trying to play with classes.
import sys
import numpy
import copy

AA = 'ACDEFGHIKLMNPQRSTVYW-'

freq_bg = { 'A':0.073, 'C':0.025, 'D':0.050, 'E':0.061, 'F':0.042, 'G':0.072, 'H':0.023, 'I':0.053, 'K':0.064, 'L':0.089, 'M':0.023, 'N':0.043, 'P':0.052, 'Q':0.040, 'R':0.052, 'S':0.073, 'T':0.056, 'V':0.063, 'W':0.013, 'Y':0.033}

def read_seq(filename, msafmt= 'fasta'):

	seq_list = []
        tag_list = []
        seqlen = 0
        if msafmt == 'fasta':
        	data = open(filename,'r')
                for line in data:
                	if line[0] == '>':
                        	try:
                                	seq_list.append(seq(sequence))
                                except:
                                	pass
                                name = line[1:-1]
                                tag_list.append(name)
                                sequence = seq('')
                        else:
                                while line.find('\n') != -1:
                        	        line = line[0:line.find('\n')] + line[line.find('\n')+2:]
                      	        #seq_dic[name]= seq_dic[name] + line
                                sequence += line
		
		seq_list.append(seq(sequence))
		L = 0
		notAligned = False
                for i in range(len(seq_list)):
                        if L == 0 or len(seq_list[i]) == L:
                                L = len(seq_list[i])
                        else:
				notAligned = True
				break
	                
		if notAligned:
			print "At this time, this class only works with Multiple Sequence Alignments and therefore your sequences needs to have the same length. If you are sure you have input a MSA, the troubled sequence could be:"
			print tag_list[i]
		else:
			return msa(tag_list, seq_list, original = {}, orgpos = [])
	else:
		print "Please provide one of the following formats: fasta"
                return -1

class msa(object):
	def __init__(self, tag_list, seq_list, original = {}, orgpos = []):
		if original == {}:
			for i in range(len(tag_list)):
				original[tag_list[i]] = seq_list[i]
		if orgpos == []:
			self.orgpos = range(len(seq_list[0]))
		else:
			self.orgpos = orgpos
		self.org = original
		self.tag = tag_list
		self.seq = seq_list
		self.num = len(tag_list)
		col_out = []
                for i in range(len(seq_list[0])):
                        col_i = ''
                        for aa in seq_list:
                                col_i += aa[i]
                        col_out.append(col(col_i))
                self.col = col_out

		
	def __str__(self):
                return "MSA with %d sequences" % len(self.tag)

	def __getitem__(self, i):
		return [ self.tag[i], self.seq[i]]

	def __len__(self):
                return len(self.tag)

	def res_num2pos(self, res_num, tag, bias = 0):
		aa_count = -1
		res_num = res_num - bias
		for i in range(len(self.org[tag])):
			if self.org[tag][i] != '-':
				aa_count += 1
			if aa_count == res_num:
				break
		try:
			pos = self.orgpos.index(i)
			return pos
		except ValueError:
			print "The " + self.org[tag][i] + str(aa_count) + " from " + tag + " was excluded from the current version of the multiple sequence alignemnt (column " + str(i) + " of the original alignment)."
			return 0

	def pos2res_num(self, pos, tag, bias = 0):
		gaps = 0
#		removed = self.rem
		res_num = bias
#		pos = [ i for i in range(len(self.org.values()[0])) if i not in self.rem ][pos]
		pos = self.orgpos[pos]
		for aa in self.org[tag][:pos+1]:
			if aa != '-':
				res_num += 1

		return res_num
			
	def remove_empty_col(self, max_percent_gap = 1):
		current = copy.deepcopy(self.orgpos)
		tags = copy.deepcopy(self.tag)
		seqs = copy.deepcopy(self.seq)
		count_out = 0
		for i in range(len(self.col)):
			prev_aa = self.col[i].prevalent(nogaps = False)
			if prev_aa[0] >= max_percent_gap and prev_aa[1] == '-':
				for j in range(len(self)):
					seqs[j] = (seqs[j][:i-count_out] + seqs[j][i+1-count_out:])
				current.pop(i-count_out)
				count_out += 1
		return msa(tags, seqs, self.org, current)
	
	def remove_hc_col(self, max_freq = 1, verbose = False):
		current = copy.deepcopy(self.orgpos)
		tags = copy.deepcopy(self.tag)
                seqs = copy.deepcopy(self.seq)
                count_out = 0
                for i in range(len(self.col)):
                        prev_aa = self.col[i].prevalent(nogaps = False)
                        if prev_aa[0] >= max_freq:
                                for j in range(len(self)):
                                        seqs[j] = (seqs[j][:i-count_out] + seqs[j][i+1-count_out:])
                                current.pop(i-count_out)
                                count_out += 1
                return msa(tags, seqs, self.org, current)
					  


	def findseq(self, tag):
		return self.seq[self.tag.index(tag)] 

	def get_column(self, i):
		col = ''
		for aa in self.seq:
			col += aa[i]
		return col

	def msa_bin(self):
		new_tag = self.tag
		msabin = numpy.zeros( [len(self.tag), len(self.col)])
		for i in range(len(self.col)):
			A = self.col[i].prevalent()
			for j in range(self.num):
				if self.seq[j][i] == A[1] and self.seq[j][i] != '-':
					msabin[j][i] = 1
		return numpy.array(msabin)

	def rel_s(self, i, aa = 'all', weight = True, nogap = True):
		if weight:
			freq_here = freq_bg
		else:
			freq_here = {}
			for j in freq_bg.keys():
				freq_here[j] = 1/float(len(freq_bg.keys()))

		if nogap:
                        AA_here = AA[:-1]
			
                else:
                        AA_here = AA
			gap_num = 0
			if weight:
				for s in self.seq:
					gap_num += s.count('-')
				freq_gap = gap_num/float(self.num*len(self.col))
				for j in freq_bg.keys():
					freq_here[j] = ( 1 - freq_gap ) * freq_bg[j]
				freq_here['-'] = freq_gap
			else:
				for j in freq_here.keys() + ['-']:
					freq_here[j] = 1/float(len(freq_bg.keys())+1)
			print "Assuming f('-') = " + str(freq_here['-'])

		


		f_i_a = self.col[i].hist()
		if aa == 'all' :
			result = {}
			for a in AA_here:
				result[a] = f_i_a[a]*numpy.log(f_i_a[a]/freq_here[a])+ (1 - f_i_a[a])*numpy.log((1 - f_i_a[a])/(1 - freq_here[a]))
			return result
		else:
			return f_i_a[aa]*numpy.log(f_i_a[aa]/freq_bg[aa])+ (1 - f_i_a[aa])*numpy.log((1 - f_i_a[aa])/(1 - freq_bg[aa]))

	def SCA(self):
		prev_array = [c.prevalent() for c in self.col]
		freq_bg_bin = numpy.reshape(numpy.array([ freq_bg[c[1]] for c in prev_array ]), [len(prev_array), 1])
		freq_bin = numpy.reshape(numpy.array([ c[0] for c in prev_array ]), [len(prev_array),1])
		#Calculating Weights
		W = numpy.log(freq_bin*(1-freq_bg_bin)/(freq_bg_bin*(1-freq_bin)))
		#Calculating Correlation SCA matrix:
		msa_bin = self.msa_bin()
		freq_pairs_bin = numpy.dot(msa_bin.T, msa_bin)/len(self) 
		C_bin = freq_pairs_bin - numpy.dot(freq_bin, freq_bin.T)
		#return weigthed Cij
		return numpy.dot(W, W.T)*abs(C_bin)

	def D_bin(self):
		prev_array = [c.prevalent() for c in self.col]
                freq_bg_bin = numpy.reshape(numpy.array([ freq_bg[c[1]] for c in prev_array ]), [len(prev_array), 1])
                freq_bin = numpy.reshape(numpy.array([ c[0] for c in prev_array ]), [len(prev_array),1])
		D_bin = freq_bin* numpy.log(freq_bin/freq_bg_bin) + (1 - freq_bin) * numpy.log((1-freq_bin)/(1-freq_bg_bin))
		return D_bin	

	def SCA_clean(self, eig_list = [1, 2, 3]):
		Cij = self.SCA()
		eigval, eigvec = numpy.linalg.eig(Cij)
		#ordering
		idx = eigval.argsort()
		eigval = eigval[idx[::-1]]
		eigvec = eigvec[:,idx[::-1]]
#		return eigval, eigvec

		C_sca_clean = numpy.zeros([len(self.col), len(self.col)])
		for i in eig_list:
			C_sca_clean += eigval[i] * eigvec[:,i][numpy.newaxis].T * eigvec[:,i]
		#protein sector
		threshold = 0.05
		sec_green = numpy.array([ i for i in range(len(eigvec[:,1])) if eigvec[i,3] > max(0.05, abs(eigvec[i,1])) ])
		sec_blue = numpy.array([ i for i in range(len(eigvec[:,1])) if -eigvec[i,1] < -max(0.05, abs(eigvec[i,3])) ])
		sec_red = numpy.array([ i for i in range(len(eigvec[:,1])) if -eigvec[i,1] > max(0.05, abs(eigvec[i,3])) ])
	
		sec_green_ord = sec_green[eigvec[sec_green,3].argsort()[::-1]]
		sec_red_ord = sec_red[eigvec[sec_red,1].argsort()]
		sec_blue_ord = sec_blue[eigvec[sec_blue,1].argsort()[::-1]]
		
		sec_all_ord = numpy.concatenate( (sec_blue_ord, sec_green_ord, sec_red_ord))
		
		C_sca_clean_sectors = C_sca_clean[sec_all_ord][:,sec_all_ord]
		return C_sca_clean_sectors

	def SCA_ICA(self, eig_list = [1,2,3]):
		Cij = self.SCA()
                eigval, eigvec = numpy.linalg.eig(Cij)
                #ordering
                idx = eigval.argsort()
                eigval = eigval[idx[::-1]]
                eigvec = eigvec[:,idx[::-1]]
		W, c = ICA(eigvec, eig_list, 20000, 0.0001)
		
		C_sca_clean = numpy.zeros([len(self.col), len(self.col)])
                for i in eig_list:
                        C_sca_clean += eigval[i] * eigvec[:,i][numpy.newaxis].T * eigvec[:,i]

		C_ICA = numpy.dot(W, eigvec[:,eig_list].T)
		return C_ICA	
		#sec_1 = numpy.array([ i for i in range(len(C_ICA[:,1])) if C_ICA[i,1] > 

#def order_sca_clean_( C_sca_clean = numpy.array([]), threshold = 0.05 ):
	

		#Gotta program ICA!
	
def ICA(eigvec =numpy.array([]), eig_list = [1,2,3,4], Niter = 1000, r = 0.0001):
	#based on the MATLAB implementation by Olivier Rivoire 2010
	x = eigvec[:,eig_list]
	w = numpy.identity(x.shape[1])
	change = []
	
	for i in range(Niter):
		u = numpy.dot( w, x.T )
		dw = r * ( numpy.identity(x.shape[1]) + numpy.dot(numpy.dot( ( 1 - ( 2 / ( 1 + numpy.exp(-u)))), u.T), w))
		w += dw
		dw = numpy.reshape(dw, [1,dw.size])
		change. append(numpy.dot(dw, dw.T))
		#print change[-1]

	return w, numpy.array(change)

def p_sectors(MSA = 0, kmax = 3):
	print MSA.__class__
	if MSA.__class__ != msa :
		print "The function needs a msa as input and the number of eigvalues to be studied"
		return False
	eig_list = range(1, kmax)
	return MSA.SCA_ICA(eig_list)
	
	
	
			

#class columns(list):
#	def __init__(self, msas):
#		col_out = []
#                for i in range(len(msas.seq[0])):
#                        col_i = ''
#                        for aa in msas.seq:
#                                col_i += aa[i]
#                        col_out.append(col(col_i))
#                self.aas = col_out
		

#	def __call__(self):
#		return "%d columns" % len(self.aas)	

#	def __getitem__(self, i):
#		return self.aas[i]

#	def __iter__(self):
#		return self.aas[
	
#	def __len__(self):
#		return len(list(self.aas))
	
class seq(str):
	def nogap(self):
		out = copy.deepcopy(self)
		while '-' in out:
			out = out.replace('-','')
		return out

class col(str):
	
	def prevalent(self, nogaps = True):
		if nogaps:
			count_list = []
                	for aa in AA[:-1]:
        	                count_list.append(self.count(aa)/float(len(self)))
	                return [ max(count_list), AA[count_list.index(max(count_list))] ]
		else:
			count_list = []
			for aa in AA:
				count_list.append(self.count(aa)/float(len(self)))
			return [ max(count_list), AA[count_list.index(max(count_list))] ]

	def hist(self):
		count_dic = {}
		for aa in AA:
			count_dic[aa] = self.count(aa)/float(len(self))
		return count_dic
	



