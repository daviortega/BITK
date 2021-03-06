import fileinput
from math import *
import sys
import numpy
import md5
import base64
import ete2
import pymongo
import time
import json


def xtractAC(tag):
	if tag.split('-')[0].count('.') >= 2:
		try:
			AC = tag.split('-')[2]
			return AC
		except:
			print "Tag " + tag + " not in format for extraction"
			sys.exit()
	else:
		print "Tag " + tag + " is not a bitktag"
	return 0

def get_mist22_client():
	print "Verifying tunnels"
	print "Mist"
	try:
		client = pymongo.MongoClient('localhost',27019)
	except:
		print "You must open a tunnel with ares.bio.utk.edu: ssh -p 32790 -f -N -L 27019:localhost:27017 ortega@ares.bio.utk.edu"
		sys.exit()
	return client.mist22

def get_seqdepot_client():
	print "Verifying tunnels"
	print "SeqDepot"
	#try:                                    
		#client = pymongo.MongoClient('localhost',27018)
        client = pymongo.MongoClient("aphrodite.bio.utk.edu",27017)
        client.the_database.authenticate('binf',open("/home/ortegad/private/mongodb_binf.txt","r").readline().strip(), mechanism='MONGODB-CR',source='admin')
	#except:
        #        print "You must open a tunnel with ares.bio.utk.edu: ssh -p 32790 -f -N -L 27018:127.0.0.1:27017 ortega@aphrodite.bio.utk.edu"
	print "Authenticated"	
	return client.seqdepot

def insertmethodinfilename(filename, method):
	filenamenew = filename.split('.')
        filenamenew.insert(len(filenamenew)-1, method)
        filenamenew = '.'.join(filenamenew)
	return filenamenew

def mistid2strain(mid_list = []):
    mid_list = [ int(i) for i in mid_list ]
    result = {}
    mist = get_mist22_client()
    for mist_id in mid_list:
         res = mist.genomes.find_one({"_id":int(mist_id)}, {"n":1})
         if res != None:
                 print str(res['_id']) + '\t' + res['n']
                 result[res['_id']] = res['n']
         else:           
                 output += mist_id + '\n'
                 print str(mist_id) + " is not a valid MIST id"
    return result

def proteinid2bitktag(proteinid_list = []):
        out_dic = {}
        mist22 = get_mist22_client()
        genes = mist22.genes.find({'_id' : { '$in' : proteinid_list }}, {'gid' : 1, 'lo' : 1, '_id' : 1, 'p.ac':1})
        genomes = []
        for card in genes:
                genomes.append(card['gid'])
		try:
	                out_dic[card['_id']] = str(card['gid']) + '-' + str(card['lo']) + '-' + str(card['p']['ac'])
		except KeyError:
			out_dic[card['_id']] = str(card['gid']) + '-NULL-' + str(card['p']['ac'])
        genomes = list(set(genomes))
        genomes_mist = mist22.genomes.find({'_id': {'$in' : genomes}}, {'sp': 1, 'g' : 1})
        gid_dic = {}
        for card in genomes_mist:
                gid_dic[card['_id']] = card['g'][:2] + '.' + card['sp'][:3] + '.'
#       print gid_dic.keys()
        for proteinid in proteinid_list:
#                print out_dic[proteinid]
                out_dic[proteinid] = gid_dic[int(out_dic[proteinid].split('-')[0])] + out_dic[proteinid]
        return out_dic

def lc2ac(lc_list = []):
        print "Locus -> Accession"
        print lc_list
        ac_list = []
        mist22 = get_mist22_client()
        print "Got the client"
        genes = mist22.genes.find({'lo' : { '$in' : lc_list }}, {'p.ac':1})
        for card in genes:
            print card
            try:
                ac_list.append(card['p']['ac'])
            except:
                print "Found exception"
                print card
                sys.exit()
        return ac_list


def accession2bitktag(accession_list = [],):
        out_dic = {}
	preout_dic = {}
        mist22 = get_mist22_client()
        genes = mist22.genes.find({'p.ac' : { '$in' : accession_list }}, {'gid' : 1, 'lo' : 1, '_id' : 1, 'p.ac':1})
        genomes = []
	print "\tProcessing cards..."
        for card in genes:
		try:
                	genomes.append(int(card['gid']))
		except:
			print "Found exception"
			print card
			sys.exit()
                try:
                        preout_dic[card['p']['ac']] = str(card['gid']) + '-' + str(card['lo']) + '-' + str(card['p']['ac'])
                except KeyError:
                        preout_dic[card['p']['ac']] = str(card['gid']) + '-NULL-' + str(card['p']['ac'])
        genomes = list(set(genomes))
        genomes_mist = mist22.genomes.find({'_id': {'$in' : genomes}}, {'sp': 1, 'g' : 1})
        gid_dic = {}
        for card in genomes_mist:
                gid_dic[card['_id']] = card['g'][:2] + '.' + card['sp'][:3] + '.'
#       print gid_dic.keys()
	errors = []

        for accession in accession_list:
#                print out_dic[proteinid]
		try:
#			print preout_dic[accession]
        	       	out_dic[accession] = gid_dic[int(preout_dic[accession].split('-')[0])] + preout_dic[accession]
        	except KeyError:
			print "This sequence has - in the name and is messing up the data collection"
			print accession
#			sleep_counter(10)
			errors.append(accession)
			pass
#		except ValueError:
#			print "Something is weird here"
#			print accession
#			print preout_dic[accession]
#			sleep_counter(1)
#			errors.append(accession)
#			pass
	return out_dic, errors

def aseq2bitktaglist(aseqs):
	print "If this is breaking your program, notice that the output value of the dictionary is a list not a string"
	mist22 = get_mist22_client()
	cards = mist22.genes.find( {'p.aid' : { '$in' :  aseqs}},)
	aseq2ac = {}
	aclist = []
	for card in cards:
		aclist.append(card['p']['ac'])
		if card['p']['aid'] not in aseq2ac.keys():
			aseq2ac[card['p']['aid']] = [ card['p']['ac'] ]
		else:
			aseq2ac[card['p']['aid']].append( card['p']['ac'] )
			
	ac2tag = accession2bitktag(aclist)[0]
	result = {}
	for aseq in aseqs:
		try:
			for ac in aseq2ac[aseq]:
				if aseq not in result.keys():
					result[aseq] = [ ac2tag[ac] ]
				else:
					result[aseq].append(ac2tag[ac])
		except KeyError:
			result[aseq] = ['None']
	return result

def accession2_id(aclist = []):
    ac2id = {}
    mist22 = get_mist22_client()
    cards = mist22.genes.find({'p.ac' : { '$in' : aclist}})
    for card in cards:
        if card['p']['ac'] not in ac2id.keys():
            ac2id[card['p']['ac']] = card['_id']
        else:
            print "same accession number in two different proteins... something is weird"
            print card['p']['ac']
            print card
    for ac in aclist:
        if ac not in ac2id.keys():
            ac2id[ac] = "not_in_mist"
    return ac2id

def accession2seq(aclist = []):
	mist22 = get_mist22_client()
	sd = get_seqdepot_client()
	cards = mist22.genes.find({'p.ac' : { '$in' : aclist }})
	aseq2ac = {}
	cards = [card for card in cards]
	for card in cards:
		if card['p']['aid'] not in aseq2ac.keys():
			aseq2ac[card['p']['aid']] = [card['p']['ac']]
		else:
			aseq2ac[card['p']['aid']].append(card['p']['ac'])
	if len(cards) != len(set(aclist)):
	        print "There are " + str(len(aclist)) + " accession numbers but only " + str(len(cards)) + " in the database"
	        print set(aclist).difference(set([ ac['p']['ac'] for ac in cards]))
		sleep_counter(10)

	cards = sd.aseqs.find({'_id': { '$in' : aseq2ac.keys() }})
	cards = [card for card in cards]
	ac2seq = {}
	for card in cards:
		for ac in aseq2ac[card['_id']]:
                    print card
                    try:      
			ac2seq[ac] = card['s']
                    except KeyError:
                        print card
                        pass
	return ac2seq

def accession2md5(aclist = []):
        out_dic = {}
        mist22 = get_mist22_client()
        sd = get_seqdepot_client()
        cards = mist22.genes.find({'p.ac' : { '$in' : aclist }})
        ac2aseq = {}
        cards = [card for card in cards]
        for card in cards:
                if card['p']['aid'] not in ac2aseq.keys():
                        ac2aseq[card['p']['ac']] = [card['p']['aid']]
                else:
                        ac2aseq[card['p']['ac']].append(card['p']['aid'])
        if len(cards) != len(aclist):
                print "There are " + str(len(aclist)) + " accession numbers but only " + str(len(cards)) + " in the database"
                print set(aclist).difference(set([ ac['p']['ac'] for ac in cards]))
                sleep_counter(10)
	return ac2aseq




def alnreader(datafile, list = 'no', just_name = 'no'):
#	"""Reads a clustalw format file and returns the information in a dictionary name:sequences"""
	data = open(datafile, 'r')
	seq_dic = {}
	name_list = []
	
	for line in data:
		if line.find('CLUSTAL') == -1 and line[0] != ' ' and line[0] != '<':
			if just_name == 'Yes':
				if line.find('/') != -1:
					name = line[:line.find('/')]
				else:
					name = line[:line.find(' ')]
			else:
				name = line[:line.find(' ')]
#			name_list.append(name)
			rest = line[line.find(' '):]
			seq  = ''
			for i in rest:
				if i != ' ' and i != '\n':
					seq = seq + i
			if name in seq_dic.keys():
				seq_dic[name] = seq_dic[name] + seq
			elif name != '':
				seq_dic[name] = seq
				name_list.append(name)
	dataout = ''
	for k, v in seq_dic.iteritems():
		dataout = dataout + '>' + k + '\n' + v + '\n'
#	output = open(datafile[:-3] + 'dav', 'w')
#	output.write(dataout)
#	output.close()
	data.close()
	if list == 'Yes':
		return name_list
	else:	
		return seq_dic

def vissareader(datafile, just_name = 'no'):
#	"""Reads the vissa format files and return the information in a dictionary name:sequences"""
	data = open(datafile, 'r')
        seq_dic = {}
        for line in data:
                if line.find('CLUSTAL') == -1 and line[0] != ' ' and line[0] != '<':
			if just_name == 'Yes':
				name = line[:line.find('/')]
			else:
                        	name = line[:line.find(' ')]
                        rest = line[line.find(' '):]
                        seq  = ''
			i = 0
			start = 0
                        while rest[i] == ' ':
                        	i = i+1
			seq = rest[i:]
			if name in seq_dic.keys():
                                seq_dic[name] = seq_dic[name] + seq
                        else:
                                seq_dic[name] = seq
#        print seq
#	print start
	dataout = ''
        for k, v in seq_dic.iteritems():
                dataout = dataout + '>' + k + '\n' + v + '\n'
#       output = open(datafile[:-3] + 'dav', 'w')
#       output.write(dataout)
#       output.close()
        data.close()
        return seq_dic


def phreader(datafile):
#	"""Reads a .ph file and return a list with the sequences in the same order as the .ph file"""
        data = open(datafile, 'r')
        seq_list = []
	for line in data:
		if line[0] != '(' and line[0] != ':':
	        	name = line[:line.find(':')]
			seq_list.append(name)
        data.close()
        return seq_list


def nwkreader2(datafile):
	seq_list = []
	tr = ete2.Tree(datafile)
	for n in tr:
		seq_list.append(n.name)
	return seq_list


def nwkreader(datafile,just_name='no'):
#	"""Reads a .nwk file and return a list with the sequences in the same order as the .nwk file"
	data = open(datafile,'r')
	seq_list = []
	forb_list = ['(', ')', ':', ',', '.', '-', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
	for line in data:
		i = 0
		while line[i] != ';':
			if line[i] in forb_list:
				i = i+1
			else:
				name = ''
				while line[i] != ':' and line[i] != ')' and line[i] != '(' and line[i] != ',' and line[i] != ';' and i <= len(line)-1:
					name = name + line[i]
					i = i+1
				if just_name == 'Yes':
					seq_list.append(name.split('/')[0])
				else:
					seq_list.append(name)
#				print name
	data.close()
	return seq_list

def aln2vissa(datafile):
#	"""If aln2ss_xml complains about the format of your MA file and you are \
#		sure it is in clustalw format, this program will solve teh problem."""
	data = open(datafile,'r')
	output = ''
	for line in data:
		if line.find('CLUSTAL') == -1 and line[0] != ' ':
			output = output + line[:line.find(' ')] + '     ' + line[line.find(' '):]
	return output

# def fastareader(datafile, just_name = 'no'):
# 	data = open(datafile,'r')
# 	seq_dic = {}
# 	list_order = []
# 	for line in data:
# 		line = line.replace('\r','')
# 		if line[0] == '>':
# 			if just_name == 'Yes':
# 				name = line.split('/')[0][1:]
# 			else:
# 				name = line[1:-1]
# 			list_order.append(name)
# 			seq_dic[name] = ''
# 		else:
# 			while line.find('\n') != -1:
# 				line = line[0:line.find('\n')] + line[line.find('\n')+2:]
# 			seq_dic[name]= seq_dic[name] + line # Let's see if this comment breaks too many things ".replace(' ','')"
# 	dataout = ''
# 	for k, v in seq_dic.iteritems():
# 	 	dataout = dataout + '>' + k + '\n' + v + '\n'
	
# 	data.close()
# 	return seq_dic, list_order

def fastareader(datafile):
	
	ListOrder = []
	SeqDic = {}

	fh = open(datafile, 'r')
	
	fastaBuffer = None

	line = fastaBuffer if fastaBuffer else fh.readline()
	while line:
		if line[0] != '>':
			print "here"
			raise Exception('Invalid FASTA file. Header line must begin with a greater than symbol\nLine: ' + line + '\n\n')
		name = line[1:-1]
		ListOrder.append(name)
		SeqDic[name] = ''
		line = fh.readline()
		while line:
			if line[0] != '>':
				SeqDic[name] += line
				line = fh.readline()
				continue

			fastaBuffer = line
			break

		SeqDic[name] = SeqDic[name].replace('\n','')
	return SeqDic, ListOrder		

def alnwriter(seq_dic, list_order = []):
	if list_order == []:
		list_order = seq_dic.keys()

	#Determinning the maximum sequence lenght and the longer name
	max_len = 0
	max_nam = 0
	for k, v in seq_dic.iteritems():
		if len(v) > max_len:
			max_len = len(v)
		if len(k) > max_nam:
			max_nam = len(k)

	#Adjusting the sequences for the same size

	for k,v in seq_dic.iteritems():
		if len(v) < max_len:
			v = v + (max_len-len(v))*'-'
			seq_dic[k] = v
			#print k + '>' + v

	groups = int(floor(max_len/60+1)) #numbers 

	output = 'CLUSTAL W (1.83) multiple sequence alignment\n\n\n'

	for i in range(groups):
		for name in list_order:
		#for k, v in seq_dic.iteritems():
			output = output + name + (max_nam+5-len(name))*' ' + seq_dic[name][i*60:(i+1)*60] + '\n'
		output = output + '\n'
	#print output
	return output

def sleep_counter(N):
	sys.stdout.write('\033[91m' + '\033[1m' + 'Waiting')
	for i in range(N):
		sys.stdout.write('...' + str(i+1))
		time.sleep(1)
		sys.stdout.flush()
	print '\033[0m'
	


def pirwriter(seq_dic):
	output = ''
	for name, seq in seq_dic.iteritems():
		output += '>P1;' + name + '\nsequence:' + name + ':::::::0.00: 0.00\n' + seq + '*\n'
	
	return output
	
def nogap(seq):
	seq = seq.replace('-','')
	return seq		

def nogaps(seq):
	if seq.__class__ == dict:
		for name, line in seq.iteritems():
			seq[name] = line.replace('-','')
		return seq
	else:
		print("No deal. bitk.nogaps only accepts dict")
		return -1

def selectseq(seq_dic, seq_list):
	# select the sequences with names on seq_list from seq_dic and return a new dictionary
	
	new_seq_dic = {}
	nothere = 'n'
	counts = 0

	for name in seq_list:
		nothere = 'y'
		for seq_name in seq_dic.keys():
			if name in seq_name:
				new_seq_dic[name] = seq_dic[seq_name]
				nothere = 'n'
		if nothere == 'y':
			print 'Sequence in list and not in alignment ---->   ' + name 
			counts = counts + 1
	print 'Total not found: ' + str(counts)
	return new_seq_dic

def ESCL(seq_dic1, seq_dic2, seq_seq = {}, norm = 'no'):

	print 'ESCL Calculations - Start'

	score_table = []
	
	### The seq_seq matrix is a matrix that correlates the seq2 tags with seq1 tags in the following fashion:
	### seq_seq[seq1_tag] = [seq2_tag1, seq2_tag2, ...] if seq1_tag and seq2_tag1,2,3, ... are predicted to interact
	
	if seq_seq == {}:
		for loci in seq_dic1.keys():
			seq_seq[loci]=[loci]
		print 'Self interaction !!!!!!!'
	
	
	#print seq_seq

	fac = lambda n:n-1 + abs(n-1) and fac(n-1)*long(n) or 1

	AA = {'A':0,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0,'K':0,'L':0,'M':0,'N':0,'P':0,'Q':0,'R':0,'S':0,'T':0,'V':0,'Y':0,'W':0,'-':0}

	AA = 'ACDEFGHILMNPQRSTVYW-'

	length1 = 0
	for seq in seq_dic1.values():
		if len(seq) > length1:
			length1 = len(seq)
	
	length2 = 0
	for seq in seq_dic2.values():
		if len(seq) > length2:
      			length2 = len(seq)
	length1 = range(length1)
	length2 = range(length2)

	#print str(length2) + '  ' + str(length1)

	for i in length1: #i position in molecule 1
		
		#Selection
		AA = {'A':0,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0,'K':0,'L':0,'M':0,'N':0,'P':0,'Q':0,'R':0,'S':0,'T':0,'V':0,'Y':0,'W':0, '-':0}
		for seq1 in seq_dic1.values():
			AA[seq1[i]] = AA[seq1[i]] + 1
		most_cons_value = max(AA.values())
		for AA1 in AA.keys():
			if AA[AA1] == most_cons_value:
				most_cons = AA1

		select_loci1 = []
		for loci1, seq1 in seq_dic1.iteritems():
			if seq1[i] == most_cons:
				select_loci1.append(loci1)
		
		#print 'Column ' +str(i) + ' -->> ' + most_cons + '/' + str(most_cons_value)
		# No need for one to one list of anything. The list are uncoupled once the selected list is built
		select_loci2 = []
		for loci1 in select_loci1:
			for loci2 in seq_seq[loci1]:
				select_loci2.append(loci2)
		
		print select_loci2

		seq_dic2_sel = {}
		for loci2, seq2 in seq_dic2.iteritems():
			if loci2 in select_loci2:
				seq_dic2_sel[loci2] = seq2
				#print 'Unselected: ' + loci2

		#print str(len(seq_dic2.keys()))
		#print '\n\n' + str(len(seq_dic2_sel.keys())) + '\n\n'
		print seq_dic2_sel.keys()
		score_col = []	
		if most_cons == '-' and AA[most_cons] > len(seq_dic2.keys())/2:
			print 'Column ' + str(i) + ' eliminated'
			score_col = [ 0 ] * len(length2)
		else:
			#Calculations
			for j in length2:
				#N
				print 'Column ' + str(j) + ' from second set'
				AA_N = {'A':0,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0,'K':0,'L':0,'M':0,'N':0,'P':0,'Q':0,'R':0,'S':0,'T':0,'V':0,'Y':0,'W':0, '-':0}
				for seq2 in seq_dic2.values():
					AA_N[seq2[j]] = AA_N[seq2[j]] + 1
					#print seq2
				
				print AA_N
				for amino in AA_N.keys():
					if AA_N[amino] > AA_N['-']:
						flag = 'ok'
						break
					else:
						flag = 'notok'
				if flag == 'notok':
					print 'Column ' + str(j) + 'eliminated'
					score_col.append(0)
				else:
					#n
					AA_n = {'A':0,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0,'K':0,'L':0,'M':0,'N':0,'P':0,'Q':0,'R':0,'S':0,'T':0,'V':0,'Y':0,'W':0, '-':0}
					for seq2 in seq_dic2_sel.values():
        	                		 AA_n[seq2[j]] = AA_n[seq2[j]] + 1
					sum_n = 0
					for count in AA_n.values():
						sum_n = sum_n + count
				
					print AA_n
	
					#N/ntotal
					AA_Nntotal = fac(len(seq_dic2.keys()))/(fac(len(seq_dic2_sel.keys()))*fac((len(seq_dic2.keys())-len(seq_dic2_sel.keys()))))


					#m
					sum_m = 0
					sum_m_old = 0
					AA_m_old = {}
					AA_m = {'A':0,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0,'K':0,'L':0,'M':0,'N':0,'P':0,'Q':0,'R':0,'S':0,'T':0,'V':0,'Y':0,'W':0,'-':0}
					for AAm in AA_m.keys():
						div = AA_N[AAm]*float(len(seq_dic2_sel.keys()))/len(seq_dic2.keys())
						AA_m[AAm] = round(div)
					for count in AA_m.values():
						sum_m = sum_m + count
					#print AA_m
					#print str(sum_n) + ' - - - ' + str(sum_m)
					while sum_m != sum_n:
						if sum_m > sum_n:
							AAm_max_value = max(AA_m.values())
							for AAm in AA_m.keys():
								if AA_m[AAm] == AAm_max_value:
									AA_m[AAm] = AA_m[AAm] - 1
									break
						else:
							 AAm_min_value = min(AA_m.values())
                		                         for AAm in AA_m.keys():
                        		                         if AA_m[AAm] == AAm_min_value:
                                		                         AA_m[AAm] = AA_m[AAm] + 1
									 break
						sum_m = 0
						for count in AA_m.values():
	        		                        sum_m = sum_m + count
						
						print str(sum_n) + ' - - - ' + str(sum_m)
	
					print AA_m
	
					#(N/n)
					AA_Nn = {'A':0,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0,'K':0,'L':0,'M':0,'N':0,'P':0,'Q':0,'R':0,'S':0,'T':0,'V':0,'Y':0,'W':0,'-':0}
				
					for AAn in AA_Nn.keys():
						AA_Nn[AAn] = fac(AA_N[AAn])/(fac(AA_n[AAn])*fac(AA_N[AAn]-AA_n[AAn]))
				
					#print AA_Nn
					#print str(fac(0))
					#(N/m)
					AA_Nm = {'A':0,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0,'K':0,'L':0,'M':0,'N':0,'P':0,'Q':0,'R':0,'S':0,'T':0,'V':0,'Y':0,'W':0,'-':0}
	
        	        	        for AAm in AA_Nm.keys():
                	        	        AA_Nm[AAm] = fac(AA_N[AAm])/(fac(AA_m[AAm])*fac(AA_N[AAm]-AA_m[AAm]))
	
					print AA_Nm	
					#Lambda
				
					L = 1
					for AAj in AA_Nn.keys():
						L = L * float(AA_Nn[AAj])/float(AA_Nm[AAj])
					print str(L)
	
					L = -1 * log(L)
	
					#print str(L)
					score_col.append(L)
				
		score_table.append(score_col)
	if norm == 'Yes':

		score_table = normtable(score_table)

#		max_score = 0
#		for i in length1:
#			if max(score_table[i]) > max_score:
#				maxi = i
#				maxj = score_table[i].index(max(score_table[i]))
#				max_score = max(score_table[i])
#		norm_score_table = []
#		for i in length1:
#			norm_line = []
#			for j in length2:
#				norm_line.append(float(score_table[i][j])/float(score_table[maxi][maxj]))
#			norm_score_table.append(norm_line)
#		print 'Norm'
#		print str(max_score) + ' ( ' + str(maxi) + ' , ' + str(maxj) + ' )' + ' ----- ' + str(score_table[maxi][maxj]) + ' ---- ' + str(norm_score_table[maxi][maxj])
#		score_table = norm_score_table


	#print score_table
	print 'ESCL Calculations - End'
	print str(len(score_table)) + ' - - - ' + str(length1)
	print str(len(score_table[1])) + ' - - - ' + str(length2)

	return score_table

def ESCL2 (seq_dic1, seq_dic2, seq_seq = {}, norm = 'no'):

	if seq_seq == {}:
                for loci in seq_dic1.keys():
                        seq_seq[loci]=[loci]
                print 'Self interaction !!!!!!!'

	empty = {}
	length1 = 0
        for seq in seq_dic1.values():
                if len(seq) > length1:
                        length1 = len(seq)
        length2 = 0
        for seq in seq_dic2.values():
                if len(seq) > length2:
                        length2 = len(seq)
	if norm == 'no':
		score_table = ESCL(seq_dic1, seq_dic2, seq_seq, norm='no')
		self_score_table = ESCL(seq_dic2, seq_dic2, empty, norm='no')
	else:
		score_table = ESCL(seq_dic1, seq_dic2, seq_seq, norm='Yes')
		self_score_table = ESCL(seq_dic2, seq_dic2, empty, norm='Yes')

	new_table = []
 	for i in range(length1):
        	new_table_line = []
                for j in range(length2):
			new_table_line.append(score_table[i][j])
		new_table.append(new_table_line)
        for i in range(length1):
	        for j in range(length2):
		        for k in range(length2):
			        if k != j:
				        new_table[i][j] = new_table[i][j] - self_score_table[j][k] * score_table[i][k]
	
	score_table = new_table
	return score_table

def normtable(score_table):

	max_score = 0
        length1 = len(score_table)
	length2 = len(score_table[1])
	for i in range(length1):
		if max(score_table[i]) > max_score:
			maxi = i
			maxj = score_table[i].index(max(score_table[i]))
	                max_score = max(score_table[i])
        norm_score_table = []
        for i in range(length1):
	        norm_line = []
	        for j in range(length2):
		        norm_line.append(float(score_table[i][j])/float(score_table[maxi][maxj]))
		norm_score_table.append(norm_line)
	print 'Norm'
	print str(max_score) + ' ( ' + str(maxi) + ' , ' + str(maxj) + ' )' + ' ----- ' + str(score_table[maxi][maxj]) + ' ---- ' + str(norm_score_table[maxi][maxj])
	score_table = norm_score_table
	return score_table


def renormtable(score_table):
	max_score = 0
        length1 = len(score_table)
        length2 = len(score_table[1])
	for i in range(length1):
                if abs(max(score_table[i])) > max_score:
			maxi = i
                        maxj = score_table[i].index(max(score_table[i]))
                        max_score = abs(max(score_table[i]))
		else:
			if abs(min(score_table[i])) > max_score:
			        maxi = i
	                	maxj = score_table[i].index(min(score_table[i]))
	                        max_score = abs(min(score_table[i]))
	renorm_score_table = []
	for i in range(length1):
                renorm_line = []
                for j in range(length2):
                        renorm_line.append(1+(float(score_table[i][j])/float(abs(score_table[maxi][maxj]))))
		renorm_score_table.append(renorm_line)
	score_table = renorm_score_table
	return score_table
	
	
def check_len(seq_dic = {}):
	for i in range(1, len(seq_dic.keys())):
		if len(seq_dic[seq_dic.keys()[i-1]]) != len(seq_dic[seq_dic.keys()[i]]):
			L = -1
			break
		else:
			L = len(seq_dic[seq_dic.keys()[i]])
	return L

def concat(seq_dic1={},seq_dic2={},seq_seq={}):
	final_dic = {}

	if seq_seq == {}:
		for loci in seq_dic1.keys():
			seq_seq[loci] = [loci]


	for loci1 in seq_seq.keys():
		seq1 = seq_dic1[loci1]
		seq2 = seq_dic2[seq_seq[loci1]]
		final_dic[loci1] = seq1 + seq2
	
	return final_dic

def check_loci(seq_list1=[],seq_list2=[]):
	"""Check both lists to find unique members. returns list of common members"""
	unique = []
	common = []
	for loci in seq_list1:
		if loci not in seq_list2:
			print 'Not in list 2 > ' + loci
			unique.append(loci)
		else:
			common.append(loci)
	
	for loci in seq_list2:
		if loci not in seq_list1:
			print 'Not in list 1 > ' + loci
			unique.append(loci)

	return common

def Identity(seq1='', seq2=''):
	""" Calculate the identity between two sequences """
	errors = 0
	total = 0
	for i in range(len(seq1)):
		if seq1[i] == seq2[i] and (seq1[i] != '-' and seq2[i]!='-'):
			total += 1
		elif (seq1[i] != '-' or seq2[i] != '-'):
			total += 1
			errors += 1
	return (float(total - errors)/float(total))

def buildIDmatrix(seq_dic={}, seq_list = [], formatliketrisup = False):
	""" Build the identity matrix between the sequences in the dictionary """
	if seq_dic != {} and seq_list == []:
		seq_list = seq_dic.keys()
	
	namelist = []
	count = 0
	
	if formatliketrisup:
		output = ''
		for i in range(len(seq_list)):
			output += seq_list[i] + ';'*i
	                for j in range(i,len(seq_list)):
				output += ';%03f' % Identity(seq_dic[seq_list[i]],seq_dic[seq_list[j]])
			output += '\n'
			print 'Done with ' + seq_list[i] + '\t\t' + str(i+1) + '/' + str(len(seq_list))
		return output
	else:
		IDMatrix = {}
		for name1, seq1 in seq_dic.iteritems():
                	count += 1
	                namelist.append(name1)
        	        if name1 not in IDMatrix.keys():
                	        IDMatrix[name1] = {}
	                for name2, seq2 in seq_dic.iteritems():
        	                if name2 not in namelist:
                	                ID = Identity(seq1,seq2)
                        	        IDMatrix[name1][name2] = ID
                                	if name2 not in IDMatrix.keys():
        	                                IDMatrix[name2] = {}
	                                IDMatrix[name2][name1] = ID	
			print 'Done with ' + name1[:-1] + '\t\t' + str(count) + '/' + str(len(seq_dic))
		return IDMatrix


def MatrixPathAnt(IDMatrix={},name1='',name2='', Nants = 100, cutoff = 0.3):
	import random as rand
	import sys
	import time
	""" Find the optimum or quase optimum sequence path between two sequences in a alignment using ant algorithm (needs to be rewrite) """
	Ants_paths = []
	rand.seed()
	for i in range(Nants):
		Ants_paths.append([name1])
	print len(Ants_paths)
	winner = 0
	firstfind = -1
	#building a pheronome triangular matrix
	pheronome = {}
	chosen_path = {}
	for seq1 in IDMatrix.keys():
		pheronome[seq1] = {}
		for seq2 in IDMatrix.keys():
			if seq1 != seq2:
				pheronome[seq1][seq2] = 1
	sim_steps = 0

	#starting the search
	while winner == 0:
		print 'Starting the step ' + str(sim_steps) + ' of the simulation'
		norm = 0
		for seq1 in IDMatrix.keys():
		 	for seq2 in IDMatrix[seq1].keys():
				pheronome[seq1][seq2] = 0.7 * pheronome[seq1][seq2] #pheronome evaporation 
		#		norm += pheronome[seq1][seq2] * IDMatrix[seq1][seq2]
		
		#print norm
	
		for ant in range(len(Ants_paths)):
			#print '\nRunning ant number ' + str(ant)
			#print 'This is its path so far: ' + str(Ants_paths[ant])
			last = Ants_paths[ant][-1]
			#print last
			#finding the possible steps
			possibles = []
			for name,id in IDMatrix[last].iteritems():
				#print name + '\t' + str(id)
				if len(Ants_paths[ant]) > 1:
					a = -2
				else:
					a = -1
				if id >= cutoff and name not in Ants_paths[ant]: 
					possibles.append(name)
			#print 'This is the possible steps: ' + str(possibles)
			if len(possibles) == 0:
				Ants_paths[ant] = [name1]
				chosen = ''
			elif name2 in possibles:
				chosen = name2
				Ants_paths[ant].append(chosen)
			else:
				#update the probabilities with the pheronome information
				prob = {}
				prev = 0
				#print 'Probabilities'
				#normalization parameter
				#norm = 0
				#for seq1 in IDMatrix.keys():
				#	for seq2 in IDMatrix[seq1].keys():
				#		norm += pheronome[seq1][seq2] * IDMatrix[seq1][seq2]
			
				for name in possibles:				
					#P = float(pheronome[last][name]*IDMatrix[last][name])/float(norm)
					prob[name] = int(float(pheronome[last][name]*IDMatrix[last][name])*100) + prev
					prev = prob[name]
					#print name + '\t' + str(prob[name])
				#time.sleep(3)
				#Pick the next step
				next = ''
				i = rand.randint(0,prev)
				for name in possibles:
					if i <= prob[name]:
						chosen = name
						break
				#print 'From those we choose :' + chosen
				Ants_paths[ant].append(chosen)
				#print 'Now the path of the ant number ' + str(ant) + ' is: ' + str(Ants_paths[ant])
				#Check if this step reaches the point
			if chosen == name2:
				#print 'The ant number ' + str(ant) + ' made it with: ' + str(Ants_paths[ant]) 
				chosen_path[ant] = []
				for seq in Ants_paths[ant]:
					chosen_path[ant].append(seq)
				Ants_paths[ant] = [name1]
				firstfind = chosen
				#time.sleep(5)
		#updating pheronome matrix:
		for ant in range(len(Ants_paths)):
			if len(Ants_paths[ant]) > 1:
				i = Ants_paths[ant][-2]
				j = Ants_paths[ant][-1]
				pheronome[i][j] += IDMatrix[i][j]

		#print chosen_path
		#check if there is a winner path 90% of the ants on it
		if firstfind != -1:
			counts = {}
			for paths in chosen_path.values():
				counts.setdefault(str(paths),0)
				counts[str(paths)] += 1
			#print counts

			max = 0
			for path in counts.keys():
				if counts[path] > max:
					max = counts[path]
					best = path
			print 'Most common path shared by ' + str(max) + ' ants has ' + str(best.count(',')+1) + ' jumps.' 
			print best
			if max >= 0.75*Nants:
				winner = eval(best)
				print " Winner: " + str(winner)
				for i in range(len(winner)-1):
					print str(i) + ' - ' + winner[i] + ' --> ' + winner[i+1] + ':\t' + str(IDMatrix[winner[i]][winner[i+1]] )
				break
		else:
			print Ants_paths[0]
			print Ants_paths[1]
		#test if cutoff need to be lower to allow at least 1 connection
		for ant in range(len(Ants_paths)):
			if len(Ants_paths[ant]) != 1:
				go = 1
				break
			else:
				go = 0
		if go == 0:
			print 'You gotta lower the cutoff'
			return 1

		sim_steps += 1
		print "\nStep number: " + str(sim_steps)
		#time.sleep(2)

	return winner

def aa_hist(seq_dic={}):
	"""Calculates the counts of each AA in each column of a given MSA. The results is a list of histograms)""" 
	results = []
	for i in range(len(seq_dic.values()[0])):
		AA_Nm = {'A':0,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0,'K':0,'L':0,'M':0,'N':0,'P':0,'Q':0,'R':0,'S':0,'T':0,'V':0,'Y':0,'W':0}
		for bug in seq_dic.keys():
			if seq_dic[bug][i] != '-':
				AA_Nm[seq_dic[bug][i]] += 1
		results.append(AA_Nm)
	return results

def trimMSA(seq_dic = {}, N1=0, N2=-1, seq_list =[] ):
    """ trim the MSA by the cordinates N1 and N2 """
    if N2 > N1:
        if seq_list == []:
            seq_list = seq_dic.keys()
        new_seq_dic = {}
        for seq in seq_list:
            new_seq_dic[seq] = seq_dic[seq][N1-1:N2]
        return new_seq_dic
    else:
        print "trimMSA doing nothing with your data because of the trim coordinates passed " + str(N1) + " " + str(N2) 
        return seq_dic

def makeconsensus(seq_dic={}):
        import random
        """Returns a sequence of the most prevalent amino-acid in each position counting "-" as 21st amino-acid"""
        results = ""
        AAlist = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','Y','W','-']
        for i in range(len(seq_dic.values()[0])):
            AAnumb = numpy.array([0]*len(AAlist))
            for seq in seq_dic.keys():
                AAnumb[AAlist.index(seq_dic[seq][i])] += 1
            M = numpy.argmax(AAnumb)
            if len(AAnumb[ AAnumb == AAnumb[M] ]) > 1:
                print "Position " + str(i) + " is ambiguous: " + " ".join( [ AAlist[i] for i in numpy.where( AAnumb == AAnumb[M] )[0] ] )
                print AAnumb
                print list([ i for i in numpy.where( AAnumb == AAnumb[M] )[0] ])
                results += AAlist[random.choice(list([ i for i in numpy.where( AAnumb == AAnumb[M] )[0] ]))]
            else:
                results += AAlist[M]
        return results

def aa_distr(seq_dic={}, pos=0):
	"""Calculates the distribution of AA in a given column of a given MSA. The result is a dictionary with the distribution"""
	AA_Nm = {'A':0,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0,'K':0,'L':0,'M':0,'N':0,'P':0,'Q':0,'R':0,'S':0,'T':0,'V':0,'Y':0,'W':0,'-':0}
	if pos > len(seq_dic[seq_dic.keys()[0]]):
		AA_Nm = {'Error':'Position out of range: ' + str(pos)}
		return AA_Nm
	for bug in seq_dic.keys():
		try:
			AA_Nm[seq_dic[bug][pos]] += 1
		except IndexError:
			print 'Error: Sequences are not aligned'
			sys.exit()
	return AA_Nm


def Xsq2sample( cont_table = [[]]):
	from scipy.stats import chi2
	"""Return the chi-square value and the p-value for a given contigency table)"""
	# calculates the degrees of freedom:
	DF = (len(cont_table) - 2) * (len(cont_table[0]) - 2)
	
	# calculates the Expected values:
	E = []
	for i in range(len(cont_table)-1):
		Ei = []
		for j in range(len(cont_table[i])-1):
			#print str(i) + '\t' + str(j) + '\n' + str(cont_table[i][:-1])
			Ei.append((cont_table[i][-1] * cont_table[-1][j]) / float(cont_table[-1][-1]))
		E.append(Ei)
	
	# calculate the Chi-square:
	X = 0
	for i in range(len(cont_table)-1):
                for j in range(len(cont_table[i])-1):
			try:
				X += ((cont_table[i][j] - E[i][j])**2)/float(E[i][j])
			except ZeroDivisionError:
				X = 'Err'
				break
		if X == 'Err':
			break


			

	#print 'DF = ' + str(DF)
	if X != 'Err':
		if X != 0:
			P = 1 - chi2.cdf(X,DF)
		else:
			P = 1.0
	else:
		P = 'Err'

	return X, P		
	
def identity_rank(seq_dic = {}):
	""" Returns the conservation level (identity) of each position)"""
	results = []
	for i in range(len(seq_dic.values()[0])):
		Total = 0
		AA_Nm = {'A':0,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0,'K':0,'L':0,'M':0,'N':0,'P':0,'Q':0,'R':0,'S':0,'T':0,'V':0,'Y':0,'W':0, 'X':0, '-':0}
		for bug in seq_dic.keys():
			AA_Nm[seq_dic[bug][i]] += 1
			Total += 1
		if AA_Nm['-'] == max(AA_Nm.values()):
			results.append(0)
		else:
			results.append(max(AA_Nm.values())/float(Total))
	return results

def consensus(seq_dic = {}, sec = [100,90,80,70,60]):
	"""It will somehow return the consensus of each position of an alignment"""
	#check for lenght:
	seq_len = 0
	for seq in seq_dic.values():
		if seq_len == 0:
			seq_len = len(seq)
		elif len(seq) != seq_len:
			print "Consider revision to the alignment: mismatched lenght"
			sys.exit()				

	#groups:
	results = {}

	for cut in sec:
		cons = ''
		cut_num = int(cut)/float(100)
		for i in range(len(seq_dic.values()[0])):
			groups_dic = {  'A':[0,['A']],
					'C':[0,['C']],
					'D':[0,['D']],
					'E':[0,['E']],
					'F':[0,['F']],
					'G':[0,['G']],
					'H':[0,['H']],
					'I':[0,['I']],
					'K':[0,['K']],
					'L':[0,['L']],
					'M':[0,['M']],
					'N':[0,['N']],
					'P':[0,['P']],
					'Q':[0,['Q']],
					'R':[0,['R']],
					'S':[0,['S']],
					'T':[0,['T']],
					'V':[0,['V']],
					'Y':[0,['Y']],
					'W':[0,['W']],
					'a':[0,['F','Y','W','H']],
					'l':[0,['I','V','L']],
					'h':[0,['F','Y','W','H','I','V','L','A','G','M','C','K','R','T']],
					'+':[0,['H','K','R']],
					'-':[0,['D','E']],
					'c':[0,['H','K','R','D','E']],
					'p':[0,['H','K','R','D','E','Q','N','S','T','C']],
					'o':[0,['S','T']],
					'u':[0,['G','A','S']],
					's':[0,['G','A','S','V','T','D','N','P','C']],
					't':[0,['G','A','S','H','K','R','D','E','Q','N','S','T','C']],
					'.':[0,['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','Y','W']] }
	
			for bug in seq_dic.keys():
				for group in groups_dic.keys():
					if seq_dic[bug][i] in groups_dic[group][1]:
						groups_dic[group][0] += 1
			for group in groups_dic.keys():
				if  groups_dic[group][0] != 0:
					groups_dic[group][0] /= float(len(seq_dic.keys()))
			winner_size = 30
			winner = ' '
                        winner_score = 0
    			for group in groups_dic.keys():
				if groups_dic[group][0] >= cut_num and len(groups_dic[group][1]) < winner_size: 
					winner = group
                                        winner_score = groups_dic[group][0]
					winner_size = len(groups_dic[group][1])
                                elif groups_dic[group][0] >= cut_num and groups_dic[group][0] > winner_score and len(groups_dic[group][1]) == winner_size:
                                        winner = group
                                        winner_size = len(groups_dic[group][1])
                                        winner_score = groups_dic[group][0]
                        cons += winner
		results[cut] = cons
	return results


def clustal_colors(seq_dic = {}, seq_list = []):
    """ Calculates what is needed and place it into a ... ... .. to output a color code like clustalX. Mainly intended to work with map_feature_ontree.v3"""
    #constants

    #two types of rules: 0 -> any amino-acids (sum), 1 -> that conservation threshold to any individual amino-acid 

    rules = [ [ "AILMFWV", [ [ 0.6 , "WLVIMAFCHP", 0 ]                                                      ], "#80a0f0" ], 
              [ "RK"     , [ [ 0.6 , "KR"        , 0 ], [ 0.80 , "KRQ"         , 1 ]                        ], "#f01505" ],
              [ "N"      , [ [ 0.5 , "N"         , 0 ], [ 0.85 , "NY"          , 1 ]                        ], "#80a0f0" ],
              [ "C"      , [ [ 0.6 , "WLVIMAFCHP", 0 ],                                                     ], "#80a0f0" ],
              [ "C"      , [ [ 1.0 , "C"         , 0 ]                                                      ], "#f08080" ],
              [ "Q"      , [ [ 0.6 , "KR"        , 0 ], [ 0.5  , "QE"          , 0 ], [ 0.85 , "QEKR" , 1 ] ], "#15c015" ],
              [ "E"      , [ [ 0.6 , "KR"        , 0 ], [ 0.5  , "QE"          , 0 ], [ 0.85 , "QED"  , 1 ] ], "#c048c0" ],
              [ "D"      , [ [ 0.6 , "KR"        , 0 ], [ 0.5  , "ED"          , 0 ], [ 0.85 , "KRQ"  , 1 ] ], "#c048c0" ],
              [ "G"      , [ [ 0.0 , "G"         , 0 ]                                                      ], "#f09048" ],
              [ "HY"     , [ [ 0.6 , "WLVIMAFCHP", 0 ], [ 0.85 , "WYACPQFHILMV", 1 ]                        ], "#15a4a4" ],
              [ "P"      , [ [ 0.0 , "P"         , 0 ]                                                      ], "#c0c000" ],
              [ "ST"     , [ [ 0.6 , "WLVIMAFCHP", 0 ], [ 0.5 , "TS"           , 0 ], [ 0.85, "TS"    , 1 ] ], "#15c015" ] ]

    rules = [ [ "A"      , [ [ 0.6  , "WLVIMAFCYHP", 0 ], [ 0.5  , "P"  , 0 ], [ 0.85 , "TSG" , 1 ] ,           ], "#80a0f0" ],
              [ "ILMFWV" , [ [ 0.6  , "WLVIMAFCYHP", 0 ], [ 0.5  , "P"  , 0 ]                                   ], "#80a0f0" ], 
              [ "RK"     , [ [ 0.6  , "KR"         , 0 ], [ 0.85 , "Q"  , 0 ]                                   ], "#f01505" ],
              [ "N"      , [ [ 0.5  , "N"          , 0 ], [ 0.85 , "D"  , 0 ]                                   ], "#15c015" ],
              [ "C"      , [ [ 0.6  , "WLVIMAFCYHP", 0 ], [ 0.5  , "P"  , 0 ], [ 0.85 , "S", 0 ]                ], "#80a0f0" ],
              [ "C"      , [ [ 0.85 , "C"          , 0 ]                                                        ], "#f08080" ],
              [ "Q"      , [ [ 0.6  , "KR"         , 0 ], [ 0.5  , "QE" , 0 ]                                   ], "#15c015" ],
              [ "E"      , [ [ 0.5  , "QE"         , 0 ], [ 0.5  , "ED" , 0 ]                                   ], "#c048c0" ],
              [ "D"      , [ [ 0.5  , "ED"         , 0 ], [ 0.5  , "N"  , 0 ]                                   ], "#c048c0" ],
              [ "G"      , [ [ 0.0  , "G"          , 0 ]                                                        ], "#f09048" ],
              [ "HY"     , [ [ 0.6  , "WLVIMAFCYHP", 0 ], [ 0.5  , "P"  , 0 ]                                   ], "#15a4a4" ],
              [ "P"      , [ [ 0.0  , "P"          , 0 ]                                                        ], "#c0c000" ],
              [ "S"      , [ [ 0.6  , "WLVIMAFCYHP", 0 ], [ 0.5 , "TS"  , 0 ]                                   ], "#15c015" ],
              [ "T"      , [ [ 0.8  , "WLVIMAFCYHP", 0 ], [ 0.5 , "TS"  , 0 ]                                   ], "#15c015" ] ]




    res_dic = {}
    
    hist = aa_hist(seq_dic)

    num_seq = len(seq_list)

    print json.dumps(hist, indent = 2 )

    MSA_len = len(hist)

    for i in range(MSA_len):
        for tag in seq_list:
            if tag not in res_dic.keys():
                res_dic[tag] = []
            AA = seq_dic[tag][i]
#            color = "#FFFFFF"
#            get_color = False
            print tag + " " + str(i) + " " + AA
            for rule in rules:
                #print rule[0]
                if AA in rule[0]:
                    for dec in rule[1]:
                        if dec[2] == 0:
                            p = 0
                            #print hist[i]
                            #print dec[1]
                            for aa in dec[1]:
                                p += float(hist[i][aa])/num_seq
                            #print p
                            if p > dec[0]:
                                print "Rule type 0 color: " + rule[2]
                                color = rule[2]
                                break
                            else:
                                color = "#FFFFFF"
                        if dec[2] == 1:
                            for aa in dec[1]:
                                if float(hist[i][aa])/num_seq > dec[0]:
                                    color = rule[2]
                                    break
                                else:
                                    color = "#FFFFFF"
                            if color != "#FFFFFF":
                                break
                    if color != "#FFFFFF":
                        print "Found rule" + rule[0]
                        break
            res_dic[tag].append(color)
    return res_dic

def sortrow(counts=[],labels=[]):
	""" To be used to sort rows by the maximum valuer in matrixes such as the ones that correlates MCP class with Che class"""
	#print len(counts)
	#print len(labels)
	for i in range(len(counts[0])-1):
		max_coord = numpy.where(counts[i:,i:]==counts[i:,i:].max())
		max_coord[0][0] = max_coord[0][0] + i
		max_coord[1][0] = max_coord[1][0] + i
		#swapping rows
		#print "Round " + str(i)
		#print "Searching Max in:"
		#print counts[i:,i:]
		#print "Max val: " + str(counts[i:,i:].max()) + " at (" + str(max_coord[0][0]) + ',' + str(max_coord[1][0]) + ')'
		#print "Swapping rows: " + str(i) + " with " + str(max_coord[0][0])
		#print "Before"
		#print counts
		temp = counts[i].copy()
		counts[i] = counts[max_coord[0][0]]
		counts[max_coord[0][0]] = temp
		#print str(len(labels)) + '\t' + str(max_coord[0][0])
		temp = labels[i]
		labels[i] = labels[max_coord[0][0]]
		labels[max_coord[0][0]] = temp
		#print "After"
		#print counts
	return counts, labels

def sortcol(counts=[], labels=[]):
	""" Same as the above but for columns """
	#print "row"
	#print len(counts[0])
	#print len(labels)
	for i in range(len(counts[0])-1):
		max_coord = numpy.where(counts[i:,i:]==counts[i:,i:].max())
		max_coord[0][0] = max_coord[0][0] + i
		max_coord[1][0] = max_coord[1][0] + i
		#print "Round " + str(i)
		#print "Searching Max in:"
		#print counts[i:,i:]
		#print "Max val: " + str(counts[i:,i:].max()) + " at (" + str(max_coord[0][0]) + ',' + str(max_coord[1][0]) + ')'
		#print "Swapping colums: " + str(i) + " with " + str(max_coord[1][0])
		#print "Before"
		#print counts
		temp = counts[:,i].copy()
		counts[:,i] = counts[:,max_coord[1][0]]
		counts[:,max_coord[1][0]] = temp
		#print str(len(labels)) + '\t' + str(max_coord[1][0])
		temp = labels[i]
		labels[i] = labels[max_coord[1][0]]
		labels[max_coord[1][0]] = temp
		#print "After"
		#print counts
	#print "fim"
	return counts, labels

def match_pairs(seq1_dic = {}, seq2_dic = {}):
	""" builds the dictionary matching the interacting sequences used in ESCL """

	seq_seq = {}
	
	for seq1_tag in seq1_dic.keys():
		orgn_ID1 = seq1_tag.split('-')[0]
		seq_seq[seq1_tag] = []
		for seq2_tag in seq2_dic.keys():
			orgn_ID2 = seq2_tag.split('-')[0]
			if orgn_ID1 == orgn_ID2:
				seq_seq[seq1_tag].append(seq2_tag)

	return seq_seq
			
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
        return Dmatrix, names

def alnpos_dic(seq1, seq2, bias1 = 0, bias2 = 0):
	"""Must input two aligned sequences. It return a dictionary with coordinates of position in seq2 aligned to certain position in seq1. It also accepts bias for seq1 and seq2"""
	if len(seq1) !=  len(seq2):
		print("Sequences must have the same length (aligned)")
		return -1
	pos1 = 1 + bias1
	pos2 = 1 + bias2
	dic = {}
	lis = []
	for i in range(len(seq1)):
		if seq1[i] != '-' and seq2[i] != '-':
			dic[pos1] = pos2
			lis.append(pos1)
			pos1 += 1
                        pos2 += 1

		elif seq1[i] != '-' and seq2[i] == '-':
			pos1 += 1
		elif seq1[i] == '-' and seq2[i] != '-':
			pos2 += 1
	return dic, lis

def threeLetter2oneLetter(seq):
	dic = { 'ALA':'A',
                'ARG':'R',
                'ASN':'N',
                'ASP':'D',
                'ASX':'B',
                'CYS':'C',
                'GLU':'E',
                'GLN':'Q',
                'GLX':'Z',
                'GLY':'G',
                'HIS':'H',
		'HSD':'H',
                'ILE':'I',
                'LEU':'L',
                'LYS':'K',
                'MET':'M',
                'PHE':'F',
                'PRO':'P',
                'SER':'S',
                'THR':'T',
                'TRP':'W',
                'TYR':'Y',
                'VAL':'V'}
        if seq.__class__ != list:
		print "Incorrect input!!!\nInput must be three letter code list"
		return -1
	else:
		return ''.join([ dic[i] for i in seq])

def getmd5(seq):
	return base64.encodestring(md5.new(seq.replace('-','')).digest()).replace('/','_').replace('=','').replace('+','-').replace('\n','')

def tree_to_phyloxml (ete_tree, che_class_list = []):
	colors_array = [ '#E41A1C' , '#377EB8' , '#4DAF4A' , '#984EA3' , '#FF7F00' , '#999999' , '#A65628' , '#F781BF', '#17BECF']
	"""
	Convert an Ete2 tree to PhyloXML.
 
	:Parameters:
	ete_tree
	An Ete2 format tree
 
	:Returns:
	PhyloXML markup as text
 
	"""
	from cStringIO import StringIO
	buffer = StringIO()
 
	def visit_node (node, buf, indent=0):
		buf.write (" " * indent)
		buf.write ("<phy:clade>\n")
		buf.write (" " * (indent+1))
		buf.write ("<phy:name>%s</phy:name>\n" % node.name )
		buf.write (" " * (indent+1))
		buf.write ("<phy:branch_length>%s</phy:branch_length>\n" % node.dist)
		buf.write (" " * (indent+1))
		buf.write ("<phy:confidence type='branch_support'>%s</phy:confidence>\n" % node.support)
		if 'uri' in list(node.features) or 'desc' in list(node.features):
			buf.write (" " * (indent+1))
			buf.write ("<phy:annotation>\n")
			if 'uri' in list(node.features):
				buf.write (" " * (indent+2))
				buf.write ("<phy:uri>%s</phy:uri>\n" % node.uri)
			if 'desc' in list(node.features):
                                buf.write (" " * (indent+2))
                                buf.write ("<phy:desc>%s</phy:desc>\n" % node.desc)
			buf.write (" " * (indent+1))
			buf.write ("</phy:annotation>\n")
		if che_class_list != [] and 'checlass' in list(node.features):
			buf.write (" " * (indent+1))
                        buf.write ("<phy:chart>\n")
			for che in che_class_list:
				if che in node.checlass:
					buf.write (" " * (indent+2))
					buf.write ("<phy:component%s>%s</phy:component%s>\n" % (che_class_list.index(che), che, che_class_list.index(che)))
				else:
					buf.write (" " * (indent+2))
                                        buf.write ("<phy:component%s>none</phy:component%s>\n" % (che_class_list.index(che),che_class_list.index(che)))
			buf.write (" " * (indent+1))
                        buf.write ("</phy:chart>\n")
 
		for c in node.get_children():
			visit_node (c, buf, indent=indent+1)
 
		buf.write (" " * indent)
		buf.write ("</phy:clade>\n")
 
	buffer.write ("<phy:Phyloxml xmlns:phy='http://www.phyloxml.org/1.10/phyloxml.xsd'>\n")
	buffer.write ("<phy:phylogeny>\n")
	buffer.write ("<phy:name>test_tree</phy:name>\n")
	buffer.write ("<phy:render>\n")
	buffer.write (" <phy:parameters>\n")
	buffer.write ("  <phy:circular>\n")
	buffer.write ("   <phy:bufferRadius>0.5</phy:bufferRadius>\n")
	buffer.write ("  </phy:circular>\n")
	buffer.write ("	 <phy:rectangular>\n")
	buffer.write ("   <phy:alignRight>1</phy:alignRight>\n")
	buffer.write ("  </phy:rectangular>\n")
	buffer.write (" </phy:parameters>\n")
	buffer.write (" <phy:charts>\n")
	for che in che_class_list:
		buffer.write ('	 <phy:component%s type="binary" thickness="20" />\n ' % che_class_list.index(che))
#	buffer.write ('	 <phy:acidity type="binary" thickness="10" disjointed="1" bufferSiblings="0.3" />\n')
	buffer.write (" </phy:charts>\n")
	buffer.write (" <phy:styles>\n")
	buffer.write ("  <phy:caffeine fill='#A93' stroke='#DDD' />\n")
	buffer.write ("	 <phy:base fill='#8b7100' stroke='#DDD' />\n")
	buffer.write ("	 <phy:other fill='#333' stroke='#DDD' />\n")
	buffer.write ("	 <phy:high fill='#666' stroke='#DDD' />\n")
	buffer.write ("  <phy:mid fill='#999' stroke='#DDD' />\n")
	buffer.write ("  <phy:low fill='#CCC' stroke='#DDD' />\n")
	buffer.write ("  <phy:none fill='#FFF' stroke='#CCC' />\n")
	for che in che_class_list:
		buffer.write ("  <phy:%s fill='%s' stroke='#DDD' />\n" % (che, colors_array[che_class_list.index(che)]))
#        buffer.write ("  <phy:ACF fill='#377EB8' stroke='#DDD' />\n")
	buffer.write (" </phy:styles>\n")
	buffer.write ("</phy:render>\n")
 
	visit_node (ete_tree.get_tree_root(), buffer)
 
 	buffer.write ("</phy:phylogeny>\n")
	buffer.write ("</phy:Phyloxml>\n")
 
	return buffer.getvalue()



