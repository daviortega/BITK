#! /usr/bin/env python 
###################################
#    Davi Ortega 3/29/2013 
###################################
import sys
import bitk
import copy
if '-h' in sys.argv:
	print 'Update tags in fasta files from the organism id provided in list\n \
	Sintaxe: updatetagwithcount.py org_id.list file.fa file2.fa ...'
	sys.exit()


orgid_list = []
orgid_dic = {}

datafile = open(sys.argv[1],'r')

start_list = [0]*(len(sys.argv[2:]))

for line in datafile:
	orgid = line.replace('\n','')
	orgid_list.append(orgid)
	orgid_dic[orgid] = copy.deepcopy(start_list)

datafile.close()

for f in range(2, len(sys.argv)):
	print sys.argv[f]
	datafile = open(sys.argv[f],'r')
	for line in datafile:
		if '>' in line:
			orgid = line.split('-')[0].split('.')[-1]
			if orgid in orgid_list:
				orgid_dic[orgid][f-2] += 1
	datafile.close()

print orgid_dic

for f in range(2, len(sys.argv)):
        print 'Working on file ' + sys.argv[f]
	seq_dic, seq_list = bitk.fastareader(sys.argv[f])
	output = ''
	for tag in seq_list:
		orgid = tag.split('-')[0].split('.')[-1]
		if orgid in orgid_list:
			output += '>' + tag + '-' + '-'.join([str(i) for i in orgid_dic[orgid]]) + '\n' + seq_dic[tag] + '\n'
	datafile = open(sys.argv[f][:-3] + '.countche.fa', 'w')
	datafile.write(output)
	datafile.close()













