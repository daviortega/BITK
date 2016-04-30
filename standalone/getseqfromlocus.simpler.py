#! /usr/bin/env python 
###################################
#    Davi Ortega 10/17/2014 
###################################
import pymongo
import sys
if '-h' in sys.argv:
	print 'Input list of locus. Output fasta with sequences'
	sys.exit()

lc_list = []
feature = {}

def get_seqdepot_client():
	print "Verifying tunnels"
	print "SeqDepot"
	try:                                    
		#client = pymongo.MongoClient('localhost',27018)
		client = pymongo.MongoClient("aphrodite.bio.utk.edu",27017)
		client.the_database.authenticate('binf',open("/Users/ortegad/private/mongodb_binf.txt","r").readline().strip(), mechanism='MONGODB-CR',source='admin')
	except:
		print "Edit get_seqdepot_client to get a SeqDepot client"
		sys.exit()
	print "Authenticated"	
	return client.seqdepot

def get_mist22_client():
	print "Verifying tunnels"
	print "Mist"
	try:
		client = pymongo.MongoClient('localhost',27019)
	except:
		print "You must open a tunnel with ares.bio.utk.edu: ssh -p 32790 -f -N -L 27019:localhost:27017 ortega@ares.bio.utk.edu"
		sys.exit()
	return client.mist22

def lo2seq_fasta(lc_list = []):
	print "Locus -> Accession"
	print lc_list
	aseq_list = []
	lo2aseq = {}
	mist22 = get_mist22_client()
	print "Got the client\nGetting the aseqs from locus. This might take a bit because MiST is not indexed by locus"
	genes = mist22.genes.find({'lo' : { '$in' : lc_list }}, {'p.aid':1, 'lo' : 1})
	for card in genes:
	    print card
	    try:
			aseq_list.append(card['p']['aid'])
			lo2aseq[card['lo']] = card['p']['aid'] 
	    except:
	        print "Found exception"
	        print card
	        sys.exit()
	print "Got them\nGetting the sequences from SeqDepot now"
	seqs = sd.aseqs.find({'_id': { '$in' : aseq_list }})
	aseq2seq = {}
	for seq in seqs:
		aseq2seq[seq['_id']] = seq['s']
	print "Got them!!!\Let's build the file"
	output = ""
	for lo in lc_list:
		output += ">" + lo + "\n" + aseq2seq[lo2aseq[lo]] + '\n'
	return output

#main
if __name__ == '__main__' :
	print "Test connections to DB before we start"
	mist22 = get_mist22_client()
	sd = get_seqdepot_client()

	with open(sys.argv[1], 'r') as f:
		for l in f:
			items = l.replace('\n','').split('-')
			lc_list.append(items[0])
			if len(items) > 1:
				feature[items[0]] = '-'.join(items[1:])

	output = lo2seq_fasta(lc_list)
	
	print "Quality control"
	if output.count(">") != len(lc_list):
		print "something went wrong"
	else:
		print " All good "
		print " Your sequences will be stored at " + sys.argv[1][:-3] + '.fa'
	
	with open(sys.argv[1][:-3] + '.fa', 'w') as f:
		f.write(output)
