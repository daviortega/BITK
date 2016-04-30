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

def lc2ac(lc_list = []):
	print "Locus -> Accession"
	print lc_list
	ac_list = []
	mist22 = get_mist22_client()
	print "Got the client\nGetting the cards now"
	genes = mist22.genes.find({'lo' : { '$in' : lc_list }}, {'p.ac':1})
	for card in genes:
	    print card
	    try:
	        ac_list.append(card['p']['ac'])
	    except:
	        print "Found exception"
	        print card
	        sys.exit()		  
	print "Got them"
	return ac_list

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

def accession2ag(accession_list = [],):
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
	errors = []
	for accession in accession_list:
		try:
			out_dic[accession] = gid_dic[int(preout_dic[accession].split('-')[0])] + preout_dic[accession]
		except KeyError:
			print "This sequence has - in the name and is messing up the data collection"
			print accession
			errors.append(accession)
			pass
	return out_dic, errors


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

	ac_list = lc2ac(lc_list)

	ac2tag, error = accession2ag(ac_list)
	ac2seq = accession2seq(ac_list)

	output = ''
	for ac in ac_list:
		if feature == {}:  
			output += '>' + ac2tag[ac] + '\n' + ac2seq[ac] + '\n'
		else:
			try:
				output += '>' + ac2tag[ac] + '-' + feature[ac] + '\n' + ac2seq[ac] + '\n'
			except KeyError:
				output += '>' + ac2tag[ac] + '\n' + ac2seq[ac] + '\n'

	with open(sys.argv[1][:-3] + '.fa', 'w') as f:
		f.write(output)

	print "Error"
	print error
