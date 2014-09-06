import sys
import time
import bitk
import os
import pymongo
import json
import multiprocessing as multip
import copy
import ast

#########################################################
cdds_path = '/home/ortegad/MIST3/TOL/cdds/'
genomes_path = '/home/ortegad/MIST3/TOL/genomes/'
TOLdb_path = '/home/ortegad/MIST3/TOL/TOLdb/'
TOLdb_name = 'CogTOL2'
TOL_cog_list = [ "gnl|CDD|223091",
        "gnl|CDD|223095",
        "gnl|CDD|223126",
        "gnl|CDD|223127",
        "gnl|CDD|223130",
        "gnl|CDD|223158",
        "gnl|CDD|223159",
        "gnl|CDD|223165",
        "gnl|CDD|223169",
        "gnl|CDD|223170",
        "gnl|CDD|223171",
        "gnl|CDD|223172",
        "gnl|CDD|223174",
        "gnl|CDD|223175",
        "gnl|CDD|223176",
        "gnl|CDD|223178",
        "gnl|CDD|223180",
        "gnl|CDD|223181",
        "gnl|CDD|223250",
        "gnl|CDD|223262",
        "gnl|CDD|223264",
        "gnl|CDD|223275",
        "gnl|CDD|223278",
        "gnl|CDD|223279",
        "gnl|CDD|223280",
        "gnl|CDD|223334",
        "gnl|CDD|223569",
        "gnl|CDD|223596",
        "gnl|CDD|223607"
	]

print "Verifying tunnels"
print "Mist"
try:
        client = pymongo.MongoClient('localhost',27019)
except:
        print "You must open a tunnel with ares.bio.utk.edu: ssh -L 27019:localhost:27017 ares.bio.utk.edu"
	sys.exit()

mist = client.mist22

print "SeqDepot"
try:
	client2 = pymongo.MongoClient('localhost', 27018)
except:
        print "You must open a tunnel with aphrodite.bio.utk.edu: ssh -L 27018:localhost:27017 aphrodite.bio.utk.edu"
	sys.exit()
seqdepot = client2.seqdepot

local = pymongo.MongoClient()

def get_mist22_client():
	print "Verifying tunnels"
	print "Mist"
	try:
		client = pymongo.MongoClient('localhost',27019)
	except:
	        print "You must open a tunnel with ares.bio.utk.edu: ssh -L 27019:localhost:27017 ares.bio.utk.edu"
	        sys.exit()
	return client.mist22

def job_rps(commands):
	command, mid = commands.split(':')
	os.system(command)
	print "Working on genome " + mid
	seq_dic = {'gid' : mid }
	cog_list_tmp = []
	rpsout = open('rps.mist22.' + mid + '.dat','r')
	for line in rpsout:
		fields = line.split('\t')
		if fields[1] not in cog_list_tmp:
			cog_list_tmp.append(fields[1])
			seq_dic[fields[1]] = [ fields[0], float(fields[10]) ]
		elif seq_dic[fields[1]][1] > float(fields[10]):
#                       print "older: " + str(seq_dic[fields[1]][1])
#                       print "new: " + str(float(fields[10]))
			seq_dic[fields[1]] = [ fields[0], float(fields[10]) ]
	rpsout.close()
	cog_list_tmp.sort()
	seq_dic['cog_list'] = cog_list_tmp
	
	#       Saving the results
	with open('rps.mist22.' + mid + '.json', 'w') as f:
		json.dump(seq_dic, f, sort_keys=True, indent = 4)



def run_rps(mid_list = [], dbname = TOLdb_name,TOLdb_path = TOLdb_path, autofix = False, NPrps = 20, cog_list = TOL_cog_list):
	mongo_card = []
	problems = []
	dbname = TOLdb_path + dbname
	#preparing list of commands
	rpscommands = []
	for mid in mid_list:
		rpscommands.append('rpsblast -query ' + genomes_path + 'mist22.' + mid + '.fa -db ' + dbname + ' -evalue 0.1 -outfmt 6 -out rps.mist22.' + mid + '.dat :' + mid )

	if NPrps <= 1:
		for c in rpscommands:
			job_rps(c)
	else:
	        if len(mid_list) < NPrps:
        	        NPrps = len(mid_list)
	        pool = multip.Pool(processes=NPrps)
        	pool.map(job_rps, rpscommands)

	problem_gid_list = []
        for mid in mid_list:
                fin = open('rps.mist22.' + mid + '.json', 'r')
                seq_dic = json.load(fin)
                fin.close()
                if len(set(cog_list)) != len(set(seq_dic['cog_list'])):
                        problem = mid + ": The genome " + mid + " does not have the following hits: " + ' '.join(list(set(cog_list) - set(seq_dic['cog_list'])))
			problem_gid_list.append(mid)
                else:
                        problem = None

                if not problem:
                        print "Genome " + mid + " is good"
                        mongo_card.append(seq_dic)
                        cog_list = seq_dic['cog_list']
                else:
                        print problem
                        problems.append(problem)

        if len(problems) != 0:
                print "List of problems: "
                print '\n'.join(problems)
                if autofix != True:
			print "There are " + str(len(problems)) + " problems with the RPS section. Check your files and if you have RPS-BLAST installed."
			for p in problems:
				print p
                        sys.exit()
                else:
                        print "AUTO-FIX: Let's ignore those..."
			for p in problems:
				print p
                        print "Saving the prelim data at data.json"
                        with open('data.json','w') as f:
                                json.dump(mongo_card, f, sort_keys=True, indent = 4)
			for p in problem_gid_list:
				mid_list.remove(p)
        else:
                print "Saving the prelim data at data.json"
                with open('data.json','w') as f:
                        json.dump(mongo_card, f, sort_keys=True, indent = 4)

	return mongo_card, mid_list


def rpsdata2mongocards(mongo_cards):
        pacs = []
        mongo_card = mongo_cards
	for card in mongo_card:
                for key, value in card.iteritems():
                        if key != 'gid' and key != 'cog_list':
                                pacs.append(value[0].split('-')[-1])
        print "Retriving information from MIST"
        all_ids = mist.genes.find({"p.ac": { "$in" : pacs }}, {"p.aid":1, "p.ac":1})
        all_ids_proc = {}
        all_p_aid = []
        print "Parsing information from MIST"
#       all_ids_json = []
        for ids in all_ids:
#               all_ids_json.append(ids)
                all_ids_proc[ids['p']['ac']] = ids['p']['aid']
                all_p_aid.append(ids['p']['aid'])
        print "Retriving information from Seqdepot"
        all_seq = seqdepot.aseqs.find( {'_id' : { '$in' : all_p_aid }}, {'s':1})
        print "Parsing information from SeqDepot"
        all_seq_proc = {}
        for seq in all_seq:
                all_seq_proc[seq['_id']] = seq['s']

        if '--debug' in sys.argv:
                with open('all_ids_proc.json', 'w') as f:
                        json.dump(all_ids_proc, f, sort_keys=True, indent = 4)
                with open('all_ids_json.json', 'w') as f:
                        json.dump(all_ids_json, f, sort_keys=True, indent = 4)



        print "Processing the data from MIST and SeqDepot"

        new_mongo = []
        for card in mongo_card:
                seq_dic = {'gid' : card['gid']}
#               print card['gid']
                for key, value in card.iteritems():
                        if key != 'gid' and key != 'cog_list':
                                pac = value[0].split('-')[-1]
                                try:
                                        _id = all_ids_proc[pac]
                                except KeyError:
                                        print "Some accesion number have not been found. Maybe the sequence tag was not correct:"
                                        print pac
                                        print value
                                        sys.exit()
                                seq = all_seq_proc[_id]
                                seq_dic[key] = {'s' : seq}
#                               p_aid = mist.genes.find_one({"p.ac":pac},{"p.aid":1})
#                               print p_aid
#                               seq = seqdepot.aseqs.find_one({"_id": p_aid["p"]["aid"]}, {"s":1})
#                               seq_dic[key] = { 't' : value[0], '_id' : p_aid['p']["aid"], 's' : seq['s']}
                new_mongo.append(seq_dic)
        print " Saving...."
        if '--build-new-db' in sys.argv:
                with open('tol.mongo.tmp.json','w') as f:
                        json.dump(new_mongo, f, sort_keys=True, indent = 4)
        else:
                with open('tol.mongo.json','w') as f:
                        json.dump(new_mongo, f, sort_keys=True, indent = 4)

	return new_mongo

def isinmongoTOL(mid, colname):
	exec( 'coll = local.tol.' + colname)
	if coll.find_one({'gid':mid}):
		return True
	else:
		return False

def update_mongoTOL(mongocards, colname = 'toltest', force = False):
	exec('coll = local.tol.' + colname)
	for card in mongocards:
		_id = card['gid']
		if isinmongoTOL(_id, colname) == False:
			print "Updating database with genome " + str(_id)
			coll.insert(card)
		elif force == True:
			print "Forced to update info on genome " + str(_id)
			coll.remove({'gid': str(_id)})
			coll.insert(card)
		else:
			print "No need to update, genome " + str(_id) + " already in the database"
	return 0

def new_mongoTOL(mongocards, colname = 'toltest'):
	exec('coll = local.tol.' + colname)
	coll.drop()
	exec('coll = local.tol.' + colname)
	coll.insert(mongocards)
	return 0

def check_mongoTOL(colname):
	exec('coll = local.tol.' + colname)
	for card in coll.find({},{'gid':1}):
		print card

def load_from_mongoTOL(mid_list, colname, update = False):
	exec('coll = local.tol.' + colname)
	data = []
	errors = []
	for mid in mid_list:
		if coll.find_one({'gid' : str(mid)}):
			data.append(coll.find_one({'gid' : str(mid)}))
		else:
			print "Genome " + mid + " not in the database (or at least not in this collection). Update the database accordingly."
			errors.append(mid)
#			return [{'Error': 'Genome not found in collection'}]
	return data, errors

def ErrorTest(datajson):
	for card in datajson:
		print card
		if 'Error' in card.keys():
			print card
			return True
	return False

def pre_aligncog(mid_list, new_mongo, ignore_problems = False):
        # building alignments fasta files
	error = False
	aln_dic = {}
	print "Loading and processing cards..."
        for card in new_mongo:
#		print card['gid']
                if card['gid'] in mid_list:
                        for cog, info in card.iteritems():
                                if cog not in ['gid', '_id']:
                                        if cog not in aln_dic.keys():
                                                aln_dic[cog] = '>' + card['gid'] + '\n' + info['s'] + '\n'
                                        else:
                                                aln_dic[cog] += '>' + card['gid'] + '\n' + info['s'] + '\n'

        print "Saving all of them and pushing a list of COGs to a new collection on the database"
        cog_list = []
        for cog in aln_dic.keys():
                print "There are " + str(aln_dic[cog].count('>')) + " in cog " + cog
                cog_list.append(cog)
                alnout = open( 'cdd.' + cog.split('|')[-1] + '.fa', 'w')
                alnout.write(aln_dic[cog])
                alnout.close()
		if aln_dic[cog].count('>') != len(mid_list):
			lines = aln_dic[cog].split('\n')
			mids = []
			for line in lines:
				if '>' in line:
					mids.append(line.replace('>',''))
			print "Genome(s) " + ' '.join(list(set(mid_list) - set(mids))) + ' is(are) missing the cog ' + cog
			print "Let's keep going to see if there is any more problems, but I will exit with error at the end"
			error = True
	if error == True and ignore_problems == False:
		print "There is a problem with the dataset. Aborting execution."
		sys.exit()
	cog_list.sort()
	return cog_list

def aligncog(cog_list):
        def submit_to_linsi(fasta = ''):
                os.system(' linsi --quiet --thread 12 ' + fasta + ' > ' + fasta[:-3] + '.linsi.fa')
        for cog in cog_list:
                print "Aligning " + cog
                submit_to_linsi( 'cdd.' + cog.split('|')[-1] + '.fa' )


def concat(cog_list, filename = 'concat.fa'):
        concat = {}
        print "Concatenating the alignments"
        for cog in cog_list:
                seq_dic, tags = bitk.fastareader('cdd.' + cog.split('|')[-1] + '.linsi.fa', 'r')
                for tag in tags:
                        if tag not in concat.keys():
                                concat[tag] = seq_dic[tag]
                        else:
                                concat[tag] += seq_dic[tag]

        print "Saving the full alignment and checking for bugs"

        output = ''
        max_len = 0

        for tag in tags:
                output += '>' + tag + '\n' + concat[tag] + '\n'
                if max_len == 0:
                        max_len = len(concat[tag])
                elif max_len != len(concat[tag]):
                        print "Something is wrong with the genome: " + tag
                        print len(concat[tag])
			fout = open(filename + '.debug.fa', 'w')
		        fout.write(output)
		        fout.close()
                        sys.exit()
		print len(concat[tag])

        fout = open(filename, 'w')
        fout.write(output)
        fout.close()
	return filename

def aligncogbykingdom(cog_list):
        def submit_to_linsi(fasta = ''):
                os.system(' linsi --quiet --thread 12 ' + fasta + ' > ' + fasta[:-3] + '.linsi.fa')
		return fasta[:-3] + '.linsi.fa'

	sets = {}
	kingdic = {}
	kinglist = []
	for cog in cog_list:
		seq_dic, tags = bitk.fastareader('cdd.' + cog.split('|')[-1] + '.fa')
		if kingdic == {}:
			tagsI = [ int(tag) for tag in tags ]
			taxinfo = mist.genomes.find({'_id' :  { '$in' : tagsI }}, { 'ta' : 1, 'n' : 1 })
			for i in taxinfo:
				kingdic[i['_id']] = i['ta'][0]
			kinglist = list(set(kingdic.values()))
			print kinglist
		for kingdom in kinglist:
			print "Aligning " + cog + ' from kingdom ' + kingdom
			if kingdom not in kinglist:
				kinglist.append(kingdom)
			output = ''
			for tag in tags:
				if kingdic[int(tag)] == kingdom:
					output += '>' + str(tag) + '\n' + seq_dic[str(tag)] + '\n'
			with open('cdd.' + cog.split('|')[-1] + '.' + kingdom + '.fa', 'w' ) as f:
				f.write(output)
			submit_to_linsi( 'cdd.' + cog.split('|')[-1] + '.' + kingdom + '.fa')

	return kinglist


def concatbykingdom(cog_list, kinglist = [], filename = 'concat.fa'):
	kingdic = {}
        print "Concatenating the alignments"
	filenames = []
	for kingdom in kinglist:
		concat = {}
	        for cog in cog_list:
        	        seq_dic, tags = bitk.fastareader('cdd.' + cog.split('|')[-1] + '.' + kingdom + '.linsi.fa', 'r')
                	for tag in tags:
                        	if tag not in concat.keys():
                                	concat[tag] = seq_dic[tag]
	                        else:
        	                        concat[tag] += seq_dic[tag]

	        print "Saving the full alignment and checking for bugs"

	        output = ''
        	max_len = 0

	        for tag in tags:
        	        output += '>' + tag + '\n' + concat[tag] + '\n'
                	if max_len == 0:
	                        max_len = len(concat[tag])
        	        elif max_len != len(concat[tag]):
                	        print "Something is wrong with the genome: " + tag
                        	print len(concat[tag])
	                        fout = open(filename + '.debug.fa', 'w')
        	                fout.write(output)
                	        fout.close()
                        	sys.exit()
	                print len(concat[tag])

		fname = filename[:-3] + '.' + kingdom + '.fa'
        	fout = open(fname, 'w')
	        fout.write(output)
        	fout.close()
		filenames.append(fname)
        return filenames

def gbbykingdom(filenames = [], par = '-b3=8 -b4=2 -b5=h') :
	profile = []
	for filename in filenames:
	        os.system('Gblocks ' + filename + ' ' + par)
        	filenamenew = filename.split('.')
	        filenamenew.insert(len(filenamenew)-1, 'gb')
        	filenamenew = '.'.join(filenamenew)
        	os.system('mv ' + filename + '-gb ' + filenamenew)
		profile.append(filenamenew)

	for i in range(len(profile)-1):
		if i == 0:
			os.system('mafft-profile ' + profile[i] + ' ' + profile[i+1] + ' > profile_concat_tmp.fa')
			os.system('mv profile_concat_tmp.fa pctmp.fa')
		else:
			os.system('mafft-profile pctmp.fa' + profile[i+1] + ' > profile_concat_tmp.fa')
			os.system('mv profile_concat_tmp.fa pctmp.fa')

	filenamenew = 'concat.gb.mafftprofile.fa'
	os.system('mv pctmp.fa ' + filenamenew)
        return filenamenew



def gb(filename = 'concat.fa', par = '-b3=8 -b4=2 -b5=h') :
        os.system('Gblocks ' + filename + ' ' + par)
	filenamenew = filename.split('.')
	filenamenew.insert(len(filenamenew)-1, 'gb')
	filenamenew = '.'.join(filenamenew)
        os.system('mv ' + filename + '-gb ' + filenamenew)
	return filenamenew

def preptree(filename = 'concat.gb.fa', nogid = False):
        seq_dic, tags = bitk.fastareader(filename)
        tags = [ int(i) for i in tags]
        names = mist.genomes.find({'_id' :  { '$in' : tags }}, {'n' : 1, '_id':1})
        names_dic = {}
        for name in names:
		if nogid == True:
			names_dic[name['_id']] = name['n']
		else:
	                names_dic[name['_id']] = str(name['_id']) + '|' + name['n']

        output = ''
        for tag in tags:
                output += '>' + names_dic[tag] + '\n' + seq_dic[str(tag)] + '\n'

	filename = filename.split('.')
	filename.insert(len(filename)-1, 'names')
	filename = '.'.join(filename)
        fout = open(filename, 'w')
        fout.write(output)
        fout.close()
	return filename


def maketree(filename = 'concat.gb.names.fa', bootstrap = 0, NPraxml = 20):
        os.system('rm RAxML*')
        os.system('fa2phy ' + filename)
	filename = filename.split('.')
	filenamenew = '.'.join(filename[:-1])
	filename = '.'.join(filename[:-1] + ['phy'])

	if bootstrap != 0:
		os.system('raxmlHPC-PTHREADS-AVX -T ' + str(NPraxml) + ' -m PROTGAMMAJTTF -p 12345 -# ' + str(bootstrap) + ' -s ' + filename + ' -n ' + filenamenew + '.besttree_raxml.nwk')
        else:
                os.system('raxmlHPC-PTHREADS-AVX -T ' + str(NPraxml) + ' -m PROTGAMMAJTTF -p 12345 -s ' + filename + ' -n ' + filenamenew + '.besttree_raxml.nwk')
	return filenamenew + '.besttree_raxml.nwk'


def maketreephyml(filename = 'concat.gb.names.fa', bootstrap = 0, NP = 20):
        os.system('fa2phy ' + filename)
        filename = filename.split('.')
        filenamenew = '.'.join(filename[:-1])
        filename = '.'.join(filename[:-1] + ['phy'])

        if bootstrap == 0:
                os.system('phyml -i ' + filename + ' -q -d aa -m JTT -c 4 -a e -v e -s SPR --no_memory_check')
        else:
                os.system('mpirun -n ' + str(NP) + 'phyml-mpi -i ' + filename + ' -q -d aa -m JTT -c 4 -a e -v e -s SPR --no_memory_check -b ' + bootstrap )
        return filenamenew + '.phy_phyml_tree.txt'



def posttree(filename = 'concat.gb.names.besttree_raxml.nwk', treefilename = 'bitkTOL.tree.nwk'):
#	E = os.system('cp RAxML_result.concat.tree ' + filename')
#       if E != 0:
#       	print "but do not worry, I got this!!!!"
        os.system('cp RAxML_bestTree.' + filename + ' ' + filename)
        os.system('rectaxontree2 ' + filename)
	filename = filename.split('.')
        filename.insert(len(filename)-1, 'rec')
        filename = '.'.join(filename)
        os.system('cp ' + filename + ' ' + treefilename)
        print "Tree is ready under the name: " + treefilename

def posttreephyml(filename = 'concat.gb.names.phy_phyml_tree.txt', treefilename = 'bitkTOL.tree.nwk'):
	os.system('rectaxontree3 ' + filename)
	filename = filename[:-4] + '.rec.nwk'
	os.system('cp ' + filename + ' ' + treefilename)
        print "Tree is ready under the name: " + treefilename


def mktreeengine(cog):
	print "Working with cog " + cog
	t = cog.split('|')[-1]
	print "here"
#	fn = gb('cdd.' + t + '.linsi.fa')
	fn = preptree('cdd.' + t + '.linsi.fa') # fn)
	fn = maketreephyml(fn)
	os.system(' rectaxontree3 ' + fn)



def maketreeofeachcogaln(cog_list):
	aligncog(cog_list)
	N = len(cog_list)
	if N > 20:
		N = 20
	print "making the trees"
	pool = multip.Pool(processes=N)
        pool.map(mktreeengine, cog_list)
	

	

