#!/usr/bin/python
# Davi Ortega - April 2016

import os
import argparse
import json
import sys
import shutil
import pymongo
import multiprocessing
import time
import datetime

#HardVariables

PIPELINE = ['init', 'fetchProtFams', 'fetchGenInfo', 'mkFastaFiles' ]
BITKTAGSEP = '|'
BITKGENSEP = '_'

class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

def FastaReader(datafile):
	
	ListOrder = []
	SeqDic = {}
	fastaBuffer = None
	
	fh = open(datafile, 'r')
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

def BuildDirectorySystem( ProjectName = "PhyloPro", force = False ):
	""" Build the local directory system and config file used by PhyloPro """

	main_path = os.getcwd()

	makedir = [ main_path + '/' + ProjectName, main_path + '/' + ProjectName + '/COGs/' , main_path + '/' + ProjectName + '/merge/', main_path + '/' + ProjectName + '/trees/' ]

	if os.path.isdir(makedir[0]) and force == False:
		print "\n" + color.RED + color.BOLD + "PhyloPro - WARNING ***************************" + color.END + "\n\nThere is a project named " + color.BOLD + ProjectName + color.END + " in this direcptry.\n\nIf you really want to proceed, please use the" + color.RED + color.BOLD + "--init-force" +color.END + " flag instead.\n\n" + color.BOLD + color.RED + "All the information in the directory and subdirectories will be lost.\n\n" + color.END
		sys.exit()
	elif force == True:
		shutil.rmtree(makedir[0])
		pass
		#REMOVE DIRECTORIES

	#Make new dir
	for localDir in makedir:
		os.mkdir(localDir)
	#Make new phylopro config file
	cfgfilename = "phylopro." + ProjectName + ".cfg.json"
	LocalConfigFile = { "ProjectName" : ProjectName ,
						"stage"		  : "",
						"history"     : { 'init' : datetime.datetime.now().isoformat() }
					  }

	with open(makedir[0] + '/' + cfgfilename, 'w') as f:
		json.dump(LocalConfigFile, f, indent = 2)

def update_stage_phylopro_cfg (ProjectName, stage):
	""" Update the PhyloPro config file to state the latest sucessful step """
	main_path = os.getcwd()
	
	with open( main_path + '/' + ProjectName + '/' + "phylopro." + ProjectName + ".cfg.json", 'r') as f:
		LocalConfigFile = json.load(f)

	LocalConfigFile['stage'] = stage
	LocalConfigFile['history'][stage] = datetime.datetime.now().isoformat()
		
	with open( main_path + '/' + ProjectName + '/' + "phylopro." + ProjectName + ".cfg.json", 'w') as f:
		json.dump(LocalConfigFile, f, indent = 2)

def addInfo2phylopro_cfg( ProjectName, Info = {}):
	LocalConfigFile = get_cfg_file(ProjectName)
	LocalConfigFile = merge_two_dicts(LocalConfigFile, Info)
	with open( main_path + '/' + ProjectName + '/' + "phylopro." + ProjectName + ".cfg.json", 'w') as f:
		json.dump(LocalConfigFile, f, indent = 2)
		
def get_cfg_file(ProjectName):
	main_path = os.getcwd()
	
	if ProjectName:
		with open( main_path + '/' + ProjectName + '/' + "phylopro." + ProjectName + ".cfg.json", 'r') as f:
			LocalConfigFile = json.load(f)
		return LocalConfigFile
	else:
		print "Something is wrong: Your config file does not seem to exist in the right place"
		print " Project Name : " + ProjectName
		sys.exit() #If you ever get this error, feel free to throw a better error.

def get_stage(ProjectName):
	LocalConfigFile = get_cfg_file(ProjectName)
	return LocalConfigFile['stage']

def get_ProjectName():
	main_path = os.getcwd()
	listdir = [d for d in os.listdir(main_path) if os.path.isdir(os.path.join(main_path, d))]
	print listdir
	if len(listdir) == 1:
		ProjectName = listdir[0]
	elif len(listdir) == 0:
		print "It seems that there is no project in this directory, please initiate PhyloPro with --init ProjectName flag."
		sys.exit()
	else:
		print "It seems that there are more than one project in this directory, please specify the project name you want in front of the flag."
		sys.exit()
	return ProjectName

def isTheRightOrder(requested_stage, LocalConfigFile):
	current_stage = get_stage(ProjectName)
	if PIPELINE.index(current_stage) + 1 == PIPELINE.index(requested_stage):
		return True
	elif PIPELINE.index(current_stage) + 1 > PIPELINE.index(requested_stage):
		print "You are about to start the pipeline from a previous stage than where it actually is. This is fine, but we will re-calculate and rewrite all the information from here on. Let's think about it for 15 seconds, shall we ?!"
		sleep_counter(15)
		print "ok... let's go!"
		return True
	else:
		print "Sorry, your pipeline cannot restart from this stage because the previous stage has not been sucessful."
		print "It seems that your project is at stage : " + LocalConfigFile['stage']
		print "You can run the PhyloProf with the --continue flag."
		print "Alternatively, you can run PhyloProf from these stages: \n " + "\t ".join(PIPELINE[:PIPELINE.index(current_stage) + 1])
		sys.exit()

def sleep_counter(N):
	sys.stdout.write('\033[91m' + '\033[1m' + 'Waiting')
	for i in range(N):
		sys.stdout.write('...' + str(i+1))
		time.sleep(1)
		sys.stdout.flush()
	print '\033[0m'

def merge_two_dicts( x, y):
	z = x.copy()
	z.update(y)
	return z

def fetchProtFams(ProjectName):
	main_path = os.getcwd()
	LocalConfigFile = get_cfg_file(ProjectName)
	if isTheRightOrder('fetchProtFams', LocalConfigFile):
		ProtFamDef = LocalConfigFile['ProtFamDef']
		# Deal with SeqDepot related Protein Families definitions
		SeqDepotQuery = mkSeqDepotQuery(ProtFamDef['SeqDepot'])
		SeqDepotResults, aseq2seq = searchSD(SeqDepotQuery)

		#Deal with custom HMM #
		# Not Implemented yet #
		#######################
		
		#Saving JSON files
		mkAseqJson(SeqDepotResults, ProjectName, main_path)
		mkAseq2seqJson(aseq2seq, ProjectName, main_path)
		
		#Update config file
		update_stage_phylopro_cfg(ProjectName, "fetchProtFams")

		#Next step in the PIPELINE
		fetchGenInfo(ProjectName)

def getProtFam( LocalConfigFile ):
	ProtFams = []
	for db in LocalConfigFile['ProtFamDef']:
		for instruction in LocalConfigFile['ProtFamDef'][db]:
			if instruction['name'] not in ProtFams:
				ProtFams.append(instruction['name'])
	return ProtFams


def separateAseq2Seq ( searchSD_results ):
	"""Not used"""
	aseq2seq = {}
	rest = []
	for data in searchSD_results:
		aseq2seq = merge_two_dicts( aseq2seq, data[1])
		rest.append(data[0])
	return rest, aseq2seq

def get_seqdepot_client(passwordfile = ""):
	if passwordfile == "":
			passwordfile = "/home/ortegad/private/mongodb_binf.txt"
	client = pymongo.MongoClient("aphrodite.bio.utk.edu",27017)
	client.the_database.authenticate('binf',open(passwordfile,"r").readline().strip(), mechanism='MONGODB-CR',source='admin')
	#print "SeqDepot connection authenticated"	
	return client.seqdepot

def processAseqsDic( results ):
	print "\n"
	merged_results = {}
	aseq2seq = {}
	for r in results:
		if r[0]['name'] not in merged_results.keys():
			merged_results[r[0]['name']] = set(r[0]['aseqs'])
		else:
			merged_results[r[0]['name']].union(set(r[0]['aseqs']))
		aseq2seq = merge_two_dicts( aseq2seq, r[1])
	for name in merged_results.keys():
		merged_results[name] = list(merged_results[name])
		print 'Number of ' + name + ' proteins =>\t' + str(len(merged_results[name]))
	return merged_results, aseq2seq

def singleSearchSD(query):
	n = query['name']
	q = query['query']
	print "Requesting information about : " + n
	sd = get_seqdepot_client("/Users/ortegad/private/mongodb_binf.txt")
	if seqLim == 0:
		cards = sd.aseqs.find( q , { '_id' : 1 , 's' : 1}, no_cursor_timeout=True)
	else:
		cards = sd.aseqs.find( q , no_cursor_timeout=True).limit(seqLim)
	results = []
	aseq2seq = {}
	allcards = []
	for card in cards:
		results.append(card['_id'])
		aseq2seq[card['_id']] = card['s']
		allcards.append(card)
	
	with open("test.json", 'w') as f:
		json.dump(allcards, f)	
	
	return [ { "name" : n , 'aseqs' : results }, aseq2seq ]
	
def searchSD( query_list ):
	NP = len(query_list)
	processes = multiprocessing.Pool(NP)
	results = processes.map( singleSearchSD, query_list)
	merged_results, aseq2seq = processAseqsDic(results)
	return merged_results, aseq2seq
	
def mkAseqJson (SeqDepotResults, ProjectName, main_path):
	for name, aseqs in SeqDepotResults.iteritems():
		filename = main_path + '/' + ProjectName + '/COGs/' + name + '.' + ProjectName + '.aseqs.json'
		print " Making file => " + filename
		with open(filename,'w') as f:
			json.dump(aseqs, f, indent = 2 )

def mkAseq2seqJson ( aseq2seq, ProjectName, main_path ):
	filename = main_path + '/' + ProjectName + '/COGs/All.' + ProjectName + '.aseqs2seq.json'
	print "\n Making file => " + filename
	with open(filename,'w') as f:
		json.dump(aseq2seq, f, indent = 2 )
		
def mkSeqDepotQuery(Instructions = {}):
	""" Builds the string to submit to SeqDepot """
	ProtFamSDQuery = []
	for ProtFamDefSD in Instructions:
		name = ProtFamDefSD['name']
		allInstructions = []
		for tool in ProtFamDefSD.keys():
			if tool not in ["name", "group"]:
				try:
					domCountIn = len(ProtFamDefSD[tool]["in"])
				except KeyError:
					domCountIn = 0
				try:
					domCountOut = len(ProtFamDefSD[tool]["out"])
				except KeyError:
					domCountOut = 0
				if domCountIn + domCountOut > 1:
					for i in range(domCountIn):
						allInstructions.append(' { "t.' + tool + '" : { "$elemMatch" : { "$elemMatch" : { "$in" : ["' + ProtFamDefSD[tool]["in"][i] + '"] }}}} ')
					if domCountOut > 0:
						allInstructions.append('{ "$nor" : [ { "t.' + tool + '" : { "$elemMatch" : { "$elemMatch" : { "$in" : [ "' + '","'.join(ProtFamDefSD[tool]["out"]) + '" ] }}}} ] } ')
				elif domCountOut == 0:
					allInstructions.append(' { "t.' + tool + '" : { "$elemMatch" : { "$elemMatch" : { "$in" : ["' + ProtFamDefSD[tool]["in"][0] + '"] }}}} ')
				elif domCountIn == 0:
					allInstructions.append(' { "t.' + tool + '" : { "$elemMatch" : { "$elemMatch" : { "$nin" : ["' + ProtFamDefSD[tool]["out"][0] + '"] }}}} ')
		if len(allInstructions) > 1:
			searchStrSD = '{ "$and" : [ ' + ','.join(allInstructions) + ']}'
		else:
			searchStrSD = allInstructions[0]
		
		#check if it is a valid json
		try:
			searchJsonSD = json.loads(searchStrSD)
		except ValueError:
			print searchStrSD
			print "There is something missing"
			sys.exit()
		
		#print json.dumps(searchJsonSD)
		
		ProtFamSDQuery.append( { "name" : name, "query" : searchJsonSD } )
	
	return ProtFamSDQuery

def fetchGenInfo(ProjectName):
	main_path = os.getcwd()
	LocalConfigFile = get_cfg_file(ProjectName)
	if isTheRightOrder('fetchGenInfo', LocalConfigFile):
		with open( main_path + '/' + ProjectName + '/COGs/All.' + ProjectName + '.aseqs2seq.json', 'r') as f:
			aseqs2seq = json.load(f)
		
		aseq2bitk_dic = aseq2bitk(aseqs2seq)
		mkAseq2bitkTagJson(aseq2bitk_dic, ProjectName, main_path)
		
		bitk2seq = {}
		notfound = []
		print " Number of Aseqs to be parsed : " + str(len(aseqs2seq.keys()))
		for i in range(len(aseqs2seq.keys())):
			
			aseq = aseqs2seq.keys()[i]
			
			if i % float(1000) == 0 and i > 999:
				print str(i) + " : " + aseq
			if aseq in aseq2bitk_dic.keys():
				for tag in aseq2bitk_dic[aseq]:
					bitk2seq[tag] = aseqs2seq[aseq]
			else:
				notfound.append(aseq) #NEED to do something with this.	
		
		print "Done"
		
		mkBitkTag2seqJson( bitk2seq, ProjectName, main_path)
		
		#Update config file
		update_stage_phylopro_cfg(ProjectName, "fetchGenInfo")
		
		#call next step:
		mkFastaFiles(ProjectName)
		

def get_mist22_client():
	print "Verifying tunnels"
	print "Mist"
	try:
		client = pymongo.MongoClient('localhost',27019)
		client.mist22.genes.find_one()
	except TypeError:
		print "You must open a tunnel with ares.bio.utk.edu: ssh -p 32790 -f -N -L 27019:localhost:27017 unsername@ares.bio.utk.edu"
		sys.exit()

	return client.mist22	
	
def aseq2bitk( aseqs2seq ):
	aseqs_list = aseqs2seq.keys()
	mist22 = get_mist22_client()
	
	if seqLim == 0:
		genes = mist22.genes.find( { 'p.aid' : { '$in' : aseqs_list }} )
	else:
		genes = mist22.genes.find( { 'p.aid' : { '$in' : aseqs_list }} )
	
	mistId_list = []
	aseq2bitk = {}
	
	for gene in genes:
		mistId = gene['gid']
		try:
			lo = gene['lo']
		except KeyError:
			lo = 'NULL'
		accession = gene['p']['ac']
		aseq = gene['p']['aid']
		bitkTag = str(mistId) + BITKTAGSEP + lo + BITKTAGSEP + accession 
		if mistId not in mistId_list:
			mistId_list.append(mistId)
		if aseq in aseq2bitk.keys():
			aseq2bitk[aseq].append(bitkTag)
		else:
			aseq2bitk[aseq] = [ bitkTag ]
	
	# getting genome Name 
	
	genomes = mist22.genomes.find( {'_id': {'$in' : mistId_list }}, {'sp': 1, 'g' : 1})
	
	gid2name = {}
	
	for genome in genomes:
		gid2name[genome['_id']] = genome['g'][:2] + BITKGENSEP + genome['sp'][:3] + BITKGENSEP
	
	for aseq in aseq2bitk.keys():
		for i in range(len(aseq2bitk[aseq])):
			mistId = int(aseq2bitk[aseq][i].split(BITKTAGSEP)[0])
			aseq2bitk[aseq][i] = gid2name[mistId] + aseq2bitk[aseq][i]
	
	return aseq2bitk

def mkAseq2bitkTagJson ( aseq2bitk_dic, ProjectName, main_path ):
	filename = main_path + '/' + ProjectName + '/COGs/All.' + ProjectName + '.aseqs2bitktag.json'
	print "\n Making file => " + filename
	with open(filename,'w') as f:
		json.dump(aseq2bitk_dic, f, indent = 2 )	
	
def mkBitkTag2seqJson ( aseq2seq, ProjectName, main_path ):
	filename = main_path + '/' + ProjectName + '/COGs/All.' + ProjectName + '.bitkseq2seq.json'
	print "\n Making file => " + filename
	with open(filename,'w') as f:
		json.dump(aseq2seq, f, indent = 2 )

def mkFastaFiles(ProjectName):
	main_path = os.getcwd()
	LocalConfigFile = get_cfg_file(ProjectName)
	if isTheRightOrder('mkFastaFiles', LocalConfigFile):
		ProtFams = getProtFam(LocalConfigFile)
		with open( main_path + '/' + ProjectName + '/COGs/All.' + ProjectName + '.bitkseq2seq.json', 'r') as f:
			bitk2seq = json.load(f)
		with open(main_path + '/' + ProjectName + '/COGs/All.' + ProjectName + '.aseqs2bitktag.json', 'r') as f:
			aseq2bitktag = json.load(f)
				
		for ProtFam in ProtFams:
			filename = main_path + '/' + ProjectName + '/COGs/' + ProtFam + '.' + ProjectName + '.aseqs.json'
			with open(filename, 'r') as f:
				infoJson = json.load(f)
			
			output = ''
			
			for aseq in infoJson:
				if aseq in aseq2bitktag.keys():
					for i in range(len(aseq2bitktag[aseq])):
						output += '>' + aseq2bitktag[aseq][i] + '\n' + bitk2seq[aseq2bitktag[aseq][i]] + '\n'
			
			filename = main_path + '/' + ProjectName + '/COGs/' + ProtFam + '.' + ProjectName + '.fa'
			print "\n Making file => " + filename
			with open(filename, 'w') as f:
				f.write(output)
	
	#Update config file
	update_stage_phylopro_cfg(ProjectName, "mkFastaFiles")
			
			
			
			
			
		
if __name__ == "__main__":
	#Parse flags.
	parser = argparse.ArgumentParser()
	group = parser.add_mutually_exclusive_group(required = True)
	group.add_argument("--init", type=str, help = 'Build the local directory system', metavar='ProjectName')
	group.add_argument("--init-force", type=str, help = 'Build the local directory system and erase any existing files. Use with caution.', metavar='ProjectName')
	group.add_argument("--continue", dest = "cont", nargs = "?", const = "", action = "store", help = 'Restart the pipeline',)
	group.add_argument("--fetchProtFams", help = 'Restart the pipeline at the fetchProtFam stage',)
	group.add_argument("--fetchGenInfo", help = 'Restart the pipeline at the fetchGenInfo stage',)
	group.add_argument("--mkFastaFiles", help = 'Restart the pipeline at the mkFastaFiles stage',)
	
	parser.add_argument("--test", action = 'store_true', help = 'Run searches with a limited number of sequences' )

	args = parser.parse_args()
	
	print args.cont

	if args.test:
		print "This will be only a test"
		seqLim = 10
	else:
		seqLim = 0

	if args.cont == "":
		ProjectName = get_ProjectName()
	elif args.cont != None:
		ProjectName = args.cont
	
	if args.cont != None:
		stage = get_stage(ProjectName)
		next_stage = PIPELINE[PIPELINE.index(stage) + 1]
		print "Continuing from stage : " + next_stage
		args = parser.parse_args([ "--" + next_stage , ProjectName])

	if args.init:
		ProjectName = args.init
		BuildDirectorySystem(ProjectName)
		print "Done with init stage. Recording changes to the local config file."
		update_stage_phylopro_cfg(ProjectName, "init")
		print "Awesome, we completed the init stage of PhyloPro. That was easy."
		print "Now, please edit the " + main_path + '/' + ProjectName + '/' + "phylopro." + ProjectName + ".cfg.json file before continue with the pipeline"
	elif args.init_force:
		ProjectName = args.init_force
		BuildDirectorySystem(ProjectName, force = True)
		print "Done with init stage. Recording changes to the local config file."
		update_stage_phylopro_cfg(ProjectName, "init")
		print "Awesome, we completed the init stage of PhyloPro. That was easy."
		print "Now, please edit the " + main_path + '/' + ProjectName + '/' + "phylopro." + ProjectName + ".cfg.json file before continue with the pipeline"
	elif args.fetchProtFams:
		print "Searching the database and making the relevant files"
		if args.fetchProtFams == "":
			ProjectName = get_ProjectName()
		else:
			ProjectName = args.fetchProtFams
		fetchProtFams(ProjectName)
	elif args.fetchGenInfo:
		print "Searching the database and making the relevant files"
		if args.fetchGenInfo == "":
			ProjectName = get_ProjectName()
		else:
			ProjectName = args.fetchGenInfo
		fetchGenInfo(ProjectName)
	elif args.mkFastaFiles:
		print "Searching the database and making the relevant files"
		if args.mkFastaFiles == "":
			ProjectName = get_ProjectName()
		else:
			ProjectName = args.mkFastaFiles
		mkFastaFiles(ProjectName)	
		