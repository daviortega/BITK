#!/usr/bin/python
# Davi Ortega - April 2016

import os
import argparse
import json
import sys
import shutil
import pymongo
#import multiprocessing
import time
import datetime
import urllib2
import distutils
import networkx

#HardVariables

PIPELINE = ['init', 'fetchProtFams', 'fetchGenInfo', 'mkFastaFiles', 'filterByGen', 'trimSeqs', 'COGFinderBLAST', 'COGFinderParser', 'COGFinderMk' ]
BITKTAGSEP = '|' #TAG separator
BITKGENSEP = '_' #GENome separator

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

def fastaReader(datafile):
	
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

	makedir = [ main_path + '/' + ProjectName, main_path + '/' + ProjectName + '/COGs/' , main_path + '/' + ProjectName + '/merge/', main_path + '/' + ProjectName + '/trees/', main_path + '/' + ProjectName + '/pfamModels/' ]

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
		print "We are going to bind to Seqdepot and get some data.\n \
			This will take a while... hang tight!"
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

def getGroupFam ( LocalConfigFile ):
	GroupFams ={}
	for db in LocalConfigFile['ProtFamDef']:
		for instruction in LocalConfigFile['ProtFamDef'][db]:
			if instruction['group'] not in GroupFams:
				GroupFams[instruction['group']] = [ instruction['name'] ]
			else:
				GroupFams[instruction['group']].append( instruction['name'] )
	return GroupFams

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
	
	sd = get_seqdepot_client("/home/ortegad/private/mongodb_binf.txt")
	if seqLim == 0:
		print "Requesting information about : " + n
		cards = sd.aseqs.find( q , { '_id' : 1 , 's' : 1}, no_cursor_timeout=True)
	else:
		print "Requesting information about : " + n + " - TEST"
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
			
			if i % float(100) == 0 and i > 99:
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
	except pymongo.errors.ConnectionFailure:
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
	
	print "Getting info about the sequences... "
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
	
	print "Getting info about the genomes..."
	genomes = mist22.genomes.find( {'_id': {'$in' : mistId_list }}, {'sp': 1, 'g' : 1})
	
	gid2name = {}
	print "Building BITK style sequence headers..."
	for genome in genomes:
		gid2name[genome['_id']] = genome['g'][:2] + BITKGENSEP + genome['sp'][:3] + BITKGENSEP
	
	print "Building dictionary of sequence headers and Aseqs..."
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
	print "\n ==== Stage mkFastaFiles ====\n"
	main_path = os.getcwd()
	LocalConfigFile = get_cfg_file(ProjectName)
	if isTheRightOrder('mkFastaFiles', LocalConfigFile):
		ProtFams = getProtFam(LocalConfigFile)
		with open( main_path + '/' + ProjectName + '/COGs/All.' + ProjectName + '.bitkseq2seq.json', 'r') as f:
			bitk2seq = json.load(f)
		with open(main_path + '/' + ProjectName + '/COGs/All.' + ProjectName + '.aseqs2bitktag.json', 'r') as f:
			aseq2bitktag = json.load(f)
				
		for ProtFam in ProtFams:
			print " Parsing JSON file for " + ProtFam
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
	
	#Next step in the PIPELINE
	filterByGen(ProjectName)
	
			
def filterByGen(ProjectName):
	print "\n ==== Stage filterByGen ====\n"
	main_path = os.getcwd()
	LocalConfigFile = get_cfg_file(ProjectName)
	if isTheRightOrder('filterByGen', LocalConfigFile):		
		ProtFams = getProtFam(LocalConfigFile)
		# getting the genomes from config file
		try:	
			genomes = LocalConfigFile['genomes']
		except KeyError:
			genomes = []
		if genomes == []:
			print 'No genomes have been passed, we will keep going with the pipeline with whole dataset. There will be only a symbolic link to the original files and the "only" files. If this is not what you want, please hit Ctrl-C, check your config file in the key "genomes" and restart the pipeline with --continue flag'
			for ProtFam in ProtFams:
				filename_all = main_path + '/' + ProjectName + '/COGs/' + ProtFam + '.' + ProjectName + '.fa'
				filename_only = main_path + '/' + ProjectName + '/COGs/' + ProtFam + '.' + ProjectName + '.only.fa'
				print "\n ==> Creating symlinks for " + ProtFam
				
				try:
					os.remove(filename_only)
				except OSError:
					pass
				os.symlink(filename_all, filename_only )
		else:
			for ProtFam in ProtFams:
				print "\n ==> Filtering sequences of " + ProtFam
				filename_all = main_path + '/' + ProjectName + '/COGs/' + ProtFam + '.' + ProjectName + '.fa'
				filename_only = main_path + '/' + ProjectName + '/COGs/' + ProtFam + '.' + ProjectName + '.only.fa'
				fasta_tmp, tags = fastaReader(filename_all)
				new_fasta = ''
				counter = 0
				for tag in tags:
					mistID = int(tag.split(BITKTAGSEP)[0].split(BITKGENSEP)[-1])
					if mistID in genomes:
						counter += 1
						new_fasta += '>' + tag + '\n' + fasta_tmp[tag] + '\n'
				
					with open(filename_only, 'w') as f:
						f.write(new_fasta)
				print str(counter) + " / " + str(len(tags))
	#Update config file
	update_stage_phylopro_cfg(ProjectName, "filterByGen")
	
	#Stop here?
	if dontStopHere(ProjectName):	
		#Next step in the PIPELINE
		trimSeqs(ProjectName)			
						
def trimSeqs(ProjectName):
	print "\n ==== Stage trimSeqs ====\n"
	main_path = os.getcwd()
	LocalConfigFile = get_cfg_file(ProjectName)
	if isTheRightOrder('trimSeqs', LocalConfigFile):
		if 'trim' in LocalConfigFile.keys():
			trimInstructs = LocalConfigFile['trim']
			ProtFams = getProtFam(LocalConfigFile)
			GroupFams = getGroupFam( LocalConfigFile )
			for trimInstruc in trimInstructs:
				if 'pfam' in trimInstruc.keys():
					hmm_model = downloadPfamModel(trimInstruc['pfam'])
					for ProtFam in GroupFams[trimInstruc['group']]:

						fasta_file = main_path + '/' + ProjectName + '/COGs/' + ProtFam + '.' + 	ProjectName + '.only.fa'
						outputfile = main_path + '/' + ProjectName + '/COGs/' + ProtFam	+ '.output.hmm.dat'
						
						runHmmsearch( hmm_model, fasta_file, outputfile)
						fasta = parseHmmData( fasta_file, outputfile )
						
						with open( main_path + '/' + ProjectName + '/COGs/' + ProtFam + '.' + 	ProjectName + '.trimmed.fa' , 'w') as f:
							f.write(fasta)					
	#Update config file
	update_stage_phylopro_cfg(ProjectName, "trimSeqs")
	
	#Stop here?
	if dontStopHere(ProjectName):
		#Next step in the PIPELINE
		COGFinderBLAST( ProjectName )
	
					
def runHmmsearch( hmm_model, fasta_file, outputfile):
	from distutils.spawn import find_executable
	if find_executable('hmmsearch'):
		os.system('hmmsearch --noali --cut_tc ' + hmm_model + ' ' + fasta_file + ' > ' + outputfile)
	else:
		print 'You must install HMMER to use this step'
		print 'Please, install HMMER or take out the "trim" key from the config file'
		print 'Aborting...'
		sys.exit()					

def downloadPfamModel( pfamNumber ):
	main_path = os.getcwd()
	
	print ' ==> Downloading ' + pfamNumber
	
	url = 'http://pfam.xfam.org/family/' + pfamNumber + '/hmm'
	response = urllib2.urlopen(url)
	pfam_model = response.read()
	
	filename = main_path + '/' + ProjectName + '/pfamModels/' + pfamNumber + '.hmm'
	
	with open( filename , 'w') as f:
		f.write(pfam_model)
	
	return filename
	
def parseHmmData ( fasta_file, filename ):
	seq_dic, seq_list = fastaReader( fasta_file )
	output = ''
	tags = []
	with open( filename , 'r') as f:
		for line in f:
			if '>>' in line:
				name = line.split(' ')[1]
			if line[60:66].replace(' ','').isdigit() and name not in tags:
				start = int(line[60:66])
				end = int(line[68:74])
				output += '>' + name + '\n' + seq_dic[name][start:end] + '\n'
				tags.append(name)		
	return output

def getTrimmedGroups ( LocalConfigFile ):
	trimmed = []
	if 'trim' in LocalConfigFile.keys():
		trimInstructs = LocalConfigFile['trim']
		for instructions in trimInstructs:
			trimmed.append(instructions['group'])
	return trimmed

def COGFinderBLAST ( ProjectName ):
	print "\n ==== Stage COGFinderBLAST ====\n"
	main_path = os.getcwd()
	LocalConfigFile = get_cfg_file(ProjectName)
	if isTheRightOrder('COGFinderBLAST', LocalConfigFile):
		GroupFams = getGroupFam( LocalConfigFile )
		trimmed = getTrimmedGroups( LocalConfigFile )
		ProtFams = []
		files = []
		for group in GroupFams.keys():
			fileNames = ''
			if group in trimmed:
				typeFasta = '.trimmed'
			else:
				typeFasta = '.only'
			for ProtFam in GroupFams[group]:
				if ProtFam not in ProtFams:
					ProtFams.append(ProtFam)
					fileNames += main_path + '/' + ProjectName + '/COGs/' + ProtFam + '.' + ProjectName + typeFasta + '.fa '
			files.append([ProjectName, group, fileNames])
		print "Running BLAST with defaults : Please refer to documentation to see how to change it\n\n"
		blastRun( files )

		#Update config file
		update_stage_phylopro_cfg(ProjectName, "COGFinderBLAST")

		if dontStopHere(ProjectName):
			#Next step in the PIPELINE
			COGFinderParser( ProjectName )

def blastRun ( query_list ):
	NP = len(query_list)
	processes = multiprocessing.Pool(NP)
	processes.map( singleBlastRun, query_list)
	
def singleBlastRun( dataBlast ):
	#making the db
	ProjectName, group, filesNames = dataBlast
	main_path = os.getcwd()
	print "Making BLAST database for group " + group
	groupedFastaName = main_path + '/' + ProjectName + '/COGs/' + group + ".g." + ProjectName + ".temp.fa"
	os.system( 'cat ' + filesNames + ' > ' + groupedFastaName )
	dbname = main_path + '/' + ProjectName + '/COGs/temp.' + group + '.g.' + ProjectName + '.blastp.db'
	os.system('formatdb -i ' + groupedFastaName + ' -n ' + dbname)

	print "Running BLAST for group " + group

	blastname = main_path + '/' + ProjectName + '/COGs/outputBlast.' + group + ".g." + ProjectName + ".dat"

	os.system('blastp -db ' + dbname + ' -query ' + groupedFastaName + ' -out ' + blastname + ' -num_threads 10 -outfmt "6 qseqid sseqid bitscore evalue qlen length qcovs slen" -evalue 0.001 -max_target_seqs 100000' )

	os.system('rm ' + groupedFastaName )
	print "Done with " + group

def COGFinderParser ( ProjectName ):
	print "\n ==== Stage COGFinderParser ====\n"
	main_path = os.getcwd()
	LocalConfigFile = get_cfg_file(ProjectName)
	if isTheRightOrder('COGFinderParser', LocalConfigFile):
		GroupFams = getGroupFam( LocalConfigFile )
		trimmed = getTrimmedGroups( LocalConfigFile )
		ProtFams = []
		files = []
		for group in GroupFams.keys():
			fileList = []
			ProtFams2Go = []
			if group in trimmed:
				typeFasta = '.trimmed'
			else:
				typeFasta = '.only'
			for ProtFam in GroupFams[group]:
				if ProtFam not in ProtFams:
					ProtFams.append(ProtFam)
					ProtFams2Go.append(ProtFam)
					fileList.append(main_path + '/' + ProjectName + '/COGs/' + ProtFam + '.' + ProjectName + typeFasta + '.fa')
			files.append([ProjectName, group, ProtFams2Go, fileList])
		blastParser ( files )
		addCOGFindeMkrCfg ( ProjectName )
		#Update config file
		update_stage_phylopro_cfg(ProjectName, "COGFinderParser")

		if dontStopHere(ProjectName):
			#Next step in the PIPELINE
			COGFinderMk( ProjectName )

def blastParser ( query_list ):
	NP = len(query_list)
	processes = multiprocessing.Pool(NP)
	processes.map( singleBlastParser, query_list)

def singleBlastParser ( dataParse ):
	ProjectName, group, ProtFams, files = dataParse

	blastname = main_path + '/' + ProjectName + '/COGs/outputBlast.' + group + ".g." + ProjectName + ".dat"

	print "\t===> Parsing BLAST results for " + group 

	seqInfo_array = []	
	allSeqDic = {}
	allTagList = []
	tagFastaFileIndex = [] # this is a trick to quickly findout where the tag came from.  

	i = 0
	for eachFile in files:
		print "\t\t Loading ==> " + eachFile
		seqDic, tagList = fastaReader( eachFile )
		seqInfo_array.append( [ seqDic, tagList ])
		if len(seqInfo_array) == 1:
			allSeqDic = seqDic
			allTagList = tagList
		else:
			allSeqDic = merge_two_dicts( allSeqDic, seqDic )
			allTagList += tagList

		tagFastaFileIndex += [i] * len(tagList)
		i += 1
	
	print "\t Parsing the BLAST file ==> " + blastname 

	data_all = {}
	output = ''

	with open( blastname, 'r') as datafile:
		for line in datafile:
			field = line.replace('\n','').split('\t')
			
			qry = field[0]
			hit = field[1]
	
			rawvalues = [ float(i) for i in field[2:] ]
			if qry not in data_all.keys():
				#print allTagList.index[qry]
				data_all[qry] = {"pF" : ProtFams[ tagFastaFileIndex[allTagList.index(qry)]], 
									"h" : {} }
			if hit not in data_all[qry]['h'].keys():
				data_all[qry]['h'][hit] = rawvalues #[ float(i) for i in field[2:] ]
			elif  data_all[qry]['h'][hit][1] > rawvalues[1]:
				data_all[qry]['h'][hit] = rawvalues

	jsonFileName = main_path + '/' + ProjectName + '/COGs/outputBlast.' + group + ".g." + ProjectName + ".parsed.json"
	
	print "\t Saving parsed data for group : " + group 

	with open( jsonFileName , 'w') as f:
		json.dump( data_all, f )

def addCOGFindeMkrCfg ( ProjectName ):
	""" Update the PhyloPro config file with default settings for COGFinderMk """
	main_path = os.getcwd()
	
	with open( main_path + '/' + ProjectName + '/' + "phylopro." + ProjectName + ".cfg.json", 'r') as f:
		LocalConfigFile = json.load(f)

	GroupFams = getGroupFam( LocalConfigFile )

	if 'COGFinderMkCfg' not in LocalConfigFile.keys():
		print "\n Adding default settings for the COG process to the config file "
		LocalConfigFile['COGFinderMkCfg'] = {}

		for group in GroupFams:
			LocalConfigFile['COGFinderMkCfg'][group] = { "eValue" : 10E-120,
														 "qCov"   : 0.9,
														 "mergeRules"  : {},
													   } 
		with open( main_path + '/' + ProjectName + '/' + "phylopro." + ProjectName + ".cfg.json", 'w') as f:
			json.dump(LocalConfigFile, f, indent = 2)
	else:
		print "\n Settings already present in the config file: Do nothing !"

def COGFinderMk( ProjectName ):
	print "\n ==== Stage COGFinderMk ====\n"
	main_path = os.getcwd()
	LocalConfigFile = get_cfg_file(ProjectName)
	if isTheRightOrder('COGFinderMk', LocalConfigFile):
		GroupFams = getGroupFam( LocalConfigFile )
	
		dataQuery = []
		for group in GroupFams.keys():
			print group
			settings = LocalConfigFile['COGFinderMkCfg'][group]
			dataQuery.append( [ ProjectName, group, GroupFams[group], settings ])

		singleCOGMaker( dataQuery[0] )
		#COGMaker ( dataQuery )

def COGMaker ( query_list ):
	NP = len(query_list)
	processes = multiprocessing.Pool(NP)
	processes.map( singleCOGMaker, query_list)

def singleCOGMaker ( dataInfo ):
	ProjectName, group, ProtFams, settings = dataInfo
	jsonFileName = main_path + '/' + ProjectName + '/COGs/outputBlast.' + group + ".g." + ProjectName + ".parsed.json"
	dataCOG = {}
	groups = []

	with open( jsonFileName , 'r' ) as f:
		print "\t Loading the data ==> " + group
		data_all = json.load(f)
	print "\n"
	print "\t Parsing the data ==> " + group
	
	bestHits = {}
	groupDict = {}
	groups = []

	for qry in data_all.keys():
		bestHits[qry] = {}
		qryOrg = qry.split(BITKTAGSEP)[0]
		for hit in data_all[qry]['h'].keys():
			hitOrg = hit.split(BITKTAGSEP)[0]
			if hitOrg != qryOrg:
				dataRawFwd = data_all[qry]['h'][hit]
				try:
					dataRawBck = data_all[hit]['h'][qry]
				except KeyError:
					break
				if dataRawFwd[3]/dataRawFwd[2] > settings['qCov'] and dataRawBck[3]/dataRawBck[2] > settings['qCov'] and dataRawFwd[1] < settings['eValue'] and dataRawBck[1] < settings['eValue'] :
					evalue = min(dataRawFwd[1], dataRawBck[1])
					if hitOrg not in bestHits[qry].keys():
						bestHits[qry][hitOrg] = [ hit, evalue ]
					elif evalue < bestHits[qry][hitOrg][1]:
						bestHits[qry][hitOrg] = [ hit, evalue ]

	print "Ps_aer_495|PA14_56010|YP_792655.1 " + str(bestHits["Ps_aer_495|PA14_56010|YP_792655.1"])
	print "Ps_aer_479|PA4307|NP_252997.1 " + str(bestHits["Ps_aer_479|PA4307|NP_252997.1"])
	print "Ps_flu_1633|PSF113_0374|YP_005205803.1 " + str(bestHits["Ps_flu_1633|PSF113_0374|YP_005205803.1"])
	

	def scanSeqs( initTag, groupId, scanned, groups, dataDic, cogNet ):
		if initTag in scanned:
			return groupId, scanned, groups, dataDic, cogNet
		else:
			scanned.append(initTag)
			cogNet.add_node( initTag, grp = data_all[initTag]['pF'] , url = "#" + initTag)
			qryOrg = initTag.split(BITKTAGSEP)[0]
			if groupId == -1:
				groups.append([initTag])
				groupId = len(groups) - 1
			
			
			for hitOrg in dataDic[initTag].keys():
				hit = dataDic[initTag][hitOrg][0]
				

				if qryOrg in dataDic[hit].keys():
					if hit not in groups[groupId] and initTag == dataDic[hit][qryOrg][0]:
						if hit == "Ps_aer_495|PA14_56010|YP_792655.1":
							print "as hit"
							print groups[groupId]
							print initTag
							print groupId
						if initTag == "Ps_aer_495|PA14_56010|YP_792655.1":
							print "as qry"
							print groups[groupId]
							print initTag
							print groupId
						groups[groupId].append(hit)
						linkID = len(cogNet.edges())
						cogNet.add_edge( initTag, hit , weigth=min(data_all[initTag]['h'][hit][1], data_all[hit]['h'][initTag][1]), id = linkID )
						linkID += 1
						groupId, scanned, groups, dataDic, cogNet = scanSeqs( hit , groupId, scanned, groups, dataDic, cogNet )


			return groupId, scanned, groups, dataDic, cogNet

	groups = []
	scanned = []
	groupId = 0
	cogNet = networkx.Graph()

	for qry in bestHits.keys():
		groupId = -1
		groupId, scanned, groups, dataDic, cogNet = scanSeqs( qry, groupId, scanned, groups, bestHits, cogNet)

	print " ==> There are " + str(len(data_all.keys())) + " " + group
	print " ==> " + group + " proteins scanned :: " + str(len(scanned))

	inGroups = []

	COGs = 0
	for grp in groups:
		for g in grp:
			if g not in inGroups:
				inGroups.append(g)
			else:
				print "Sequence " + g + " present in 2 or more groups"
				for grp2 in groups:
					if g in grp2:
						print str(groups.index(grp2)) + '\t' + str(grp2)
				
		if len(grp) > 1:
			COGs += 1
	
	print " ==> Number of COGs with at least 2 " + group + " proteins :: " + str(COGs) 
	

	print [ cogNet.degree(k) for k in data_all.keys()] 
	print [ len(grp) for grp in groups ]
	print len(groups)
	print groups[:2]
	print cogNet.edges()[:5]

	def sorting_groups(groups, G):
		print "number of groups: " + str(len(groups))
		print "sorting..."
		grp_size = [ len(i) for i in groups ]
		grp_size = list(set(grp_size))
		grp_size.sort()
		grp_size.reverse()
		groups_order = []
		for i in grp_size:
			for group in groups:
				if len(group) == i:
					grp_deg = [ G.degree(k) for k in group]
	#               print grp_deg
					grp_deg = list(set(grp_deg))
					grp_deg.sort()
					grp_deg.reverse()
					group_sorted = []
					for j in grp_deg:
						for tag in group:
							if j == G.degree(tag):
								group_sorted.append(tag)
					groups_order.append(group_sorted)
	#   for i in range(len(groups)):
	#       index[i] = groups_order.index(groups[i])
		return groups_order

	newgroup = sorting_groups(groups, cogNet)

	print [ len(grp) for grp in newgroup ]

	print newgroup[0]
	return 0


def dontStopHere( ProjectName ):
	main_path = os.getcwd()
	LocalConfigFile = get_cfg_file(ProjectName)
	if 'stop' in LocalConfigFile.keys():
		if LocalConfigFile['stage'] == LocalConfigFile['stop']:
			return None
	return True	
	 
			
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
	group.add_argument("--filterByGen", help = 'Restart the pipeline at the filterByGen stage',)
	group.add_argument("--trimSeqs", help = 'Restart the pipeline at the trimSeqs stage',)
	group.add_argument("--COGFinderBLAST", help = ' Restart pipeline at the COGFinderBLAST stage')
	group.add_argument("--COGFinderParser", help = ' Restart pipeline at the COGFinderParser stage')
	group.add_argument("--COGFinderMk", help = ' Restart pipeline at the COGFinderMk stage')
	
	parser.add_argument("--test", action = 'store_true', help = 'Run searches with a limited number of sequences' )

	args = parser.parse_args()
	main_path = os.getcwd()
	
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
		if args.mkFastaFiles == "":
			ProjectName = get_ProjectName()
		else:
			ProjectName = args.mkFastaFiles
		mkFastaFiles(ProjectName)	
	elif args.filterByGen:
		if args.filterByGen == "":
			ProjectName = get_ProjectName()
		else:
			ProjectName = args.filterByGen
		filterByGen(ProjectName)
	elif args.trimSeqs:
		if args.trimSeqs == "":
			ProjectName = get_ProjectName()
		else:
			ProjectName = args.trimSeqs
		trimSeqs(ProjectName)
	elif args.COGFinderBLAST:
		if args.COGFinderBLAST == "":
			ProjectName = get_ProjectName()
		else:
			ProjectName = args.COGFinderBLAST
		COGFinderBLAST(ProjectName)
	elif args.COGFinderParser:
		if args.COGFinderParser == "":
			ProjectName = get_ProjectName()
		else:
			ProjectName = args.COGFinderParser
		COGFinderParser(ProjectName)
	elif args.COGFinderMk:
		if args.COGFinderMk == "":
			ProjectName = get_ProjectName()
		else:
			ProjectName = args.COGFinderMk
		COGFinderMk(ProjectName)