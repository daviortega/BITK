#!/usr/bin/python
# Davi Ortega - April 2016

import os
import argparse
import json
import sys
import shutil

#HardVariables

PIPELINE = ['init', 'mkfiles', ]

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
						"stage"		  : ""}

	with open(makedir[0] + '/' + cfgfilename, 'w') as f:
		json.dump(LocalConfigFile, f, indent = 2)

def update_phylopro_cfg (ProjectName, stage):
	""" Update the PhyloPro config file to state the latest sucessful step """
	main_path = os.getcwd()
	
	with open( main_path + '/' + ProjectName + '/' + "phylopro." + ProjectName + ".cfg.json", 'r') as f:
		LocalConfigFile = json.load(f)

	LocalConfigFile['stage'] = stage

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
	else:
		print "Sorry, your pipeline cannot restart from this stage because the previous stage has not been sucessful."
		print "It seems that your project is at stage : " + LocalConfigFile['stage']
		print "You can run the PhyloProf with the --continue flag."
		print "Alternatively, you can run PhyloProf from these stages: " + " ".join(PIPELINE[:PIPELINE.index(current_stage)])
		sys.exit()
	

def mkfiles(ProjectName):
	LocalConfigFile = get_cfg_file(ProjectName)
	if isTheRightOrder('mkfiles', LocalConfigFile):
		LocalConfigFile = get_cfg_file(ProjectName)
		ProtFamDef = LocalConfigFile['ProtFamDef']
		SeqDepotQuery = mkSeqDepotQuery(ProtFamDef['SeqDepot'])
		print SeqDepotQuery
		
		
def mkSeqDepotQuery(Instructions = {}):
	""" Builds the string to submit to SeqDepot """
	ProtFamSDQuery = {}
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
					for i in range(domCountOut):
						allInstructions.append(' { "t.' + tool + '" : { "$elemMatch" : { "$elemMatch" : { "$nin" : ["' + ProtFamDefSD[tool]["out"][i] + '"] }}}} ')
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
		if name in ProtFamSDQuery.keys():
			ProtFamSDQuery[name].append(searchStrSD)
		else:
			ProtFamSDQuery[name] = [ searchStrSD ]
	
	return ProtFamSDQuery

if __name__ == "__main__":
	#Parse flags.
	parser = argparse.ArgumentParser()
	group = parser.add_mutually_exclusive_group(required = True)
	group.add_argument("--init", type=str, help = 'Build the local directory system', metavar='ProjectName')
	group.add_argument("--init-force", type=str, help = 'Build the local directory system and erase any existing files. Use with caution.', metavar='ProjectName')
	group.add_argument("--continue", dest = "cont", nargs = "?", const = "", action = "store", help = 'Restart the pipeline',)
	group.add_argument("--mkfiles", help = 'Restart the pipeline ath the mkfiles stage',)
	
	args = parser.parse_args()

	
	print args.cont

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
		update_phylopro_cfg(ProjectName, "init")
	elif args.init_force:
		ProjectName = args.init_force
		BuildDirectorySystem(ProjectName, force = True)
		print "Done with init stage. Recording changes to the local config file."
		update_phylopro_cfg(ProjectName, "init")
	elif args.mkfiles:
		print "Searching the database and making the relevant files"
		if args.mkfiles == "":
			ProjectName = get_ProjectName()
		else:
			ProjectName = args.mkfiles
		mkfiles(ProjectName)
		
		
			
			








