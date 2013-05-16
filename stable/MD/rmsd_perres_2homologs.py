#! /usr/bin/env python 
###################################
#    Davi Ortega 5/15/2013 
###################################
import bitk
import sys
import MDAnalysis as MD
import MDAnalysis.analysis.align as MDA
import MDAnalysis.selections.base as MDS

if '-h' in sys.argv:
	print 'Calculates the rmsd of aligned residues between two homologs. Output format is compatible to data import on VMD\n \
	Sintax: homolog1.pdb homolog2.pdb -aln homolog1and2aligned.fa \n '
	sys.exit()

pdb1 = sys.argv[1]
pdb2 = sys.argv[2]

ref = MD.Universe(pdb1)
trg = MD.Universe(pdb2)

if ref.residues.resnames() == trg.residues.resnames():
	if '-aln' in sys.argv:
		seq_dic, seq_list = bitk.fastareader(sys.argv[sys.argv.index('-aln')+1])
		if bitk.threeLetter2oneLetter(ref.residues.resnames()) == seq_dic[seq_list[0]].replace('-',''):
			alndic, alnlis = bitk.alnpos_dic(seq_dic[seq_list[0]], seq_dic[seq_list[1]])
			pick = 1
		elif bitk.threeLetter2oneLetter(ref.residues.resnames()) == seq_dic[seq_list[1]].replace('-',''):
			alndic, alnlis = bitk.alnpos_dic(seq_dic[seq_list[1]], seq_dic[seq_list[0]])
			pick = 2
		else:
			print("fasta file irrelevant for requested operation... ignoring")
			alnlis = ref.residues.resnames()
	else:
		alnlis = ref.residues.resnames()

	for res in alnlis:
		try:
			refcrd = ref.selectAtoms("backbone and resid " + str(res)).coordinates()
		        trgcrd = trg.selectAtoms("backbone and resid " + str(res)).coordinates()
		except IndexError:
			print res
			print alndic[res]
			sys.exit()
		rmsd = MDA.rmsd(refcrd, trgcrd)
		if alnlis == ref.residues.resnames() or pick == 1:
			print(str(res) + ' * ' + str(rmsd))
                        output = str(res) + ' * ' + str(rmsd) + '\n'
		elif pick == 2:
			print(str(alndic[res]) + ' * ' + str(rmsd))
                        output = str(alndic[res]) + ' * ' + str(rmsd) + '\n'
else:
	if '-aln' not in sys.argv:
		print("Structure from different homologs or truncated. Please insert alignment")
		sys.exit()
	else:
		seq_dic, seq_list = bitk.fastareader(sys.argv[sys.argv.index('-aln')+1]) 
		if bitk.threeLetter2oneLetter(ref.residues.resnames()) != seq_dic[seq_list[0]].replace('-',''):
			print("input order is important. Please, make the input order coherent with the alignment file")
			sys.exit()
		else:
			alndic, alnlis = bitk.alnpos_dic(seq_dic[seq_list[0]], seq_dic[seq_list[1]])
	for res in alnlis:
                try:
                        refcrd = ref.selectAtoms("backbone and resid " + str(res)).coordinates()
                        trgcrd = trg.selectAtoms("backbone and resid " + str(alndic[res])).coordinates()
                except IndexError:
                        print res
                        print alndic[res]
                        sys.exit()
                rmsd = MDA.rmsd(refcrd, trgcrd)
		print(str(res) + ' * ' + str(rmsd))
                output = str(res) + ' * ' + str(rmsd) + '\n'

sys.exit()


if '-aln' in sys.argv:
	aln = sys.argv[sys.argv.index('-aln')+1]
	seq_dic, seq_list = bitk.fastareader(aln)
	if len(seq_list) > 2:
		print "Alignment must contain 2 sequences only. Order is important... Ref, Target"
		sys.exit()
	else:
		alndic, alnlis = bitk.alnpos_dic(seq_dic[seq_list[0]], seq_dic[seq_list[1]])
else:
        aln = 0
        alndic = {}

		
if '-useressel' in sys.argv:
	aln = sys.argv[sys.argv.index('-useressel')+1]
	try:
		pick = sys.argv[sys.argv.index('-useressel')+2]
	except:
		print "-useressel fastafile (ref or trg)"
		sys.exit()
	seq_dic, seq_list = bitk.fastareader(aln)
	if len(seq_list) > 2:
		print "Alignment must contain 2 sequences only. Order is important... Ref, Target"
		sys.exit()
	else:
		alndic, alnlis = bitk.alnpos_dic(seq_dic[seq_list[0]], seq_dic[seq_list[1]])
else:
	pick = ''

ref = MD.Universe(pdb1)
trg = MD.Universe(pdb2)

output = ''

if aln == 0 or '-useressel' in sys.argv:
	if '-skipaln' not in sys.argv:
		MDA.alignto(trg, ref, mass_weighted=True)
	
	if '-useressel' in sys.argv:
		ressel = alnlis
	else:
		ressel = ref.residues.resids()

	for res in ressel:
                refcrd = ref.selectAtoms("backbone and resid " + str(res)).coordinates()
                trgcrd = trg.selectAtoms("backbone and resid " + str(res)).coordinates()
                rmsd = MDA.rmsd(refcrd, trgcrd)
        	if pick == 'ref':
		        print(str(res) + ' * ' + str(rmsd))
        	        output = str(res) + ' * ' + str(rmsd) + '\n'
		else:
			print(str(alndic[res]) + ' * ' + str(rmsd))
                        output = str(alndic[res]) + ' * ' + str(rmsd) + '\n'

	
	


#	if pick == 'ref':
#		for res in ressel:
#			refcrd = ref.selectAtoms("backbone and resid " + str(res)).coordinates()
#			trgcrd = trg.selectAtoms("backbone and resid " + str(res)).coordinates()
#			rmsd = MDA.rmsd(refcrd, trgcrd)
#	                print(str(res) + ' * ' + str(rmsd))
 #       	        output = str(res) + ' * ' + str(rmsd) + '\n'
#	else:
#		for res in ressel:
#			refcrd = ref.selectAtoms("backbone and resid " + str(alndic[res])).coordinates()
#			trgcrd = trg.selectAtoms("backbone and resid " + str(alndic[res])).coordinates()
#			rmsd = MDA.rmsd(refcrd, trgcrd)
#			print(str(alndic[res]) + ' * ' + str(rmsd))
#			output = str(alndic[res]) + ' * ' + str(rmsd) + '\n'
else:
#	print("Assuming structures aligned")
	
	for res in alnlis:
		try:
			refcrd = ref.selectAtoms("backbone and resid " + str(res)).coordinates()
			trgcrd = trg.selectAtoms("backbone and resid " + str(alndic[res])).coordinates()
		except:
			print res
			print alndic[res]
			sys.exit()
		rmsd = MDA.rmsd(refcrd, trgcrd)
		print(str(res) + ' * ' + str(rmsd))
		output = str(res) + ' * ' + str(rmsd) + '\n'

