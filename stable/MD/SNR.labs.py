#! /usr/bin/env python 
###################################
#    Davi Ortega 2/11/2013 
###################################
import sys
import os

import numpy
import scipy.integrate
import MDAnalysis as MD
import MDAnalysis.analysis.distances as MDdist
import matplotlib.pylab as mp


if '-h' in sys.argv:
	print("To be used with the MCP (Anton) data. Input: pdbfile1.pdb, dcdfile1.pdb number_of_layers cutoff_for_the_fourier_filter tip_coordinate\n\n  \
Calculates: \n \
1) C - N', C'- N, C - N, C' - N', C - C' and N - N' distances (alpha carbon) \n \
2) Apply fourier filter\n \
3) Calculates Amplitude\n \
4) Calculate SNR.\n\n \
Options: -plot to plot FFT spectrum and filtered data.")
	sys.exit()


#Main
pdb_file = sys.argv[1]
dcd_file = sys.argv[2]
num_layers = int(sys.argv[3])
cutoff = int(sys.argv[4])
tip = int(sys.argv[5])

print "Loading trajectory"
mol = MD.Universe(pdb_file, dcd_file)

name_file_dcd = sys.argv[2][:-3]
while "/" in name_file_dcd:
        name_file_dcd = name_file_dcd[name_file_dcd.index("/")+1:]

#layer_list = []

#for i in range(num_layers):
#	layer_list.append([391 - i, 391 + i])

F_data_header = []
F_data = []

output = 'Layer\tL\tMS\tSNR_CANB\tSNR_CBNA\tSNR_CACB\tSNR_CANA\tSNR_CBNB\tSNR_NANB\n'
hdisp_R = 'Frame\tLayer\tMS\tCA2NB\tCB2NA\tCA2CB\tCA2NA\tCB2NB\tNA2NB\n' 


if '-test' in sys.argv:
	LIST = range(1,5)
else:
	LIST = range(1,num_layers)



for L in LIST:
#for L in [66]:

	res_C = tip - L
	res_N = tip + L

	C_A = mol.selectAtoms("resid " + str(res_C) + " and name CA and segid A")
        C_B = mol.selectAtoms("resid " + str(res_C) + " and name CA and segid B")
	N_A = mol.selectAtoms("resid " + str(res_N) + " and name CA and segid A")
        N_B = mol.selectAtoms("resid " + str(res_N) + " and name CA and segid B")

	CA2NB = []
	CB2NA = []
	CA2CB = []
	CA2NA = []
	CB2NB = []
	NA2NB = []
	

	f = 0
	for ts1 in mol.trajectory:

		CA2NB.append(MDdist.dist(C_A, N_B)[2])
		CB2NA.append(MDdist.dist(C_B, N_A)[2])
		CA2CB.append(MDdist.dist(C_A, C_B)[2])
                CA2NA.append(MDdist.dist(C_A, N_A)[2])
		CB2NB.append(MDdist.dist(C_B, N_B)[2])
                NA2NB.append(MDdist.dist(N_A, N_B)[2])


		hdisp_R	+= str(f) + '\t' + str(res_C) + '_' + str(res_N) + '\t' + name_file_dcd.split('.')[0] + '\t%.3f' % MDdist.dist(C_A, N_B)[2] + '\t%.3f' % MDdist.dist(C_B, N_A)[2] + '\t%.3f' % MDdist.dist(C_A, C_B)[2] + '\t%.3f' % MDdist.dist(C_A, N_A)[2] + '\t%.3f' % MDdist.dist(C_B, N_B)[2] + '\t%.3f\n' % MDdist.dist(N_A, N_B)[2]
		f += 1
				

	CA2NB = numpy.array(CA2NB)
	CB2NA = numpy.array(CB2NA)
	CA2CB = numpy.array(CA2CB)
        CA2NA = numpy.array(CA2NA)
	CB2NB = numpy.array(CB2NB)
        NA2NB = numpy.array(NA2NB)

	#take DC out:
#	CA2NB -= numpy.mean(CA2NB)
#       CB2NA -= numpy.mean(CB2NA)

	FFT_CA2NB = numpy.fft.rfft(CA2NB, axis = 0)
	FFT_CB2NA = numpy.fft.rfft(CB2NA, axis = 0)
	FFT_CA2CB = numpy.fft.rfft(CA2CB, axis = 0)
        FFT_CA2NA = numpy.fft.rfft(CA2NA, axis = 0)
	FFT_CB2NB = numpy.fft.rfft(CB2NB, axis = 0)
        FFT_NA2NB = numpy.fft.rfft(NA2NB, axis = 0)



	if '-plot' in sys.argv:
		mp.plot(abs(FFT_CB2NA), '-b')
        	mp.grid(True)
	#       mp.xscale('log')
		mp.yscale('log')
	        mp.show()

	#apply filter:
	
	LF_CA2NB = numpy.zeros(len(FFT_CA2NB), dtype=complex)
	LF_CA2NB[range(cutoff)] = FFT_CA2NB[range(cutoff)]
	LF_CA2NB[0] = 0
	HF_CA2NB = numpy.zeros(len(FFT_CA2NB), dtype=complex)
	HF_CA2NB[range(cutoff,len(HF_CA2NB))] = FFT_CA2NB[range(cutoff,len(HF_CA2NB))]
	

	LF_CB2NA = numpy.zeros(len(FFT_CB2NA), dtype=complex)
	LF_CB2NA[range(cutoff)] = FFT_CB2NA[range(cutoff)]
	LF_CB2NA[0] = 0
        HF_CB2NA = numpy.zeros(len(FFT_CB2NA), dtype=complex)
	HF_CB2NA[range(cutoff,len(HF_CB2NA))] = FFT_CB2NA[range(cutoff,len(HF_CB2NA))]

	LF_CA2CB = numpy.zeros(len(FFT_CA2CB), dtype=complex)
        LF_CA2CB[range(cutoff)] = FFT_CA2CB[range(cutoff)]
	LF_CA2CB[0] = 0
        HF_CA2CB = numpy.zeros(len(FFT_CA2CB), dtype=complex)
        HF_CA2CB[range(cutoff,len(HF_CA2CB))] = FFT_CA2CB[range(cutoff,len(HF_CA2CB))]

        LF_CA2NA = numpy.zeros(len(FFT_CA2NA), dtype=complex)
        LF_CA2NA[range(cutoff)] = FFT_CA2NA[range(cutoff)]
	LF_CA2NA[0] = 0
        HF_CA2NA = numpy.zeros(len(FFT_CA2NA), dtype=complex)
        HF_CA2NA[range(cutoff,len(HF_CA2NA))] = FFT_CA2NA[range(cutoff,len(HF_CA2NA))]

	LF_CB2NB = numpy.zeros(len(FFT_CB2NB), dtype=complex)
        LF_CB2NB[range(cutoff)] = FFT_CB2NB[range(cutoff)]
	LF_CB2NB[0] = 0
        HF_CB2NB = numpy.zeros(len(FFT_CB2NB), dtype=complex)
        HF_CB2NB[range(cutoff,len(HF_CB2NB))] = FFT_CB2NB[range(cutoff,len(HF_CB2NB))]

        LF_NA2NB = numpy.zeros(len(FFT_NA2NB), dtype=complex)
        LF_NA2NB[range(cutoff)] = FFT_NA2NB[range(cutoff)]
	LF_NA2NB[0] = 0
        HF_NA2NB = numpy.zeros(len(FFT_NA2NB), dtype=complex)
        HF_NA2NB[range(cutoff,len(HF_NA2NB))] = FFT_NA2NB[range(cutoff,len(HF_NA2NB))]


	SIG_CA2NB = numpy.fft.irfft(LF_CA2NB, axis = 0)
	NOI_CA2NB = numpy.fft.irfft(HF_CA2NB, axis = 0)
	
	SIG_CB2NA = numpy.fft.irfft(LF_CB2NA, axis = 0)
        NOI_CB2NA = numpy.fft.irfft(HF_CB2NA, axis = 0)

	SIG_CA2CB = numpy.fft.irfft(LF_CA2CB, axis = 0)
        NOI_CA2CB = numpy.fft.irfft(HF_CA2CB, axis = 0)

        SIG_CA2NA = numpy.fft.irfft(LF_CA2NA, axis = 0)
        NOI_CA2NA = numpy.fft.irfft(HF_CA2NA, axis = 0)

	SIG_CB2NB = numpy.fft.irfft(LF_CB2NB, axis = 0)
        NOI_CB2NB = numpy.fft.irfft(HF_CB2NB, axis = 0)

        SIG_NA2NB = numpy.fft.irfft(LF_NA2NB, axis = 0)
        NOI_NA2NB = numpy.fft.irfft(HF_NA2NB, axis = 0)

	if "-plot" in sys.argv:
		CB2NA = numpy.reshape(CB2NA, len(CB2NA))
#		mp.plot(SIG_CA2NB,'-r')	
		mp.plot(CB2NA, '-b', NOI_CB2NA, '-g', SIG_CB2NA,'-r')
		mp.grid(True)
		mp.show()

	Samp_CA2NB = numpy.sqrt(scipy.integrate.trapz(SIG_CA2NB*SIG_CA2NB)/len(SIG_CA2NB))
	Namp_CA2NB = numpy.sqrt(scipy.integrate.trapz(NOI_CA2NB*NOI_CA2NB)/len(NOI_CA2NB))
	
	Samp_CB2NA = numpy.sqrt(scipy.integrate.trapz(SIG_CB2NA*SIG_CB2NA)/len(SIG_CB2NA))
        Namp_CB2NA = numpy.sqrt(scipy.integrate.trapz(NOI_CB2NA*NOI_CB2NA)/len(NOI_CB2NA))

	Samp_CA2CB = numpy.sqrt(scipy.integrate.trapz(SIG_CA2CB*SIG_CA2CB)/len(SIG_CA2CB))
        Namp_CA2CB = numpy.sqrt(scipy.integrate.trapz(NOI_CA2CB*NOI_CA2CB)/len(NOI_CA2CB))

        Samp_CA2NA = numpy.sqrt(scipy.integrate.trapz(SIG_CA2NA*SIG_CA2NA)/len(SIG_CA2NA))
        Namp_CA2NA = numpy.sqrt(scipy.integrate.trapz(NOI_CA2NA*NOI_CA2NA)/len(NOI_CA2NA))

	Samp_CB2NB = numpy.sqrt(scipy.integrate.trapz(SIG_CB2NB*SIG_CB2NB)/len(SIG_CB2NB))
        Namp_CB2NB = numpy.sqrt(scipy.integrate.trapz(NOI_CB2NB*NOI_CB2NB)/len(NOI_CB2NB))

        Samp_NA2NB = numpy.sqrt(scipy.integrate.trapz(SIG_NA2NB*SIG_NA2NB)/len(SIG_NA2NB))
        Namp_NA2NB = numpy.sqrt(scipy.integrate.trapz(NOI_NA2NB*NOI_NA2NB)/len(NOI_NA2NB))



	print( str(res_C) + '-' + str(res_N) + '\t' + str(L) + '\t' + str(Samp_CA2NB/Namp_CA2NB) + '\t' + str(Samp_CB2NA/Namp_CB2NA) + '\t' + str(Samp_CA2CB/Namp_CA2CB) + '\t' + str(Samp_CA2NA/Namp_CA2NA) + '\t' + str(Samp_CB2NB/Namp_CB2NB) + '\t' + str(Samp_NA2NB/Namp_NA2NB))

	output += str(res_C) + '-' + str(res_N) + '\t' + str(L) + '\t' + name_file_dcd.split('.')[0] + '\t' +  str(Samp_CA2NB/Namp_CA2NB) + '\t' + str(Samp_CB2NA/Namp_CB2NA) + '\t' +  str(Samp_CA2CB/Namp_CA2CB) + '\t' + str(Samp_CA2NA/Namp_CA2NA) + '\t' +  str(Samp_CB2NB/Namp_CB2NB) + '\t' + str(Samp_NA2NB/Namp_NA2NB) + '\n'

	F_data.append(SIG_CA2NB)
	F_data_header.append(str(res_C) + '-' + str(res_N) + '_SIG-CA2NB')   
        F_data.append(SIG_CB2NA)
	F_data_header.append(str(res_C) + '-' + str(res_N) + '_SIG-CB2NA')
        F_data.append(SIG_CA2CB)
	F_data_header.append(str(res_C) + '-' + str(res_N) + '_SIG-CA2CB')
        F_data.append(SIG_CA2NA)
	F_data_header.append(str(res_C) + '-' + str(res_N) + '_SIG-CA2NA')
        F_data.append(SIG_CB2NB)
	F_data_header.append(str(res_C) + '-' + str(res_N) + '_SIG-CB2NB')
        F_data.append(SIG_NA2NB)
	F_data_header.append(str(res_C) + '-' + str(res_N) + '_SIG-NA2NB')

dataout = open(name_file_dcd.split('.')[0] + '.SNR.dat', 'w')
dataout.write(output)
dataout.close()

dataout =  open(name_file_dcd.split('.')[0] + '.hdisp_R.dat', 'w')
dataout.write(hdisp_R)
dataout.close()

F_data = numpy.array(F_data).T

numpy.savetxt(name_file_dcd.split('.')[0] + '.F_data.dat', F_data, fmt = '%.3f', delimiter = '\t')
Fout = open(name_file_dcd.split('.')[0] + '.F_data_header.dat', 'w')
Fout.write('\t'.join(F_data_header) + '\n')
Fout.close()

os.system('cat ' + name_file_dcd.split('.')[0] + '.F_data_header.dat ' + name_file_dcd.split('.')[0] + '.F_data.dat > ' + name_file_dcd.split('.')[0] + '.filtered.dat')
os.system('rm ' + name_file_dcd.split('.')[0] + '.F_data_header.dat ' + name_file_dcd.split('.')[0] + '.F_data.dat')

sys.exit()




