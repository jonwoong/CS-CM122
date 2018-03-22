########## HEADERS ##########

import numpy as np
import pandas as pd
from helper_functions import *
import itertools
import operator
import random
from collections import defaultdict, namedtuple

########## MAIN ##########

def main(argv=None):
	# store strains into an array known_strains[]
	with open('hw4_W_0/hw4_W_0_strains.txt') as strains_file:
		strains_file_lines = strains_file.read().splitlines()
	strains_file.close()

	known_strains = []
	for strain in range(1,len(strains_file_lines)):
		known_strains.append(strains_file_lines[strain])

	number_of_strains = len(known_strains)
	
	# store reads into an array reads[]
	with open('hw4_W_0/hw4_W_0_reads.txt') as reads_file:
		reads_file_lines = reads_file.read().splitlines()
	reads_file.close()

	reads = []
	for read in range(1,len(reads_file_lines)):
		reads.append(reads_file_lines[read])

	# calculate coverage
	coverage = (len(reads)*50)/(4*len(known_strains[1]))

	# build reference index for each strain
	strain_index = []
	for strain in range(number_of_strains):
		strain_index.append(build_index(known_strains[strain]))

	# align reads to all strains
	strain_matrix = [defaultdict(list) for i in range(number_of_strains)] # for each strain, store position of snp and whether a read matched that strain
	bases = [defaultdict() for i in range(number_of_strains)] # stores the real base letter for all strains
	for read_number in range(len(reads)):
		collected_SNPs = defaultdict(list) # a list of (position,read base,strain base) of a snp
		for strain_number in range(number_of_strains):
			SNPs = align(reads[read_number],known_strains[strain_number],strain_index[strain_number]) # a list of (position,read base,strain base) of a snp
			if len(SNPs) > 0:
				collected_SNPs[strain_number].append(SNPs) 
		unique_snps = set() # all snps found for this iteration of 50-base long read in reads[]
		for strain_number in range(number_of_strains):
			for snp_list in collected_SNPs[strain_number]:
				if len(snp_list) > 0:
					for snp in snp_list:
						if snp not in unique_snps:
							unique_snps.add(snp)
		# for every unique (position,read base,strain base) snp that this 50-base read found when aigning to all 4 strains,
		# check which strain it actually comes from
		for snp in unique_snps:
			for strain_number in range(number_of_strains):
				strain = known_strains[strain_number]
				strain_base = known_strains[strain_number][snp[0]]
				read_base = snp[1]
				if strain_base == read_base: # the snp came from this strain_number
					strain_matrix[strain_number][snp[0]].append(1)
					bases[strain_number][snp[0]] = snp[2]
				else: # the snp did not come from this strain_number
					strain_matrix[strain_number][snp[0]].append(0)
				#if strain_base != read_base: # the snp came from this strain_number
				#	strain_matrix[strain_number][snp[0]].append(1)
					#bases[strain_number][snp[0]] = snp[2]
	#print strain_matrix[2][8288]
	#print_list(bases)

	number_of_SNPs = len(strain_matrix[0].keys())
	for strain_number in range(len(strain_matrix)):
		for snp_position in strain_matrix[strain_number].keys():
			snp_index = strain_matrix[strain_number].keys().index(snp_position)
			strain_matrix[strain_number][snp_index] = strain_matrix[strain_number].pop(snp_position)

	# count snps
	snp_counts = defaultdict(list)
	for strain_number in range(number_of_strains):
		for snp_number in strain_matrix[strain_number].keys():
			snp_counts[snp_number] += strain_matrix[strain_number][snp_number]
	#print snp_counts

	# calculate snp/variant frequencies
	snp_freqs = calculate_snp_freq(snp_counts)
	print snp_freqs

	strain_snp_freq = [[0 for i in range(number_of_SNPs)] for j in range(number_of_strains)]
	for strain_number in range(number_of_strains):
		for snp_number in strain_matrix[strain_number]:
			strain_snp_freq[strain_number][snp_number] = float(sum(strain_matrix[strain_number][snp_number]))/len(strain_matrix[strain_number][snp_number])
			#strain_snp_freq[strain_number][snp_number] = float(sum(strain_matrix[strain_number][snp_number]))/len(snp_counts[snp_number])
			#strain_snp_freq[strain_number][snp_number] = sum(strain_matrix[strain_number][snp_number])
	#print_list(strain_snp_freq)

	#for strain_number in range(number_of_strains):
	#	print sum(strain_snp_freq[strain_number])

	snp_freq_in_order = sorted(snp_freqs,reverse=True)
	#print snp_freq_in_order

	unique_snp_freqs = []
	for freq in snp_freq_in_order:
		if freq not in unique_snp_freqs:
			unique_snp_freqs.append(freq)
	#print unique_snp_freqs

	# build strain matrix
	strains = []
	test_strains = generate_strain_matrix(n_strains=4,n_snps=7)
	f1 = np.linalg.lstsq(np.array(test_strains).T,snp_freqs)
	print f1[0]

	f2 = get_freqs_from_strains_and_counts(snp_counts,test_strains)
	print f2


	
	#hw_0
	#strains.append([0,1,1,1,1,1,0])
	#strains.append([1,0,0,0,0,1,1])
	#strains.append([0,0,1,1,1,1,1])
	#strains.append([0,0,0,0,1,0,0])
	#test_freq = [0.56,0.78,0.7,0.64,0.76,0.94,0.61]
	#x = np.linalg.lstsq(np.array(strains).T,np.array(test_freq))
	#print x[0]
	#hw_1
	#strains.append([1,0,0,1,1,0,0])
	#strains.append([1,1,1,1,0,1,1])
	#strains.append([0,1,0,0,1,0,0])
	#strains.append([0,0,0,0,1,0,0])
	# how did i get this strain matrix?
	# i looked at strain_snp_freq for columns that pairs which added up to 1.0
	# the 3rd column satisfies this
	# strain 0 column 3: 0.64
	# strain 1 column 3: 0.64
	# strain 2 column 3: 0.36
	# strain 3 column 3: 0.36
	# I know that one of the strains must have frequency 0.5, since this is one of the variant frequencies in snp_freqs
	# So I let strain 1 have a strain frequency of 0.5 and strain 2 have a frequency of 0.14 (they sum to 0.64)
	# And I do the same for strain 2 and strain 3: 
	# I split 0.36 into 0.3 and 0.06
	# so let strain 3 have frequency 0.3 and strain 2 have frequency 0.06
	# I then fill the strain matrix by hand like so
	# 
	#			SNP
	# STRAIN 0 (0.14): 	1		0		0		1 		1 		0 		0
	# STRAIN 1 (0.5):	1		1		1		1		0		1		1
	# STRAIN 2 (0.06):	0		1		0		0		1		0		0
	# STRAIN 3 (0.3):	0		0		0		0		1		0		0
	# SUM:				0.64 	0.56	0.5		0.36	0.5		0.5		0.5
	# VARIANT FREQUENCY:0.61	0.54	0.5		0.61	0.5		0.51	0.5
	#
	# I wasn't sure how to code this process, but I worked it out by hand for hw_W_0 and hw_W_1
	# and it seems to work?
	# Im confident about the strain frequency numbers but not about which strain gets what frequency

	#final_strain_freqs = get_freqs_from_strains_and_counts(snp_counts, np.array(strains))
	#print final_strain_freqs

	#with open("hw4_w_1_ans.txt","w+") as answer_file:
	#	answer_file.write(">hw4_W_1\n")
	#	for strain_number in range(number_of_strains):
	#		entry = str(final_strain_freqs[strain_number]) + ',' + str(known_strains[strain_number]) + '\n'
	#		answer_file.write(entry)
	#answer_file.close()

if __name__=="__main__":
	main()
