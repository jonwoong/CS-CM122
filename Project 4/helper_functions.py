import numpy as np
import pandas as pd
import itertools
import operator
import random
from collections import defaultdict, namedtuple

def generate_test_SNPs(coverage, set_seed, strain_freqs, n_snps):
    ## We're going to assume there is overlap between
    ## Strains for each SNP.
    
    if set_seed: np.random.seed(96000)
    
    n_strains = len(strain_freqs)
    assert sum(strain_freqs) == 1
    
    powerset_iterable = itertools.chain.from_iterable(itertools.combinations(range(n_strains), r)
                                             for r in range(2, n_strains))
    
    frozen_powerset = [index_tuple for index_tuple in powerset_iterable]
    
    while True:
        these_indices = np.random.choice(frozen_powerset, n_snps)
        output_strains = [[] for i in range(n_strains)]
        for index_list in these_indices:
            for i in range(len(output_strains)):
                if i in index_list:
                    output_strains[i].append(1)
                else:
                    output_strains[i].append(0)
                    
        output_strains = [tuple(strain) for strain in output_strains]
        n_unique_strains = len(set(output_strains))
        # test if it satisfies a few conditions 
        # if it does, we return
        #otherwise; loop back and try again
        if n_unique_strains == n_strains:
            break
    
    #for freq, strain in zip(strain_freqs, output_strains):
    #     print freq, strain
        
    
    sampled_strains = np.random.choice(range(n_strains),
                                       size=coverage*n_strains,
                                       p=strain_freqs)
    #print sampled_strains[:10]
    
    output_snp_counts = defaultdict(list)
    
    for strain_index in sampled_strains:
        snp_index = np.random.choice(range(n_snps))
        snp_value = output_strains[strain_index][snp_index]
        output_snp_counts[snp_index].append(snp_value)
    
    print sampled_strains
    
    output_snp_freqs = {k: float(sum(v))/len(v) for k,v in 
                       output_snp_counts.items()}
    
#     for k, v in output_snp_counts.items(): print k, v[:10]
#     for k, v in output_snp_freqs.items(): print k, v
    return output_snp_counts, output_strains

def get_freqs_from_strains_and_counts(input_snp_counts, input_strains):
 #   print input_snp_counts
    snp_freqs = []
    for k in range(len(input_snp_counts)):
        raw_snps = input_snp_counts[k]
        snp_freq = float(sum(raw_snps))/len(raw_snps)
        snp_freqs.append(snp_freq)
    
    snp_freqs = np.array(snp_freqs)
    #print snp_freqs
    
    strain_matrix = np.array(input_strains)
    #print strain_matrix
    
    strain_freqs = np.linalg.lstsq(strain_matrix.T, snp_freqs)
    #for freq, strain in zip(strain_freqs[0], input_strains):
    #    print freq, strain
    return strain_freqs[0]

# builds index for a strain, keys are length 10
def build_index(strain):
    fragments = 5
    k = 10
    dict = {}
    for i in range(len(strain)-10):
        kmer = strain[i:i+10]
        if kmer in dict:
            dict[kmer].append(i)
        else:
            dict[kmer] = [i]
    return dict

# align a read + collect SNPs 
def align(read,strain,strain_index):
    kmers = [read[:10]] # split read into 10-mers
    strain_length = len(strain)
    SNPs = [] 
    for kmer_number in range(len(kmers)): 
        key = kmers[kmer_number]
        if strain_index.has_key(key): 
            for pos in strain_index[key]:
                possible_SNPs = []
                index = read.find(key) # position where 10-mer starts in read
                mismatches = 0
                if (pos-index<0 or pos-index+len(read)>strain_length): # out of bounds check
                    continue
                for i in range(len(read)): # align read to strain
                    if (read[i]!=strain[pos-index+i]):
                        snp_position = pos-index+i # SNP position in strain
                        read_base = read[i] # read base
                        strain_base = strain[pos-index+i] # strain base
                        snp = (snp_position,read_base,strain_base) # (position,read,strain)
                        possible_SNPs.append(snp)
                        mismatches += 1
                if mismatches > 2:
                    continue
                else:
                    if len(possible_SNPs) > 0:
                        for snp in possible_SNPs:
                            SNPs.append(snp)
    return SNPs

# calculate snp frequencies 
def calculate_snp_freq(input_snp_counts):
    snp_freqs = []
    for k in range(len(input_snp_counts)):
        raw_snps = input_snp_counts[k]
        snp_freq = float(sum(raw_snps))/len(raw_snps)
        snp_freqs.append(snp_freq)
    return snp_freqs

# print list
def print_list(list1):
    for item in list1:
        print item
