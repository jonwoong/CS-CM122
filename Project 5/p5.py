########## HEADERS ##########

import numpy as np
import pandas as pd
from helper import *

########## MAIN ##########

if __name__ == "__main__":

    ##### I/O SETUP #####

    input_folder = 'hw5_W_0'
    input_fn_start = '{0}/{0}'.format(input_folder)
    exons_fn = '{}_exons.txt'.format(input_fn_start)
    isoforms_fn = '{}_isoforms.txt'.format(input_fn_start)
    exon_counts_fn = '{}_exon_counts.txt'.format(input_fn_start)
    reads_fn = '{}_reads.txt'.format(input_fn_start)
    reference_fn = 'ref_hw5.txt'

    ##### READ DATA #####

    # store refrence into a string reference
    with open(reference_fn) as reference_file:
        reference_file_lines = reference_file.read().splitlines()
    reference_file.close()

    reference = ""
    for line in range(1,len(reference_file_lines)):
        reference += reference_file_lines[line]

    # calculate length of reference
    length_reference = len(reference)

    # store reads into an array reads[]
    with open(reads_fn) as reads_file:
        reads_file_lines = reads_file.read().splitlines()
    reads_file.close()

    reads = []
    for read in range(1,len(reads_file_lines)):
        reads.append(reads_file_lines[read])

    # calculate total length of all reads
    total_length_reads = len(reads) * len(reads[0])

    # calculate coverage
    coverage = total_length_reads/float(length_reference)
    #print coverage

    ##### ALIGN READS #####

    # build hash table for reference
    hash_table = build_index(reference)

    # count of how many reads matched in each position of the reference
    matched_reference = [0 for i in range(len(reference))]

    # align reads to reference + collect possible snps
    for read in reads:
        match_positions = align(read,reference,hash_table)
        for position in match_positions: # +1 for every position the read matched to
            for i in range(len(read)):
                matched_reference[position+i] += 1

    ##### GET EXON READ COUNTS #####

    # get the length and count of each exon
    exon_length_dict = {}
    exon_read_count_dict = {}
    with open(exons_fn) as exons_file:
        exons_file.readline()
        for line in exons_file:
            exon, start_end = line.strip().split(':')
            start, end = (int(x) for x in start_end.split(','))
            exon_length = end - start

            ##########
            #
            read_count = max(matched_reference[start:end]) # the number of times the exon appears after aligning all reads
            exon_read_count_dict[exon] = read_count
            #print exon_read_count_dict
            #
            ##########

            exon_length_dict[exon] = exon_length
            #print exon_length_dict
    
    exon_length_df = pd.DataFrame.from_dict(exon_length_dict, orient='index').reset_index().rename(columns={'index':'exon', 0:'length'})
    #print exon_length_df.head()
    exon_read_count_df = pd.DataFrame.from_dict(exon_read_count_dict, orient='index').reset_index().rename(columns={'index':'exon', 0:'read_count'})
    
    ##########
    #
    exon_complete_df = pd.merge(exon_length_df, exon_read_count_df, on='exon')
    #print exon_complete_df
    #
    ##########
    
    # get the exons of each isoform
    isoform_exons = []
    with open(isoforms_fn) as isoforms_file:
        isoforms_file.readline()
        for line in isoforms_file:
            isoform, exons = line.strip().split(':')
            exons = exons.split(',')
            gene_id = isoform.split('_ISO')[0]
            full_exons = ['{}_{}'.format(gene_id, exon) for exon in exons]
            these_exons = [(isoform, exon) for exon in full_exons]
            isoform_exons += these_exons
    #print isoform_exons

    isoforms_df = pd.DataFrame(isoform_exons).rename(columns={0:'isoform', 1:'exon'})
    #print isoforms_df.head()

    isoform_exon_df = pd.merge(isoforms_df, exon_complete_df, on='exon')

    ##########
    #
    isoform_exon_df['total_length'] = isoform_exon_df.apply(lambda row:(row['length']*row['read_count']),axis=1)
    #print isoform_exon_df.head()
    #
    ##########

    combined_df = pd.pivot_table(isoform_exon_df, columns='isoform', index='exon', values='total_length').fillna(0)

    isoform_columns = combined_df.columns
    combined_df.reset_index(inplace=True)
    #print combined_df.head()
    #print combined_df.columns

    exon_counts = []
    with open(exon_counts_fn) as exon_counts_file:
        exon_counts_file.readline()
        for line in exon_counts_file:
            exon, count = line.strip().split(':')
            exon_counts.append((exon, int(count)))
    exon_counts_df = pd.DataFrame(exon_counts).rename(columns={0:'exon', 1:'exon_count'})
    #print exon_counts_df.head()
    complete_df = pd.merge(combined_df, exon_counts_df, on='exon')
    #print complete_df#.head()

    exon_counts_array = complete_df.exon_count.values
    #print exon_counts_array.shape
    exon_isoform_matrix = complete_df.loc[:, isoform_columns].values

    ##########
    #
    # divide total exon length by length of read
    for row_index in range(len(exon_isoform_matrix)):
        for value_index in range(len(exon_isoform_matrix[row_index])):
            exon_isoform_matrix[row_index][value_index] /= float(len(read[0]))
    #
    ##########

    #print exon_isoform_matrix
    #print exon_isoform_matrix.shape

    implied_frequencies = np.linalg.lstsq(exon_isoform_matrix, exon_counts_array)[0]
    implied_frequencies /= implied_frequencies.sum()

    ##########
    #
    # remove negative frequencies
    for frequency_index in range(len(implied_frequencies)):
        if implied_frequencies[frequency_index] < 0:
            implied_frequencies[frequency_index]=0
    print implied_frequencies
    #
    ##########

    output_fn = '{}_output.txt'.format(input_fn_start)
    with open(output_fn, 'w') as output_file:
        output_file.write('>{}\n'.format(input_folder))
        for isoform, freq in sorted([_ for _ in zip(isoform_columns, implied_frequencies)],
                                    key=lambda x: (int(x[0].split('_')[1]), int(x[0].split('_')[3]))):
            output_file.write('{}:{}\n'.format(isoform, freq))