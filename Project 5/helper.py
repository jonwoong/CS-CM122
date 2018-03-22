from collections import defaultdict

# builds index for a reference, keys are length 10
def build_index(reference):
    fragments = 5
    k = 10
    dict = {}
    for i in range(len(reference)-10):
        kmer = reference[i:i+10]
        if kmer in dict:
            dict[kmer].append(i)
        else:
            dict[kmer] = [i]
    return dict

# align a read + collect SNPs 
def align(read,reference,hash_table):
    kmers = [read[:10]] # split read into 10-mers
    reference_length = len(reference)
    matched_positions = []
    for kmer_number in range(len(kmers)): 
        key = kmers[kmer_number]
        if hash_table.has_key(key): 
            for pos in hash_table[key]:
                index = read.find(key) # position where 10-mer starts in read
                if (pos-index<0 or pos-index+len(read)>reference_length): # out of bounds check
                    continue

                mismatches = 0
                for i in range(len(read)): # align read to strain
                    if (read[i]!=reference[pos-index+i]):
                        snp_position = pos-index+i # SNP position in reference
                        read_base = read[i] # read base
                        reference_base = reference[pos-index+i] # reference base
                        mismatches += 1
                if mismatches > 4:
                    continue
                else:
                    matched_positions.append(pos)
    return matched_positions