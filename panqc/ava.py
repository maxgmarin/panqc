# pqgc/ava.py

from .kmerlib import all_vs_all_kmer_MaxJC, read_kmers_from_file_ToHashesDict
import time 

#import logging
# Set the logging level to INFO
#logging.basicConfig(level=logging.INFO)

def ava(input_PG_Ref_FA, kmer_size):

    ## Parse and hash all k-mers for each representative nucleotide sequence
    print("Beginning parsing of input FASTA")

    start = time.time()

    Ref_DictOf_Hashes, Ref_DictOf_SeqLen = read_kmers_from_file_ToHashesDict(input_PG_Ref_FA, kmer_size)             

    All_SeqIDs = list(Ref_DictOf_Hashes.keys())

    end = time.time()
    time_diff = end - start
    print(f"Time to parse and hash all k-mers: {round(time_diff, 2)} seconds")


    ## Calculate the maximum Jaccard Containment (JC) between all pairs of sequences.
    ### NOTE: The maximum JC between sets a and b will always be symetrical, while JC is not
    print(f"Beginning all vs all comparison of k-mer profiles:")

    start = time.time()

    PG_AvA_DF = all_vs_all_kmer_MaxJC(All_SeqIDs, Ref_DictOf_Hashes, Ref_DictOf_SeqLen) 

    end = time.time()
    time_diff = end - start
    print(f"Time for all vs all comparison of k-mer profiles: {round(time_diff, 2)} seconds")
    
    return PG_AvA_DF
