# pgqc/module1.py


import pandas as pd
import numpy as np

import time

import screed
import mmh3


##### k-mer functions #####

def build_kmers(sequence, ksize):
    kmers = []
    n_kmers = len(sequence) - ksize + 1
    
    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        kmers.append(kmer)
        
    return kmers


def read_kmers_from_file(filename, ksize):
    all_kmers = []
    for record in screed.open(filename):
        sequence = record.sequence
        
        kmers = build_kmers(sequence, ksize)
        all_kmers += kmers

    return all_kmers


def hash_kmer(kmer):
    # calculate the reverse complement
    rc_kmer = screed.rc(kmer)
    
    # determine whether original k-mer or reverse complement is lesser
    if kmer < rc_kmer:
        canonical_kmer = kmer
    else:
        canonical_kmer = rc_kmer
        
    # calculate murmurhash using a hash seed of 42
    hash = mmh3.hash64(canonical_kmer, 42)[0]
    if hash < 0: hash += 2**64
        
    # done
    return hash

def hash_kmers_ToSet(kmers):
    hashes = set()
    for kmer in kmers:
        hashes.add(hash_kmer(kmer))
    return hashes



def read_kmers_from_file_ToHashesDict(filename, ksize):

    all_hashes_Set_Dict = {}
    seqLen_Dict = {}
    
    NumParsedRecords = 0
    
    for record in screed.open(filename):
        
        ShortName = record.name.split(" ")[-1]

        NumParsedRecords += 1
        sequence = record.sequence

        kmers = build_kmers(sequence, ksize)
        hashes_Set = hash_kmers_ToSet(kmers)
        
        all_hashes_Set_Dict[ShortName] = hashes_Set
        seqLen_Dict[ShortName] = len(sequence)

    print(NumParsedRecords, " total records were parsed")
    
    return all_hashes_Set_Dict, seqLen_Dict



def jaccard_containment_FromSets(a, b):
    '''
    This function returns the Jaccard Containment between sets a and b.
    '''
    
    intersection = len(a.intersection(b))
    
    return intersection / len(a)

def jaccard_similarity_FromSets(a, b):
    '''
    This function returns the Jaccard Similarity between sets a and b.
    '''
    intersection = len(a.intersection(b))
    union = len(a.union(b))
    
    return intersection / union

def jaccard_containment_MaxVal_FromSets(a, b):
    '''
    This function returns the maximum possible Jaccard Containment between sets a and b.
    '''
    
    intersection = len(a.intersection(b))

    min_Len = min(len(a), len(b) )

    return intersection / min_Len





######## Define functions for performing the all vs all comparison of k-mer profiles ########

def all_vs_all_kmer_JC(all_SeqIDs, dictOf_Hashes_Set, dictOf_SeqLen):

    listOfTuples = []
    
    for record_Name_1 in all_SeqIDs:
        for record_Name_2 in all_SeqIDs: 
            if record_Name_1 != record_Name_2:
        
                record_1and2_JC = jaccard_containment_FromSets(dictOf_Hashes_Set[record_Name_1],
                                                                  dictOf_Hashes_Set[record_Name_2] )


                if record_1and2_JC != 0:
                    seqlen_1 = dictOf_SeqLen[record_Name_1]
                    seqlen_2 = dictOf_SeqLen[record_Name_2]

                    listOfTuples.append((record_Name_1, record_Name_2, seqlen_1, seqlen_2, record_1and2_JC) )
    if listOfTuples != []:
        PG_AvA_DF = pd.DataFrame(listOfTuples)
        PG_AvA_DF.columns = ["RecordID_1", "RecordID_2", "Record1_Len", "Record2_Len", "JaccardContain"]
        PG_AvA_DF = PG_AvA_DF.sort_values(["JaccardContain", "RecordID_1", "RecordID_2"], ascending=False)
    else:
        PG_AvA_DF = pd.DataFrame(columns=["RecordID_1", "RecordID_2", "Record1_Len", "Record2_Len", "JaccardContain"])

    return PG_AvA_DF

def all_vs_all_kmer_MaxJC(all_SeqIDs, dictOf_Hashes_Set, dictOf_SeqLen):

    listOfTuples = []
    
    for i, record_Name_1 in enumerate(all_SeqIDs) :
        for j, record_Name_2 in enumerate(all_SeqIDs) : 
            if i < j: # Check the seqID index so that the same pair of sequences is not compared twice

                record_1and2_JC = jaccard_containment_MaxVal_FromSets(dictOf_Hashes_Set[record_Name_1],
                                                                      dictOf_Hashes_Set[record_Name_2] )

                if record_1and2_JC != 0:
                    seqlen_1 = dictOf_SeqLen[record_Name_1]
                    seqlen_2 = dictOf_SeqLen[record_Name_2]    

                    listOfTuples.append((record_Name_1, record_Name_2, seqlen_1, seqlen_2, record_1and2_JC) )
    
    if listOfTuples != []:
        PG_AvA_DF = pd.DataFrame(listOfTuples)
        PG_AvA_DF.columns = ["RecordID_1", "RecordID_2", "Record1_Len", "Record2_Len", "MaxJaccardContain"]
        PG_AvA_DF = PG_AvA_DF.sort_values(["MaxJaccardContain", "RecordID_1", "RecordID_2"], ascending=False)
    else:
        PG_AvA_DF = pd.DataFrame(columns=["RecordID_1", "RecordID_2", "Record1_Len", "Record2_Len", "MaxJaccardContain"])

    return PG_AvA_DF


########################################################################################################################



