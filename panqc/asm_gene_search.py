# pgqc/asm_gene_search.py

import pandas as pd
import tqdm as tqdm

import mappy as mp
from tqdm import tqdm

#import logging
# Set the logging level to INFO
#logging.basicConfig(level=logging.INFO)

from pgqc.utils import parse_PresAbs_CSV_General, get_columns_excluding, parse_PG_Ref_FA


#### Define function for parsing Mappy alignment hits

def parse_AlnHits_To_DF(i_AsmAlner, QuerySeq):

    hits = list(i_AsmAlner.map(QuerySeq, cs = True))
    
    i_QueryLen = len(QuerySeq)
    
    listOfAlnRows = []

    if hits != []:
        
        for hit in hits:
        
            SeqID = hit.mlen / hit.blen
            LenQueryAligned = hit.q_en - hit.q_st
            
            QueryCov = LenQueryAligned / i_QueryLen
            
            i_AlnRow = [hit.ctg, hit.r_st, hit.r_en, i_QueryLen, LenQueryAligned,
                        QueryCov,  hit.mlen, hit.blen, SeqID, hit.is_primary]
            listOfAlnRows.append(i_AlnRow)
            
    Aln_DF = pd.DataFrame(listOfAlnRows, columns = ["RefSeq", "Ref_Start", "Ref_End", "QueryLen", "LenQueryAligned", "QueryCoverage", "AlnMatches", "AlnLen", "AlnSeqID", "IsPriAln"])

    return Aln_DF

######################################################








# General function for searching for gene DNA seqs in each assembly

def PresAbsQC_CheckAsmForGeneSeq(i_Gene_PresAbs_DF, i_PG_Ref_NucSeqs,
                                 i_AsmFA_Dict, i_SampleIDs,
                                 MinQueryCov = 0.9, MinQuerySeqID = 0.9):
    """
    This function takes in i) a gene presence/absence dataframe,
                           ii) a dictionary of protein-coding gene reference sequences,
                           iii) a dictionary of sample-specific genome assemblies,
                           iv) a list of sample IDs,
                           and optional parameters
                           v) for minimum query coverage
                           vi) and minimum query sequence identity. 
                           
    It searches the genome assemblies for each absent gene in each sample,
    and updates the gene presence/absence dataframe accordingly. 
    If a gene sequence is found with high confidence, the gene is marked as "not present" at the protein level but "present" at the DNA level.
    The function returns the updated gene presence/absence dataframe.

    Args:
    - i_Gene_PresAbs_DF: pandas DataFrame containing gene presence/absence information for each sample. The DataFrame should have the following columns: "Gene" (gene name), "NumTotalGenomes" (total number of genomes), and one column for each sample ID containing 0 (absent) or 1 (present) to indicate gene presence/absence.
    - i_PG_Ref_NucSeqs: dictionary of protein-coding gene reference sequences. The keys are gene names and the values are nucleotide sequences.
    - i_AsmFA_Dict: dictionary of sample-specific genome assemblies. The keys are sample IDs and the values are file paths to the genome assembly FASTA files.
    - i_SampleIDs: list of sample IDs to process.
    - MinQueryCov: optional float specifying the minimum query coverage required for a gene sequence alignment to be considered a hit. Default is 0.9.
    - MinQuerySeqID: optional float specifying the minimum query sequence identity required for a gene sequence alignment to be considered a hit. Default is 0.9.

    Returns:
    - i_Gene_PresAbs_DF_Updated: pandas DataFrame containing updated gene presence/absence information for each sample. The DataFrame has the same columns as i_Gene_PresAbs_DF, plus an additional column "NumAsm_WiGene_DNASeq" containing the number of samples in which the gene sequence was found with high confidence.

    """

    # Assert that sequence IDs (gene names) are identical between the gene presence CSV & the pan-enome reference fasta
    PG_Ref_SeqIDs_Set = set(sorted(list(i_PG_Ref_NucSeqs.keys())))
    PresMatrix_SeqIDs_Set = set(list(i_Gene_PresAbs_DF.index))
    
    ErrorMessage = "ERROR: The geneIDs in the gene presence/absence matrix does not match the provided nucleotide reference FASTA."
    assert PG_Ref_SeqIDs_Set == PresMatrix_SeqIDs_Set, ErrorMessage

    i_Gene_PresAbs_DF_Updated = i_Gene_PresAbs_DF.copy().set_index("Gene", drop=False)
    
    for i_SampleID in tqdm(i_SampleIDs) :

        i_Alner_Asm = mp.Aligner(i_AsmFA_Dict[i_SampleID], preset="asm10")  # load or build index
        if not i_Alner_Asm: raise Exception(f"ERROR: failed to load/build index for SR Asm - {i_SampleID}")
        
        i_SampleOnly_GenePres = i_Gene_PresAbs_DF.set_index("Gene")[i_SampleID]
        
        i_AbsentGenes = list( i_SampleOnly_GenePres[i_SampleOnly_GenePres == 0].index)
        
        for gene in i_AbsentGenes:
        
            gene_seq = i_PG_Ref_NucSeqs[gene]
        
            i_AlnToAsm_DF = parse_AlnHits_To_DF(i_Alner_Asm, gene_seq)

            i_AlnToAsm_Filt = i_AlnToAsm_DF[(i_AlnToAsm_DF["QueryCoverage"] > MinQueryCov) & (i_AlnToAsm_DF["AlnSeqID"] > MinQuerySeqID)]
            
            Num_Hits_Pass = i_AlnToAsm_Filt.shape[0]
    
            if Num_Hits_Pass > 0: # Update Pres/Abs matrix to have 2, meaning that the Protein-level annotation is Not present, BUT the gene sequence can be found with high confidence.
                i_Gene_PresAbs_DF_Updated.loc[gene, i_SampleID] = 2

    i_Gene_PresAbs_DF_Updated["NumAsm_WiGene_DNASeq"] = i_Gene_PresAbs_DF_Updated[i_SampleIDs].applymap(lambda x: 1 if x > 0 else 0).sum(axis = 1)

    return i_Gene_PresAbs_DF_Updated


def create_AsmFA_PATH_Dict(i_AsmFA_TSV):
    AsmFA_DF = pd.read_csv(i_AsmFA_TSV, sep="\t")
    AsmFA_Dict = dict(AsmFA_DF[['SampleID', 'Genome_ASM_PATH']].values)
    return AsmFA_Dict

def create_LRandSR_AsmFA_PATH_Dict(i_AsmFA_LRandSR_TSV):
    AsmFA_LRandSR_DF = pd.read_csv(i_AsmFA_LRandSR_TSV, sep="\t")

    LR_AsmFA_Dict = dict(AsmFA_LRandSR_DF[['SampleID', 'Genome_LR_ASM_PATH']].values)
    SR_AsmFA_Dict = dict(AsmFA_LRandSR_DF[['SampleID', 'Genome_SR_ASM_PATH']].values)

    return LR_AsmFA_Dict, SR_AsmFA_Dict

PresAbs_NonSampleID_ColNames = ['Gene', 'NumAsm_WiGene', 'NumAsm_WiGene_DNASeq',
                                'Non-unique Gene name', 'Annotation', 'No. isolates',
                                'No. sequences', 'Avg sequences per isolate', 'Genome Fragment',
                                'Order within Fragment', 'Accessory Fragment', 'Accessory Order with Fragment', 'QC',
                                'Min group size nuc', 'Max group size nuc', 'Avg group size nuc']

def asmseqcheck_frompaths(i_Gene_PresAbs_CSV_PATH,
                            i_PG_Ref_FA_PATH,
                            i_AsmFA_TSV,
                            MinQueryCov, MinQuerySeqID):
    """
    Searches for gene sequences in genome assemblies and returns a DataFrame with presence/absence information.

    Args:
        i_Gene_PresAbs_CSV_PATH (str): Path to the gene presence/absence matrix CSV file.
        i_PG_Ref_FA_PATH (str): Path to the pan-genome (nucleotide) reference fasta file.
        i_AsmFA_TSV (str): Path to the TSV file containing the genome assembly fasta file paths.
        MinQueryCov (float): Minimum query coverage required for a gene sequence to be considered present.
        MinQuerySeqID (float): Minimum query sequence identity required for a gene sequence to be considered present.

    Returns:
        pandas.DataFrame: A DataFrame with presence/absence information for each gene in each sample.
    """

    # 1) Parse in Gene Pres/Abs matrix and select sampleIDs

    print("Parsing input gene presence/absence matrix...")
    Gene_PresAbs_DF = parse_PresAbs_CSV_General(i_Gene_PresAbs_CSV_PATH)
    
    i_SampleIDs = get_columns_excluding(Gene_PresAbs_DF, PresAbs_NonSampleID_ColNames)

    print("Finished parsing gene presence/absence matrix\n")

    # 2) Read in pan-genome (nucleotide) reference fasta
    print("Parsing input gene reference fasta")
    PG_Ref_NucSeqs = parse_PG_Ref_FA(i_PG_Ref_FA_PATH)
    print("Finished parsing gene reference fasta\n")

    # 3) Read in genome assembly fasta PATH TSV, convert to dict
    print("Parsing TSV of input assembly PATHs")
    AsmFA_Dict = create_AsmFA_PATH_Dict(i_AsmFA_TSV)
    print("Finished parsing input assembly PATHs \n")
    
    # 4) Run alignments to check for gene sequences in assemblies
    print("Running alignments to check for gene sequences in assemblies...")

    Gene_PresAbs_WiAsmSeqCheck_DF = PresAbsQC_CheckAsmForGeneSeq(Gene_PresAbs_DF, PG_Ref_NucSeqs,
                                                                AsmFA_Dict, i_SampleIDs,
                                                                MinQueryCov, MinQuerySeqID)

    print("Finished searching for gene sequences in assemblies \n")


    return Gene_PresAbs_WiAsmSeqCheck_DF





def SRandLRQC_asmseqcheck_frompaths(i_SR_Gene_PresAbs_CSV_PATH,
                                    i_SR_PG_Ref_FA_PATH,
                                    i_AsmFA_LRandSR_TSV,
                                    MinQueryCov, MinQuerySeqID):

    # 1) Parse in Gene Pres/Abs matrix and select sampleIDs

    print("Parsing input gene presence/absence matrix...")
    SR_Gene_PresAbs_DF = parse_PresAbs_CSV_General(i_SR_Gene_PresAbs_CSV_PATH)
    
    i_SampleIDs = get_columns_excluding(SR_Gene_PresAbs_DF, PresAbs_NonSampleID_ColNames)

    print("Finished parsing gene presence/absence matrix\n")

    # 2) Read in pan-genome (nucleotide) reference fasta (For SR based analysis)
    print("Parsing input gene reference fasta")
    SR_PG_Ref_NucSeqs = parse_PG_Ref_FA(i_SR_PG_Ref_FA_PATH)
    print("Finished parsing gene reference fasta\n")

    # 3) Read in LR and SR assemblies PATH TSV, convert to dictionaries (For LR and SR based analysis)
    print("Parsing TSV of input assembly PATHs")
    LR_AsmFA_Dict, SR_AsmFA_Dict = create_LRandSR_AsmFA_PATH_Dict(i_AsmFA_LRandSR_TSV)
    print("Finished parsing input assembly PATHs \n")

    # 4) Run alignments to check for gene sequences in assemblies
    print("Running alignments to check for gene sequences in assemblies...")

    Gene_PresAbs_WiAsmSeqCheck_DF = SRAsm_PresAbsQC_CheckInLRAsm(SR_Gene_PresAbs_DF, SR_PG_Ref_NucSeqs,
                                                                 LR_AsmFA_Dict, SR_AsmFA_Dict, i_SampleIDs,
                                                                 MinQueryCov, MinQuerySeqID)

    print("Finished searching for gene sequences in both LR and SR assemblies \n")

    return Gene_PresAbs_WiAsmSeqCheck_DF





def SRAsm_PresAbsQC_CheckInLRAsm(i_SR_Gene_PresAbs_DF, i_SR_PG_Ref_NucSeqs,
                                 i_SR_AsmFA_Dict, i_LR_AsmFA_Dict, i_SampleIDs,
                                 MinQueryCov = 0.9, MinQuerySeqID = 0.9):
    """
    Check the presence/absence of genes in short-read (SR) and long-read (LR) genome assemblies,
    and update a gene presence/absence matrix accordingly.

    Args:
    - i_SR_Gene_PresAbs_DF: pandas DataFrame containing the gene presence/absence matrix for the short-read assembly analysis.
    - i_SR_PG_Ref_NucSeqs: dictionary containing the reference nucleotide sequences for each gene.
    - i_SR_AsmFA_Dict: dictionary containing the file paths of the SR genome assemblies for each sample.
    - i_LR_AsmFA_Dict: dictionary containing the file paths of the LR genome assemblies for each sample.
    - i_SampleIDs: list of sample IDs to process.
    - MinQueryCov: minimum query coverage required for a sequence alignment to be considered valid (default: 0.9).
    - MinQuerySeqID: minimum sequence identity required for a sequence alignment to be considered valid (default: 0.9).

    Returns:
    - i_SR_Gene_PresAbs_DF_Updated: pandas DataFrame containing the updated gene presence/absence matrix for each sample,
    with additional column indicating the number of genome assemblies where the gene was found.

    Note on output values in Presence/Absence matrix:
    - 0 means CDS AA seq & DNA seq not found in SR
    - 1 means CDS AA seq present in SR
    - 3 Means DNA Seq "Not in SR, but In LR Asm", ("Incomplete Assembly")
    - 4 Means DNA Seq "In SR Asm, NOT In LR Asm" 
    - 5 Means DNA Seq "In SR Asm, In LR Asm" ("Annotation Discrepancy")

    """

    TotalNum_All_SR_MissingGenes = 0
    TotalNum_All_Abs_InSRandLR = 0
    TotalNum_All_Abs_InSR_NotInLR = 0
    
    i_SR_Gene_PresAbs_DF_Updated = i_SR_Gene_PresAbs_DF.copy().set_index("Gene", drop=False)
    
    PresAbs_Updates = []
    for i_SampleID in tqdm(i_SampleIDs) :
    
        i_SRAsm_GenePres = i_SR_Gene_PresAbs_DF.set_index("Gene")[i_SampleID]
        
        i_Alner_SRAsm = mp.Aligner(i_SR_AsmFA_Dict[i_SampleID], preset="asm10")  # load or build index
        if not i_Alner_SRAsm: raise Exception(f"ERROR: failed to load/build index for SR Asm - {i_SampleID}")
        
        i_Alner_LRAsm = mp.Aligner(i_LR_AsmFA_Dict[i_SampleID], preset="asm10")  # load or build index
        if not i_Alner_LRAsm: raise Exception(f"ERROR: failed to load/build index for LR Asm - {i_SampleID}")
        
        i_SRAsm_AbsentGenes = list( i_SRAsm_GenePres[i_SRAsm_GenePres == 0].index)
        
        NumSRAsm_NotInSR_InLR = 0
        NumAbsent_ButInSRandLR = 0
        NumAbsent_InSR_NotInLR = 0

        
        for gene in i_SRAsm_AbsentGenes:
        
            gene_seq = i_SR_PG_Ref_NucSeqs[gene]
        
            i_AlnToSRAsm_DF = parse_AlnHits_To_DF(i_Alner_SRAsm, gene_seq)
            i_AlnToLRAsm_DF = parse_AlnHits_To_DF(i_Alner_LRAsm, gene_seq)
        
            #i_AlnToSRAsm_Filt = i_AlnToSRAsm_DF[(i_AlnToSRAsm_DF["QueryCoverage"] > MinQueryCov) & (i_AlnToSRAsm_DF["AlnSeqID"] > MinQuerySeqID)]
            #i_AlnToLRAsm_Filt = i_AlnToLRAsm_DF[(i_AlnToLRAsm_DF["QueryCoverage"] > MinQueryCov) & (i_AlnToLRAsm_DF["AlnSeqID"] > MinQuerySeqID)]

            Num_SR_Hits = i_AlnToSRAsm_DF[(i_AlnToSRAsm_DF["QueryCoverage"] > MinQueryCov) & (i_AlnToSRAsm_DF["AlnSeqID"] > MinQuerySeqID)].shape[0]
            Num_LR_Hits = i_AlnToLRAsm_DF[(i_AlnToLRAsm_DF["QueryCoverage"] > MinQueryCov) & (i_AlnToLRAsm_DF["AlnSeqID"] > MinQuerySeqID)].shape[0]
    
            #if Num_SR_Hits > 0: # Update Pres/Abs matrix to have 2, meaning that the Protein-level annotation is Not present, BUT the gene sequence can be found with high confidence.
            #    i_SR_Gene_PresAbs_DF_Updated.loc[gene, i_SampleID] = 2
    
            
            if (Num_SR_Hits == 0) & (Num_LR_Hits > 0):
                #i_SR_Gene_PresAbs_DF_Updated.loc[gene, i_SampleID] = 3 # Means Not in SR, In LR Asm

                PresAbs_Updates.append((gene, i_SampleID, 3))
                
                NumSRAsm_NotInSR_InLR += 1
                
            if (Num_SR_Hits > 0) & (Num_LR_Hits == 0):  
                #i_SR_Gene_PresAbs_DF_Updated.loc[gene, i_SampleID] = 4 # Means In SR Asm, NOT In LR Asm
                PresAbs_Updates.append((gene, i_SampleID, 4))
                
                NumAbsent_InSR_NotInLR += 1
    
            if (Num_SR_Hits > 0) & (Num_LR_Hits > 0):
                #i_SR_Gene_PresAbs_DF_Updated.loc[gene, i_SampleID] = 5 # Means In SR Asm, In LR Asm
                PresAbs_Updates.append((gene, i_SampleID, 5))
                
                NumAbsent_ButInSRandLR += 1
            
        TotalNum_All_SR_MissingGenes += NumSRAsm_NotInSR_InLR
    
        TotalNum_All_Abs_InSR_NotInLR += NumAbsent_InSR_NotInLR
        
        TotalNum_All_Abs_InSRandLR += NumAbsent_ButInSRandLR
        
        #print("------------")

    # Apply all updates
    for gene, sample, value in PresAbs_Updates:
        i_SR_Gene_PresAbs_DF_Updated.loc[gene, sample] = value
            
    print("Across all samples, total missing genes - Not in SR Asm, but in LR Asm:", TotalNum_All_SR_MissingGenes)
    print("Across all samples, total missing genes - In SR Asm, Not in LR Asm:", TotalNum_All_Abs_InSR_NotInLR)
    print("Across all samples, total missing genes - In BOTH SR Asm and LR Asm:", TotalNum_All_Abs_InSRandLR)
    
    i_SR_Gene_PresAbs_DF_Updated["NumAsm_WiGene_DNASeq"] = i_SR_Gene_PresAbs_DF_Updated[i_SampleIDs].applymap(lambda x: 1 if x in [1, 4, 5] else 0).sum(axis = 1)
    
    return i_SR_Gene_PresAbs_DF_Updated






def get_SRAsm_Vs_LRAsm_QCStats(i_Gene_PresAbs_DF, i_SampleIDs, print_stats = True):

    N_AbsentCDS = (i_Gene_PresAbs_DF[i_SampleIDs] == 0).sum().sum()  
    N_PresentCDS = (i_Gene_PresAbs_DF[i_SampleIDs] == 1).sum().sum()  
    
    N_AbsentCDS_DNASeq_NotInSR_InLR = (i_Gene_PresAbs_DF[i_SampleIDs] == 3).sum().sum()  
    
    N_AbsentCDS_DNASeq_InSR_NotInLR = (i_Gene_PresAbs_DF[i_SampleIDs] == 4).sum().sum()  
    
    N_AbsentCDS_DNASeq_InSR_InLR = (i_Gene_PresAbs_DF[i_SampleIDs] == 5).sum().sum()  

    N_AnyType_AbsentCDS = (i_Gene_PresAbs_DF[i_SampleIDs] != 1).sum().sum()  

    
    if print_stats:
        print("# of absent genes (CDS and DNA level):", N_AbsentCDS, round(N_AbsentCDS/ N_AnyType_AbsentCDS, 4))
        print("# of present genes (CDS):", N_PresentCDS)
    
        print("# of absent genes (DNA only found in LR Asm):", 
              N_AbsentCDS_DNASeq_NotInSR_InLR,
              round(N_AbsentCDS_DNASeq_NotInSR_InLR/ N_AnyType_AbsentCDS, 4))
    
        print("# of absent genes (DNA only found in SR Asm):", N_AbsentCDS_DNASeq_InSR_NotInLR, round(N_AbsentCDS_DNASeq_InSR_NotInLR/ N_AnyType_AbsentCDS, 4))
    
        print("# of absent genes (DNA only found in BOTH SR & LR Asm):",
              N_AbsentCDS_DNASeq_InSR_InLR,
              round(N_AbsentCDS_DNASeq_InSR_InLR / N_AnyType_AbsentCDS, 4))

        print("# of TOTAL absent genes (At Any Level):", N_AnyType_AbsentCDS)

    
    return N_AbsentCDS, N_PresentCDS, N_AbsentCDS_DNASeq_NotInSR_InLR, N_AbsentCDS_DNASeq_InSR_NotInLR, N_AbsentCDS_DNASeq_InSR_InLR



def get_AsmSeqCheck_QCStats(i_Gene_PresAbs_DF, print_stats = True): # i_SampleIDs

    i_SampleIDs = get_columns_excluding(i_Gene_PresAbs_DF, PresAbs_NonSampleID_ColNames)

    N_AbsentCDS = (i_Gene_PresAbs_DF[i_SampleIDs] == 0).sum().sum()  
    N_PresentCDS = (i_Gene_PresAbs_DF[i_SampleIDs] == 1).sum().sum()  
    
    N_AbsentCDS_DNASeq_InAsm = (i_Gene_PresAbs_DF[i_SampleIDs] == 2).sum().sum()  

    N_AnyType_AbsentCDS = (i_Gene_PresAbs_DF[i_SampleIDs] != 1).sum().sum()  

    if print_stats:
        print("# of absent genes (CDS and DNA level):", N_AbsentCDS, round(N_AbsentCDS/ N_AnyType_AbsentCDS, 4))
        print("# of present genes (CDS):", N_PresentCDS)
    
        print("# of absent genes (but DNA seq found in Asm):", 
              N_AbsentCDS_DNASeq_InAsm,
              round(N_AbsentCDS_DNASeq_InAsm / N_AnyType_AbsentCDS, 4))
    
        print("# of TOTAL absent genes (At Any Level):", N_AnyType_AbsentCDS)

    
    return N_AbsentCDS, N_PresentCDS, N_AbsentCDS_DNASeq_InAsm, N_AnyType_AbsentCDS






































