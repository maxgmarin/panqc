# pgqc/utils.py

import pandas as pd
import screed 



def parse_PresAbs_Rtab(PresAbs_Rtab_PATH):
    '''
    This function parsesthe `gene_presence_absence.csv` file output by Panaroo '''

    i_Gene_PresAbs_DF = pd.read_csv(PresAbs_Rtab_PATH, sep = "\t")

    ### Relabel Columns for presence/absence tracking
    #i_Gene_PresAbs_DF.columns = [ x.split(".Bakta")[0] for x in i_Gene_PresAbs_DF.columns ]


    ListOf_SampleID_Cols = list(i_Gene_PresAbs_DF.drop(["Gene"], axis=1).columns)
    
    i_Gene_PresAbs_DF["NumAsm_WiGene"] = i_Gene_PresAbs_DF[ListOf_SampleID_Cols].sum(axis = 1)

    i_Gene_PresAbs_DF = i_Gene_PresAbs_DF.sort_values(by='NumAsm_WiGene', ascending=False)
    i_Gene_PresAbs_DF = i_Gene_PresAbs_DF.set_index("Gene", drop=False)

    return i_Gene_PresAbs_DF





def get_columns_excluding(df, exclude_columns):
    """
    Get all column names from a dataframe excluding the defined columns.

    Parameters:
    df (pd.DataFrame): The input dataframe.
    exclude_columns (list): A list of column names to exclude.

    Returns:
    list: A list of column names excluding the defined columns.
    """
    return [col for col in df.columns if col not in exclude_columns]


PresAbs_NonSampleID_ColNames = ['Gene', 'NumAsm_WiGene', 'NumAsm_WiGene_DNASeq',
                                'Non-unique Gene name', 'Annotation', 'No. isolates',
                                'No. sequences', 'Avg sequences per isolate', 'Genome Fragment',
                                'Order within Fragment', 'Accessory Fragment', 'Accessory Order with Fragment', 'QC',
                                'Min group size nuc', 'Max group size nuc', 'Avg group size nuc']



def parse_PresAbs_CSV_General(PresAbs_CSV_PATH):
    '''
    This function parses the `gene_presence_absence.csv` file output by Panaroo OR Roary '''


    ColNames_ToRemove = ['Non-unique Gene name', 'Annotation', 'No. isolates',
                        'No. sequences', 'Avg sequences per isolate', 'Genome Fragment',
                        'Order within Fragment', 'Accessory Fragment', 'Accessory Order with Fragment', 'QC',
                        'Min group size nuc', 'Max group size nuc', 'Avg group size nuc']


    i_Gene_PresAbs_DF = pd.read_csv(PresAbs_CSV_PATH, low_memory=False)

    ### Relabel Columns for presence/absence tracking
    i_Gene_PresAbs_DF.columns = [ x.split(".Bakta")[0] for x in i_Gene_PresAbs_DF.columns ]

    ListOf_SampleID_Cols = get_columns_excluding(i_Gene_PresAbs_DF, PresAbs_NonSampleID_ColNames)
    #print(ListOf_SampleID_Cols)
    ColsToRemove = [col for col in i_Gene_PresAbs_DF.columns if col in ColNames_ToRemove]

    i_Gene_PresAbs_DF = i_Gene_PresAbs_DF.drop(ColsToRemove, axis=1)

    # https://stackoverflow.com/questions/12741092/pandas-dataframe-apply-function-to-all-columns
    i_Gene_PresAbs_DF[ListOf_SampleID_Cols] = i_Gene_PresAbs_DF[ListOf_SampleID_Cols].applymap(lambda x: 1 if isinstance(x, str) else 0)         
    i_Gene_PresAbs_DF["NumAsm_WiGene"] = i_Gene_PresAbs_DF[ListOf_SampleID_Cols].sum(axis = 1)

    i_Gene_PresAbs_DF = i_Gene_PresAbs_DF.sort_values(by='NumAsm_WiGene', ascending=False)
    #print(i_Gene_PresAbs_DF.head(1))

    i_Gene_PresAbs_DF = i_Gene_PresAbs_DF.set_index("Gene", drop=False)
    #print(i_Gene_PresAbs_DF.head(1))
    return i_Gene_PresAbs_DF  


def PresAbs_InferSampleColOnly(i_Gene_PresAbs_DF):
    PresAbs_NonSampleID_ColNames = ['Gene', 'NumAsm_WiGene', 'NumAsm_WiGene_DNASeq',
                                'Non-unique Gene name', 'Annotation', 'No. isolates',
                                'No. sequences', 'Avg sequences per isolate', 'Genome Fragment',
                                'Order within Fragment', 'Accessory Fragment', 'Accessory Order with Fragment', 'QC',
                                'Min group size nuc', 'Max group size nuc', 'Avg group size nuc']


    ListOf_SampleID_Cols = get_columns_excluding(i_Gene_PresAbs_DF, PresAbs_NonSampleID_ColNames)
    #print(ListOf_SampleID_Cols)
    ColsToRemove = [col for col in i_Gene_PresAbs_DF.columns if col in PresAbs_NonSampleID_ColNames]

    i_Gene_PresAbs_DF = i_Gene_PresAbs_DF.drop(ColsToRemove, axis=1)

    return i_Gene_PresAbs_DF, ListOf_SampleID_Cols






def parse_PresAbs_CSV_Panaroo(PresAbs_CSV_PATH):
    '''
    This function parses the `gene_presence_absence.csv` file output by Panaroo '''

    'Non-unique Gene name', 'Annotation',

    Panaroo_ExtraCol_ToRemove = ['Non-unique Gene name', 'Annotation']

    i_Gene_PresAbs_DF = pd.read_csv(PresAbs_CSV_PATH)

    ### Relabel Columns for presence/absence tracking
    i_Gene_PresAbs_DF.columns = [ x.split(".Bakta")[0] for x in i_Gene_PresAbs_DF.columns ]

    ListOf_SampleID_Cols = list(i_Gene_PresAbs_DF.columns[3:].values)
    
    i_Gene_PresAbs_DF = i_Gene_PresAbs_DF.drop(Panaroo_ExtraCol_ToRemove, axis=1)

    # https://stackoverflow.com/questions/12741092/pandas-dataframe-apply-function-to-all-columns
    i_Gene_PresAbs_DF[ListOf_SampleID_Cols] = i_Gene_PresAbs_DF[ListOf_SampleID_Cols].applymap(lambda x: 1 if isinstance(x, str) else 0)         
    i_Gene_PresAbs_DF["NumAsm_WiGene"] = i_Gene_PresAbs_DF[ListOf_SampleID_Cols].sum(axis = 1)

    i_Gene_PresAbs_DF = i_Gene_PresAbs_DF.sort_values(by='NumAsm_WiGene', ascending=False)
    i_Gene_PresAbs_DF = i_Gene_PresAbs_DF.set_index("Gene", drop=False)

    return i_Gene_PresAbs_DF



def parse_PresAbs_CSV_Roary(PresAbs_CSV_PATH):
    '''
    This function parsesthe `gene_presence_absence.csv` file output by Roary
    '''

    Roary_ExtraCol_ToRemove = ['Non-unique Gene name', 'Annotation', 'No. isolates',
                                'No. sequences', 'Avg sequences per isolate', 'Genome Fragment',
                                'Order within Fragment', 'Accessory Fragment', 'Accessory Order with Fragment', 'QC',
                                'Min group size nuc', 'Max group size nuc', 'Avg group size nuc']

    i_Gene_PresAbs_DF = pd.read_csv(PresAbs_CSV_PATH)

    ### Relabel Columns for presence/absence tracking
    i_Gene_PresAbs_DF.columns = [ x.split(".Bakta")[0] for x in i_Gene_PresAbs_DF.columns ]

    ListOf_SampleID_Cols = list(i_Gene_PresAbs_DF.columns[14:].values)
    i_Gene_PresAbs_DF = i_Gene_PresAbs_DF.drop(Roary_ExtraCol_ToRemove, axis=1)

    # https://stackoverflow.com/questions/12741092/pandas-dataframe-apply-function-to-all-columns
    i_Gene_PresAbs_DF[ListOf_SampleID_Cols] = i_Gene_PresAbs_DF[ListOf_SampleID_Cols].applymap(lambda x: 1 if isinstance(x, str) else 0)         
    i_Gene_PresAbs_DF["NumAsm_WiGene"] = i_Gene_PresAbs_DF[ListOf_SampleID_Cols].sum(axis = 1)

    i_Gene_PresAbs_DF = i_Gene_PresAbs_DF.sort_values(by='NumAsm_WiGene', ascending=False)
    
    i_Gene_PresAbs_DF = i_Gene_PresAbs_DF.set_index("Gene", drop=False)

    return i_Gene_PresAbs_DF




def get_PG_Stats_FromPresAbs(i_Gene_PresAbs_DF, NumSamples, Verbose=False):

    # Total number of assemblies
    
    total_assemblies = NumSamples
    
    # Defining the threshold for accessory genes (less than 99% of all assemblies)
    accessory_threshold = total_assemblies * 0.99
    print("Accessory Thresh:", accessory_threshold)
    
    Gene_PresAbs_Core_DF = i_Gene_PresAbs_DF.query(f"NumAsm_WiGene >= {accessory_threshold}")
    
    # Calculating the number of core genes (present in all assemblies)
    num_core_genes = Gene_PresAbs_Core_DF.shape[0]
    
    # Calculating the number of accessory genes (present in less than 99% of all assemblies)
    
    Gene_PresAbs_Acc_DF = i_Gene_PresAbs_DF.query(f"NumAsm_WiGene < {accessory_threshold}")
    
    num_accessory_genes = Gene_PresAbs_Acc_DF.shape[0]

    num_total_genes = i_Gene_PresAbs_DF.shape[0]
    if Verbose:
        print("# of core genes:", num_core_genes)
        print("# of accessory genes:", num_accessory_genes)

    return num_total_genes, num_core_genes, num_accessory_genes


def get_PG_Stats_FromDNASeqPresAbs(i_Gene_PresAbs_DF, NumSamples, Verbose=False):

    # Total number of assemblies
    
    total_assemblies = NumSamples
    
    # Defining the threshold for accessory genes (less than 99% of all assemblies)
    accessory_threshold = total_assemblies * 0.99
    print("Accessory Thresh:", accessory_threshold)
    
    Gene_PresAbs_Core_DF = i_Gene_PresAbs_DF.query(f"NumAsm_WiGene_DNASeq >= {accessory_threshold}")
    
    # Calculating the number of core genes (present in all assemblies)
    num_core_genes = Gene_PresAbs_Core_DF.shape[0]
    
    # Calculating the number of accessory genes (present in less than 99% of all assemblies)
    
    Gene_PresAbs_Acc_DF = i_Gene_PresAbs_DF.query(f"NumAsm_WiGene_DNASeq < {accessory_threshold}")
    
    num_accessory_genes = Gene_PresAbs_Acc_DF.shape[0]

    num_total_genes = i_Gene_PresAbs_DF.shape[0]
    if Verbose:
        print("# of core genes:", num_core_genes)
        print("# of accessory genes:", num_accessory_genes)

    return num_total_genes, num_core_genes, num_accessory_genes


def get_PG_Stats_AdjByIncompCDSAsm(i_Gene_PresAbs_DF, NumSamples):

    # Total number of assemblies
    
    total_assemblies = NumSamples
    
    # Defining the threshold for accessory genes (less than 99% of all assemblies)
    accessory_threshold = total_assemblies * 0.99
    print("Accessory Thresh:", accessory_threshold)
    
    Gene_PresAbs_Core_DF = i_Gene_PresAbs_DF.query(f"NumAsm_WiGene_AdjByIncompCDSAsm >= {accessory_threshold}")
    
    # Calculating the number of core genes (present in all assemblies)
    num_core_genes = Gene_PresAbs_Core_DF.shape[0]
    
    # Calculating the number of accessory genes (present in less than 99% of all assemblies)
    
    Gene_PresAbs_Acc_DF = i_Gene_PresAbs_DF.query(f"NumAsm_WiGene_AdjByIncompCDSAsm < {accessory_threshold}")
    
    num_accessory_genes = Gene_PresAbs_Acc_DF.shape[0]

    num_total_genes = i_Gene_PresAbs_DF.shape[0]
    
    print("# of core genes:", num_core_genes)
    print("# of accessory genes:", num_accessory_genes)

    return num_total_genes, num_core_genes, num_accessory_genes



def get_PG_Stats_FromPresAbs_V2(i_Gene_PresAbs_DF, NumSamples, genefreq_col = "NumAsm_WiGene"):

    # Total number of assemblies
    
    total_assemblies = NumSamples
    
    # Defining the threshold for accessory genes (less than 99% of all assemblies)
    accessory_threshold = total_assemblies * 0.99
    print("Accessory Thresh:", accessory_threshold)
    
    Gene_PresAbs_Core_DF = i_Gene_PresAbs_DF.query(f"{genefreq_col} >= {accessory_threshold}")
    
    # Calculating the number of core genes (present in all assemblies)
    num_core_genes = Gene_PresAbs_Core_DF.shape[0]
    
    # Calculating the number of accessory genes (present in less than 99% of all assemblies)
    
    Gene_PresAbs_Acc_DF = i_Gene_PresAbs_DF.query(f"{genefreq_col} < {accessory_threshold}")
    
    num_accessory_genes = Gene_PresAbs_Acc_DF.shape[0]

    num_total_genes = i_Gene_PresAbs_DF.shape[0]
    
    print("# of core genes:", num_core_genes)
    print("# of accessory genes:", num_accessory_genes)

    return num_total_genes, num_core_genes, num_accessory_genes




def parse_PG_Ref_FA(in_FA_PATH):
    """
    Parses a FASTA file containing pan-genome nucleotide reference sequences (`pan_genome_reference.fa`)
    and returns a dictionary containing the names of the sequences as keys and the sequences themselves as values.
    
    Args:
    - in_FA_PATH (str): The path to the input FASTA file
    
    Returns:
    - dictOf_PG_Ref_Seqs (dict): A dictionary containing the short names of the sequences as keys and the sequences 
    themselves as values.
    """
    
    dictOf_PG_Ref_Seqs = {}
    NumParsedRecords = 0 
    
    for record in screed.open(in_FA_PATH):
        
        NumParsedRecords += 1
        sequence = record.sequence
        name = record.name
        
        if len(name.split(" ")) > 1:
            short_name = " ".join(name.split(" ")[1:])
        else:
            short_name = name

        dictOf_PG_Ref_Seqs[short_name] = sequence

    return dictOf_PG_Ref_Seqs





