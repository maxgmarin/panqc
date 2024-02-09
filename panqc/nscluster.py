# pqgc/nscluster.py

import networkx as nx
import pandas as pd
import numpy as np

import time 

from pgqc.utils import get_columns_excluding


import logging
# Set the logging level to INFO
logging.basicConfig(level=logging.INFO)


def create_MaxKmerSim_JC_Dict(in_AvA_DF):

    i_MaxSim_Dict = {}
    
    # Iterating through the dataframe to find the maximum similarity value for each pair of genes
    for index, row in in_AvA_DF.iterrows():
        
        pair = tuple(sorted([row['RecordID_1'], row['RecordID_2']])) # Sorting to ensure consistency
        
        current_similarity = row['MaxJaccardContain']
        if pair not in i_MaxSim_Dict or current_similarity > i_MaxSim_Dict[pair]:
            
            if row['RecordID_1'] != row['RecordID_2']:
                i_MaxSim_Dict[pair] = current_similarity


    # Converting the maximum similarity values into a DataFrame for easier manipulation
    i_MaxSim_DF = pd.DataFrame.from_dict(i_MaxSim_Dict, orient='index', columns=['MaxJaccardContain'])
    i_MaxSim_DF.reset_index(inplace=True)
    i_MaxSim_DF[['RecordID_1', 'RecordID_2']] = pd.DataFrame(i_MaxSim_DF['index'].tolist(),
                                                             index = i_MaxSim_DF.index)
    
    i_MaxSim_DF.drop(columns=['index'], inplace=True)
    
    return i_MaxSim_Dict



def get_cluster_label(cluster):
    return '---'.join(cluster)


def make_ClusterID_Maps(i_clusters):
    
    dict_ClusterToGenes = {}
    dict_ClusterNumToGenes = {}
    
    dict_GeneToCluster = {}
    dict_GeneToClusterNum = {}
    
    listOf_Gene_ClusterInfo = []
    
    for n, cluster in enumerate(i_clusters):
        
        cluster_label = get_cluster_label(cluster)
        cluster_n = n
            
        dict_ClusterToGenes[cluster_label] = cluster
        dict_ClusterNumToGenes[cluster_n] = cluster
    
        for node in cluster:
            dict_GeneToCluster[node] = cluster_label
            dict_GeneToClusterNum[node] = cluster_n
    
            i_row = (node, cluster_label, cluster_n)
    
            listOf_Gene_ClusterInfo.append(i_row)


    ClusterInfo_Dict = {}
    ClusterInfo_Dict["ClusterToGenes"] = dict_ClusterToGenes
    ClusterInfo_Dict["ClusterNumToGenes"] = dict_ClusterNumToGenes 
    ClusterInfo_Dict["GeneToCluster"] = dict_GeneToCluster
    ClusterInfo_Dict["GeneToClusterNum"] = dict_GeneToClusterNum


    Gene_ClusterInfo = pd.DataFrame(listOf_Gene_ClusterInfo)
    
    Gene_ClusterInfo.columns = ["GeneID", "NS_ClusterName", "NS_ClusterNum"]

    Gene_ClusterInfo["NS_ClusterID"] = "NS_" + Gene_ClusterInfo["NS_ClusterNum"].astype(str) + "_" + Gene_ClusterInfo["NS_ClusterName"]

    
    return Gene_ClusterInfo, ClusterInfo_Dict




def create_MST_FiltByJC(in_AvA_DF, i_MaxSim_Dict, threshold ):

    # Creating the filtered graph again
    filtered_graph_max = nx.Graph()
    for index, row in in_AvA_DF.iterrows():
        
        pair = tuple(sorted([row['RecordID_1'], row['RecordID_2']])) # Sorting to ensure consistency
        
        # Adding an edge only if the RecordIDs are different (avoiding self-loops)
        if row['RecordID_1'] != row['RecordID_2']:
            if i_MaxSim_Dict[pair] >= threshold:
                
                filtered_graph_max.add_edge(row['RecordID_1'], row['RecordID_2'], weight = i_MaxSim_Dict[pair] )

    # Finding the minimum spanning tree
    MST_Filt = nx.minimum_spanning_tree(filtered_graph_max)


    # Extracting connected components (clusters)
    clusters = [list(component) for component in nx.connected_components(MST_Filt)]

    
    return MST_Filt, clusters



PresAbs_NonSampleID_ColNames = ['Gene', 'NumAsm_WiGene', 'NumAsm_WiGene_DNASeq',
                                'Non-unique Gene name', 'Annotation', 'No. isolates',
                                'No. sequences', 'Avg sequences per isolate', 'Genome Fragment',
                                'Order within Fragment', 'Accessory Fragment', 'Accessory Order with Fragment', 'QC',
                                'Min group size nuc', 'Max group size nuc', 'Avg group size nuc']


def make_NS_ClusterMerged_Pres_DF(i_Gene_PresAbs_DF, i_Filt_Cluster_DF, ResetToBinary=True):
    
    i_SampleIDs = get_columns_excluding(i_Gene_PresAbs_DF, PresAbs_NonSampleID_ColNames)
    # Extracting the original presence/absence information for genes (excluding the extra rows)
    original_gene_presence_absence = i_Gene_PresAbs_DF.set_index('Gene')[i_SampleIDs]
    
    # Checking the unique genes in the cluster information
    clustered_genes = set(i_Filt_Cluster_DF['GeneID'])

    # Creating a dictionary to map each gene to its corresponding cluster ID
    gene_to_clusterID_map = dict(zip(i_Filt_Cluster_DF["GeneID"], i_Filt_Cluster_DF["NS_ClusterID"]))
    
    # Extracting the gene presence/absence information for clustered genes
    clustered_gene_presence_absence = original_gene_presence_absence.loc[original_gene_presence_absence.index.isin(clustered_genes)]
    
    # Initializing a dictionary to hold the merged presence/absence information for clusters
    merged_cluster_presence_absence = {}
    
    # Iterating through the clustered gene presence/absence information
    for gene, row in clustered_gene_presence_absence.iterrows():
        cluster_id = gene_to_clusterID_map.get(gene)
        # Merging the presence/absence information for genes within a cluster
        merged_cluster_presence_absence[cluster_id] = [max(val, merged_cluster_presence_absence.get(cluster_id, [0] * len(i_SampleIDs))[idx]) for idx, val in enumerate(row)]
    
    # Converting the merged cluster presence/absence information into a DataFrame
    merged_cluster_presence_absence_DF = pd.DataFrame.from_dict(merged_cluster_presence_absence, orient='index', columns=i_SampleIDs).reset_index().rename(columns={'index': 'Gene'})
    
    # Concatenating with the original presence/absence matrix to include non-clustered genes and extra rows
    final_presence_absence_DF_with_clusters_and_all_info = pd.concat([i_Gene_PresAbs_DF[~i_Gene_PresAbs_DF['Gene'].isin(clustered_genes)], merged_cluster_presence_absence_DF], ignore_index=True)
    
    # Displaying the first few rows of the final presence/absence matrix
    # final_presence_absence_DF_with_clusters_and_all_info.head()
    
    PG_Pres_WiNSC_DF = final_presence_absence_DF_with_clusters_and_all_info.copy()
    
    #PG_Pres_WiNSC_DF["NumAsm_WiGene"] = PG_Pres_WiNSC_DF[i_SampleIDs].sum(axis = 1)

    PG_Pres_WiNSC_DF["NumAsm_WiGene"] = PG_Pres_WiNSC_DF[i_SampleIDs].applymap(lambda x: 1 if x > 0 else 0).sum(axis = 1)

    PG_Pres_WiNSC_DF = PG_Pres_WiNSC_DF.sort_values(by='NumAsm_WiGene', ascending=False)

    if ResetToBinary: # This step will reset all NON zero values to 1. (Default Behavior)
        PG_Pres_WiNSC_DF[i_SampleIDs] = PG_Pres_WiNSC_DF[i_SampleIDs].applymap(lambda x: 1 if x > 0 else 0)

    return PG_Pres_WiNSC_DF

    
def summarize_NSClusters(i_Cluster_DF, i_CoreGenes_List, printStats = True):

    i_Cluster_DF["GeneType"] = i_Cluster_DF["GeneID"].isin(i_CoreGenes_List).replace(True, "Core").replace(False, "Accessory")

    i_Cluster_AccAndCore_Class = i_Cluster_DF.groupby("NS_ClusterID").apply(lambda x: "-".join(list(np.sort(x["GeneType"].unique()))))
    i_CountOf_ClusterType = i_Cluster_AccAndCore_Class.value_counts().to_dict()
    i_Cluster_NumGenes = i_Cluster_DF.groupby("NS_ClusterID")["GeneID"].nunique()

    # Create a new summary DF 
    i_ClusterSumm_DF = pd.concat([i_Cluster_AccAndCore_Class,
                                  i_Cluster_NumGenes,], axis = 1)
    i_ClusterSumm_DF.columns = ["Type", "NumGenes"]
    
    i_ClusterSumm_DF = i_ClusterSumm_DF.reset_index()

    if printStats:
        Num_NSClusters = i_ClusterSumm_DF.shape[0]
        Num_ClusteredGenes = i_Cluster_DF.shape[0]
    
        Num_CoreOnly = i_CountOf_ClusterType.get("Core", 0)
        Num_AccOnly = i_CountOf_ClusterType.get("Accessory", 0)
        Num_AccAndCore = i_CountOf_ClusterType.get("Accessory-Core", 0)
    
        print("# NucSim Clusters:", Num_NSClusters)
        print("# of total clustered genes:", Num_ClusteredGenes)
        print("# of NS clusters w/ only Core genes:", Num_CoreOnly)
        print("# of NS clusters w/ only Accessory genes:", Num_AccOnly)
        print("# of NS clusters w/ BOTH Core & Accessory genes:", Num_AccAndCore)

    return i_ClusterSumm_DF, i_CountOf_ClusterType

def summarize_NSClusters_Simple(i_Cluster_DF, printStats = True):

    i_Cluster_NumGenes = i_Cluster_DF.groupby("NS_ClusterID")["GeneID"].nunique()

    # Create a new summary DF 
    i_ClusterSumm_DF = pd.concat([i_Cluster_NumGenes,], axis = 1)
    i_ClusterSumm_DF.columns = ["NumGenes"]
    
    i_ClusterSumm_DF = i_ClusterSumm_DF.reset_index()

    if printStats:
        Num_NSClusters = i_ClusterSumm_DF.shape[0]
        Num_ClusteredGenes = i_ClusterSumm_DF["NumGenes"].sum()

        print("# NucSim Clusters:", Num_NSClusters)
        print("# of total clustered genes:", Num_ClusteredGenes)

    return i_ClusterSumm_DF





## Define complete function for clustering and merging of PG matrix

def clusterBy_KmerJC(in_AvA_DF, i_Gene_PresAbs_DF, JC_threshold):
    
    # Creat a Dict of the maximum similarity between all pairs
    AvA_MaxSim_Dict = create_MaxKmerSim_JC_Dict(in_AvA_DF)
    
    MST_Filt, Filt_Clusters = create_MST_FiltByJC(in_AvA_DF, AvA_MaxSim_Dict, JC_threshold)
    
    Filt_Cluster_DF, Filt_ClusterDict = make_ClusterID_Maps(Filt_Clusters)

    ClusterInfoGraph = {}
    ClusterInfoGraph["MST_Filt_Graph"] = MST_Filt
    ClusterInfoGraph["Filt_Clusters"] = Filt_Clusters
    ClusterInfoGraph["Filt_Cluster_DF"] = Filt_Cluster_DF
    ClusterInfoGraph["Filt_ClusterDict"] = Filt_ClusterDict

    # Filt_GenesInACluster = Filt_Cluster_DF["GeneID"]
    
    i_Gene_PresAbs_NSC_Filt_DF = make_NS_ClusterMerged_Pres_DF(i_Gene_PresAbs_DF, Filt_Cluster_DF)

    return i_Gene_PresAbs_NSC_Filt_DF, ClusterInfoGraph




def run_nscluster(in_AvA_DF, i_Gene_PresAbs_DF, JC_threshold = 0.8):

    PresAbs_NSC_Filt08_DF, ClusterInfoGraphDict  = clusterBy_KmerJC(in_AvA_DF, i_Gene_PresAbs_DF, 0.8)

    Filt08_Cluster_DF = ClusterInfoGraphDict["Filt_Cluster_DF"]

    print("not yet implemented")
    return None














