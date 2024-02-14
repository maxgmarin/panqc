#!/usr/bin/env python3

### Authors: Max Marin (maximillian_marin@hms.harvard.edu)
# Pan-genome QC toolkit (PQGC)


import sys
import argparse
import os 
from ._version import __version__

from .ava import ava
from .nscluster import clusterBy_KmerJC
from .asm_gene_search import asmseqcheck_frompaths, get_AsmSeqCheck_QCStats, get_AsmSeqCheck_QCStatsDF

import pandas as pd


#from colored import fg, bg, attr

def _nrc_cli(args):
    ## 1) Set input parameters and PATHs ####
    input_Assemblies_TSV = args.asms
    input_PG_Ref_FA = args.pg_ref
    input_PresAbs_CSV = args.gene_matrix
    results_dir = args.results_dir
    min_query_cov = args.min_query_cov
    min_seq_id = args.min_seq_id
    kmer_size = args.kmer_size
    ksim_cluster_thresh = args.min_ksim
    prefix = args.prefix

    ## 2) Run the assembly sequence check function ####
    
    Gene_PresAbs_WiAsmSeqCheck_DF = asmseqcheck_frompaths(input_PresAbs_CSV,
                                                        input_PG_Ref_FA,
                                                        input_Assemblies_TSV,
                                                        min_query_cov,
                                                        min_seq_id)
    
    # 3) Print the general QC Stats 
    ASC_Stats_DF = get_AsmSeqCheck_QCStatsDF(Gene_PresAbs_WiAsmSeqCheck_DF)

    # 4) Run all vs all comparison of sequence k-mer profiles
    PG_AvA_DF = ava(input_PG_Ref_FA, kmer_size)

    # 3) Run the nucleotide similarity clustering and produce updated Presence/Absence matrix + cluster info
    PresAbs_NSC_DF, ClusterInfoGraphDict  = clusterBy_KmerJC(PG_AvA_DF,
                                                            Gene_PresAbs_WiAsmSeqCheck_DF,
                                                            ksim_cluster_thresh)

    ClusterInfo_DF = ClusterInfoGraphDict["Filt_Cluster_DF"]


    # 4) Output files to the results directory

    Step1_OutDir = f"{results_dir}/Step1_AsmSeqCheck"
    Step2_OutDir = f"{results_dir}/Step2_SeqClustering"

    # Create the target directories
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(Step1_OutDir, exist_ok=True)
    os.makedirs(Step2_OutDir, exist_ok=True)


    # Define path to all output directories
    if prefix != "": prefix = prefix + "."

    out_nsc_gene_matrix_TSV = f"{results_dir}/{prefix}gene_presence_absence.NRCUpdated.tsv"

    out_asc_gene_matrix_TSV = f"{Step1_OutDir}/{prefix}gene_presence_absence.AsmSeqCheck.tsv"
    out_ASC_Summ_TSV = f"{Step1_OutDir}/{prefix}AsmSeqCheck.Stats.tsv"

    output_AvA_TSV = f"{Step2_OutDir}/{prefix}AllVsAll.KmerSimilarity.tsv"
    out_cluster_tsv = f"{Step2_OutDir}/{prefix}NSC.ClusterInfo.tsv"


    # Output the NSC updated Presence/Absence matrix (Step #1 & #2)
    print(f"Saving the adjusted gene presence/absence matrix to: {out_nsc_gene_matrix_TSV} \n")
    PresAbs_NSC_DF.to_csv(out_nsc_gene_matrix_TSV,
                          sep = "\t", index = False)

    # Output the ASC updated Presence/Absence matrix (Step #1)
    #print(f" Saving the intermediate (After Assembly Seq Check Step) adjusted gene presence/absence matrix to: {out_asc_gene_matrix_TSV}")
    Gene_PresAbs_WiAsmSeqCheck_DF.to_csv(out_asc_gene_matrix_TSV,
                          sep = "\t", index = False)

    # output the AsmSeqCheck Summary
    ASC_Stats_DF.to_csv(out_ASC_Summ_TSV, sep = "\t", index = False)

    ## Output the AvA comparison table
    print(f"Saving the All vs All comparison table (TSV) to: {output_AvA_TSV} \n")
    PG_AvA_DF.to_csv(output_AvA_TSV, sep = "\t", index = False)

    # Output the cluster info table
    print(f"Saving info for all identified nucleotide similarity clusters (TSV) to: {out_cluster_tsv} \n")

    ClusterInfo_DF.to_csv(out_cluster_tsv,
                          sep = "\t", index = False)






def _ava_cli(args):
    ## 1) Set input parameters and PATHs ####
    input_PG_Ref_FA = args.in_pg_ref

    output_AvA_TSV = args.out_ava_tsv

    kmer_size = args.kmer_size

    ## Run the All vs All comparison function 
    PG_AvA_DF = ava(input_PG_Ref_FA, kmer_size)

    ## Output the AvA comparison table
    print(f" Saving the All vs All comparison table (TSV) to: {output_AvA_TSV}")
    PG_AvA_DF.to_csv(output_AvA_TSV, sep = "\t", index = False)


def _asmseqcheck_cli(args):

    ## 1) Set input parameters and PATHs ####
    
    # Define input paths
    input_PG_Ref_FA = args.in_pg_ref

    input_AsmFA_TSV = args.in_assemblies

    input_PresAbs_CSV = args.in_gene_matrix

    # Define output path
    output_PresAbs_WiDNASeqCheck = args.out_gene_matrix_wi_geneseqcheck
    
    # Set the minimum query coverage and sequence identity based on user input 
    min_query_cov = args.min_query_cov
    min_seq_id = args.min_seq_id

    ## 2) Run the assembly sequence check function ####

    Gene_PresAbs_WiAsmSeqCheck_DF = asmseqcheck_frompaths(input_PresAbs_CSV,
                                                          input_PG_Ref_FA,
                                                          input_AsmFA_TSV,
                                                          min_query_cov,
                                                          min_seq_id)

    # 3) Print the general QC Stats

    _ = get_AsmSeqCheck_QCStats(Gene_PresAbs_WiAsmSeqCheck_DF)


    # 4) Output TSV 
    Gene_PresAbs_WiAsmSeqCheck_DF.to_csv(output_PresAbs_WiDNASeqCheck,
                                         sep = "\t",
                                         index = False)


def _nscluster_cli(args):

    ## 1) Set input parameters and PATHs ####
    in_AvA_TSV = args.in_ava_tsv

    gene_matrix_ASC_TSV = args.in_gene_matrix_tsv

    ksim_cluster_thresh = args.min_ksim

    out_nsc_gene_matrix_TSV = args.out_nsc_gene_matrix

    out_cluster_tsv = args.out_clusterinfo_tsv


    # 2) Read in the input files as dataframes

    in_AvA_DF = pd.read_csv(in_AvA_TSV, sep = "\t")

    i_Gene_PresAbs_DF = pd.read_csv(gene_matrix_ASC_TSV, sep = "\t")

    # 3) Run the nucleotide similarity clustering and produce updated Presence/Absence matrix + cluster info
    PresAbs_NSC_DF, ClusterInfoGraphDict  = clusterBy_KmerJC(in_AvA_DF,
                                                             i_Gene_PresAbs_DF,
                                                             ksim_cluster_thresh)

    ClusterInfo_DF = ClusterInfoGraphDict["Filt_Cluster_DF"]

    # 4) Output matrix and cluster info

    # Output the NSC updated Presence/Absence matrix
    PresAbs_NSC_DF.to_csv(out_nsc_gene_matrix_TSV,
                          sep = "\t", index = False)

    # Output the cluster info table
    ClusterInfo_DF.to_csv(out_cluster_tsv,
                          sep = "\t", index = False)



# fg("blue") 
def main():
    parser = argparse.ArgumentParser(description = "Toolkit for focused on augmenting common CDS based pan-genome analysis with nucleotide sequence comparison.")
    sub_parser_1 = parser.add_subparsers(required=True, help='Please select one of the pipelines of the PanQC toolkit.')


    nrc_parser = sub_parser_1.add_parser("nrc", help="(Nucleotide Redundancy Correction)   Adjusts for nucleotide redundancy in estimated pan-genomes")
    nrc_parser.add_argument('-a', '--asms', type=str, required=True, metavar="PathToAsms.tsv",
                                                    help="Table with SampleID & Paths to each input assemblies. (TSV)")

    nrc_parser.add_argument('-r', '--pg-ref', type=str, required=True,  metavar="pan_genome_reference.fasta",
                                                    help="Input pan-genome nucleotide reference. Typically output as `pan_genome_reference.fasta` by Panaroo/Roary (FASTA)")

    nrc_parser.add_argument('-m', '--gene_matrix', type=str, required=True, metavar="gene_presence_absence.csv",
                                                    help="Input pan-genome gene presence/absence matrix. Typically output as `gene_presence_absence.csv` by Panaroo/Roary (CSV)")

    nrc_parser.add_argument('-o', '--results_dir', type=str, required=True,
                                                    help="Output directory for analysis results.")
    
    nrc_parser.add_argument('-p', '--prefix', type=str, required=False, default="",
                            help="Prefix to append to output files")

    nrc_parser.add_argument('-c', '--min-query-cov', type=float, default=0.9,
                            help="Minimum query coverage (ranging from 0 to 1) to classify a gene as present within an assembly (Default: 0.9)")

    nrc_parser.add_argument('-i', '--min-seq-id', type=float, default=0.9, 
                            help="Minimum sequence identity (ranging from 0 to 1) to classify a gene as present within an assembly  (Default: 0.9)")
    
    nrc_parser.add_argument('-k', '--kmer_size',type=int, default=31,
                            help="k-mer size (bp) to use for generating profile of unique k-mers for each sequence (Default: 31))")
    
    nrc_parser.add_argument('-t', '--min-ksim',type=float, default=0.8,
                            help='Minimum k-mer similarity (maximum jaccard containment of k-mers between pair of sequences) to cluster sequences into the same "nucleotide similarity cluster" (Default: 0.8))')

    nrc_parser.set_defaults(func=_nrc_cli)


    utils_parser = sub_parser_1.add_parser("utils", help="(Utilities)   Sub-pipelines related to NRC approach or other types of pan-genome workflows")
    #utils_parser.set_defaults(func=_utils_help)
    # Store the parser in the args for access in the default function
    #utils_parser.set_defaults(parser=utils_parser)

    utils_sub_parser = utils_parser.add_subparsers(required=True, help='Please select one of the utilility pipelines of the PanQC toolkit.')

    asmseqcheck_parser = utils_sub_parser.add_parser("asmseqcheck", help="")
    asmseqcheck_parser.add_argument('-a', '--in_assemblies', type=str, required=True,
                                    help="Paths to input assemblies. (TSV)")

    asmseqcheck_parser.add_argument('-r', '--in_pg_ref', type=str, required=True,
                                    help="Input pan-genome nucleotide reference. Typically output as `pan_genome_reference.fasta` by Panaroo/Roary (FASTA)")

    asmseqcheck_parser.add_argument('-m', '--in_gene_matrix', type=str, required=True,
                                    help="Input pan-genome gene presence/absence matrix. Typically output as `gene_presence_absence.csv` by Panaroo/Roary (CSV)")

    asmseqcheck_parser.add_argument('-o', '--out_gene_matrix_wi_geneseqcheck',type=str, required=True,
                                    help="Output pan-genome gene presence/absence matrix with updated gene presence/absence calls. (CSV). \n NOTE: 2 reflects that similar gene sequence is present at the nucleotide level (CSV)")

    asmseqcheck_parser.add_argument('-c', '--min_query_cov', type=float, default=0.9,
                            help="Minimum query coverage to classify a gene as present within an assembly (0-1)")

    asmseqcheck_parser.add_argument('-i', '--min_seq_id', type=float, default=0.9,
                            help="Minimum sequence identity to classify a gene as present within an assembly (0-1)")

    asmseqcheck_parser.set_defaults(func=_asmseqcheck_cli)


    ava_parser = utils_sub_parser.add_parser("ava", help="")
    ava_parser.add_argument('-i', '--in_pg_ref', type=str, required=True, metavar="pan_genome_reference.fasta",
                            help="Input Panaroo pan-genome nucleotide reference. Typically output as `pan_genome_reference.fasta` by Panaroo (FASTA)")
    ava_parser.add_argument('-o', '--out_ava_tsv',type=str, required=True, metavar="AllVsAll.KmerSimilarity.tsv",
                            help="All vs all comparison of sequence k-mer profiles. (TSV)")
    ava_parser.add_argument('-k', '--kmer_size',type=int, default=31,  metavar="k-mer size (bp)",
                            help="k-mer size (bp) to use for generating profile of unique k-mers for each sequence (Default: 31))")
    ava_parser.set_defaults(func=_ava_cli)


    cluster_parser = utils_sub_parser.add_parser("nscluster", help="")
    cluster_parser.add_argument('-a', '--in_ava_tsv',type=str, required=True,
                            help="Input table with all vs all comparison of sequence k-mer profiles. (TSV)")

    cluster_parser.add_argument('-m', '--in_gene_matrix_tsv',type=str, required=True,
                                    help="Input pan-genome gene presence/absence matrix (CSV). \n NOTE: 2 reflects that similar gene sequence is present at the nucleotide level (CSV)")

    cluster_parser.add_argument('--min_ksim',type=float, default=0.8,
                            help='Minimum k-mer similarity (maximum jaccard containment of k-mers between pair of sequences) to cluster sequences into the same "nucleotide similarity cluster" (Default: 0.8))')

    cluster_parser.add_argument('-o', '--out_nsc_gene_matrix',type=str, required=True,
                            help="Nucleotide Similarity Cluster adjusted Gene Presence Matrix (TSV/Rtab)")

    cluster_parser.add_argument('-c', '--out_clusterinfo_tsv',type=str, required=True,
                            help='Summary table with all genes that belong to a NSC (Nucleotide Similarity Cluster) (TSV)')
    cluster_parser.set_defaults(func=_nscluster_cli)


    # Check if no arguments were provided, if so exit with help message
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    # Output utils help message if no subcommand is provided
    if "utils" in sys.argv and len(sys.argv) == 2:
        utils_parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    sys.exit(main())
