<img width="500" src="https://github.com/maxgmarin/panqc/raw/main/Images/panqc.logo.png" alt="panqc logo">

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!---[![Build Status]()]()
[![github release version]()]()
[![DOI]()]()
--->

A **pan**-genome **q**uality **c**ontrol toolkit for evaluating nucleotide redundancy in pan-genome analyses. 

<!---
> TBD Reference
--->

[TOC]: #
## Table of Contents
- [Motivation](#motivation)
- [Installation](#installation)
  - [Install locally](#install-locally)
  - [`pip`](#pip)
  - [`conda`](#conda)
- [Basic usage](#basic-usage)
  - [Analyzing included test data set](#analyzing-included-test-data-set)
- [Full usage](#full-usage)
  - [`nrc`](#panqc-nrc)
  - [`utils`](#panqc-utils)
- [Contributing and Issues](#contributing-and-issues)
- [Citing](#citing)


## Motivation
![PanQC_NRC_Diagram](https://github.com/maxgmarin/panqc/raw/main/Images/PanQC_NRC_Diagram.png)

The Nucleotide Redundancy Correction (PGQC-NRC) pipeline adjusts for redundancy at the DNA level in two steps (Methods). In step one, all genes predicted to be absent at the Amino Acid (AA) level are compared to their corresponding assembly at the nucleotide level. In cases where the nucleotide sequence is found with high coverage and sequence identity (Query Coverage & Sequence Identity > 90%), the gene is marked as “present at the DNA level”. Next, all genes are clustered and merged using a k-mer based metric of nucleotide similarity. Cases where two or more genes are divergent at the AA level but highly similar at the nucleotide level will be merged into a single “nucleotide similarity gene cluster”. After applying this method the pan-genome gene presence matrix is readjusted according to these results.

<!---
**When to use this software**:
**When you should probably NOT use this software:** If you care about... 
--->

## Installation
### Install locally
Currently, `panqc` can be installed by cloning this repository and installing with `pip`.

```
git clone git@github.com:maxgmarin/pgqc.git

cd pgqc

pip install . 
```

### `pip`
```
pip install panqc
```

### `conda`
>🚧 Check back soon 🚧

## Basic usage

```
panqc nrc -a InputAsmPaths.tsv -r pan_genome_reference.fa -m gene_presence_absence.csv -o NRC_results/
```

The above command will output an adjusted gene presence absence matrix along with additional statistics to the specified output directory (`NRC_results/`).

### Analyzing included test data set

If you wish to run an `panqc nrc` on an artifical (and abridged) test data set, you can run the following commands:

```
cd tests/data

# Define path to the 3 needed input files:

# 1) Gene presence absence matrix (As output by Panaroo or Roary)
PG_Matrix_CSV="TestSet1.gene_presence_absence.csv"

# 2) Pan-genome nucleotide reference (As output by Panaroo or Roary)
PG_Ref_FA="TestSet1.pan_genome_reference.fa.gz"

# 3) SampleID + Path for all assemblies used in analysis
Asm_TSV="TestSet1.InputAsmPaths.tsv"

time panqc nrc -a ${Asm_TSV} -r ${PG_Ref_FA} -m ${PG_Matrix_CSV} -o test_results/
```
NOTE: Make sure that your current working directory (CWD) is `tests/data` within the repository. The `TestSet1.InputAsmPaths.tsv` file describes the path to each genome assembly relative to your CWD.


## Full usage

`panqc` has 2 sub-commands:
- `nrc` - Run the full **N**ucleotide **R**edundancy **C**orrection pipeline on a pan-genome analyses.
- `utils` - Run utlity scripts and sub-pipelines of the full NRC pipeline

---

### `panqc nrc`

Run the complete Nucleotide Redundancy Correction pipeline

```
$ panqc nrc --help

usage: panqc nrc [-h] -a PathToAsms.tsv -r pan_genome_reference.fasta -m gene_presence_absence.csv -o RESULTS_DIR [-p PREFIX] [-c MIN_QUERY_COV] [-i MIN_SEQ_ID] [-k KMER_SIZE] [-t MIN_KSIM]

optional arguments:
  -h, --help            show this help message and exit

  -a, --asms PathToAsms.tsv
                        Table with SampleID & Paths to each input assemblies. (TSV)

  -r, --pg-ref pan_genome_reference.fasta
                        Input pan-genome nucleotide reference. Typically output as `pan_genome_reference.fasta` by Panaroo/Roary (FASTA)

  -m, --gene_matrix gene_presence_absence.csv
                        Input pan-genome gene presence/absence matrix. Typically output as `gene_presence_absence.csv` by Panaroo/Roary (CSV)

  -o, --results_dir RESULTS_DIR
                        Output directory for analysis results.

  -p, --prefix PREFIX
                        Prefix to append to output files

  -c, --min-query-cov MIN_QUERY_COV
                        Minimum query coverage (ranging from 0 to 1) to classify a gene as present within an assembly (Default: 0.9)

  -i, --min-seq-id MIN_SEQ_ID
                        Minimum sequence identity (ranging from 0 to 1) to classify a gene as present within an assembly (Default: 0.9)

  -k, --kmer_size KMER_SIZE
                        k-mer size (bp) to use for generating profile of unique k-mers for each sequence (Default: 31))

  -t, --min-ksim MIN_KSIM
                        Minimum k-mer similarity (maximum jaccard containment of k-mers between pair of sequences) to cluster sequences into the same "nucleotide similarity cluster" (Default: 0.8))
```


### `panqc utils`

Within `utils` there are 3 sub-commands that run specific components of the NRC pipeline:
- `utils asmseqcheck` - Perform alignment of all genes classified as absent to their respective assemblies.
- `utils ava` - Perform all vs all comparison of k-mer profiles of input sequences. 
- `utils nscluster` - Perform nucleotide similarity clustering and readjust pan-genome estimates.

```
$ panqc utils --help

usage: panqc utils [-h] {asmseqcheck,ava,nscluster} ...

positional arguments:
  {asmseqcheck,ava,nscluster}
                        Please select one of the utilility pipelines of the panqc toolkit.
    asmseqcheck
    ava
    nscluster

optional arguments:
  -h, --help            show this help message and exit

```

>🚧 Check back soon for full usage for each of the utility sub-pipelines of the NRC pipeline 🚧

## Contributing and Issues
>🚧 Check back soon 🚧

## Citing
>🚧 Check back soon 🚧

<!---
If you use `panqc` in your work, please cite:
> TBD
--->
