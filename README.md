# snp_tools
A set of tools and add-ons for SNP analysis 

##Introduction 
Tools and add-ons primarily for working with the output of [kSNP](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0081760), an alignment-free SNP identifying software package. 

##Requirements
All packages are written in Perl and require the bash shell. Prerequistes below must be available in your PATH to function. 

Required for all:
* bash
* perl

Required for some (ksnp_snp_confirm.pl, kmer_core.pl, kmer_compare.pl):
* [jellyfish](http://www.cbcb.umd.edu/software/jellyfish/) (>= v1.1.2).  This is likely already installed and in your path if you have kSNP.

Optional for some:
* gzip

##Installation
Simply download the folder containing the scripts and place it somewhere you'll be able to find it.
**IMPORTANT**: The _lib_ directory must stay in the same directory as the parent scripts. Do not move or rename this directory or any of its contents.

##Programs
###ksnp_snp_confirm.pl
```
ksnp_snp_confirm.pl [options] -f file_of_files.txt -a SNPs_all
```
Searches sequence reads for SNP kmers output by kSNP, confirms the SNP assignments, and outputs new SNPs_all and SNPs_all_matrix.fasta files with corrected positions.
Does the following corrections:

1. If a kSNP SNP position is found with multiple midddle bases in the reads (exceeding a given error rate, -e), it will be replaced with a "-" in the matrix.fasta file and omitted from the SNPs_all file.
2. If a kSNP SNP position was listed as absent in the genome, but a minimum number of reads contain the SNP position (as given by -l) and are not too mixed (maximum error given by -e), the base will be included in the matrix.fasta file and the SNPs_all file.

For simplicity, corrected genome positions in the SNPs_all file will be listed as "0" and all kmer directions will be listed as "F".
Will also output, to STDOUT, a listing of modified SNP positions. The categories are:

| Header | Description |
| --- | --- |
| `#_positions` | Total number of SNP positions identified by kSNP |
| `#_confirmed` | Total number of SNP positions where the base call was confirmed by the reads |
| `#_mixed` | Total number of SNP positions where a mixture of bases were found in the reads with an error rate above the threshold set by -e. |
| `#_omitted` | Total number of SNP positions that were found to be present and not mixed in the reads, but were not listed for this genome in the SNPs_all file. |

**REQUIRED:**

`-f`    file of files.  Contains genome names (as provided to kSNP) and paths to one or more fastq read files. Gzipped or uncompressed fastq read files may be used. If no read file is listed, the SNP positions will be included in the output files as given in the SNPs_all input file.

Format of example file:
```
genome_A<tab>/path/to/genome_A_reads_1.fastq.gz<tab>/path/to/genome_A_reads_2.fastq.gz
genome_B<tab>/path/to/genome_B_reads.fastq.gz
genome_C
genome_D<tab>/path/to/genome_C_reads.fastq
etc..
```

`-a`    SNPs_all file (as output by kSNP)

**OPTIONS:**

`-l`    minimum number of times a kmer must be present to be counted. Higher numbers should run faster (default 5)

`-e`    maximum allowable error rate between number of reads with reference base and number of reads with discrepant base. This is calcuated as: `1 - (Nref / (Na + Nc + Ng + Nt))` where Nref is the number of reads with the reference base and Na is the number of reads with "A" at the SNP position, Nc is the number of reads with "C" at the SNP position, etc. If the `-d` option is given and the error rate is greater than this value, the position will be output. (default: 0.1)

`-p`    output file prefix (default: "corrected")

`-t`    number of threads (default: 15)

###kmer_core.pl

###kmer_compare.pl

###kmer_compare_groups.pl

###ksnp_matrix_filter.pl

###ksnp_matrix_to_diff_matrix.pl
