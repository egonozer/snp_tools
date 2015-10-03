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

`-e`    maximum allowable error rate between number of reads with reference base and number of reads with discrepant base. This is calcuated as: `1 - (Nref / (Na + Nc + Ng + Nt))` where Nref is the number of reads with the reference base and Na is the number of reads with "A" at the SNP position, Nc is the number of reads with "C" at the SNP position, etc. (default: 0.1)

`-p`    output file prefix (default: "corrected")

`-t`    number of threads (default: 15)

###kmer_core.pl
```
kmer_core.pl
```
**DESCRIPTION:** Determines the set of kmers that are present in a given subset of a group of input genomes. This can be used as input to kSNP.

**PARAMETERS:**

Required:

`-f`    text file with names and paths to input genomes, separated by any white space. Each file listed should be in fasta format. Multi-contig genome files are permitted, but each sequence file will be assumed to contain sequence for only one genome.

*OPTIONAL:* After the genome name and the path to the genome sequence file, paths to one or more sequencing read files (in fastq or fasta format) can be included. These will be used to find kmers with conflicting alleles that were collapsed in assembly and/or kmers present in reads, but for which assembly failed in one or more strains. New kmers will not be identified from read sets. They will only be used to confirm kmers found in at least one of the input genomic sequence files. Gzipped read files (with .gz extension) can be given.

*FORMAT:*

`<genome_id> </path/to/sequence.fasta> </path/to/reads.fastq OPTIONAL> </path/to/reads2.fastq.gz OPTIONAL> etc.`

*EXAMPLE:*
```        
PAO1 seqs/PAO1_contigs.fasta
PA14 seqs/PA14_contigs.fasta reads/PA14_1.fastq.gz reads/PA14_2.fastq.gz
PA7 seqs/PA7_contigs.fasta reads/PA7_1.fastq.gz
```
---
Optional:

`-a`    minimum percentage of the total input genomes in which a kmer must be found to be considered core (default: 100)

`-k`    kmer size. Maximum value is 31. (default: 31)

`-n`    number of N's with which to separate kmers in the core kmers sequence output file (default: 1)

`-t`    number of threads (default: 15)

`-o`    output files prefix (default: output)

`-m`    [in reads file(s)] minimum number of times a kmer must be present to be counted. (default: 5)

`-e`    [in reads file(s)] maximum allowable error rate between number of reads with reference base and number of reads with discrepant base. This is calcuated as: `1 - (Nref / (Na + Nc + Ng + Nt))` where `Nref` is the number of reads with the reference base and `Na` is the number of reads with "A" at the SNP position, `Nc` is the number of reads with "C" at the SNP position, etc. (default: 0.1)


###kmer_compare.pl
```
kmer_compare.pl
```
Determines SNP differences among pairs of sequences.

**Required:**

`-f`    file containing sequences or kmers of input strains
        
Strains may be provided assembled as fasta files or raw read files or both. To save time (if performing multiple comparisons), can also pre-process sequence files using 'jellyfish count' and 'jellyfish dump' and provide the resulting non-binary kmer quantity files. Just make sure all the kmer files have the same length kmers or else the script will exit with an error.
        
*VERY IMPORTANT:* If providing pre-processed jellyfish files, you must ensure that jellyfish count was run with the '-C' (both-strands) option or this program's output will almost certainly be incorrect.
        
Each line in the file should have three columns separated by tabs.

1. 1st column: Unique strain ID
2. 2nd column: Path to sequence or kmer file
3. 3rd column: File type. The possible file types are as follows
  * `fs`    = assembled sequence file in fasta format
  * `fr`    = sequence reads in fasta or fastq format
  * `js`    = kmer file
  * `jr`    = kmer file from reads

If the third column entry does not match these types or is left blank, the file will be assumed to be of type 'fs'.
        
The script accepts gzipped or uncompressed files.
        
If multiple files are provided for a strain (i.e. forward and reverse reads), their paths must be provided on separate lines with the same strain ID (column 1). If both assembly and read sequences are given for a strain, only kmers found in the assembly will be compared and reads will only be used to identify conflicting kmers within a genome.
        
Optional comments may follow any line as long as they are preceded by a `#`.
        
Example file:
```
genA    /path/to/genA.fasta fs          #genA is a sequence assembly
genB    /path/to/genB.fasta fs          #genB is a sequence assembly
genB    /path/to/genB_1.fastq.gz    fr  #including sequence reads for genB
genC    /path/to/genC_1.fastq.gz    fr  #only providing reads for genC
genC    /path/to/genC_2.fastq.gz    fr
genD    /path/to/genD_kmers.txt     js  #file of kmers from the genD assembly
genE    /path/to/genE_kmers.txt.gz  jr  #file of kmers from the genE reads (gzipped)
```
        
*OPTIONAL:* If you would only like to compare one or more reference genomes to a group of other genomes, and not every possible all-vs-all pairwise comparison, then add a `*` character at the end of the desired reference genomes'line(s).

For example:

```
genA    /path/to/genA.fasta fs  *
genB    /path/to/genB.fasta fs          
genB    /path/to/genB_1.fastq.gz    fr
genC    /path/to/genC_1.fastq.gz    fr
genC    /path/to/genC_2.fastq.gz    fr
genD    /path/to/genD_kmers.txt     js
genE    /path/to/genE_kmers.txt.gz  jr
```

In the above example, instead of all possible pairwise comparisons of the 5 genomes, the comparisons performed will be:

* genA to genB
* genA to genC
* genA to genD
* genA to genE
* done.
        
**Optional:**

`-c`    File of kmers to which to limit the analysis, i.e. the file of core kmers output by 'kmer_core.pl'. This should be a fasta file of kmers with one or more N's separating the kmers. If this file is given, it will provide the kmer size and any entry to (-k) will be ignored. (default: all SNP kmers in common between two strains will be counted)

`-k`    kmer size.  This setting will be ignored if jellyfish kmer files are provided or if a list of kmers is given using (`-c`) (default: 31)

`-m`    [in reads file(s)] minimum number of times a kmer must be present to be counted. (default: 5)

`-e`    [in reads file(s)] maximum allowable error rate between number of reads with reference base and number of reads with discrepant base. This is calcuated as: `1 - (Nref / (Na + Nc + Ng + Nt))` where `Nref` is the number of reads with the reference base and `Na` is the number of reads with "A" at the SNP position, `Nc` is the number of reads with "C" at the SNP position, etc. (default: 0.1)

`-o`    output prefix (default: "output")

`-t`    threads (default: 15)

###kmer_compare_groups.pl

###ksnp_matrix_filter.pl

###ksnp_matrix_to_diff_matrix.pl
