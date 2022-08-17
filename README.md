# snp_tools

A set of tools and add-ons for SNP analysis 

## Introduction 

Tools and add-ons primarily for working with the output of [kSNP](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0081760), an alignment-free SNP identifying software package. The basic functions of ksnp_matrix_filter.pl and ksnp_matrix_to_diff_matrix.pl will work on any aligned sequences in fasta format. 

## Requirements

All packages are written in Perl and require the bash shell. Prerequistes below must be available in your PATH to function. 

Required for all:
* bash
* perl

Required for some (ksnp_snp_confirm.pl, kmer_core.pl, kmer_compare.pl):
* [jellyfish](http://www.cbcb.umd.edu/software/jellyfish/) (>= v1.1.2).  This is likely already installed and in your path if you have kSNP.

Optional for some:
* gzip

## Installation

Simply download the folder containing the scripts and place it somewhere you'll be able to find it.

**IMPORTANT**: The `lib` directory must stay in the same directory as the parent scripts. Do not move or rename this directory or any of its contents.

## Programs

### ksnp_snp_confirm.pl

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

### kmer_core.pl

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


### kmer_compare.pl

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

### kmer_compare_groups.pl

```
kmer_compare_groups.pl <genome list file> <kmer_compare.raw.txt>
```

Outputs values and statistics of SNP differences based on different groupings.

**Genome list format:**  
First line should contain a listing of the different categories separated by commas. It should start with a '#' symbol.  
For example, if the categories your genomes can belong to inlcude species, clade, and outbreak, then your first line should look like this:

```
#species,clade,outbreak
```

All subsequent lines should start with the name of the genome (same names as in kmer_compare.raw.txt) and an integer for each group listed above to indicate which group the genome belongs to, all separated by commas. You can also give a '0' to indicate no group membership or unknown group membership.

For example, using the gropus outlined above, the next lines could be:

```
genA,1,1,1     #genome A belongs to species 1, clade 1, and outbreak 1
genB,1,2,0     #genome B belongs to species 1, clade 2, and no outbreak
genC,2,0,0     #genome C belongs to species 2, unknown clade or outbreak
```

Group memberships are independent of each other, so in the example above, if genomes B and C were given 1 and 2 respectively to indicate different strains, but both were given 2 for clade, they will still be grouped together as clade-mates.

**Options:**  
`-x`    if given, will only count intergroup comparisons as different. For example, a comparison between a group 1 and a group 2 strain would be counted as "different". A comparison between a group 1 strain and a group 0 strain would be counted as "different".  A comparison between a group 0 strain and another group 0 strain would not be counted as "same" or "different." Default: Comparisons to or between group 0 strains will not be counted at all.

`-p`    Output prefix. Default: 'output'

### ksnp_matrix_filter.pl

```
ksnp_matrix_filter.pl [options] <SNPs_all_matrix.fasta> > <filtered_matrix.fasta>
```

Filters the SNP matrix output by kSNP.

**Options:**  
`-m`    minimum fraction with locus, between 0 and 1.  1 means only loci present in all of the input genomes will be output.  0.5 means only loci present in at least half of the input genomes will be output, etc.  
(default: all loci will be output)

`-k`    keep locus order.  If given, any loci that do not meet the minimum fraction requirement given above will be replaced with "-" in every genome. Use this option if you have annotation information from kSNP.  
(default: non-matching loci will be removed)

`-i`    file with list of genomes to include in the ouput. List should be genome names, one per line. Fraction limits given by `-m` will be based on these genomes (i.e. if `-m` is 0.5 and you enter list of 4 genomes for `-i`, output will be loci present in at least 2 of these 4 genomes, regardless of how many total genomes were present in the `SNPs_all_matrix.fasta` file)  
(default: all genomes will be included)

`-e`    file with list of genomes to exclude from the output. List should be genome names, one per line. Fraction limits given by `-m` will based on those genomes not listed in this file.  
(default: no genomes will be excluded)

If lists are given to both `-i` and `-e`, fraction calculations will be based on and output will only include those genomes present in `-i`, excluding any genomes that are also present in `-e`.

`-a`    If selected, will calculate minimum locus fraction based on genomes given by `-i` and/or `-e`, but will output all input genomes, either omitting the filtered SNP loci or, if `-k` is given, substituting a "-" at each of the filtered positions.

`-g`    If selected, will keep all positions that match the criteria above, even if a locus has gaps in one or more genomes. This is automatically selected if the `-k` option is given above.  
(default behavior is to only output true SNP loci, omitting loci with the same base in every genome they were found in, but not found in every genome, i.e. loci with gaps '-')

*Below options need to be given only if you want to filter SNPs based on coordinates of one of the genomes input into kSNP*
        
`-r`    Coordinate restriction genome name

`-I`    coordinates within the restriction genome to include in the output. Start and stop coordinates must be in the last two columns, respectively

`-E`    coordinates within the restriction genome to exclude from the output. Start and stop coordinates must be in the last two columns, respectively

If both `-I` and `-E` are given, exclusions listed in `-E` will trump inclusions listed in `-I`.

`-s`    "SNPs_all" file from kSNP. This is REQUIRED if you want to filter by coordinates.

### ksnp_matrix_to_diff_matrix.pl

```
ksnp_matrix_to_diff_matrix.pl [options] <SNPs_all_matrix.fasta> 
```

Takes matrix.fasta file produced by kSNP and outputs a distance matrix of all SNP differences.

**Options:**

`-g`    If given, will count a difference between two genomes at a position where no base was present (i.e. genome 1 has "A", genome 2 has "-")  
(default: only base differences will be counted)

`-m`    minimum fraction of genomes in which a SNP position must be found in order to be counted (i.e. "1" for SNP position found in all genomes, "0.7" for SNP position found in at least 70% of the genomes)  
(default: any SNP position found in at least 2 of the input genomes will be included, i.e. "0")

`-i`    file with list of genomes to include in the ouput. List should be genome names, one per line. Fraction limits given by `-m` will be based on these genomes (i.e. if `-m` is 0.5 and you enter list of 4 genomes for `-i`, output will be loci present in at least 2 of these 4 genomes, regardless of how many total genomes were present in the `SNPs_all_matrix.fasta` file)  
(default: all genomes will be included)

`-e`    file with list of genomes to exclude from the output. List should be genome names, one per line. Fraction limits given by -m will based on those genomes not listed in this file.  
(default: no genomes will be excluded)
        
If lists are given to both `-i` and `-e`, fraction calculations will be based on and output will only include those genomes present in `-i`, excluding any genomes that are also present in `-e`.

`-a`    If kSNP was run with a reference genome and information about SNP annotation is desired, give the `SNP_annotations` file here.

*Below options apply only if a `SNP_annotations` file was given to option `-a`*

`-b`    `SNPs_all` file output by kSNP. If given, will output genomic position(s) of each SNP in the genomes being compared

`-p`    Only output annotations on proteins

`-s`    Only output non-synonymous SNPs (forces `-p`, obviously)
 
`-u`    If `-p` or `-s` are given, also output SNPs that were not found in the annotated genome
 
`-o`    Output prefix for annotation information (default "output")
  
## Contact

For questions or concerns, contact me: e-ozer@northwestern.edu
