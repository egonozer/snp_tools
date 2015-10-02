# snp_tools
A set of tools and add-ons for SNP analysis 

##Introduction 
Tools and add-ons primarily for working with the output of [kSNP](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0081760), an alignment-free SNP identifying software package. 

##Requirements
All packages are written in Perl and require the bash shell. Prerequistes below must be available in your PATH to function. 

Required for all:
* bash
* perl

Required for some (see documentation for individual programs for more information):
* [jellyfish](http://www.cbcb.umd.edu/software/jellyfish/) (>= v1.1.2).  This is likely already installed and in your path if you have kSNP.

##Installation
Simply download the folder containing the scripts and place it somewhere you'll be able to find it.
**IMPORTANT**: The _lib_ directory must stay in the same directory as the parent scripts. Do not move or rename this directory or any of its contents.

##ksnp_snp_confirm.pl

##kmer_core.pl

##kmer_compare.pl

##kmer_compare_groups.pl

##ksnp_matrix_filter.pl

##ksnp_matrix_to_diff_matrix.pl
