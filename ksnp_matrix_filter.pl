#!/usr/bin/perl

#    Copyright (C) 2015  Egon A. Ozer

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

my $version = "0.1";

use strict;
use warnings;

$|++;

my $usage = "
ksnp_matrix_filter.pl [options] <SNPs_all_matrix.fasta> > <filtered_matrix.fasta>

Filters the SNP matrix output by kSNP.

Options:
  -m    minimum fraction with locus, between 0 and 1.  1 means only loci present
        in all of the input genomes will be output.  0.5 means only loci present
        in at least half of the input genomes will be output, etc.
        (default: all loci will be output)
  -k    keep locus order.  If given, any loci that do not meet the minimum fraction
        requirement given above will be replaced with \"-\" in every genome.
        Use this option if you have annotation information from kSNP.
        (default: non-matching loci will be removed)
  -i    file with list of genomes to include in the ouput. List should be genome
        names, one per line.
        Fraction limits given by -m will be based on these genomes (i.e. if
        -m is 0.5 and you enter list of 4 genomes for -i, output will be loci
        present in at least 2 of these 4 genomes, regardless of how many total
        genomes were present in the SNPs_all_matrix.fasta file)
        (default: all genomes will be included)
  -e    file with list of genomes to exclude from the output. List should be
        genome names, one per line.
        Fraction limits given by -m will based on those genomes not listed
        in this file.
        (default: no genomes will be excluded)
        If lists are given to both -i and -e, fraction calculations will be
        based on and output will only include those genomes present in -i,
        excluding any genomes that are also present in -e.
  -a    If selected, will calculate minimum locus fraction based on genomes
        given by -i and/or -e, but will output all input genomes, either
        omiting the filtered SNP loci or, if -k is given, substituting a \"-\"
        at each of the filtered positions.
  -g    If selected, will keep all positions that match the criteria above,
        even if a locus has gaps in one or more genomes. This is automatically
        selected if the -k option is given above.
        (default behavior is to only output true SNP loci, omitting loci with
        the same base in every genome they were found in, but not found in
        every genome, i.e. loci with gaps '-')
        
        ** Below options need to be given only if you want to filter SNPs based
        on coordinates of one of the genomes input into kSNP **
        
  -r    Coordinate restriction genome name        
  -I    coordinates within the restriction genome to include in the output.
        Start and stop coordinates must be in the last two columns, respectively
  -E    coordinates within the restriction genome to exclude from the output
        Start and stop coordinates must be in the last two columns, respectively
            **If both -I and -E are given, exclusions listed in -E will trump
            inclusions listed in -I.
  -s    \"SNPs_all\" file from kSNP. This is REQUIRED if you want to filter by
        coordinates.

";

die $usage unless @ARGV;
print STDERR "\n";
use Getopt::Std;
our ($opt_m, $opt_k, $opt_i, $opt_e, $opt_a, $opt_r, $opt_I, $opt_E, $opt_s, $opt_g);
getopts('m:ki:e:ar:I:E:s:g');
my $frac    = $opt_m ? $opt_m : 0;
die "ERROR: -m must be a number between 0 and 1, inclusive\n" if $frac =~ m/D/ or $frac < 0 or $frac > 1;
my $in_file = $opt_i if $opt_i;
my $ex_file = $opt_e if $opt_e;

my $ref     = $opt_r if $opt_r;
my $in_crd  = $opt_I if $opt_I;
my $ex_crd  = $opt_E if $opt_E;
my $snp_all = $opt_s if $opt_s;
$opt_g = 1 if $opt_k;

if ($ref or $in_crd or $ex_crd or $snp_all){
    if ($in_crd or $ex_crd){
        die "ERROR: Need to give genome name and SNPs_all file if you give filtering coordinates\n" if (!$ref or !$snp_all);
    }
    if ($ref or $snp_all){
        die "ERROR: Need to give filtering coordinates if you give genome name and SNPs_all file\n" if (!$in_crd and !$ex_crd);
        die "ERROR: Need to give genome name, SNPs_all file, and filtering coordiates\n" if (!$ref or !$snp_all);
    }
}

my %incl;
my %excl;
if ($ex_file){
    open (my $in, "<$ex_file") or die "Can't open $ex_file: $!\n";
    while (my $line = <$in>){
        chomp $line;
        $excl{$line}++;
    }
    close ($in);
}
if ($in_file){
    open (my $in, "<$in_file") or die "Can't open $in_file: $!\n";
    while (my $line = <$in>){
        chomp $line;
        $incl{$line}++;
    }
    close ($in);
}

my %icrd;
my %ecrd;
if ($in_crd){
    open (my $in, "<$in_crd") or die "Can't open $in_crd: $!\n";
    while (my $line = <$in>){
        chomp $line;
        my @tmp = split(" ", $line);
        my ($start, $stop) = ($tmp[$#tmp - 1], $tmp[$#tmp]);
        next if $start =~ m/\D/ or $stop =~ m/\D/;
        ($start, $stop) = ($stop, $start) if ($start > $stop);
        for my $i ($start .. $stop){
            $icrd{$i}++;
        }
    }
    close ($in);
    print STDERR "WARNING: No coordinates found in $in_crd\n" if !%icrd;
}
if ($ex_crd){
    open (my $in, "<$ex_crd") or die "Can't open $ex_crd: $!\n";
    while (my $line = <$in>){
        chomp $line;
        my @tmp = split(" ", $line);
        my ($start, $stop) = ($tmp[$#tmp - 1], $tmp[$#tmp]);
        next if $start =~ m/\D/ or $stop =~ m/\D/;
        ($start, $stop) = ($stop, $start) if ($start > $stop);
        for my $i ($start .. $stop){
            $ecrd{$i}++;
        }
    }
    close ($in);
    print STDERR "WARNING: No coordinates found in $ex_crd\n" if !%ecrd;
}

my %binfo;
if ($snp_all){
    open (my $in, "<$snp_all") or die "ERROR: Can't open $snp_all: $!\n";
    while (my $line = <$in>){
        chomp $line;
        next if $line =~ m/^\s*$/; #skip blank lines
        $line =~ s/\s+$//;
        next unless $line =~ m/$ref$/;
        my ($locus, $kmer, $base, $pos, $gen) = split("\t", $line);
        $pos =~ s/\s+\D+\s*$//;
        next if $ecrd{$pos};
        if (%icrd){
            $binfo{$locus}++ if $icrd{$pos};
        } else {
            $binfo{$locus}++
        }
    }
    close ($in);
}

open (my $in, "<$ARGV[0]") or die "Can't open matrix file: $!\n";
my @gens;
my @snparrays;
my ($id, $seq);
my $count = 0;
my $laststring = 1;
my $excluded = 0;
my @gens_per_pos;
my %bases_per_pos;
my $dont_count;
my $kept_gens = 0;
while (my $line = <$in>){
    chomp $line;
    if ($line =~ m/^>/){
        if ($id){
            my @tmp = split('', $seq);
            my $tleng = scalar @tmp;
            if ($count > 0 and $tleng != $count){
                die "ERROR: number of positions not identical between records.\ntleng = $tleng, count = $count\n";
            }
            $count = $tleng;
            for my $k (0 .. $#tmp){
                $gens_per_pos[$k] = 0 if !$gens_per_pos[$k];
                unless ($dont_count){
                    $gens_per_pos[$k]++ if $tmp[$k] ne "-";
                    $bases_per_pos{$k}{uc($tmp[$k])}++;
                }
                #$bases_per_pos{$k}{uc($tmp[$k])}++;
            }
            push @snparrays, [@tmp];
        }
        $seq = "";
        $id = substr($line, 1);
        my $skip;
        if ($in_file){
            $skip = 1 if (!$incl{$id});
        }
        $skip = 1 if $excl{$id};
        my $string = "Reading $id"; #status update
        $string = "SKIPPING $id" if $skip; #status update
        my @blanks = " " x $laststring; #status update
        print STDERR "\r", join("", @blanks), "\r$string"; #status update
        $laststring = length $string; #status update
        if ($skip){
            $excluded ++;
            $dont_count = 1 if $opt_a;
            if ($opt_a){
                push @gens, $id;
            } else {
                $id = "";
            }
        } else {
            push @gens, $id;
            $dont_count = "";
            $kept_gens++;
        }
        next;
    }
    $line =~ s/\s*$//;
    $seq .= $line;
    next;
}
if ($id){
    my @tmp = split('', $seq);
    my $tleng = scalar @tmp;
    if ($count > 0 and $tleng != $count){
        die "ERROR: number of positions not identical between records.";
    }
    $count = $tleng;
    for my $k (0 .. $#tmp){
        $gens_per_pos[$k] = 0 if !$gens_per_pos[$k];
        unless ($dont_count){
            $gens_per_pos[$k]++ if $tmp[$k] ne "-";
            $bases_per_pos{$k}{uc($tmp[$k])}++;
        }
    }
    push @snparrays, [@tmp];
} else {
    die "ERROR: No fasta records found in file.\n" if !@gens;
}
$seq = "";

close ($in);
print STDERR "\nDone!\nTotal SNPs: $count\n";
print STDERR "Excluded $excluded genome(s)\n";

die "ERROR: Need to keep at least 2 genomes (only $kept_gens genomes kept).\n" if $kept_gens < 2;

my $gencount = $kept_gens;
my $mingen = roundup($gencount * $frac);
$mingen = 2 if $mingen < 2;
my $frac_count = 0;
my %filter_pos;
for my $i (0 .. $#gens_per_pos){
    if ($gens_per_pos[$i]){
        if ($gens_per_pos[$i] >= $mingen){
            $frac_count++;
        } else {
            $filter_pos{$i}++;
        }
    } else {
        $filter_pos{$i}++;
    }
}
print STDERR "Total loci in at least ",100 * $frac,"% of the input genomes ($mingen / $gencount): $frac_count\n";

my $true_snps = 0;
my $out_snps = 0;
for my $i (0 .. $#snparrays){
    my $id = $gens[$i];
    print ">$id\n";
    my @bases = @{$snparrays[$i]};
    for my $j (0 .. $#bases){
        my $skip;
        if ($ref and !$binfo{$j}){ #determines whether SNP positions are in included coordinate regions  
            $skip = 1;
        }
        if ($skip){
            if ($opt_k){
                print "-";
                $out_snps++ if $i == 0;
            }
            next;
        } else {
            if ($filter_pos{$j}){
                if ($opt_k){
                    print "-";
                    $out_snps++ if $i == 0;
                }
                next;
            } else {
                my %hash = %{$bases_per_pos{$j}};
                my $num_keys = keys %hash;
                if ($num_keys > 1){ #only print SNP positions, i.e. exclude loci with all the same base
                    if ($num_keys == 2 and $hash{"-"}){
                        if ($opt_g){
                            print "$bases[$j]";
                            $out_snps++ if $i == 0;
                        }
                    } else {
                        print "$bases[$j]";
                        $true_snps++ if $i == 0;
                        $out_snps++ if $i == 0;
                    }
                } else {
                    if ($opt_k){
                        print "-";
                        $out_snps++ if $i == 0;
                    }
                }
            }
        }
    }
    print "\n";
}
print STDERR "Output $true_snps loci with at least 1 SNP difference (total $out_snps loci output)\n\n";

#--------------------------------------------------------------------------------
sub roundup {
    my $n = shift;
    return(($n == int($n)) ? $n : int($n + 1));
}
