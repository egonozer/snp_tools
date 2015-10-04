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
ksnp_matrix_to_diff_matrix.pl [options] <SNPs_all_matrix.fasta> 

Takes matrix.fasta file produced by kSNP and outputs a distance matrix of all
SNP differences.

Options:
  -g    If given, will count a difference between two genomes at a position where
         no base was present (i.e. genome 1 has \"A\", genome two has \"-\")
        (default: only base differences will be counted)
  -m    minimum fraction of genomes in which a SNP position must be found in
        order to be counted (i.e. \"1\" for SNP position found in all genomes,
        \"0.7\" for SNP position found in at least 70% of the genomes)
        (default: any SNP position found in at least 2 of the input genomes
        will be included, i.e. \"0\")
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
  

  -a    If kSNP was run with a reference genome and information about SNP
        annotation is desired, give the \"SNP_annotations\" file here.
  -b    \"SNPs_all\" file output by kSNP. If given, will output genomic
        position(s) of each SNP in the genomes being compared
  -p    Only output annotations on proteins
  -s    Only output non-synonymous SNPs (forces -p, obviously)
  -u    If -p or -s are given, also output SNPs that were not found in the
        annotated genome
  -o    Output prefix for annotation information (default \"output\")

";

die $usage unless @ARGV;

# command line processing
use Getopt::Std;
our ($opt_g, $opt_a, $opt_p, $opt_s, $opt_u, $opt_o, $opt_m, $opt_i, $opt_e, $opt_b);
getopts('ga:psuo:m:i:e:b:');
my $pref = $opt_o ? $opt_o: "output";
my $count_gaps;
$count_gaps = 1 if $opt_g;
my $afile = $opt_a if $opt_a;
my $bfile = $opt_b if $opt_b;
my $frac = $opt_m ? $opt_m : 0;
die "ERROR: -f must be a number between 0 and 1, inclusive\n" if $frac =~ m/D/ or $frac < 0 or $frac > 1;
my $in_file = $opt_i if $opt_i;
my $ex_file = $opt_e if $opt_e;

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

open (my $in, "<$ARGV[0]") or die "Can't open matrix file: $!\n";
my %gen_ids;
my @gens;
my @snparrays;
my %diffs;
my ($id, $seq);
my $count = 0;
my $excluded = 0;
my $laststring = 1;
my @gens_per_pos;
my %bases_per_pos;
while (my $line = <$in>){
    chomp $line;
    if ($line =~ m/^>/){
        if ($id){
            my @tmp = split('', $seq);
            my $tleng = scalar @tmp;
            if ($count > 0 and $tleng != $count){
                die "ERROR: number of positions not identical between records.\n";
            }
            $count = $tleng;
            for my $k (0 .. $#tmp){
                $gens_per_pos[$k] = 0 if !$gens_per_pos[$k];
                $gens_per_pos[$k]++ if $tmp[$k] ne "-";
                $bases_per_pos{$k}{uc($tmp[$k])}++;
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
            $id = "";
        } else {
            push @gens, $id;
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
        $gens_per_pos[$k]++ if $tmp[$k] ne "-";
        $bases_per_pos{$k}{uc($tmp[$k])}++;
    }
    push @snparrays, [@tmp];
} else {
    die "ERROR: No fasta records found in file.\n" if !@gens;
}
$seq = "";
close ($in);
print STDERR "\nDone!\nTotal SNPs: $count\n";
print STDERR "Excluded $excluded genome(s)\n";

my %ainfo;
if ($afile){
    open (my $in, "<", $afile) or die "ERROR: Can't open $afile: $!\n";
    while (my $line = <$in>){
        chomp $line;
        next if $line =~ m/^LocusNum/;
        my @tmp = split("\t", $line);
        my $id = $tmp[0];
        my @tmp2;
        for my $i (2,3,4,5,6,8,11,14){
            push @tmp2, $tmp[$i];
            #if ($i == 2 and !$tmp[$i]){
            #    push @tmp2, 0;
            #}
        }
        push @{$ainfo{$id}}, [@tmp2];
    }
    close ($in);
}
my %binfo;
if ($bfile){
    open (my $in, "<$bfile") or die "ERROR: Can't open $bfile: $!\n";
    while (my $line = <$in>){
        chomp $line;
        next if $line =~ m/^\s*$/; #skip blank lines
        my ($locus, $kmer, $base, $pos, $gen) = split("\t", $line);
        $pos =~ s/\s+\D+\s*$//;
        push @{$binfo{$locus}{$gen}}, $pos;
    }
    close ($in);
}

my $gencount = scalar @gens;
my $mingen = roundup($gencount * $frac);
$mingen = 2 if $mingen < 2;
my $frac_count = 0;
for my $i (0 .. $#gens_per_pos){
    if ($gens_per_pos[$i]){
        $frac_count++ if $gens_per_pos[$i] >= $mingen;
    }
}

print STDERR "Total loci in at least ",100 * $frac,"% of the input genomes ($mingen / $gencount): $frac_count\n";

my %poscount;
my %diffpos;
for my $i (0 .. ($#gens - 1)){
    my $refgen = $gens[$i];
    my @refsnps = @{$snparrays[$i]};
    for my $j ($i+1 .. $#gens){
        my $qrygen = $gens[$j];
        my @qrysnps = @{$snparrays[$j]};
        my $string = "Comparing $refgen - $qrygen"; #status update
        my @blanks = " " x $laststring; #status update
        print STDERR "\r", join("", @blanks), "\r$string"; #status update
        $laststring = length $string; #status update
        for my $k (0 .. ($count - 1)){
            next if $gens_per_pos[$k] < $mingen;
            my $rsnp = $refsnps[$k];
            my $qsnp = $qrysnps[$k];
            if (!$count_gaps){
                next if ($rsnp eq "-" or $qsnp eq "-");
            }
            $poscount{$refgen}{$qrygen}++ unless keys %{$bases_per_pos{$k}} == 1;
            if ($rsnp ne $qsnp){
                $diffs{$refgen}{$qrygen}++;
                push @{$diffpos{$refgen}{$qrygen}}, $k if $afile;
            }
        }
    }
}
print STDERR "\nDone!\nOutputting...\n";

@gens = sort{$a cmp $b} @gens;

#open (my $raw, "> $pref.matrix.raw.tab");
#open (my $avg, "> $pref.matrix.norm_to_tot_snp.tab");


print "totalSNPs\t$count\n";
print "totalSNPs_$frac\t$frac_count\n\n";
print " \t", join("\t", @gens), "\n";
my $tpos_string = " \t" . join("\t", @gens) . "\n";
my $ratio_string = " \t" . join("\t", @gens) . "\n";
#print $avg " \t", join("\t", @gens), "\n";
my $aout;
if ($afile){
    open ($aout, "> $pref\_annotations.txt");
    print $aout "Genome1\tGenome2\tGenome1Base\tGenome2Base";
    print $aout "\tGenome1Pos\tGenome2Pos" if $bfile;
    print $aout "\tSNPLocusNum\tNonSynonymous\tAnnotationType\tAminoAcids\tCodons\tProteinPos\tRefGenPos\tLocusID\tProduct\n";
}
for my $i (0 .. $#gens){
    my $refgen = $gens[$i];
    my @refsnps = @{$snparrays[$i]};
    print "$refgen";
    $tpos_string .= "$refgen";
    $ratio_string .= "$refgen";
    #print $avg "$refgen";
    for my $j (0 .. $#gens){
        my $qrygen = $gens[$j];
        my $val = 0;
        $val = $diffs{$refgen}{$qrygen} if $diffs{$refgen}{$qrygen};
        $val = $diffs{$qrygen}{$refgen} if $diffs{$qrygen}{$refgen};
        my $tpos = 0;
        $tpos = $poscount{$refgen}{$qrygen} if $poscount{$refgen}{$qrygen};
        $tpos = $poscount{$qrygen}{$refgen} if $poscount{$qrygen}{$refgen};
        #print STDERR "$refgen\t$val\t$tpos\n";
        if ($afile){ #get and output annotation information
            my @dpos;
            @dpos = @{$diffpos{$refgen}{$qrygen}} if $diffpos{$refgen}{$qrygen};
            @dpos = @{$diffpos{$qrygen}{$refgen}} if $diffpos{$qrygen}{$refgen};
            my $printedsomething;
            if (@dpos){
                my @qrysnps = @{$snparrays[$j]};
                for my $k (0 .. $#dpos){
                    my $pos = $dpos[$k];
                    my $rsnp = $refsnps[$pos];
                    my $qsnp = $qrysnps[$pos];
                    my @deets = @{$ainfo{$pos}};
                    for my $l (0 .. $#deets){
                        my @tmp = @{$deets[$l]};
                        #print STDERR "***$k***", join("***", @tmp), "\n";
                        unless ($tmp[1] eq "NotInAnnotatedGenome" and $opt_u){
                            next if ($tmp[1] ne "OnProtein" and $opt_p);
                            next if ($tmp[0] != 1 and $opt_s);
                        }
                        my ($rpos, $qpos) = ("?") x 2;
                        if ($bfile){
                            if ($binfo{$pos}{$refgen}){
                                $rpos = join(",", @{$binfo{$pos}{$refgen}})
                            }
                            if ($binfo{$pos}{$qrygen}){
                                $qpos = join(",", @{$binfo{$pos}{$qrygen}})
                            }
                        }
                        print $aout "$refgen\t$qrygen\t$rsnp\t$qsnp\t$rpos\t$qpos\t$pos";
                        my $first;
                        while (@tmp){
                            my $val = shift @tmp;
                            if (!$first){
                                print $aout defined $val ? "\t$val" : 0;
                                $first = 1;
                            } else {
                                print $aout "\t$val" if $val;
                            }
                        }
                        print $aout "\n";
                        $printedsomething = 1;
                    }
                }
            }
            print $aout "\n" if $printedsomething;
        }
        my $ratio = 0;
        if ($tpos > 0){
            $ratio = sprintf("%.04f", $val/$tpos);
        }
        print "\t$val";
        $tpos_string .= "\t$tpos";
        $ratio_string .= "\t$ratio";
        #my $average = sprintf("%.5f", ($val / $count));
        #print $avg "\t$average";
        if ($qrygen eq $refgen){
            print "\n";
            $tpos_string .= "\n";
            $ratio_string .= "\n";
            #print $avg "\n";
            last;
        }
    }
}

print "\nTotal shared positions\n$tpos_string";
print "\nRatio SNPs/total\n$ratio_string";

if ($afile){
    close ($aout);
}

#--------------------------------------------------------------------------------
sub roundup {
    my $n = shift;
    return(($n == int($n)) ? $n : int($n + 1));
}
