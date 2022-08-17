#!/usr/bin/perl

my $version = "0.3.3";

#changes from v0.3.2 to v0.3.3
#Add option to output positions not matching filters as "N" instead of "-"

#changes from v0.3.1 to v0.3.2
#count any ambiguous base as a gap, i.e. only the bases A,C,T, or G will be counted as bases. N's or other ambiguity codes will be counted as gaps (same as "-"). This reflects the use of the program for thinning alignments other than just those produced by kSNP.
#option to output non-core sites containing SNPs (for ClonalFrameML)

#changes from v0.2 to v0.3
#speed improvements in outputting sequences by not cycling through the entire sequence just to output filtered SNP positions.
#also improve speed by not splitting the sequence into an array of characters, then cycling through the array. Instead just remove the first base of the seq string until there is no more seq string left.
#more speed improvement: changed hash of hashes %bases_per_pos to array of hashes @bases_per_pos_a. This almost halves the amount of time required to add a base count to the internal hash.
#change with v0.3.1: Added -x option to replace filtered sites with a reference base instead of gap character -
#change with v0.3.1: Option to output phylip formatted sequence as well

#changes from v0.1 to v0.2
#stop storing all sequences in memory while determining which columns to keep. That's a big 'ol memory hog. Run through the input file twice instead.

use strict;
use warnings;

$|++;

my $usage = "
ksnp_matrix_filter.pl [options] <SNPs_all_matrix.fasta> > <filtered_matrix.fasta>
version: $version

Filters the SNP matrix output by kSNP.

Options:
  -m    minimum fraction with locus, between 0 and 1.  1 means only loci present
        in all of the input genomes will be output.  0.5 means only loci present
        in at least half of the input genomes will be output, etc.
        (default: all loci will be output)
  -k    keep locus order.  If given, any loci that do not meet the minimum fraction
        requirement given above will be replaced with \"-\" in every genome. Set
        the -l option below if you would rather replace with another character.
        Use this option if you have annotation information from kSNP.
        (default: non-matching loci will be removed)
  -l    Character with which to replace filtered positions if -k is given. Can be
        set if you want to use \"N\" instead of a gap character.
        (default: \"-\")
  -x    reference sequence ID. If this is value is given with -k above, then
        instead of replacing filtered sites with \"-\", the sites will instead
        be replaced with the base in the reference sequence. Setting -x will
        automatically trigger -k.
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
  -p    file name of phylip formatted matrix to output.
        (default: no phylip matrix will be output)
  -n    file name of list of non-core positions to be output.
        (default: no list of non-core positions will be output)
        
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
our ($opt_m, $opt_k, $opt_l, $opt_i, $opt_e, $opt_a, $opt_r, $opt_I, $opt_E, $opt_s, $opt_g, $opt_x, $opt_p, $opt_n);
getopts('m:kl:i:e:ar:I:E:s:gx:p:n:');
my $frac    = $opt_m ? $opt_m : 0;
die "ERROR: -m must be a number between 0 and 1, inclusive\n" if $frac =~ m/D/ or $frac < 0 or $frac > 1;
my $in_file = $opt_i if $opt_i;
my $ex_file = $opt_e if $opt_e;

my $gap_c   = $opt_l ? $opt_l : "-";
die "The length of the replacement character (-l) can only be one. $opt_l is too long\n" unless (length($gap_c) == 1);

my $phy     = $opt_p if $opt_p;
my $nc      = $opt_n if $opt_n;
my $ref     = $opt_r if $opt_r;
my $in_crd  = $opt_I if $opt_I;
my $ex_crd  = $opt_E if $opt_E;
my $snp_all = $opt_s if $opt_s;
my $refgen  = $opt_x if $opt_x;
$opt_k = 1 if $opt_x;
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

my $file = $ARGV[0];

open (my $in, "<$file") or die "Can't open matrix file: $!\n";
my %gens_to_keep;
my $gen_count = 0;
my ($id, $seq);
my $count = 0;
my $laststring = 1;
my $excluded = 0;
my @gens_per_pos;
my %bases_per_pos;
my @bases_per_pos_a;
my $dont_count;
my $kept_gens = 0;
my @ref_bases;
while (my $line = <$in>){
    chomp $line;
    if ($line =~ m/^>/){
        if ($id){
            my $tleng = length($seq);
            if ($count > 0 and $tleng != $count){
                die "ERROR: number of positions not identical between records.\ntleng = $tleng, count = $count\n";
            }
            $count = $tleng;
            $seq = uc($seq);
            my $scount = $tleng - 1;
            while ($seq){
                my $base = chop $seq;
                $gens_per_pos[$scount] = 0 if !$gens_per_pos[$scount];
                unless ($dont_count){
                    my $tmp_base = $base;
                    if ($base =~ m/[ACGT]/){
                        $gens_per_pos[$scount]++
                    } else {
                        $tmp_base = "-"; #classify any non-ACTG bases as gaps
                    }
                    $bases_per_pos_a[$scount]{$tmp_base}++;
                }
                $scount--;
                if ($refgen and $id eq $refgen){
                    unshift @ref_bases, $base;
                }
            }
        }
        $seq = "";
        $id = substr($line, 1);
        $gen_count++;
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
                $gens_to_keep{$gen_count}++;
            } else {
                $id = "";
            }
        } else {
            $gens_to_keep{$gen_count}++;
            $dont_count = "";
            $kept_gens++;
        }
        next;
    }
    $line =~ s/\s//g;
    $seq .= $line;
    next;
}
if ($id){
    my $tleng = length($seq);
    if ($count > 0 and $tleng != $count){
        die "ERROR: number of positions not identical between records.";
    }
    $count = $tleng;
    
    #speed test
    $seq = uc($seq);
    my $scount = $tleng - 1;
    while ($seq){
        my $base = chop $seq;
        $gens_per_pos[$scount] = 0 if !$gens_per_pos[$scount];
        unless ($dont_count){
            my $tmp_base = $base;
            if ($base =~ m/[ACGT]/){
                $gens_per_pos[$scount]++
            } else {
                $tmp_base = "-"; #classify any non-ACTG bases as gaps
            }
            $bases_per_pos_a[$scount]{$tmp_base}++;
        }
        $scount--;
        if ($refgen and $id eq $refgen){
            unshift @ref_bases, $base;
        }
    }
} else {
    die "ERROR: No fasta records found in file.\n" if !%gens_to_keep;
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
my @kpos;
my $non_snp_loci = 0;
my $nc_out;
if ($nc){
    open ($nc_out, ">$nc") or die "ERROR: Can't open $nc for writing: $!\n";
}

for my $i (0 .. $#gens_per_pos){
    if ($gens_per_pos[$i]){
        if ($gens_per_pos[$i] >= $mingen){
            $frac_count++;
            next if $ref and !$binfo{$i};
            my %hash = %{$bases_per_pos_a[$i]};
            my $num_keys = keys %hash;
            if ($num_keys > 1){ #this will exclude loci that have the same base in all genomes
                if ($num_keys == 2 and $hash{"-"}){
                    if ($opt_g){
                        $filter_pos{$i}++;
                        push @kpos, $i;
                        $non_snp_loci++;
                        print $nc_out $i + 1, "\n" if $nc;
                    }
                } else {
                    $filter_pos{$i}++;
                    push @kpos, $i;
                    print $nc_out $i + 1, "\n" if $nc;
                }
            }
        }
    }
}
close ($nc_out) if $nc_out;
print STDERR "Total loci in at least ",100 * $frac,"% of the input genomes ($mingen / $gencount): $frac_count\n";

my $phy_out;
if ($phy){
    open ($phy_out, ">$phy") or die "ERROR: Can't open $phy: $!\n";
    #kept_gens
    my $num_pos = scalar @kpos;
    $num_pos = $count if $opt_k;
    print $phy_out "$kept_gens   $num_pos\n";
}

open (my $in2, "<$file"); #read through the input file again instead of saving the sequences in memory from the first read-through.
my $true_snps = 0;
($id, $seq) = ("") x 2;
$gen_count = 0;
my $out_gens = 0;
my $keep;
while (my $line = <$in2>){
    chomp $line;
    if ($line =~ m/^>/){
        if ($id and $keep){
            print ">$id\n";
            printf $phy_out "%-90.90s",$id if $phy_out;
            $out_gens++;
            print STDERR "\rOutputting SNPs from genome $out_gens/$kept_gens";
            output_seq($seq);
        }
        $keep = "";
        $seq = "";
        $id = substr($line, 1);
        $gen_count++;
        $keep = 1 if ($gens_to_keep{$gen_count});
        next;
    }
    $line =~ s/\s//g;
    $seq .= $line;
}
close ($in2);
if ($id and $keep){
    print ">$id\n";
    printf $phy_out "%-90.90s",$id if $phy_out;
    $out_gens++;
    print STDERR "\rOutputting SNPs from genome $out_gens/$kept_gens";
    output_seq($seq);
}
print STDERR "\n";

close ($phy_out) if $phy_out;

my $out_snps = $true_snps + $non_snp_loci;
print STDERR "Output $true_snps loci with at least 1 SNP difference (total $out_snps loci output)\n\n";

#--------------------------------------------------------------------------------
sub roundup {
    my $n = shift;
    return(($n == int($n)) ? $n : int($n + 1));
}

sub output_seq {
    my $seq = shift;
    my @bases = split('', $seq);
    my $get_nums;
    $get_nums = 1 if ($true_snps == 0); #only need to calculate for the first run-through.
    
    if ($opt_k){ #with -k, all positions will be output so we'll just cycle through all of the @bases array
        for my $j (0 .. $#bases){
            if ($ref and !$binfo{$j}){ #determines whether SNP positions are in included coordinate regions
                my $val = $gap_c;
                $val = $ref_bases[$j] if @ref_bases;
                print "$val";
                print $phy_out "$val" if $phy_out;
                $non_snp_loci++ if $get_nums;
                next;
            } else {
                if (!$filter_pos{$j}){
                    my $val = $gap_c;
                    $val = $ref_bases[$j] if @ref_bases;
                    print "$val";
                    print $phy_out "$val" if $phy_out;
                    $non_snp_loci++ if $get_nums;
                    next;
                } else {
                    my %hash = %{$bases_per_pos_a[$j]};
                    my $num_keys = keys %hash;
                    if ($num_keys > 1){ #only print SNP positions, i.e. exclude loci with all the same base
                        if ($num_keys == 2 and $hash{$gap_c}){
                            if ($opt_g){ ## This if statement is unnecessary. opt_g is automatically 1 if opt_k is selected
                                print "$bases[$j]";
                                print $phy_out "$bases[$j]" if $phy_out;
                            }
                        } else {
                            print "$bases[$j]";
                            print $phy_out "$bases[$j]" if $phy_out;
                            $true_snps++ if $get_nums;
                        }
                    } else {
                        my $val = $gap_c;
                        $val = $ref_bases[$j] if @ref_bases;
                        print "$val";
                        print $phy_out "$val" if $phy_out;
                        $non_snp_loci++ if $get_nums;
                    }
                }
            }
        }
    } else { #otherwise we'll just cycle through the @kpos array of SNP positions that match the filter criteria which is likely shorter and should go faster
        foreach my $j (@kpos){
            print "$bases[$j]";
            print $phy_out "$bases[$j]" if $phy_out;
            $true_snps++ if $get_nums;
        }
    }
    print "\n";
    print $phy_out "\n" if $phy_out;
    return();
}
