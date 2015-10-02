#!/usr/bin/perl

use warnings;
use strict;

$|++;

my $version = 0.3;
#this version (which is supposed to work with kmer_compare_v0.5) doesn't add both files to a hash, then loop through the hash to calculate snps and shared. Instead it reads in the kmers for the first genome and does the comparison during the reading of the second genome.
#also outputs results to a file rather than adding them to a hash. Should ultimately save memory and time.

my $usage = "$0 <total permutations> <path to temporary folder> <threads>\n";

die $usage unless @ARGV >= 2;

my ($tot_perm, $temp_path, $threads) = @ARGV;

my @bases = qw(A C G T);
#load the mers array (4-mer)
my @mers;
for my $i (0 .. $#bases){
    for my $j (0 .. $#bases){
        for my $k (0 .. $#bases){
            for my $l (0 .. $#bases){
                my $string = "$bases[$i]$bases[$j]$bases[$k]$bases[$l]";
                push @mers, $string;
            }
        }
    }
}

open (my $in, "< $temp_path/key.txt") or die "ERROR: Can't open $temp_path/key.txt: $!\n";
my @order;
my (@refs, @non_refs);
my $ref_run;
my %id_nums;
while (my $line = <$in>){
    chomp $line;
    my ($id, $num, $count, $types, $is_ref) = split("\t", $line);
    $id_nums{$num} = $id;
    if ($is_ref eq "-"){
        push @order, ([$id, $num, $count, $types]);
    } else {
        $ref_run = 1;
        push @refs, ([$id, $num, $count, $types]) if $is_ref eq "Y";
        push @non_refs, ([$id, $num, $count, $types]) if $is_ref eq "N";
    }
}
close ($in);
unlink("$temp_path/key.txt");

my ($ref_stop, $non_ref_start) = (0) x 2;
if ($ref_run){
    my $num_ref = scalar @refs;
    $ref_stop = $num_ref;
    $non_ref_start = $num_ref;
    @order = (@refs, @non_refs);
}

#my %results; #array order: $tot_snps, $tot_shared, $tot_a, $tot_b
open (my $results_out, "> $temp_path/temp_results.txt");
my $perm_count = 0;
my $last_line_leng = 1;
my $num_threads_running = 0;
for my $i (0 .. ($#order - 1)){
    last if ($ref_run and $i == $ref_stop);
    my ($a_id, $a_num, $a_count, $a_types) = @{$order[$i]};
    my @a_typ = split(",", $a_types);
    my $a_assm = 1 if $a_types =~ m/assm/;
    my $j_start = ($i + 1);
    $j_start = $non_ref_start if $ref_run;
    for my $j ($j_start .. $#order){
        my ($b_id, $b_num, $b_count, $b_types) = @{$order[$j]};
        my @b_typ = split(",", $b_types);
        my $b_assm = 1 if $b_types =~ m/assm/;
        $perm_count++;
        my @blank_line = (" ") x $last_line_leng;
        print STDERR "\r\t", join("", @blank_line);
        my $stat_line = "Comparing $a_id - $b_id ($perm_count / $tot_perm)";
        print STDERR "\r\t$stat_line";
        $last_line_leng = length($stat_line);
        print $results_out "$a_num $b_num 0 0 0 0\n"; #for the unlikely situation that two empty files will be compared
        #$results{$a_id}{$b_id} = ([0, 0, 0, 0]); #initialize the result array
        foreach my $mer (@mers){
            if ($num_threads_running == $threads){
                my $pid = wait;
                open (my $in, "< $temp_path/$pid.tmp.txt") or die "ERROR: Can't open $temp_path/$pid.tmp.txt";
                my $line = <$in>;
                close ($in);
                print $results_out "$line";
                unlink ("$temp_path/$pid.tmp.txt");
                $num_threads_running--;
            }
            if ($num_threads_running < $threads){
                my $pid = fork;
                if ($pid == 0){
                    my %compare;
                    my ($a_count, $b_count) = (0) x 2;
                    my ($shared_kmers, $snps) = (0) x 2;
                    foreach my $type (@a_typ){
                        next unless -e "$temp_path/$type\_kmers_$a_num\_$mer.txt";
                        open (my $f_in, "< $temp_path/$type\_kmers_$a_num\_$mer.txt") or die "ERROR: Can't open temporary file \"$temp_path/$type\_kmers_$a_num\_$mer.txt\": $!\n";
                        while (my $line = <$f_in>){
                            chomp $line;
                            my ($seq, $base) = split("\t", $line);
                            next if exists $compare{$seq}{"assm"}; #skip kmers in both reads and assemblies
                            $compare{$seq}{$type} = $base;
                            $a_count++;
                        }
                        close ($f_in);
                    }
                    my %skip;
                    foreach my $type (@b_typ){
                        next unless -e "$temp_path/$type\_kmers_$b_num\_$mer.txt";
                        open (my $f_in, "< $temp_path/$type\_kmers_$b_num\_$mer.txt") or die "ERROR: Can't open temporary file \"$temp_path/$type\_kmers_$b_num\_$mer.txt\": $!\n";
                        while (my $line = <$f_in>){
                            chomp $line;
                            my ($seq, $base) = split("\t", $line);
                            next if $skip{$seq};
                            $b_count++;
                            if ($compare{$seq}){
                                my %tmp = %{$compare{$seq}};
                                if ($tmp{"assm"}){
                                    my $a_base = $tmp{"assm"};
                                    $shared_kmers++;
                                    $snps++ if $base ne $a_base;
                                    $skip{$seq}++;
                                    next;
                                } elsif ($tmp{"read"}) { #if the sequence was only found in reads in the first genome
                                    my $a_base = $tmp{"read"};
                                    if ($type eq "assm"){ #if genome b has the sequence in an assembly
                                        $shared_kmers++;
                                        $snps++ if $base ne $a_base;
                                        $skip{$seq}++;
                                    } else {
                                        next if ($a_assm and $b_assm); #skip kmers that are only in reads when assemblies are available for both genomes
                                        $shared_kmers++;
                                        $snps++ if $base ne $a_base;
                                        #$skip{$seq}++;
                                    }
                                }
                            }
                        }
                        close ($f_in);
                    }
                    open (my $out, "> $temp_path/$$.tmp.txt") or die "ERROR: Can't open $temp_path/$$.tmp.txt: $!\n";
                    print $out "$a_num $b_num $snps $shared_kmers $a_count $b_count\n";
                    close ($out);
                    exit(0);
                }
                $num_threads_running++;
                next;
            }
        }
    }
}
while (my $pid = wait){
    last if $pid == -1;
    open (my $in, "< $temp_path/$pid.tmp.txt") or die "ERROR: Can't open $temp_path/$pid.tmp.txt";
    my $line = <$in>;
    close ($in);
    print $results_out "$line";
    unlink ("$temp_path/$pid.tmp.txt");
    $num_threads_running--;
}
print STDERR "\n";
close ($results_out);

#sort and output the results
print STDERR "\nSorting and dumping pairwise results...\n";
open (my $results_in, "sort -n -k1,1 -k2,2 $temp_path/temp_results.txt |") or die "ERROR: Could not sort and open $results_out: $!\n";
open (my $out, "> $temp_path/snp_counts.txt");
my ($last_a, $last_b);
my ($tot_snps, $tot_shared, $tot_a, $tot_b) = (0) x 4;
while (my $line = <$results_in>){
    chomp $line;
    my ($a_num, $b_num, $snps, $shared, $acount, $bcount) = split(" ", $line);
    if ($last_a){
        if ($last_a != $a_num or $last_b != $b_num){
            my @blank_line = (" ") x $last_line_leng;
            print STDERR "\r\t", join("", @blank_line);
            my ($a_id, $b_id) = ($id_nums{$last_a}, $id_nums{$last_b});
            my $stat_line = "$a_id - $b_id";
            print STDERR "\r\t$stat_line";
            $last_line_leng = length($stat_line);
            my $dist = $tot_snps;
            $dist += ((($tot_a + $tot_b) - ($tot_shared * 2)) * 2); #indel mismatches (where a kmer is present in one strain but not the other), are given a distance of 2
            print $out "$a_id\t$b_id\t$tot_snps\t$tot_shared\t$dist\n";
            ($tot_snps, $tot_shared, $tot_a, $tot_b) = (0) x 4;
        }
    }
    ($last_a, $last_b) = ($a_num, $b_num);
    $tot_snps += $snps;
    $tot_shared += $shared;
    $tot_a += $acount;
    $tot_b += $bcount;
}
if ($last_a){
    my @blank_line = (" ") x $last_line_leng;
    print STDERR "\r\t", join("", @blank_line);
    my ($a_id, $b_id) = ($id_nums{$last_a}, $id_nums{$last_b});
    my $stat_line = "$a_id - $b_id";
    print STDERR "\r\t$stat_line";
    $last_line_leng = length($stat_line);
    my $dist = $tot_snps;
    $dist += ((($tot_a + $tot_b) - ($tot_shared * 2)) * 2); #indel mismatches (where a kmer is present in one strain but not the other), are given a distance of 2
    print $out "$a_id\t$b_id\t$tot_snps\t$tot_shared\t$dist\n";
    ($tot_snps, $tot_shared, $tot_a, $tot_b) = (0) x 4;
} else {
    die "ERROR: kmer_compare_pairwise_mt produced no results\n";
}
print STDERR "\n";
close ($results_in);
#unlink("$temp_path/temp_results.txt");
