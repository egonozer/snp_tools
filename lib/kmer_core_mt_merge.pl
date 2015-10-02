#!/usr/bin/perl

use warnings;
use strict;

$|++;

my $version = 0.2;
#changes from v0.1:
#updated to work with kmer_core_mt_v0.3.pl
#stop checking reverse-complemented kmers (not needed with -C setting in jellyfish)

my $usage = "perl $0 threads a_num count\n";

die $usage unless @ARGV >= 3;

my $threads = $ARGV[0];
my $a_num = $ARGV[1];
my $count = $ARGV[2];

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

my $num_threads_running = 0;
foreach my $group (@mers){
    print STDERR "\r\tGrouping $group-mers";
    unless (-e "tmp_kmers_$group.txt"){
        print STDERR "\n\t\t $group-mers kmers file not found, skipping\n";
        next;
    }
    if ($num_threads_running == $threads){
        my $pid = wait;
        my $status = $?;
        die "\nERROR: Process $pid returned a status of $status\n" if ($status != 0);
        my $sys_status = system("cat tmp_out.$pid.txt >> tmp_merge_out.txt");
        die "\nERROR: \"cat tmp_out.$pid.txt >> tmp_merge_out.txt\" returned a status of $sys_status\n" if ($sys_status != 0);
        unlink ("tmp_out.$pid.txt");
        $num_threads_running--;
    }
    my $pid = fork;
    if ($pid == 0){
        my $grp = $group;
        open (my $in, "< tmp_kmers_$grp.txt") or exit("ERROR: Temporary file tmp_kmers_$group.txt not opened: $!\n");
        my %counts;
        while (my $line = <$in>){
            chomp $line;
            my ($seq, $base, $gen) = split("\t", $line);
            $counts{$seq}{$base}++;
            push @{$counts{$seq}{"gen"}}, $gen;
        }
        close ($in);
        my %read_hash;
        if (-e "tmp_kmers_reads_$group.txt"){
            open (my $in, "< tmp_kmers_reads_$grp.txt") or exit("ERROR: Temporary file tmp_kmers_reads_$grp.txt not opened: $!\n");
            while (my $line = <$in>){
                chomp $line;
                my ($seq, $base, $gen) = split("\t", $line);
                $read_hash{$seq}{$gen} = $base;
            }
            close ($in);
        }
        open (my $out, "> tmp_out.$$.txt");
        foreach my $seq (sort {$a cmp $b} keys %counts){
            if (%read_hash){ #check to see if reads were given
                my $added; #avoid resorting the genome array if no changes were made. Saves time.
                if (exists $read_hash{$seq}){
                    my %lhash;
                    foreach my $gen (@{$counts{$seq}{"gen"}}){
                        $lhash{$gen}++;
                    }
                    for my $i (0 .. ($count - 1)){
                        next if $lhash{$i}; # skip genomes that had the kmer in the assembly
                        if (exists $read_hash{$seq}{$i}){
                            $added = 1;
                            my $base;
                            $base = $read_hash{$seq}{$i};
                            $counts{$seq}{$base}++;
                            push @{$counts{$seq}{"gen"}}, $i;
                        }
                    }
                }
                @{$counts{$seq}{"gen"}} = sort{$a <=> $b} @{$counts{$seq}{"gen"}} if $added;
            }
            ## determine the core kmers
            my $sum = 0;
            my @b_counts;
            for my $i (0 .. $#bases){
                my $base = $bases[$i];
                my $val = 0;
                $val = $counts{$seq}{$base} if $counts{$seq}{$base};
                $sum += $val;
                push @b_counts, ([$val, $base]);
            }
            my $out_seq = "x";
            if ($sum >= $a_num){
                @b_counts = sort{$b->[0] <=> $a->[0]} @b_counts;
                my $base = $b_counts[0][1];
                $out_seq = $seq;
                $out_seq =~ s/\./$base/;
            }
            my $gens = join(",", @{$counts{$seq}{"gen"}});
            my $gen_num = scalar @{$counts{$seq}{"gen"}};
            print $out "$out_seq\t$gens\t$gen_num\n";
        }
        close ($out);
        exit (0);
    }
    $num_threads_running++;
}
while (my $pid = wait){
    last if $pid < 0;
    my $status = $?;
    die "\nERROR: Process $pid returned a status of $status\n" if ($status != 0);
    my $sys_status = system("cat tmp_out.$pid.txt >> tmp_merge_out.txt");
    die "\nERROR: \"cat tmp_out.$pid.txt >> tmp_merge_out.txt\" returned a status of $sys_status\n" if ($sys_status != 0);
    unlink ("tmp_out.$pid.txt");
    $num_threads_running--;
}
