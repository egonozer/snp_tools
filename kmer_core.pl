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

use warnings;
use strict;
use Cwd 'abs_path';
use File::Path 'rmtree';
use File::Basename;

$|++;

my $version = "0.3";
## Changes from v0.2
## Removed reverse kmer checking (not needed if -C option is given to jellyfish)
## Do jellyfishing of individual genomes first, then process kmers in parallel to try to speed thing sup.


## Changes from v0.1
## Cleaned up code
## Added some debugging tests
## Added some secret support for pre-processed, gzipped jellyfish files (.jelly.txt.gz). 



my $dirname = dirname(__FILE__);

my $usage = "
kmer_core.pl

DESCRIPTION: Determines the set of kmers that are present in a given subset of
a group of input genomes. This can be used as input to kSNP.

PARAMETERS:
Required:
  -f    text file with names and paths to input genomes, separated by any white
        space. Each file listed should be in fasta format. Multi-contig genome
        files are permitted, but each sequence file will be assumed to contain
        sequence for only one genome.
        OPTIONAL:
        After the genome name and the path to the genome sequence file, paths
        to one or more sequencing read files (in fastq or fasta format) can be
        included. These will be used to find kmers with conflicting alleles
        that were collapsed in assembly and/or kmers present in reads, but for
        which assembly failed in one or more strains. New kmers will not be
        identified from read sets. They will only be used to confirm kmers
        found in at least one of the input genomic sequence files. Gzipped
        read files (with .gz extension) can be given.
        FORMAT:
        <genome_id> </path/to/sequence.fasta> </path/to/reads.fastq OPTIONAL> </path/to/reads2.fastq.gz OPTIONAL> etc.
        EXAMPLE:
        PAO1 seqs/PAO1_contigs.fasta
        PA14 seqs/PA14_contigs.fasta reads/PA14_1.fastq.gz reads/PA14_2.fastq.gz
        PA7 seqs/PA7_contigs.fasta reads/PA7_1.fastq.gz
        
Optional:
  -a    minimum percentage of the total input genomes in which a kmer must be
        found to be considered core
        (default: 100)
  -k    kmer size. Maximum value is 31.
        (default: 31)
  -n    number of N's with which to separate kmers in the core kmers sequence
        output file
        (default: 1)
  -t    number of CPUs / threads
        (default: 15)
  -j    path to jellyfish
        (default: searches for jellyfish in PATH)
  -o    output files prefix
        (default: output)
  -m    [in reads file(s)] minimum number of times a kmer must be present to be
        counted.
        (default: 5)
  -e    [in reads file(s)] maximum allowable error rate between number of reads
        with reference base and number of reads with discrepant base. This is
        calcuated as:
        1 - (Nref / (Na + Nc + Ng + Nt))
        where Nref is the number of reads with the reference base and Na is the
        number of reads with \"A\" at the SNP position, Nc is the number of reads
        with \"C\" at the SNP position, etc.
        (default: 0.1)

";

use Getopt::Std;
use vars qw( $opt_f $opt_a $opt_k $opt_n $opt_t $opt_j $opt_o $opt_m $opt_e );
getopts('f:a:k:n:t:j:o:m:e:');

die $usage unless $opt_f;

my $a_pct   = $opt_a ? $opt_a : 100;
my $kmer    = $opt_k ? $opt_k : 31;
my $n_num   = $opt_n ? $opt_n : 1;
my $threads = $opt_t ? $opt_t : 15;
my $j_path  = $opt_j if $opt_j;
my $pref    = $opt_o ? $opt_o : "output";
my $min_k   = $opt_m ? $opt_m : 5;
my $err     = $opt_e ? $opt_e : .1;

die "ERROR: kmer size (-k) must be an odd number\n" unless ($kmer % 2 == 1);

my $armleng = int($kmer / 2);

my $start_time = time;

#read in file of files
my @files;
open (my $f_in, "<$opt_f") or die "ERROR: Can't open $opt_f: $!\n";
while (my $line = <$f_in>){
    chomp $line;
    next if $line =~ m/^\s*$/;
    my @tmp = split(" ", $line);
    push @files, [@tmp];
}
close ($f_in);
my $nog = scalar @files;

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

open (my $stats, "> $pref.statistics.txt");
print $stats "version: $version\n";
print $stats "-f $opt_f\n";
print $stats "-a $a_pct\n";
print $stats "-k $kmer\n";
print $stats "-n $n_num\n";
print $stats "-t $threads\n";
print $stats "-m $min_k\n";
print $stats "-e $err\n\n";
print $stats "genome\tgenome_num\tkmers_found\treads?\tconflicts_removed\n";

#set up Temporary Files folder
mkdir "TemporaryFiles" or die "Can't make TemporaryFiles folder: $!\n";
my $temp_path = abs_path("TemporaryFiles");
print STDERR "temp_path: $temp_path\n";

my $script_path = abs_path($0);
$script_path = dirname($script_path);
print STDERR "script_path: $script_path\n";

#produce kmer files from genomic sequences and reads
my $count = 0;
for my $i (0 .. $#files){
    my @tmp = @{$files[$i]};
    my $id = shift @tmp;
    my $genfile = shift @tmp;
    my $is_jelly_gz = 1 if $genfile =~ m/.jelly.txt.gz$/;
    my $has_reads = "N";
    $has_reads = "Y" if @tmp;
    (my $shortfile) = $genfile =~ m/\/([^\/]+)$/;
    $shortfile = $genfile if !$shortfile;
    print STDERR "Counting kmers in $id ($shortfile)... ";
    my $gen_start = time;
    my $gen_jelly_start = time;
    if ($is_jelly_gz){
        my $status = system("gzip -cd $genfile >> TemporaryFiles/assm_$count.kmers.txt");
        die "ERROR: Command 'gzip -cd $genfile >> TemporaryFiles/assm_$count.kmers.txt' exited with status $status\n" if $status;
    } else {
        `jellyfish count -o Jellytmp -m $kmer -s 1000000000 -t $threads -C $genfile`; #run jellyfish count (count kmers in infput file) with kmer size -m of $kmer, hash size -s, output file prefix -o of Jellytmp, using number of threads -t and only outputting cannonical kmers (-C)
        die "ERROR: jellyfish count died with status $?\n" if $?;
        #print STDERR "\tDumping counts...\n";
        my @jfiles = glob("Jellytmp*");
        foreach my $jfile (@jfiles){
            `jellyfish dump -tc $jfile >> TemporaryFiles/assm_$count.kmers.txt`;
            die "ERROR: jellyfish dump died with status $?\n" if $?;
        }
        unlink <Jellytmp*>;
    }
    #sort the kmers for splitting later. Will use GNU sort since modern versions will automatically run in parallel
    print STDERR "Sorting... ";
    `sort TemporaryFiles/assm_$count.kmers.txt > TemporaryFiles/assm_$count.kmers_sorted.txt`;
    die "ERROR: Command 'sort TemporaryFiles/assm_$count.kmers.txt > TemporaryFiles/assm_$count.kmers_sorted.txt' died with status $?\n" if $?;
    unlink("TemporaryFiles/assm_$count.kmers.txt");
    my $end_gen_jelly_time = time - $gen_jelly_start;
    print STDERR "Done. ($end_gen_jelly_time secs)\n";
    if (@tmp){
        my $jf_instring;
        my $read_time = time;
        my $is_jelly_read_gz;
        my $num_zips = 0;
        print STDERR "\tCounting kmers in reads... ";
        my $j_read_time = time;
        foreach my $fastq (@tmp){
            if ($fastq =~ m/.jelly.txt.gz$/){
                `gzip -cd $fastq > TemporaryFiles/read_$count.kmers.txt`;
                die "ERROR: Command 'gzip -cd $fastq > TemporaryFiles/read_$count.kmers.txt' died with status $?\n" if $?;
                $is_jelly_read_gz = $fastq;
                last; #only expects one jellyfish file of kmers from both forward and reverse reads.
            }
            if ($fastq =~ m/\.gz$/){
                $jf_instring .= "<(gzip -cd $fastq) ";
                $num_zips++;
            } else {
                $jf_instring .= "<$fastq ";
            }
            die "ERROR parsing read file $fastq: status $?\n" if $?;
        }
        if (!$is_jelly_read_gz) {
            my $j_threads = $threads;
            $j_threads -- if $num_zips > 0;
            my $command = "bash -c \"jellyfish count -o Jellytmp -m $kmer -s 1000000000 -t $j_threads -L $min_k -C $jf_instring\""; #same as above, but only kmers with read depths of at least $min_k will be output (-L)
            `$command`;
            die "ERROR: jellyfish count died with status $?\n" if $?;
            print STDERR "\tFinding conflicting kmers in reads...";
            my @jfiles = glob("Jellytmp*");
            foreach my $jfile (@jfiles){
                `jellyfish dump -tc $jfile >> TemporaryFiles/read_$count.kmers.txt`;
                die "ERROR: jellyfish dump died with status $?\n" if $?;
            }
            unlink <Jellytmp*>;
        }
        #sort the kmers for splitting later. Will use GNU sort since modern versions will automatically run in parallel
        print STDERR "Sorting... ";
        `sort TemporaryFiles/read_$count.kmers.txt > TemporaryFiles/read_$count.kmers_sorted.txt`;
        die "ERROR: Command 'sort TemporaryFiles/read_$count.kmers.txt > TemporaryFiles/read_$count.kmers_sorted.txt' died with status $?\n" if $?;
        unlink("TemporaryFiles/read_$count.kmers.txt");
        my $end_read_jelly_time = time - $j_read_time;
        print STDERR "Done. ($end_read_jelly_time secs)\n";
    }
    $count++;
}



my $reads_given;
for my $i (0 .. $#files){
    my @tmp = @{$files[$i]};
    my $id = shift @tmp;
    my $genfile = shift @tmp;
    my $is_jelly_gz = 1 if $genfile =~ m/.jelly.txt.gz$/;
    my $has_reads = "N";
    $has_reads = "Y" if @tmp;
    (my $shortfile) = $genfile =~ m/\/([^\/]+)$/;
    $shortfile = $genfile if !$shortfile;
    print STDERR "Counting kmers in $id ($shortfile)... ";
    my $gen_start = time;
    my $gen_jelly_start = time;
    my %results;
    my $k_count = 0;
    if ($is_jelly_gz){
        open (my $dump_in, "gzip -cd $genfile | ") or die "ERROR: Can't open gzipped jellyfish file $genfile: $!\n";
        while (my $line = <$dump_in>){
            chomp $line;
            my ($seq, $cnt) = split(" ", $line);
            $k_count++;
            my $base = substr($seq, $armleng, 1, ".");
            $results{$seq}{$base} += $cnt;
        }
        close ($dump_in);
    } else {
        `jellyfish count -o Jellytmp -m $kmer -s 1000000000 -t $threads -C $genfile`; #run jellyfish count (count kmers in infput file) with kmer size -m of $kmer, hash size -s, output file prefix -o of Jellytmp, using number of threads -t and only outputting cannonical kmers (-C)
        die "ERROR: jellyfish count died with status $?\n" if $?;
        #print STDERR "\tDumping counts...\n";
        my @jfiles = glob("Jellytmp*");
        foreach my $jfile (@jfiles){
            open (my $dump_in, "jellyfish dump -tc $jfile |") or die "ERROR: Can't run jellyfish dump: $!\n";
            while (my $line = <$dump_in>){
                chomp $line;
                my ($seq, $cnt) = split(" ", $line);
                $k_count++;
                my $base = substr($seq, $armleng, 1, ".");
                $results{$seq}{$base} += $cnt;
            }
            close ($dump_in);
        }
        unlink <Jellytmp*>;
    }
    my $end_gen_jelly_time = time - $gen_jelly_start;
    print STDERR "Done. ($end_gen_jelly_time secs, $k_count unique kmers)\n";
    #if reads were given, screen them for conflicting kmers
    my %read_conflicts;
    if (@tmp){
        $reads_given = 1;
        my $jf_instring;
        my %read_results;
        my ($r_kcount, $r_knum) = (0) x 2;
        my $read_time = time;
        my $is_jelly_read_gz;
        while (@tmp){
            my $fastq = shift @tmp;
            if ($fastq =~ m/.jelly.txt.gz$/){
                $is_jelly_read_gz = $fastq;
                last; #only expects one jellyfish file of kmers from both forward and reverse reads.
            }
            if ($fastq =~ m/\.gz$/){
                `gzip -cd $fastq >> tmp_reads.fastq`;
            } else {
                `cat $fastq >> tmp_reads.fastq`;
            }
            die "ERROR parsing read file $fastq: status $?\n" if $?;
        }
        my $conflict_time;
        if ($is_jelly_read_gz){
            open (my $dump_in, "gzip -cd $is_jelly_read_gz | ") or die "Can't open gzipped jellyfish read file $is_jelly_read_gz: $!\n";
            print STDERR "\tFinding conflicting kmers in reads...";
            $conflict_time = time;
            while (my $line = <$dump_in>){
                chomp $line;
                my ($seq, $cnt) = split(" ", $line);
                next if $cnt < $min_k;
                $r_kcount ++;
                $r_knum += $cnt;
                my $base = substr($seq, $armleng, 1, ".");
                $read_results{$seq}{$base} += $cnt;
            }
            close ($dump_in);
        } else {
            print STDERR "\tRunning jellyfish count on reads...\n";
            my $command = "bash -c \"jellyfish count -o Jellytmp -m $kmer -s 1000000000 -t $threads -L $min_k -C tmp_reads.fastq\""; #same as above, but only kmers with read depths of at least $min_k will be output (-L)
            `$command`;
            die "ERROR: jellyfish count died with status $?\n" if $?;
            unlink ("tmp_reads.fastq");
            print STDERR "\tFinding conflicting kmers in reads...";
            $conflict_time = time;
            my @jfiles = glob("Jellytmp*");
            foreach my $jfile (@jfiles){
                open (my $dump_in, "jellyfish dump -tc $jfile |") or die "Can't run jellyfish dump: $!\n";
                while (my $line = <$dump_in>){
                    chomp $line;
                    my ($seq, $cnt) = split(" ", $line);
                    $r_kcount ++;
                    $r_knum += $cnt;
                    my $base = substr($seq, $armleng, 1, ".");
                    $read_results{$seq}{$base} += $cnt;
                }
                close ($dump_in);
            }
            unlink <Jellytmp*>;
        }
        my $end_conflict_time = time - $conflict_time;
        print STDERR " Done ($end_conflict_time seconds, $r_kcount unique kmers, $r_knum total kmers)\n";
        
        #parse the read kmers both to output the results and identify kmers with conflicts
        my $last_group = "X";
        my $out;
        my $cluster_time = time;
        foreach my $seq (sort {$a cmp $b} keys %read_results){
            my $group = substr($seq, 0, 4);
            if ($group ne $last_group){
                close ($out) if $out;
                print STDERR "\r\t\tSplitting $group-mers";
                open ($out, "> tmp_cluster_$group.txt");
            }
            my $tmp = "$seq";
            foreach my $base (@bases){
                my $val = 0;
                $val = $read_results{$seq}{$base} if $read_results{$seq}{$base};
                $tmp .= ",$val";
            }
            print $out "$tmp\n";
            $last_group = $group;
        }
        close ($out) if $out;
        undef %read_results;
        my $end_cluster_time = time - $cluster_time;
        print STDERR "\r\t\tSplitting $last_group-mers ... Done! ($end_cluster_time secs)\n";
        
        my $process_time = time;
        my @r_seqs;
        my $sys_status = system("perl $dirname/lib/kmer_core_mt_conflict.pl $count $threads $err");
        die "\nERROR: perl $dirname/lib/kmer_core_mt_conflict.pl failed with status $sys_status\n" if $sys_status != 0;
        open (my $in, "< tmp_read_conflicts.txt") or die "\nERROR: Can't open tmp_read_conflicts.txt: $!\n";
        my $rc_count = 0;
        while (my $line = <$in>){
            chomp $line;
            $read_conflicts{$line}++;
            $rc_count++;
        }
        close ($in);
        unlink ("tmp_read_conflicts.txt");
        my $end_process_time = time - $process_time;
        print STDERR " ... Done! $rc_count conflicts identified in reads ($end_process_time secs)\n";
        
        my $fin = time - $read_time;
        my $fin_form = tform($fin);
        print STDERR " \t\tTotal read process time: $fin secs ($fin_form)\n";
        
    }
    ## output non-conflicting kmers
    ## for multithreading purposes, will output into separate 4-mer files
    my $removed = 0;
    my $k_output = 0;
    print STDERR "\tSorting kmers...";
    my $k_out;
    my $last_group;
    my $first;
    foreach my $seq (sort {$a cmp $b} keys %results) {
        if (!$first){
            print STDERR "\n\r\tConflicting kmers removed: $removed";
            $first = 1;
        }
        my $group = substr($seq, 0, 4);
        if (!$last_group){
            open ($k_out, ">> tmp_kmers_$group.txt");
        } elsif ($group ne $last_group){
            close ($k_out);
            open ($k_out, ">> tmp_kmers_$group.txt");
        }
        $last_group = $group;
        my %tmp = %{$results{$seq}};
        if ($read_conflicts{$seq} or scalar keys %tmp > 1){
            $removed++;
            print STDERR "\r\tConflicting kmers removed: $removed";
            next;
        }
        my $base;
        foreach my $key (keys %tmp){
            $base = $key;
        }
        print $k_out "$seq\t$base\t$count\n";
        $k_output++;
    }
    close ($k_out);
    print STDERR "\r\tConflicting kmers removed: $removed\n";
    print STDERR "\tUnique kmers output: $k_output\n";
    print $stats "$id\t$count\t$k_output\t$has_reads\t$removed\n";
    $count++;
    my $fin = time - $gen_start;
    my $fin_form = tform($fin);
    my $tot_time = time - $start_time;
    my $tot_form = tform($tot_time);
    print STDERR "\tRun time: $fin_form, Total run time: $tot_form\n";
}

print STDERR "\nMerging kmers from $count genomes\n";
my $a_num = roundup(($a_pct / 100) * $nog);
my @core_kmers;
my %gen_kmers;
my $merge_time = time;
my $sys_status = system("perl $dirname/lib/kmer_core_mt_merge.pl $threads $a_num $count");
die "\nERROR: perl $dirname/lib/kmer_core_mt_merge.pl failed with status $sys_status\n" if $sys_status != 0;
open (my $in, "< tmp_merge_out.txt") or die "\nERROR: Can't open tmp_merge_out.txt: $!\n";
while (my $line = <$in>){
    chomp $line;
    my ($seq, $gens, $gen_count) = split("\t", $line);
    push @core_kmers, $seq if $seq ne "x";
    $gen_kmers{$gen_count}{$gens}++;
}
close ($in);
unlink ("tmp_merge_out.txt");
unlink <tmp_kmers_*>;
my $end_merge_time = time - $merge_time;
my $f_emt = tform($end_merge_time);
print STDERR " ... Done! ($end_merge_time secs, $f_emt)\n";
@core_kmers = sort{$a cmp $b} @core_kmers;
my $num_core_kmers = scalar @core_kmers;
print STDERR "$num_core_kmers core kmers found\n";
print $stats "\n\n$num_core_kmers core kmers found\n";
close ($stats);

## output the core kmer sequence file
my @ns = ("N") x $n_num;
my $nstring = join ("", @ns);
open (my $out, "> $pref\_core_kmers.fasta");
print $out ">$pref\_core_kmers\n";
print $out join($nstring, @core_kmers), "\n";
close ($out);

## output the kmer counts file
open (my $c_out, "> $pref\_kmer_counts.txt");
foreach my $num (sort {$a <=> $b} keys %gen_kmers){
    my %tmp = %{$gen_kmers{$num}};
    foreach my $gens (sort {$a cmp $b} keys %tmp){
        my $val = $tmp{$gens};
        print $c_out "$num\t$val\t$gens\n";
    }
}
close ($c_out);

my $end_time = time - $start_time;
my $f_end_time = tform($end_time);
print STDERR "\n\n***Finished running kmer_core. Total time: $end_time secs [$f_end_time]\nThanks!\n";

#-----------------------------------------------------
sub roundup {
    my $n = shift;
    return(($n == int($n)) ? $n : int($n + 1));
}

sub tform {
    my $T = shift;
    my @out = reverse($T % 60, ($T/=60) % 60, ($T/=60) % 24, ($T/=24));
    my $out = sprintf "%03d:%02d:%02d:%02d", @out;
    return "$out";
}
