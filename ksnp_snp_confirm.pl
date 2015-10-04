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

my $version = "0.2";

##changes from v0.1
## removed search for reverse complement of kmers.  If -C option is included with jellyfish count, reverse-complement kmers should not be present
## multithread analysis using external script

#based on ksnp_snp_confirm_v0.2.pl

my $usage = "
ksnp_snp_confirm_batch.pl [options] -f file_of_files.txt -a SNPs_all

Searches sequence reads for SNP kmers output by kSNP, confirms the SNP
assignments, and outputs new SNPs_all and SNPs_all_matrix.fasta files with
corrected positions.

Does the following corrections:
1) If a Ksnp SNP position is found with multiple midddle bases in the reads
   (exceeding a given error rate, -e), it will be replaced with a \"-\" in the
   matrix.fasta file and omitted from the SNPs_all file
2) If a Ksnp SNP position was listed as absent in the genome, but a minimum
   number of reads contain the SNP position (as given by -l) and are not too
   mixed (maximum error given by -e), the base will be included in the
   matrix.fasta file and the SNPs_all file
   
For simplicity, corrected genome positions in the SNPs_all file will be listed
as \"0\" and all kmer directions will be listed as \"F\".

Will also output, to STDOUT, a listing of modified SNP positions. The
categories are:
  #_positions:  Total number of SNP positions identified by kSNP
  #_confirmed:  Total number of SNP positions where the base call was confirmed
                by the reads
  #_mixed:      Total number of SNP positions where a mixture of bases were
                found in the reads with an error rate above the threshold set
                by -e.
  #_omitted:    Total number of SNP positions that were found to be present and
                not mixed in the reads, but were not listed for this genome
                in the SNPs_all file.

REQUIRED:
  -f    file of files.  Contains genome names (as provided to Ksnp) and paths
        to one or more fastq read files. Gzipped or uncompressed fastq read
        files may be used. If no read file is listed, the SNP positions
        will be included in the output files as given in the SNPs_all input
        file.
            Format of example file:
            genome_A<tab>/path/to/genome_A_reads_1.fastq.gz<tab>/path/to/genome_A_reads_2.fastq.gz
            genome_B<tab>/path/to/genome_B_reads.fastq.gz
            genome_C
            genome_D<tab>/path/to/genome_C_reads.fastq
            etc..
  -a    SNPs_all file (as output by Ksnp)

OPTIONS:
  -l    minimum number of times a kmer must be present to be counted. Higher
        numbers should run faster
        (default 5)
  -e    maximum allowable error rate between number of reads with reference
        base and number of reads with discrepant base. This is calcuated as:
        1 - (Nref / (Na + Nc + Ng + Nt))
        where Nref is the number of reads with the reference base and Na is the
        number of reads with \"A\" at the SNP position, Nc is the number of reads
        with \"C\" at the SNP position, etc.
        If the -d option is given and the error rate is greater than this value,
        the position will be output.
        (default: 0.1)
  -p    output file prefix
        (default: \"corrected\")
  -t    number of threads
        (default: 15)

";

use Getopt::Std;
use vars qw( $opt_t $opt_l $opt_e $opt_f $opt_a $opt_p );
getopts('l:e:t:f:a:p:');

my $threads = $opt_t ? $opt_t : 15;
my $min_k   = $opt_l ? $opt_l : 5;
my $err     = $opt_e ? $opt_e : .1;
my $pref    = $opt_p ? $opt_p : "corrected";

die $usage unless $opt_f and $opt_a;

my $start_time = time;

mkdir "TemporaryFiles" or die "Can't make TemporaryFiles folder: $!\n";
my $temp_path = abs_path("TemporaryFiles");
print STDERR "temp_path: $temp_path\n";

my $script_path = abs_path($0);
$script_path = dirname($script_path);
print STDERR "script_path: $script_path\n";


open (my $snp_in, "<$opt_a") or die "Can't open $opt_a: $!\n";
print STDERR "Reading in SNPs_all file\n";
## SNPs_all file should have kmers already sorted alphabetically
my @kmers;
my %gens_per_snp;
my %p_hash;
my $armleng;
my $k_size;
my %ref_bases;
my %gen_ids;
my $gen_count = 0;
my $out_refs;
my $last_group = "x";
while (my $line = <$snp_in>){
    chomp $line;
    next if $line =~ m/^\s*$/; #skip blank lines
    my ($id, $kmer, $base, $pos, $dir, $gen) = split(" ", $line);
    my $tempid = $id + 1;
    unless (defined $gen_ids{$gen}){
        $gen_count++;
        $gen_ids{$gen} = $gen_count;
    }
    my $group = substr($kmer, 0, 4);
    if ($group ne $last_group){
        if ($out_refs){
            print $out_refs "\n";
            close $out_refs;
        }
        print STDERR "\r\tLoading $group-mers...";
        open ($out_refs, "> $temp_path/refs_$group.txt");
    }
    $last_group = $group;
    unless ($armleng){
        $k_size = length($kmer);
        $armleng = int($k_size / 2);
    }
    unless ($p_hash{$kmer}){
        $p_hash{$kmer} = $tempid;
        push @kmers, ([$kmer, $tempid]);
        print $out_refs "\n$tempid,$kmer\t";
    }
    $gens_per_snp{$tempid}{$gen}++;
    push @{$ref_bases{$gen}{$tempid}}, ([$base, $pos, $dir]);
    my $gen_id = $gen_ids{$gen};
    print $out_refs ",$gen_id-$base";
}
print STDERR "\n";
close ($snp_in);
my $num_pos = scalar @kmers;
if ($out_refs){
    print $out_refs "\n";
    close $out_refs;
}

#jellyfish the files
open (my $fof, "<$opt_f") or die "Can't open file of files $opt_f: $!\n";
my @output_order;
while (my $line = <$fof>){
    chomp $line;
    my @tmp = split ("\t", $line);
    my $gen = shift @tmp;
    push @output_order, $gen;
    if (!$ref_bases{$gen}){
        print STDERR "WARNING: no data for genome $gen in input SNPs_all file. Skipping.\n";
        next;
    }
    if (@tmp){
        print STDERR "Running jellyfish on $gen...";
        my $jellytime = time;
        my $gen_id = $gen_ids{$gen};
        my $num_zips = 0;
        my $jf_instring;
        foreach my $fastq (@tmp){
            if ($fastq =~ m/\.gz$/){
                $jf_instring .= "<(gzip -cd $fastq) ";
                $num_zips++;
                #`gzip -cd $fastq >>TemporaryFiles/temp.fastq`;
                #die "\nERROR: Command 'gzip -cd $fastq >>TemporaryFiles/temp.fastq' exited with status $?\n" if $?;
            } else {
                $jf_instring .= "<$fastq ";
                #my $status = system("bash -c\"cat $fastq >>TemporaryFiles/temp.fastq\"");
                #die "\nERROR: Command 'cat $fastq >>TemporaryFiles/temp.fastq' exited with status $status\n" if $status;
            }
        }
        #my $j_threads = $threads - $num_zips;
        my $j_threads = $threads - 1 if $num_zips > 0;
        my $command = "bash -c \"jellyfish count -o TemporaryFiles/Jellytmp -m $k_size -s 1000000000 -t $j_threads -L $min_k -C $jf_instring\"";
        `$command`;
        die "\nERROR: jellyfish count died with status $?\n" if $?;
        #unlink ("TemporaryFiles/temp.fastq") if -e "TemporaryFiles/temp.fastq";
        my @jfiles = glob("TemporaryFiles/Jellytmp*");
        foreach my $jfile (@jfiles){
            my $status = system("jellyfish dump -tc $jfile >>TemporaryFiles/kmers.txt");
            die "\nERROR: command 'jellyfish dump -tc $jfile >>TemporaryFiles/kmers.txt' exited with status $status\n" if $status;
            unlink ($jfile);
        }
        `sort TemporaryFiles/kmers.txt > TemporaryFiles/gen$gen_id.kmers_sorted.txt`;
        die "\nERROR: command 'sort TemporaryFiles/kmers.txt >TemporaryFiles/gen$gen_id.kmers_sorted.txt' exited with status $?\n" if $?;
        unlink ("TemporaryFiles/kmers.txt") if -e "TemporaryFiles/kmers.txt";
        my $j_finish = time - $jellytime;
        print STDERR " Done ($j_finish seconds)\n";
    }
}

#correct the kmers, using external program to facilitate multithreading
my %corrected;
my @stats;
foreach my $gen (@output_order){
    my $gen_id = $gen_ids{$gen};
    if (-e "TemporaryFiles/gen$gen_id.kmers_sorted.txt"){
        print STDERR "Correcting $gen...\n";
        my $status = system("perl", "$script_path/lib/ksnp_snp_confirm_batch.correct_mt.pl", "$gen_id", "$armleng", "$err", "$threads", "$temp_path");
        die "ERROR: ksnp_snp_confirm_batch.correct_mt.pl exited with status: $status. Message '$?'\n" if $status;
        my ($num_mixed, $num_omitted) = (0) x 2;
        open (my $in, "< TemporaryFiles/gen$gen_id.corr_results.txt") or die "ERROR: Can't open TemporaryFiles/gen$gen_id.corr_results.txt: $!\n";
        while (my $line = <$in>){
            chomp $line;
            my ($id, $base, $corr) = split(",", $line);
            $num_mixed++ if $corr eq "m";
            $num_omitted++ if $corr eq "o";
            $corrected{$gen}{$id} = $base;
        }
        close ($in);
        unlink ("TemporaryFiles/gen$gen_id.corr_results.txt");
        print STDERR "\t#snp positions: $num_pos, #mixed: $num_mixed, #omitted: $num_omitted\n";
        my $num_correct = $num_pos - ($num_mixed + $num_omitted);
        push @stats, ([$gen, $num_pos, $num_correct, $num_mixed, $num_omitted]);
    } else {
        print STDERR "Outputting $gen, no correction\n"
    }
}

#output SNP matrix
print STDERR "Outputting $pref\_SNPs_all_matrix.fasta...\n";
open (my $mat_out, "> $pref\_SNPs_all_matrix.fasta");
foreach my $gen (@output_order){
    print $mat_out ">$gen\n";
    for my $i (0 .. $#kmers){
        my $id = $kmers[$i][1];
        if ($corrected{$gen}{$id}){
            print $mat_out "$corrected{$gen}{$id}";
        } elsif ($ref_bases{$gen}{$id}) {
            print $mat_out "$ref_bases{$gen}{$id}[0][0]"; #even with multiple instances of a kmer in a genome, they should all have the same middle base
        } else {
            print $mat_out "-";
        }
    }
    print $mat_out "\n";
}

#output corrected_SNPs_all
print STDERR "Outputting $pref\_SNPs_all...\n";
open (my $all_out, ">$pref\_SNPs_all");
print $all_out "\n";
for my $i (0 .. $#kmers){
    my ($kmer, $id) = @{$kmers[$i]};
    my $act_id = $id - 1;
    my $printed_something;
    foreach my $gen (@output_order){
        my $base;
        my @recs = @{$ref_bases{$gen}{$id}} if $ref_bases{$gen}{$id};
        my $corr_base = $corrected{$gen}{$id} if $corrected{$gen}{$id};
        next unless (@recs or $corr_base);
        if ($corr_base){
            next if $corr_base eq "-";
            $base = $corr_base;
        }
        if (@recs){
            if (!$base){
                $base = $recs[0][0];
            }
            for my $i (0 .. $#recs){
                print $all_out "$act_id\t$kmer $base\t$recs[$i][1] $recs[$i][2]\t$gen\n";
            }
        } else {
            print $all_out "$act_id\t$kmer $base\t0 F\t$gen\n";
        }
        $printed_something = 1;
    }
    print $all_out "\n" if $printed_something;
}
close ($all_out);

open (my $stat_out, ">$pref\_correction_stats.txt");
print $stat_out "genome\t#_positions\t#_confirmed\t#_mixed\t#_omitted\n";
for my $i (0 .. $#stats){
    my @tmp = @{$stats[$i]};
    print $stat_out join("\t", @tmp), "\n";
}
close ($stat_out);

rmtree(["TemporaryFiles"]);

my $end_time = time - $start_time;
my $t_time = tform($end_time);
print STDERR "Finished. $end_time seconds ($t_time)\n";


#--------------------
sub tform {
    my $T = shift;
    my @out = reverse($T % 60, ($T/=60) % 60, ($T/=60) % 24, ($T/=24));
    my $out = sprintf "%03d:%02d:%02d:%02d", @out;
    return "$out";
}
