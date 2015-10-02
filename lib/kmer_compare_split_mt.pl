#!/usr/bin/perl

use warnings;
use strict;
use File::Path 'rmtree';

$|++;

my $version = 0.2;

my $usage = "$0 <file count> <armleng> <path to temporary folder> <threads>";

die $usage unless @ARGV >= 4;
my ($file_count, $armleng, $temppath, $threads) = @ARGV;

my @bases = qw(A C G T);

#load the file array and results hash;
my %results; #results fields: id_number, has_assm, assm_mers, has_read, read_mers
my @file_list;
for my $i (1 .. $file_count){
    @{$results{$i}} = (0, 0, 0, 0);
    if (-e "$temppath/kmers$i.s.txt"){
        push @file_list, ([$i, "s"]);
        $results{$i}[0] = 1;
    }
    if (-e "$temppath/kmers$i.r.txt"){
        push @file_list, ([$i, "r"]);
        $results{$i}[2] = 1;
    }
}
die "ERROR: No temporary kmer files found at $temppath\n" if !@file_list;
my $total_count = scalar @file_list;

my $num_threads_running = 0;
my $done_count = 0;
while (@file_list){
    if ($num_threads_running == $threads){
        my $pid = wait;
        my $status = $?;
        die "ERROR: Process $pid exited with status $status\n" if $status;
        open (my $in, "< $temppath/split_results.$pid.txt") or die "ERROR: Can't open '$temppath/split_results.$pid.txt': $!\n";
        my $line = <$in>;
        close ($in);
        unlink ("$temppath/split_results.$pid.txt");
        my ($id_num, $type, $count) = split(",", $line);
        if ($type eq "s"){
            $results{$id_num}[1] = $count;
        } else {
            $results{$id_num}[3] = $count;
        }
        unlink ("$temppath/kmers$id_num.$type.txt");
        $done_count++;
        print STDERR "\r\tFinished splitting $done_count/$total_count files";
        $num_threads_running--
    }
    my ($id_num, $type) = @{shift @file_list};
    my $pid = fork;
    if ($pid == 0){
        #sort the kmers
        my ($last_seq, $last_val, $last_mer);
        my $mer_out;
        my $count = 0;
        my $out_type = "assm";
        $out_type = "read" if $type eq "r";
        open (my $in, "sort $temppath/kmers$id_num.$type.txt |") or die "command 'sort $temppath/kmers$id_num.$type.txt' exited with status $!\n";
        while (my $line = <$in>){
            chomp $line;
            my ($seq, $val) = split(" ", $line);
            my $mer = substr($seq, 0, 4);
            if (!$last_seq){
                open ($mer_out, "> $temppath/kmers$id_num.$out_type\_$mer.txt");
                ($last_seq, $last_val, $last_mer) = ($seq, $val, $mer);
                next;
            }
            if ($seq eq $last_seq){
                $last_val += $val;
                next;
            }
            if ($seq ne $last_seq){
                my $base = substr($last_seq, $armleng, 1, ".");
                print $mer_out "$last_seq,$base,$last_val\n";
                $count++;
            }
            if ($mer ne $last_mer){
                close ($mer_out);
                open ($mer_out, "> $temppath/kmers$id_num.$out_type\_$mer.txt");
            }
            ($last_seq, $last_val, $last_mer) = ($seq, $val, $mer);
        }
        my $base = substr($last_seq, $armleng, 1, ".");
        print $mer_out "$last_seq,$base,$last_val\n";
        $count++;
        close ($mer_out);
        close ($in);
        unlink ("$temppath/kmers$id_num.$type.txt");
        open (my $out, "> $temppath/split_results.$$.txt");
        print $out "$id_num,$type,$count";
        exit(0);
    }
    $num_threads_running++;
}
while (my $pid = wait){
    last if $pid == -1;
    my $status = $?;
    die "ERROR: Process $pid exited with status $status\n" if $status;
    open (my $in, "< $temppath/split_results.$pid.txt") or die "ERROR: Can't open '$temppath/split_results.$pid.txt': $!\n";
    my $line = <$in>;
    close ($in);
    unlink ("$temppath/split_results.$pid.txt");
    my ($id_num, $type, $count) = split(",", $line);
    if ($type eq "s"){
        $results{$id_num}[1] = $count;
    } else {
        $results{$id_num}[3] = $count;
    }
    unlink ("$temppath/kmers$id_num.$type.txt");
    $done_count++;
    print STDERR "\r\tFinished splitting $done_count/$total_count files";
    $num_threads_running--
}
print STDERR "\n";

open (my $out, "> $temppath/split_results.txt");
for my $i (1 .. $file_count){
    my @array = @{$results{$i}};
    print $out "$i,", join(",", @array), "\n";
}
close ($out);
