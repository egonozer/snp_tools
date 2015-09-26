#!/usr/bin/perl

use warnings;
use strict;

my $version = "0.1";

## processing subscript for ksnp_snp_confirm_batch_v0.2

my $usage = "
$0 <gen_id> <armleng> <err> <threads> <temp_path>

";

die $usage unless @ARGV >= 5;

my @bases = qw(A C G T);

my ($gen_id, $armleng, $err, $threads, $temp_path) = @ARGV;

my %hash;
my $last_group;
my $num_threads_running = 0;

open (my $k_in, "< $temp_path/gen$gen_id.kmers_sorted.txt") or die "ERROR: Can't open $temp_path/gen$gen_id.kmers_sorted.txt: $!\n";
open (my $all_out, "> $temp_path/gen$gen_id.corr_results.txt") or die "ERROR: Can't open $temp_path/gen$gen_id.corr_results.txt: $!\n";
while (my $line = <$k_in>){
    chomp $line;
    my ($seq, $cnt) = split(" ", $line);
    my $base = substr($seq, $armleng, 1, ".");
    my $group = substr($seq, 0, 4);
    if ($last_group){
        if ($last_group ne $group and keys %hash >= 1){
            if ($num_threads_running == $threads){
                my $pid = wait;
                open (my $in, "< $temp_path/$pid.results.txt") or die "ERROR: Can't open $temp_path/$pid.results.txt: $!\n";
                while (my $line = <$in>){
                    print $all_out "$line";
                }
                close ($in);
                unlink ("$temp_path/$pid.results.txt");
                $num_threads_running--;
            }
            print STDERR "\r\tCorrecting $last_group-mers...";
            fork_off();
            $num_threads_running++;
            %hash = ();
        }
    }
    push @{$hash{$seq}}, ([$cnt, $base]);
    $last_group = $group;
}
close ($k_in);
if ($last_group){
    if (keys %hash >= 1){
        if ($num_threads_running == $threads){
            my $pid = wait;
            open (my $in, "< $temp_path/$pid.results.txt") or die "ERROR: Can't open $temp_path/$pid.results.txt: $!\n";
            while (my $line = <$in>){
                print $all_out "$line";
            }
            close ($in);
            unlink ("$temp_path/$pid.results.txt");
            $num_threads_running--;
        }
        print STDERR "\r\tCorrecting $last_group-mers...";
        fork_off();
        $num_threads_running++;
        %hash = ();
    }
}
print STDERR "\n";
while (my $pid = wait){
    last if $pid == -1;
    open (my $in, "< $temp_path/$pid.results.txt") or die "ERROR: Can't open $temp_path/$pid.results.txt: $!\n";
    while (my $line = <$in>){
        print $all_out "$line";
    }
    close ($in);
    unlink ("$temp_path/$pid.results.txt");
    $num_threads_running--;
}
close ($all_out);

#---------------------------
sub fork_off {
    my $pid = fork;
    if ($pid == 0){
        open (my $out, "> $temp_path/$$.results.txt");
        #print STDERR "\rProcessing $last_group-mers";
        if (-e "$temp_path/refs_$last_group.txt"){
            open (my $in, "< $temp_path/refs_$last_group.txt");
            while (my $line = <$in>){
                chomp $line;
                next if $line =~ m/^\s*$/;
                my ($one, $two) = split("\t", $line);
                my ($id, $kmer) = split(",", $one);
                my $ref_base = "-";
                if ($two =~ m/,$gen_id-(\w)/){
                    $ref_base = $1;
                }
                unless ($hash{$kmer}){
                    print $out "$id,$ref_base,u\n";
                    next;
                }
                my @array = @{$hash{$kmer}};
                my $sum = 0;
                my $ref_count = 0;
                foreach my $rec (@array){
                    my ($cnt, $base) = @{$rec};
                    $sum += $cnt;
                    $ref_count = $cnt if $base eq $ref_base;
                }
                my $error = 0;
                $error = 1 - ($ref_count / $sum) if $sum > 0;
                if ($ref_base eq "-" and $sum > 0){
                    @array = sort{$b->[0] <=> $a->[0]} @array;
                    ($ref_count, $ref_base) = @{$array[0]};
                    my $new_error = 1 - ($ref_count / $sum);
                    if ($new_error > $err){
                        print $out "$id,-,u\n";
                        next;
                    }
                    print $out "$id,$ref_base,o\n";
                    next;
                }
                if ($error > $err){
                    print $out "$id,-,m\n";
                } else {
                    print $out "$id,$ref_base,u\n";
                }
            }
            close ($in);
        }
        close $out;
        exit(0);
    }
    return;
}
