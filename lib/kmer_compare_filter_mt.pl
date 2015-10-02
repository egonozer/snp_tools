#!/usr/bin/perl

use warnings;
use strict;

$|++;

my $version = 0.3;

my $usage = "$0 <id number> <assm (y or n)> <read (y or n)> <core (y or n)> <min_k> <err> <path to temporary folder> <threads>";

die $usage unless @ARGV >= 8;

my $id_number = $ARGV[0];
my $assm = lc($ARGV[1]);
my $read = lc($ARGV[2]);
my $core = lc($ARGV[3]);
my $min_k = $ARGV[4];
my $err = $ARGV[5];
my $temp_path = $ARGV[6];
my $threads = $ARGV[7];

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
my ($filt_a, $filt_r) = (-1) x 2;
my @long_file_string;
foreach my $mer (@mers){
    if ($num_threads_running == $threads){
        my $pid = wait;
        open (my $in, "< $temp_path/$pid.tmp.txt") or die "ERROR: Can't open $temp_path/$pid.tmp.txt: $!\n";
        my $line = <$in>;
        chomp $line;
        my @tmp = split(",", $line);
        close $in;
        unlink ("$temp_path/$pid.tmp.txt");
        my $p_filt_a = shift @tmp;
        my $p_filt_r = shift @tmp;
        $filt_a += $p_filt_a;
        $filt_r += $p_filt_r;
        if (@tmp){
            if (@long_file_string){
                @long_file_string = @tmp if scalar @tmp > scalar @long_file_string;
            } else {
                @long_file_string = @tmp;
            }
        }
        $num_threads_running--;
    }
    print STDERR "\r\tFiltering $mer-mers";
    my $pid = fork;
    if ($pid == 0){
        my ($p_filt_a, $p_filt_r) = (0) x 2;
        my @file_string;
        my (%core_kmers, %k_assm, %k_read);
        if ($core eq "y" and -e "$temp_path/core_$mer.txt"){
            open (my $in, "< $temp_path/core_$mer.txt");
            while (my $line = <$in>){
                chomp $line;
                $core_kmers{$line}++;
            }
            close ($in);
        }
        if ($assm eq "y" and -e "$temp_path/kmers$id_number.assm_$mer.txt"){
            open (my $in, "< $temp_path/kmers$id_number.assm_$mer.txt");
            while (my $line = <$in>){
                chomp $line;
                my ($seq, $base, $val) = split(",", $line);
                $k_assm{$seq}{$base} = $val;
            }
            close ($in);
            unlink ("$temp_path/kmers$id_number.assm_$mer.txt");
        }
        if ($read eq "y" and -e "$temp_path/kmers$id_number.read_$mer.txt"){
            open (my $in, "< $temp_path/kmers$id_number.read_$mer.txt");
            while (my $line = <$in>){
                chomp $line;
                my ($seq, $base, $val) = split(",", $line);
                $k_read{$seq}{$base} = $val;
            }
            close ($in);
            unlink ("$temp_path/kmers$id_number.read_$mer.txt");
        }
        #filtering begins
        my %filt_k_assm;
        if (%k_assm){
            foreach my $seq (keys %k_assm){
                my $rev = reverse($seq);
                $rev =~ tr/ACTG/TGAC/;
                #skip kmers not found in core
                if ($core eq "y"){
                    next unless ($core_kmers{$seq} or $core_kmers{$rev});
                }
                #remove conflicting kmers
                next if scalar keys %{$k_assm{$seq}} > 1;
                my $a_base;
                foreach my $base (@bases){
                    $a_base = $base if $k_assm{$seq}{$base};
                }
                #if reads were given, remove any kmers with conflicts in the reads
                if (%k_read){
                    my %tmp;
                    %tmp = %{$k_read{$seq}} if $k_read{$seq};
                    %tmp = %{$k_read{$rev}} if $k_read{$rev};
                    if (%tmp){
                        my $sum = 0;
                        my @counts;
                        foreach my $base (@bases){
                            my $val = 0;
                            $val = $tmp{$base} if $tmp{$base};
                            $sum += $val;
                            push @counts, ([$val, $base]);
                        }
                        @counts = sort{$b->[0] <=> $a->[0]} @counts;
                        my ($ref_count, $base) = @{$counts[0]};
                        next if $ref_count < $min_k;
                        my $error = 1 - ($ref_count / $sum);
                        next if ($error > $err);
                    }
                }
                #add filtered kmers to a hash
                $filt_k_assm{$seq} = $a_base;
            }
            undef %k_assm; #not sure this will save memory, but probably can't hurt?
        }  
        if (%k_read and !%filt_k_assm){ #only look for new kmers in reads if no assembly was given
            my $k_sum = 0;
            my @array;
            foreach my $seq (sort {$a cmp $b} keys %k_read){
                my $rev = reverse($seq);
                $rev =~ tr/ACTG/TGAC/;
                #skip kmers not found in core
                if ($core eq "y"){
                    next unless ($core_kmers{$seq} or $core_kmers{$rev});
                }
                #remove conflicting kmers
                my $sum = 0;
                my @counts;
                foreach my $base (@bases){
                    my $val = 0;
                    $val = $k_read{$seq}{$base} if $k_read{$seq}{$base};
                    $sum += $val;
                    push @counts, ([$val, $base]);
                }
                @counts = sort{$b->[0] <=> $a->[0]} @counts;
                my ($ref_count, $base) = @{$counts[0]};
                next if $ref_count < $min_k;
                my $error = 1 - ($ref_count / $sum);
                next if ($error > $err);
                if ($core eq "y"){
                    if ($assm eq "y"){ #if both a core file and an assembly file were given, add any non-conflicting core kmers not found in the assembly to that list of kmers, then just output one list of kmers
                        next unless (!$filt_k_assm{$seq} and !$filt_k_assm{$rev});
                    }
                    # if a core file and only read file(s) were given without an assembly, add the filtered reads to the filt_k_assm hash and output filtered read kmers as if they were assembly kmers
                    $filt_k_assm{$seq} = $base;
                    next;
                }
                #if either a core file or an assembly file for this genome (or both) was not given, then will output the read kmers for pairwise comparisons
                $p_filt_r++;
                next if $filt_k_assm{$seq} or $filt_k_assm{$rev}; #to save space and time, only output read kmers if they aren't already present in the assembly
                push @array, "$seq\t$base\n";
                #print $r_out "$seq\t$base\n";
            }
            if (@array){
                push @file_string, "read";
                open (my $r_out, "> $temp_path/read_kmers_$id_number\_$mer.txt") or die "ERROR: Can't write to $temp_path/read_kmers_$id_number\_$mer.txt: $!\n";
                while (@array){
                    my $line = shift @array;
                    print $r_out "$line";
                }
                close ($r_out);
            }
            undef %k_read;
        }
        #output the assembly kmers, if given
        if (%filt_k_assm){
            push @file_string, "assm";
            open (my $a_out, "> $temp_path/assm_kmers_$id_number\_$mer.txt") or die "ERROR: Can't write to $temp_path/assm_kmers_$id_number\_$mer.txt: $!\n";
            foreach my $seq (sort {$a cmp $b} keys %filt_k_assm){
                my $base = $filt_k_assm{$seq};
                print $a_out "$seq\t$base\n";
            }
            my $val = scalar keys %filt_k_assm;
            $p_filt_a = $val;
            undef %filt_k_assm;
        }
        open (my $out, "> $temp_path/$$.tmp.txt") or die "ERROR: Can't write to $temp_path/$$.tmp.txt: $!\n";
        print $out "$p_filt_a,$p_filt_r";
        if (@file_string){
            print $out ",", join(",", @file_string);
        }
        print $out "\n";
        close ($out);
        exit(0);
    }
    $num_threads_running++;
    
}
print STDERR "\n";
while (my $pid = wait){
    last if $pid == -1;
    open (my $in, "< $temp_path/$pid.tmp.txt") or die "ERROR: Can't open $temp_path/$pid.tmp.txt: $!\n";
    my $line = <$in>;
    chomp $line;
    my @tmp = split(",", $line);
    close $in;
    unlink ("$temp_path/$pid.tmp.txt");
    my $p_filt_a = shift @tmp;
    my $p_filt_r = shift @tmp;
    $filt_a += $p_filt_a;
    $filt_r += $p_filt_r;
    if (@tmp){
        if (@long_file_string){
            @long_file_string = @tmp if scalar @tmp > scalar @long_file_string;
        } else {
            @long_file_string = @tmp;
        }
    }
    $num_threads_running--;
}
open (my $out, "> $temp_path/counts.txt") or die "ERROR: Can't write to $temp_path/counts.txt: $!\n";
print $out "$filt_a,$filt_r";
if (@long_file_string){
    print $out ",", join(",", @long_file_string);
}
print $out "\n";
close ($out);
