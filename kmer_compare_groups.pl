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

$|++;

my $usage = "
kmer_compare_groups.pl <genome list file> <kmer_compare.raw.txt>

Outputs values and statistics of SNP differences based on different groupings.

Genome list format:
First line should contain a listing of the different categories separated by
commas. It should start with a '#' symbol.
For example, if the categories your genomes can belong to inlcude species,
clade, and outbreak, then your first line should look like this:

#species,clade,outbreak

All subsequent lines should start with the name of the genome (same names as in
kmer_compare.raw.txt) and an integer for each group listed above to indicate
which group the genome belongs to, all separated by commas. You can also give a
'0' to indicate no group membership or unknown group membership.
For example, using the gropus outlined above, the next lines could be:
genA,1,1,1     #genome A belongs to species 1, clade 1, and outbreak 1
genB,1,2,0     #genome B belongs to species 1, clade 2, and no outbreak
genC,2,0,0     #genome C belongs to species 2, unknown clade or outbreak

Group memberships are independent of each other, so in the example above, if
genomes B and C were given 1 and 2 respectively to indicate different strains,
but both were given 2 for clade, they will still be grouped together as clade-
mates.  

Options:
  -x    if given, will only count intergroup comparisons as different. For
        example, a comparison between a group 1 and a group 2 strain would be
        counted as \"different\". A comparison between a group 1 strain and a
        group 0 strain would be counted as \"different\".  A comparison between
        a group 0 strain and another group 0 strain would not be counted as
        \"same\" or \"different.\"
        Default: Comparisons to or between group 0 strains will not be counted
        at all
  -p    Output prefix
        Default: 'output'

";

use Getopt::Std;
use vars qw( $opt_x $opt_p );
getopts('xp:');

die $usage unless @ARGV >= 2;

my ($list, $raw) = @ARGV;
my $pref = $opt_p ? $opt_p : "output";

open (my $lin, "<$list") or die "Can't open $list: $!\n";
my @groups;
my %group_hash;
while (my $line = <$lin>){
    chomp $line;
    if (!@groups){
        next unless $line =~ m/^#/;
        $line = substr($line, 1);
        @groups = split(",", $line);
        next;
    }
    $line =~ s/\s*#.*$//; #delete any comments
    next if $line =~ m/^\s*$/; #skip any blank lines
    my @tmp = split(",", $line);
    my $gen = shift @tmp;
    die "ERROR: Not enough groups in this line: `$line`\n" if scalar @tmp != scalar @groups;
    for my $i (0 .. $#tmp){
        $group_hash{$gen}{$groups[$i]} = $tmp[$i];
    }
}
close ($lin);
die "ERROR: List file $list was either empty or missing the groups definition line.\n" unless @groups and %group_hash;

open (my $rin, "<$raw") or die "Can't open $raw: $!\n";
my %results;
while (my $line = <$rin>){
    chomp $line;
    next if $line =~ m/^#/;
    my ($genA, $genA_pos, $genB, $genB_pos, $snps, $shared, $frac) = split("\t", $line);
    next unless ($group_hash{$genA} and $group_hash{$genB});
    for my $i (0 .. $#groups){
        my $group = $groups[$i];
        next if ($group_hash{$genA}{$group} == 0 or $group_hash{$genB}{$group} == 0) and !$opt_x; #don't compare any genomes that were given a value of zero with anything else.
        next if ($group_hash{$genA}{$group} == 0 and $group_hash{$genB}{$group} == 0);
        my $group_mem = "diff";
        if ($group_hash{$genA}{$group} == $group_hash{$genB}{$group}){
            $group_mem = "same";
        }
        push @{$results{$group}{$group_mem}{"snps"}}, $snps;
        push @{$results{$group}{$group_mem}{"shared"}}, $shared;
        push @{$results{$group}{$group_mem}{"frac"}}, $snps/$shared;
        push @{$results{$group}{$group_mem}{"frac_per_100k"}}, $snps/($shared/100000); #snps per 100,000 shared
        push @{$results{$group}{$group_mem}{"gens"}}, "$genA:$genB";
        
        #debug
        if ($group eq "species" and $group_mem eq "diff" and $shared > 1100000){
            print STDERR "$line\n";
        }
    }
}
close ($rin);

open (my $out_stats, "> $pref.groups.stats.txt");
print $out_stats "group\tsame/diff\tmeasure\tsum_obs\tnum_obs\tminimum\tmaximum\tmean\tmedian\tmode\tmode_freq\n";
my @measures = qw(gens snps shared frac frac_per_100k);
my @opps = qw(same diff);
my @for_raw;
my ($for_raw_col, $for_raw_max) = (0) x 2;
foreach my $group (@groups){
    foreach my $opp (@opps){
        foreach my $meas (@measures){
            my @stat_results = (0) x 8;
            my @tmp;
            if ($results{$group}{$opp}{$meas}){
                @tmp = @{$results{$group}{$opp}{$meas}};
                @stat_results = stats(\@tmp) unless $meas eq "gens";
            }
            print $out_stats "$group\t$opp\t$meas\t", join("\t", @stat_results), "\n" unless $meas eq "gens";
            push @{$for_raw[$for_raw_col]}, ($group, $opp, $meas);
            my $tmp_size = 3;
            if (@tmp){
                $tmp_size += scalar @tmp;
                push @{$for_raw[$for_raw_col]}, @tmp;
            }
            if ($tmp_size > $for_raw_max){
                $for_raw_max = $tmp_size;
            }
            $for_raw_col ++;
        }   
    }    
}
close ($out_stats);

open (my $out_raw, "> $pref.groups.raw.txt");
for my $i (0 .. ($for_raw_max - 1)){
    for my $j (0 .. ($for_raw_col - 1)){
        if (exists $for_raw[$j][$i]){
            print $out_raw "$for_raw[$j][$i]";
        }
        if ($j == ($for_raw_col - 1)){
            print $out_raw "\n";
        } else {
            print $out_raw "\t";
        }
    }
}
close ($out_raw);

#-----------------------------------
sub stats{
    my @lengths = @{$_[0]};
    my ($sum, $num, $min, $maxi, $mean, $median, $mode, $mode_freq);
    my %seen;
    my @sorted_leng = sort {$a <=> $b} @lengths;
    for my $i (0 .. $#sorted_leng){
        $sum += $sorted_leng[$i];
        $seen{$sorted_leng[$i]}++;
    }
    $num = $#sorted_leng + 1;
    $min = $sorted_leng[0];
    $maxi = $sorted_leng[$#sorted_leng];
    $mean = $sum/$num;
    my $rounded_mean = sprintf("%.2f", $mean);
    my @modes;
    foreach my $leng (sort {$seen{$b} <=> $seen{$a}} keys %seen) {
        push @modes, ([$leng, $seen{$leng}]);
    }
    $mode = $modes[0][0];
    $mode_freq = $modes[0][1];
    my $mid = int @sorted_leng/2;
    if (@sorted_leng % 2){
        $median = $sorted_leng[$mid];
    } else {
        $median = ($sorted_leng[$mid-1] + $sorted_leng[$mid])/2;
    }
    return ($sum, $num, $min, $maxi, $rounded_mean, $median, $mode, $mode_freq);
}
