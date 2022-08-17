#!/usr/bin/perl

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
  
  -t    threads
        (default: 80)

";

die $usage unless @ARGV;

# command line processing
use Getopt::Std;
our ($opt_g, $opt_a, $opt_p, $opt_s, $opt_u, $opt_o, $opt_m, $opt_i, $opt_e, $opt_b, $opt_t);
getopts('ga:psuo:m:i:e:b:t:');
my $pref = $opt_o ? $opt_o: "output";
my $count_gaps;
$count_gaps = 1 if $opt_g;
my $afile = $opt_a if $opt_a;
my $bfile = $opt_b if $opt_b;
my $frac = $opt_m ? $opt_m : 0;
die "ERROR: -f must be a number between 0 and 1, inclusive\n" if $frac =~ m/D/ or $frac < 0 or $frac > 1;
my $in_file = $opt_i if $opt_i;
my $ex_file = $opt_e if $opt_e;
my $threads = $opt_t ? $opt_t : 80;
my $num_threads_running = 0;

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

my $mfile = $ARGV[0];
open (my $in, "<$mfile") or die "ERROR: Can't open matrix file: $!\n";
my %gen_ids;
my @gens;
my @snparrays;
#my ($id, $seq);
my $count = 0;
my $excluded = 0;
my $laststring = 1;
#my @gens_per_pos;
my @gaps_per_gen;
my @gaps_per_pos;
#my @bases_per_pos;


my ($id, $seq, $idtell, $seqtell);
my ($seqleng, $lasttell) = (0) x 2;

while (my $faline = <$in>){
    my $tell = tell $in;
    my @lines;
    my $lastpos = -1;
    while ($faline =~ /\R/gi){
        my $crpos = pos($faline)-1;
        my $sub = substr($faline, $lastpos + 1, $crpos - ($lastpos + 1));
        my $newtell = $lasttell + $crpos + 1;
        push @lines, ([$sub, $newtell]); #$newtell is the file position of the CR character representing the end of a line represented by $sub
        $lastpos = $crpos;
    }
    unless (($lastpos + 1) + $lasttell == $tell){ #in case the line didn't end with a CR character (which should be the case nearly 100% of the time)
        my $sub = substr($faline, $lastpos + 1, $tell - (($lastpos + 1) + $lasttell));
        push @lines, ([$sub, $tell]);
    }
    foreach my $slice (@lines){
        my ($sub, $subtell) = @{$slice};
        if ($sub =~ m/^\s*>/){
            if ($id){
                #sequence positions
                my $dist = $seqtell - $idtell; #total length of the sequence (including spaces and non-printing characters)
                my $start = $idtell; #starting position in the file of the sequence, i.e. right after the id line
                push @snparrays, ([$start, $dist]);
                #utf8::downgrade($seq); #try with and without for speed
                my $tleng = length($seq);
                if ($count == 0){
                    @gaps_per_pos = (0) x $tleng;
                    $count = $tleng;
                }
                if ($count > 0 and $tleng != $count){
                    die "ERROR: number of positions not identical between records.\n";
                }
                #$seq = uc($seq);
                #my %gaps;
                #my @garray;
                #my $gapcount = 0;
                while ($seq =~ /[-NX]/gi){
                    $gaps_per_pos[$-[0]]++;
                    #$gapcount++;
                    #$gaps{$-[0]} = 1;
                    #push @garray, $-[0];
                }
                #$gaps{'total'} = $gapcount;
                #push @gaps_per_gen, \%gaps;
                #push @gaps_per_gen, [@garray];
            }
            $seq = "";
            $sub =~ s/^\s*|\s*$//g; #remove any leading or trailing spaces / non-printing characters from the id
            $id = substr($sub, 1);
            $idtell = $subtell; #set idtell as the end of the id line
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
        unless (!$id){
            $sub =~ s/\s//g;
            $seq .= $sub;
            $seqtell = $subtell; #set seqtell as the end of the last sequence line in the sequence record
        }
        next;
    }
    $lasttell = $tell;
}
close ($in);
if ($id){
    #sequence positions
    my $dist = $seqtell - $idtell; #total length of the sequence (including spaces and non-printing characters)
    my $start = $idtell; #starting position in the file of the sequence, i.e. right after the id line
    push @snparrays, ([$start, $dist]);
    #utf8::downgrade($seq); #try with and without for speed
    my $tleng = length($seq);
    if ($count == 0){
        @gaps_per_pos = (0) x $tleng;
    }
    if ($count > 0 and $tleng != $count){
        die "ERROR: number of positions not identical between records.\n";
    }
    $count = $tleng;
    #$seq = uc($seq);
    my %gaps;
    #my @garray;
    #my $gapcount = 0;
    while ($seq =~ /[-NX]/gi){
        $gaps_per_pos[$-[0]]++;
        #$gapcount++;
        #$gaps{$-[0]} = 1;
        #push @garray, $-[0];
    }
    #$gaps{'total'} = $gapcount;
    #push @gaps_per_gen, \%gaps;
    #push @gaps_per_gen, [@garray];
} else {
    die "ERROR: No fasta records found in file.\n" if !@gens;
}
$seq = "";

print STDERR "\nDone!\nTotal sequence length: $count\n";
print STDERR "Excluded $excluded genome(s)\n";

my @gens_per_pos = (scalar @snparrays) x $count;
for my $i (0 .. $#gaps_per_pos){
    $gens_per_pos[$i] -= $gaps_per_pos[$i];
}

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
my @core = (0) x $count;
for my $i (0 .. $#gens_per_pos){
    if ($gens_per_pos[$i] >= $mingen){
        $frac_count++;
        $core[$i] = 1;
    }
}

print STDERR "Total loci in at least ",100 * $frac,"% of the input genomes ($mingen / $gencount): $frac_count\n";

my %diffs;
my %poscount;
my %diffpos;
for my $i (0 .. ($#gens - 1)){
    my $refgen = $gens[$i];
    my ($rstart, $rdist) = @{$snparrays[$i]};
    open (my $rin, "$mfile") or die "ERROR: Can't open $mfile: $!\n";
    my $refseq;
    seek $rin, $rstart, 0; #move to the sequence start site in the file
    read $rin, $refseq, $rdist; #read the sequence from the file;
    close ($rin);
    $refseq =~ s/\s//g;
    $refseq = uc($refseq);
    utf8::downgrade($refseq);
    # count gaps in the reference
    #my %refgaps = %{$gaps_per_gen[$i]};
    #my $refnumgaps = $refgaps{'total'};
    my %refgaps;
    my $refnumgaps = 0;
    #foreach (@{$gaps_per_gen[$i]}){
    #    $refgaps{$_} = 1;
    #    $refnumgaps++;
    #}
    #while ($refseq =~ /[-NX]/gi){
    #    $refgaps{$-[0]} = 1;
    #    $refnumgaps++;
    #}
    
    $refnumgaps = $refseq =~ tr/\-NX/-/;
    
    #while ($refseq =~ /[-NX]/g){
    #    $refgaps{$-[0]} = 1;
    #    $refnumgaps++;
    #}
    for my $j ($i+1 .. $#gens){
        my $qrygen = $gens[$j];
        my $string = "Comparing $refgen(". ($i+1) . ") - $qrygen(" . ($j+1) . ")"; #status update
        my @blanks = " " x $laststring; #status update
        print STDERR "\r", join("", @blanks), "\r$string"; #status update
        $laststring = length $string; #status update
        
        if ($num_threads_running == $threads){
            my $pid = wait;
            die "ERROR: pid $pid exited with status $?\n" if $?;
            open (my $tin, "<temp.kmtdm.$pid.txt");
            my ($lref, $lqry);
            while (my $line = <$tin>){
                chomp $line;
                my @tmp = split("\t", $line);
                my $type = shift @tmp;
                if ($type eq "gens"){
                    ($lref, $lqry) = @tmp;
                } elsif ($type eq "d"){
                    $diffs{$lref}{$lqry} = $tmp[0];
                } elsif ($type eq "s"){
                    $poscount{$lref}{$lqry} = $tmp[0];
                } elsif ($type eq "dp"){
                    my @vals = split(",", $tmp[0]);
                    @{$diffpos{$lref}{$lqry}} = @vals;
                }
            }
            close ($tin);
            unlink "temp.kmtdm.$pid.txt";
            $num_threads_running --;
        }
        my $pid = fork();
        if ($pid == 0){
            my ($qstart, $qdist) = @{$snparrays[$j]};
            open (my $qin, "$mfile") or die "ERROR: Can't open $mfile: $!\n";
            my $qryseq;
            seek $qin, $qstart, 0; #move to the sequence start site in the file
            read $qin, $qryseq, $qdist; #read the sequence from the file;
            close ($qin);
            $qryseq =~ s/\s//g;
            $qryseq = uc($qryseq);
            utf8::downgrade($qryseq);
            # count gaps in the query
            #my %qrygaps = %{$gaps_per_gen[$j]};
            #my $qrynumgaps = $qrygaps{'total'};
            my %qrygaps;
            my $qrynumgaps = 0;
            #foreach (@{$gaps_per_gen[$j]}){
            #    $qrygaps{$_} = 1;
            #    $qrynumgaps++;
            #}
            
            #my $sharedgaps = 0;
            #while ($qryseq =~ /[-NX]/gi){
            #    $qrygaps{$-[0]} = 1;
            #    $sharedgaps++ if exists $refgaps{$-[0]};
            #    $qrynumgaps++;
            #}
            
            $qrynumgaps = $qryseq =~ tr/\-NX/-/;
            
            #foreach my $gpos (keys %refgaps){
            #    next if $gpos eq "total";
            #    $sharedgaps++ if exists $qrygaps{$gpos};
            #}
            
            #while ($qryseq =~ /[-NX]/g){
            #    $qrygaps{$-[0]} = 1;
            #    $qrynumgaps++;
            #    if ($refgaps{$-[0]}){
            #        $sharedgaps++;
            #    }
            #}
            # calculate shared positions
            
            my $unsharedgaps = 0;
            my $runsharedgaps = 0;
            
            #$diffs{$refgen}{$qrygen} = 0;
            my $fdiffs = 0;
            my @fdiffpos;
            # xor the two strings
            my $mask = $refseq ^ $qryseq;
            while ($mask =~ m/[^\0]/g){
                my $rgap = 1 if substr($refseq, $-[0], 1) eq "-";
                my $qgap = 1 if substr($qryseq, $-[0], 1) eq "-";
                if (substr($refseq, $-[0], 1) eq "-" or substr($qryseq, $-[0], 1) eq "-"){
                    $unsharedgaps++;
                    $runsharedgaps++ if $rgap;
                    next unless $count_gaps;
                }
                #if ($refgaps{$-[0]} or $qrygaps{$-[0]}){
                #    next unless $count_gaps;
                #}
                next unless $core[$-[0]];
                #$diffs{$refgen}{$qrygen}++;
                $fdiffs++;
                #push @{$diffpos{$refgen}{$qrygen}}, $-[0] if $afile
                push @fdiffpos, $-[0] if $afile;
            }
            
            my $shared = $count - $unsharedgaps - ($refnumgaps - $runsharedgaps);
            $shared = $count - ($refnumgaps - $runsharedgaps) if $count_gaps;
            #my $shared = $count - $refnumgaps - $qrynumgaps + $sharedgaps;
            $poscount{$refgen}{$qrygen} = $shared;
            open (my $out, ">temp.kmtdm.$$.txt") or exit(1);
            print $out "gens\t$refgen\t$qrygen\n";
            print $out "d\t$fdiffs\n";
            print $out "s\t$shared\n";
            print $out "dp\t", join(",", @fdiffpos), "\n" if @fdiffpos;
            close ($out);
            exit;
        }
        $num_threads_running++;
    }
}
while (my $pid = wait){
    last if $pid == -1;
    die "ERROR: pid $pid exited with status $?\n" if $?;
    open (my $tin, "<temp.kmtdm.$pid.txt");
    my ($lref, $lqry);
    while (my $line = <$tin>){
        chomp $line;
        my @tmp = split("\t", $line);
        my $type = shift @tmp;
        if ($type eq "gens"){
            ($lref, $lqry) = @tmp;
        } elsif ($type eq "d"){
            $diffs{$lref}{$lqry} = $tmp[0];
        } elsif ($type eq "s"){
            $poscount{$lref}{$lqry} = $tmp[0];
        } elsif ($type eq "dp"){
            my @vals = split(",", $tmp[0]);
            @{$diffpos{$lref}{$lqry}} = @vals;
        }
    }
    close ($tin);
    unlink "temp.kmtdm.$pid.txt";
    $num_threads_running --;
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
    my ($rstart, $rdist) = @{$snparrays[$i]};
    open (my $rin, "$mfile") or die "ERROR: Can't open $mfile: $!\n";
    my $refseq;
    seek $rin, $rstart, 0; #move to the sequence start site in the file
    read $rin, $refseq, $rdist; #read the sequence from the file;
    close ($rin);
    $refseq =~ s/\s//g;
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
                my ($qstart, $qdist) = @{$snparrays[$j]};
                open (my $qin, "$mfile") or die "ERROR: Can't open $mfile: $!\n";
                my $qryseq;
                seek $qin, $qstart, 0; #move to the sequence start site in the file
                read $qin, $qryseq, $qdist; #read the sequence from the file;
                close ($qin);
                $qryseq =~ s/\s//g;
                for my $k (0 .. $#dpos){
                    my $pos = $dpos[$k];
                    my $rsnp = substr($refseq, $pos, 1);
                    my $qsnp = substr($qryseq, $pos, 1);
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
