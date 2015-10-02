#!/usr/bin/perl

use warnings;
use strict;
use Cwd 'abs_path';
use File::Path 'rmtree';
use File::Basename;

$|++;

my $version = 0.6;

my $usage = "
kmer_compare.pl

Determines SNP differences among pairs of sequences.

Required:
  -f    file containing sequences or kmers of input strains
        
        Strains may be provided assembled as fasta files or raw read files or
        both. To save time (if performing multiple comparisons), can also
        pre-process sequence files using 'jellyfish count' and 'jellyfish dump'
        and provide the resulting non-binary kmer quantity files. Just make sure
        all the kmer files have the same length kmers or else the script will
        exit with an error.
        
        VERY IMPORTANT: If providing pre-processed jellyfish files, you must
        ensure that jellyfish count was run with the '-C' (both-strands) option
        or this program's output will almost certainly be incorrect.
        
        Each line in the file should have three columns separated by tabs.
        1st column: Unique strain ID
        2nd column: Path to sequence or kmer file
        3rd column: File type. The possible file types are as follows
            fs    = assembled sequence file in fasta format
            fr    = sequence reads in fasta or fastq format
            js    = kmer file
            jr    = kmer file from reads
            If the third column entry does not match these types or is left
            blank, the file will be assumed to be of type 'fs'.
        
        The script accepts gzipped or uncompressed files.
        
        If multiple files are provided for a strain (i.e. forward and reverse
        reads), their paths must be provided on separate lines with the same
        strain ID (column 1). If both assembly and read sequences are given
        for a strain, only kmers found in the assembly will be compared and
        reads will only be used to identify conflicting kmers within a genome.
        
        Optional comments may follow any line as long as they are preceded by
        a '#'.
        
        Example file:
        genA    /path/to/genA.fasta fs          #genA is a sequence assembly
        genB    /path/to/genB.fasta fs          #genB is a sequence assembly
        genB    /path/to/genB_1.fastq.gz    fr  #including sequence reads for genB
        genC    /path/to/genC_1.fastq.gz    fr  #only providing reads for genC
        genC    /path/to/genC_2.fastq.gz    fr
        genD    /path/to/genD_kmers.txt     js  #file of kmers from the genD assembly
        genE    /path/to/genE_kmers.txt.gz  jr  #file of kmers from the genE reads (gzipped)
        
        OPTIONAL:
        If you would only like to compare one or more reference genomes to a
        group of other genomes, and not every possible all-vs-all pairwise
        comparison, then add a `*` character at the end of the desired reference
        genomes'line(s).
        
        For example:
        genA    /path/to/genA.fasta fs  *
        genB    /path/to/genB.fasta fs          
        genB    /path/to/genB_1.fastq.gz    fr
        genC    /path/to/genC_1.fastq.gz    fr
        genC    /path/to/genC_2.fastq.gz    fr
        genD    /path/to/genD_kmers.txt     js
        genE    /path/to/genE_kmers.txt.gz  jr
        
        In the above example, instead of all possible pairwise comparisons of
        the 5 genomes, the comparisons performed will be:
        genA - genB
        genA - genC
        genA - genD
        genA - genE
        done.
        
Optional:
  -c    File of kmers to which to limit the analysis, i.e. the file of core
        kmers output by 'kmer_core.pl'. This should be a fasta file of kmers
        with one or more N's separating the kmers. If this file is given, it
        will provide the kmer size and any entry to (-k) will be ignored.
        (default: all SNP kmers in common between two strains will be counted)
  -k    kmer size.  This setting will be ignored if jellyfish kmer files are
        provided or if a list of kmers is given using (-c)
        (default: 31)
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
  -o    output prefix
        (default: output)
  -t    threads
        (default: 15)
  -j    path to jellyfish
        (default: will look for jellyfish in PATH)

";

use Getopt::Std;
use vars qw( $opt_f $opt_c $opt_k $opt_t $opt_j $opt_o $opt_m $opt_e );
getopts('f:c:k:t:j:o:m:e:');

die $usage unless $opt_f;

my $core    = $opt_c if $opt_c;
my $ukmer   = $opt_k ? $opt_k : 31;
my $threads = $opt_t ? $opt_t : 15;
my $j_path  = $opt_j if $opt_j;
my $pref    = $opt_o ? $opt_o : "output";
my $min_k   = $opt_m ? $opt_m : 5;
my $err     = $opt_e ? $opt_e : .1;

die "ERROR: kmer size (-k) must be an odd number\n" unless ($ukmer % 2 == 1);

my $armleng;
my $kmer;
my @bases = qw(A C G T);

#load the mers array (4-mer)
#do I need this in this script or just the pairwise comparison subscript?
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

my $start_time = time;

mkdir "TemporaryFiles" or die "Can't make TemporaryFiles folder: $!\n";
my $temp_path = abs_path("TemporaryFiles");
print STDERR "temp_path: $temp_path\n";

my $script_path = abs_path($0);
$script_path = dirname($script_path);
print STDERR "script_path: $script_path\n";

#read in core, if given
my $is_core = "n";
my %core_kmers;
if ($opt_c){
    print STDERR "Reading kmer filter file ... ";
    open (my $in, "<", $core) or die "ERROR: Can't open $core: $!\n";
    my $seq;
    while (my $line = <$in>){
        chomp $line;
        $line =~ s/\s.*$//;
        if ($line =~ m/^>/){
            if ($seq){
                my @tmp = split(/N+/, $seq);
                foreach my $mer (@tmp){
                    if (!$kmer){
                        $kmer = length($mer);
                        $armleng = int($kmer / 2);
                    }
                    $mer = uc($mer);
                    my $base = substr($mer, $armleng, 1, ".");
                    my $rev = reverse($mer);
                    $rev =~ tr/ACTG/TGAC/;
                    if ($core_kmers{$rev}){
                        my $revbase = $base;
                        $revbase =~ tr/ACTG/TGAC/;
                        push @{$core_kmers{$rev}}, $revbase;
                        next;
                    }
                    push @{$core_kmers{$mer}}, $base;
                }
            }
            next;
        }
        $seq .= $line;
    }
    close ($in);
    my @tmp = split(/N+/, $seq);
    foreach my $mer (@tmp){
        if (!$kmer){
            $kmer = length($mer);
            $armleng = int($kmer / 2);
        }
        $mer = uc($mer);
        my $base = substr($mer, $armleng, 1, ".");
        my $rev = reverse($mer);
        $rev =~ tr/ACTG/TGAC/;
        if ($core_kmers{$rev}){
            my $revbase = $base;
            $revbase =~ tr/ACTG/TGAC/;
            push @{$core_kmers{$rev}}, $revbase;
            next;
        }
        push @{$core_kmers{$mer}}, $base;
    }
    my $num = scalar keys %core_kmers;
    print STDERR "$num kmers\n";
   
    my $out;
    my $last_mer;
    foreach my $c_kmer (sort {$a cmp $b} keys %core_kmers){
        my $mer = substr($c_kmer, 0, 4);
        if (!$last_mer or $mer ne $last_mer){
            close ($out) if $out;
            open ($out, "> TemporaryFiles/core_$mer.txt") or die "ERROR: Can't open TemporaryFiles/core_$mer.txt: $!\n";
        }
        print $out "$c_kmer\n";
        $last_mer = $mer;
    }
    close ($out);
    
    undef %core_kmers;
    $is_core = "y";
    #end multithreading output (get rid of all of the above if we don't want to multithread kmer filtering)
}

#read in sequence files
my %files;
my @order;
my $file_count = 0;
my (%id_num, %num_id);
my %refs;
open (my $f_in, "<", $opt_f) or die "ERROR: Can't open $opt_f: $!\n";
while (my $line = <$f_in>){
    chomp $line;
    next if $line =~ m/^\s*$/; #skip blank lines
    $line =~ s/\s*#.*$//; #remove comments
    my ($id, $file, $type) = split("\t", $line);
    if ($line =~ m/\*\s*$/){
        $refs{$id}++;
    }
    $type = "fs" unless ($type);
    $type = lc($type);
    $type = "fs" unless ($type =~ m/(?:fs|fr|js|jr)/);
    if ($type =~ m/^j/){
        #check if the kmer lengths match
        my $in;
        if ($file =~ m/\.gz$/){
            open ($in, "gzip -cd $file |") or die "ERROR: Can't open $file: $!\n";
        } else {
            open ($in, "< $file") or die "ERROR: Can't open $file: $!\n";
        }
        my $line = <$in>; #read just the first line
        chomp $line;
        close ($in);
        (my $seq) = $line =~ m/^([ACTG]+)\s+\d+$/;
        die "ERROR: file '$file' does not appear to be a jellyfish output file\n" if !$seq;
        my $test = length($seq);
        if ($kmer){
            die "ERROR: kmer lengths in file '$file' ($test) do not match kmer lengths in core file and/or other jellyfish kmer files ($kmer)\n" if ($test != $kmer);
        } else {
            $kmer = $test;
        }
    }
    unless ($files{$id}){ #don't put duplicate IDs in the order list
        push @order, $id;
        $file_count++;
        $id_num{$id} = $file_count;
        $num_id{$file_count} = $id;
    }
    push @{$files{$id}}, ([$file, $type]);
}
close ($f_in);

if (!$kmer){
    $kmer = $ukmer;
}
$armleng = int($kmer / 2);

my (@ref_order, @non_ref_order);
if (%refs){
    print STDERR "Reference sequences given:\n";
    foreach my $gen (@order){
        if ($refs{$gen}){
            push @ref_order, $gen;
            print STDERR "$gen\n";
        } else {
            push @non_ref_order, $gen;
        }
    }
    die "ERROR: Can't mark all sequences as reference\n" if scalar @ref_order == scalar @order;
} else {
    print STDERR "No reference sequences given. Will perform all-vs-all comparisons\n";
}

#process the input files into kmers and filter

### Let's further breakup and multithread the steps to increase speed
## First, we'll do all our jellyfishing
my $status_line_length = 1;
my %kmer_files;
print STDERR "Separating sequence(s) into $kmer-mers...\n";
foreach my $id (sort keys %files){
    my $id_number = $id_num{$id};
    my @members = @{$files{$id}};
    foreach my $rec (@members){
        my ($file, $type) = @{$rec};
        my $rs = substr($type, 1, 1);
        (my $shortfile) = $file =~ m/\/([^\/]+)$/;
        $shortfile = $file if !$shortfile;
        my @blanks = (" ") x $status_line_length;
        print STDERR "\r\t", join("", @blanks);
        my $outfile = "TemporaryFiles/kmers$id_number.$rs.txt";
        if ($type =~ m/^f/){
            my $status_line = "Running jellyfish on $shortfile";
            print STDERR "\r\t$status_line";
            $status_line_length = length($status_line);
            my $infile = $file;
            if ($file =~ m/\.gz$/){
                my $status = system("gzip -cd $file > TemporaryFiles/seq.txt");
                die "ERROR: Unzipping file $file with gzip exited with status $status" unless $status == 0;
                $infile = "TemporaryFiles/seq.txt";
            }
            my $status = system("jellyfish count -o TemporaryFiles/Jellytmp -m $kmer -s 1000000000 -t $threads -C $infile"); #run jellyfish count (count kmers in infput file) with kmer size -m of $kmer, hash size -s, output file prefix -o of Jellytmp, using number of threads -t and only outputting cannonical kmers (-C)
            die "ERROR: jellyfish count died with status $status\n" unless $status == 0;
            unlink ("TemporaryFiles/seq.txt") if -e "TemporaryFiles/seq.txt";
            my @jfiles = glob("TemporaryFiles/Jellytmp*");
            foreach my $jfile (@jfiles){
                my $status = system("jellyfish dump -tc $jfile >> $outfile");
                die "ERROR: command 'jellyfish dump -tc $jfile >> $outfile' exited with status $status\n" unless $status == 0;
                unlink ($jfile);
            }
        } else { #if a pre-processed file is given
            my $status_line = "Adding kmer file $shortfile";
            print STDERR "\r\t$status_line";
            $status_line_length = length($status_line);
            my $status;
            if ($file =~ m/\.gz$/){ #unzip gzipped files
                $status = system("gzip -cd $file >> $outfile");
                die "ERROR: command 'gzip -cd $file >> $outfile' exited with status $status\n" unless $status == 0;
            } else {
                $status = system("cat $file >> $outfile");
                die "ERROR: command 'cat $file >> $outfile' exited with status $status\n" unless $status == 0;
            }
        }
        push @{$kmer_files{$id_number}}, $outfile;
    }
}
print STDERR "\n";

## Second, we'll combine kmer files and output mers using a separate script that can be multithreaded
print STDERR "Processing and splitting kmer files...\n";
my $status = system("/usr/bin/perl","$script_path/lib/kmer_compare_split_mt.pl","$file_count","$armleng","$temp_path","$threads");
die "ERROR: kmer_compare_split_mt.pl exited with status: $status\n" if $status;

## Third, we'll filter the kmers one genome at a time (since that's how kmer_compare_filter_mt.pl is already written)
print STDERR "Filtering kmers...\n";
open (my $split_results_in, "< TemporaryFiles/split_results.txt") or die "ERROR: Can't open TemporaryFiles/split_results.txt: $!\n";
my %k_nums;
my %k_files;
while (my $line = <$split_results_in>){
    chomp $line;
    my ($id_num, $assm, $assm_count, $read, $read_count) = split(",", $line);
    my $id = $num_id{$id_num};
    print STDERR "**$id";
    print STDERR " ($assm_count unique kmers in sequence(s))" if $assm > 0;
    print STDERR " ($read_count unique kmers in reads)" if $read > 0;
    print STDERR "\n";
    my ($is_assm, $is_read) = ("n") x 2;
    $is_assm = "y" if $assm > 0;
    $is_read = "y" if $read > 0;
    my $status = system("/usr/bin/perl","$script_path/lib/kmer_compare_filter_mt.pl","$id_num","$is_assm","$is_read","$is_core","$min_k","$err","$temp_path","$threads");
    die "ERROR: kmer_compare_filter_mt.pl exited with status: $status\n" if $status;
    open (my $results_in, "< TemporaryFiles/counts.txt") or die "ERROR: Can't open TemporaryFiles/counts.txt: $!\n";
    my $line = <$results_in>;
    chomp $line;
    my @tmp = split(",", $line);
    close ($results_in);
    unlink("TemporaryFiles/gen$id_num/counts.txt");
    my $filt_a = shift @tmp;
    my $filt_r = shift @tmp;
    if (@tmp){
        @{$k_files{$id}} = @tmp;
    }
    print STDERR "\t", $filt_a + 1 ," kmers remain after filtering\n" if $filt_a >= 0;
    print STDERR "\t", $filt_r + 1 ," kmers remain after filtering reads\n" if $filt_r >= 0;
    $filt_a++;
    $filt_r++;
    $k_nums{$id}{'filt'}{'assm'} = $filt_a;
    $k_nums{$id}{'filt'}{'read'} = $filt_r;
}
close ($split_results_in);
unlink ("TemporaryFiles/split_results.txt");

# perform pairwise comparisons
# output key file for kmer_compare_pairwise_mt.pl
open (my $key_out, "> TemporaryFiles/key.txt") or die "ERROR: Can't open TemporaryFiles/key.txt: $!\n";
for my $i (0 .. $#order){
    my $id = $order[$i];
    my $num = $id_num{$id};
    my $count = 0;
    $count = $k_nums{$id}{'filt'}{'read'} if $k_nums{$id}{'filt'}{'read'};
    $count = $k_nums{$id}{'filt'}{'assm'} if $k_nums{$id}{'filt'}{'assm'};
    my $types = join(",", sort{$a cmp $b} @{$k_files{$id}}); #sort to make sure assemblies are parsed before reads
    my $is_ref = "-";
    if (%refs){
        if ($refs{$id}){
            $is_ref = "Y";
        } else {
            $is_ref = "N";
        }
    }
    print $key_out "$id\t$num\t$count\t$types\t$is_ref\n";
}
close ($key_out);

#calculate the number of permutations. Will do this externally since the bignum function is so slow, I don't want it in this script or the pairwise script
my $tot_perm;
if (%refs){
    $tot_perm = scalar @ref_order * scalar @non_ref_order;
} else {
    $file_count = scalar @order;
    $tot_perm = `/usr/bin/perl $script_path/kmer_compare_permute.pl $file_count`;
    die "ERROR: Command '/usr/bin/perl $script_path/kmer_compare_permute.pl $file_count' died with error $?\n" if $?;
    chomp $tot_perm;
}

#run pairwise comparisons
print STDERR "Performing pairwise comparisons...\n";
$status = system("perl $script_path/lib/kmer_compare_pairwise_mt.pl $tot_perm $temp_path $threads");
die "\nERROR: kmer_compare_pairwise_mt.pl exited with status: $status\n" if $status;
my %results;
open (my $results_in, "< TemporaryFiles/snp_counts.txt") or die "\nERROR: Can't open TemporaryFiles/snp_counts.txt: $!\n";
while (my $line = <$results_in>){
    chomp $line;
    my ($a_id, $b_id, $snps, $shared_kmers, $dist) = split("\t", $line);
    $results{$a_id}{$b_id}{"shared"} = $shared_kmers;
    $results{$a_id}{$b_id}{"snps"} = $snps;
    $results{$a_id}{$b_id}{"dist"} = $dist;
}
close ($results_in);
unlink("TemporaryFiles/snp_counts.txt");

rmtree(["TemporaryFiles"]);

open (my $out, "> $pref.kmer_compare.table.txt");
open (my $r_out, "> $pref.kmer_compare.raw.txt");
print $r_out "#genome_A\tgenome_A_kmers\tgenome_B\tgenome_B_kmers\tSNPs\tShared_kmers\tSNPs_per_Shared\n";
#output results
#snp counts
print $out "**SNPs\n";
if (%refs){
    print $out "\t", join("\t", @non_ref_order), "\n";
    for my $i (0 .. $#ref_order){
        my $a_id = $ref_order[$i];
        print $out "$a_id";
        for my $j (0 .. $#non_ref_order){
            my $b_id = $non_ref_order[$j];
            my $snps = 0;
            $snps = $results{$a_id}{$b_id}{"snps"} if $results{$a_id}{$b_id}{"snps"};
            $snps = $results{$b_id}{$a_id}{"snps"} if $results{$b_id}{$a_id}{"snps"};
            print $out "\t$snps";
        }
        print $out "\n";
    }
} else {
    print $out "\t", join("\t", @order), "\n";
    for my $i (0 .. $#order){
        my $a_id = $order[$i];
        print $out "$a_id";
        for my $j (0 .. $i){
            if ($i == $j){
                print $out "\t0\n";
                next;
            }
            my $b_id = $order[$j];
            my $snps = 0;
            $snps = $results{$a_id}{$b_id}{"snps"} if $results{$a_id}{$b_id}{"snps"};
            $snps = $results{$b_id}{$a_id}{"snps"} if $results{$b_id}{$a_id}{"snps"};
            print $out "\t$snps";
        }
    }
}
print $out "\n";

#shared counts
print $out "**Shared kmers\n";
if (%refs){
    print $out "\t", join("\t", @non_ref_order), "\n";
    for my $i (0 .. $#ref_order){
        my $a_id = $ref_order[$i];
        print $out "$a_id";
        for my $j (0 .. $#non_ref_order){
            my $b_id = $non_ref_order[$j];
            my $snps = 0;
            $snps = $results{$a_id}{$b_id}{"shared"} if $results{$a_id}{$b_id}{"shared"};
            $snps = $results{$b_id}{$a_id}{"shared"} if $results{$b_id}{$a_id}{"shared"};
            print $out "\t$snps";
        }
        print $out "\n";
    }
} else {
    print $out "\t", join("\t", @order), "\n";
    for my $i (0 .. $#order){
        my $a_id = $order[$i];
        print $out "$a_id";
        for my $j (0 .. $i){
            if ($i == $j){
                my ($a_num, $r_num) = (0) x 2;
                $a_num = $k_nums{$a_id}{'filt'}{'assm'} if $k_nums{$a_id}{'filt'}{'assm'};
                $r_num = $k_nums{$a_id}{'filt'}{'read'} if $k_nums{$a_id}{'filt'}{'read'};
                my $sum = $a_num + $r_num;
                print $out "\t$sum\n";
                next;
            }
            my $b_id = $order[$j];
            my $shared = 0;
            $shared = $results{$a_id}{$b_id}{"shared"} if $results{$a_id}{$b_id}{"shared"};
            $shared = $results{$b_id}{$a_id}{"shared"} if $results{$b_id}{$a_id}{"shared"};
            print $out "\t$shared";
        }
    }
}
print $out "\n";

#snps per shared
print $out "**SNPs per shared kmers\n";
if (%refs){
    print $out "\t", join("\t", @non_ref_order), "\n";
    for my $i (0 .. $#ref_order){
        my $a_id = $ref_order[$i];
        print $out "$a_id";
        for my $j (0 .. $#non_ref_order){
            my $b_id = $non_ref_order[$j];
            my $snps = 0;
            $snps = $results{$a_id}{$b_id}{"snps"} if $results{$a_id}{$b_id}{"snps"};
            $snps = $results{$b_id}{$a_id}{"snps"} if $results{$b_id}{$a_id}{"snps"};
            my $shared = 0;
            $shared = $results{$a_id}{$b_id}{"shared"} if $results{$a_id}{$b_id}{"shared"};
            $shared = $results{$b_id}{$a_id}{"shared"} if $results{$b_id}{$a_id}{"shared"};
            my $frac = 0;
            if ($shared > 0){
                $frac = $snps / $shared;
            }
            $frac = sprintf("%.4f", $frac);
            print $out "\t$frac";
            my ($a_count, $b_count) = (0) x 2;
            $a_count = $k_nums{$a_id}{'filt'}{'read'} if $k_nums{$a_id}{'filt'}{'read'};
            $a_count = $k_nums{$a_id}{'filt'}{'assm'} if $k_nums{$a_id}{'filt'}{'assm'};
            $b_count = $k_nums{$b_id}{'filt'}{'read'} if $k_nums{$b_id}{'filt'}{'read'};
            $b_count = $k_nums{$b_id}{'filt'}{'assm'} if $k_nums{$b_id}{'filt'}{'assm'};
            print $r_out "$a_id\t$a_count\t$b_id\t$b_count\t$snps\t$shared\t$frac\n";
        }
        print $out "\n";
    }
} else {
    print $out "\t", join("\t", @order), "\n";
    for my $i (0 .. $#order){
        my $a_id = $order[$i];
        print $out "$a_id";
        for my $j (0 .. $i){
            if ($i == $j){
                print $out "\t0.0000\n";
                next;
            }
            my $b_id = $order[$j];
            my $snps = 0;
            $snps = $results{$a_id}{$b_id}{"snps"} if $results{$a_id}{$b_id}{"snps"};
            $snps = $results{$b_id}{$a_id}{"snps"} if $results{$b_id}{$a_id}{"snps"};
            my $shared = 0;
            $shared = $results{$a_id}{$b_id}{"shared"} if $results{$a_id}{$b_id}{"shared"};
            $shared = $results{$b_id}{$a_id}{"shared"} if $results{$b_id}{$a_id}{"shared"};
            my $frac = 0;
            if ($shared > 0){
                $frac = $snps / $shared;
            }
            $frac = sprintf("%.4f", $frac);
            print $out "\t$frac";
            my ($a_count, $b_count) = (0) x 2;
            $a_count = $k_nums{$a_id}{'filt'}{'read'} if $k_nums{$a_id}{'filt'}{'read'};
            $a_count = $k_nums{$a_id}{'filt'}{'assm'} if $k_nums{$a_id}{'filt'}{'assm'};
            $b_count = $k_nums{$b_id}{'filt'}{'read'} if $k_nums{$b_id}{'filt'}{'read'};
            $b_count = $k_nums{$b_id}{'filt'}{'assm'} if $k_nums{$b_id}{'filt'}{'assm'};
            print $r_out "$a_id\t$a_count\t$b_id\t$b_count\t$snps\t$shared\t$frac\n";
        }
    }
}
print $out "\n";

close ($out);
close ($r_out);

unless (%refs){
    open (my $p_out, "> $pref.kmer_compare.dist.txt");
    my $num_seqs = scalar @order;
    print $p_out " $num_seqs\n";
    for my $i (0 .. $#order){
        my $a_id = $order[$i];
        my $id10char = $a_id;
        $id10char =~ s/^(\S+)\s.*/$1/;
        $id10char=sprintf("%-90.90s",$id10char);
        print $p_out "$id10char";
        for my $j (0 .. $#order){
            my $b_id = $order[$j];
            my $dist = 0;
            $dist = $results{$a_id}{$b_id}{"dist"} if $results{$a_id}{$b_id}{"dist"};
            $dist = $results{$b_id}{$a_id}{"dist"} if $results{$b_id}{$a_id}{"dist"};
            print $p_out "$dist ";
        }
        print $p_out "\n";
    }
    close ($p_out);
}
