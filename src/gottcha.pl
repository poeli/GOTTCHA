#!/usr/bin/env perl
# [COPYRIGHT]
#
# Copyright (2014).  Los Alamos National Security, LLC. This material was produced
#  under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Labora
# tory (LANL), which is operated by Los Alamos National Security, LLC for the U.S.
#  Department of Energy. The U.S. Government has rights to use, reproduce, and dis
# tribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY,
#  LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE US
# E OF THIS SOFTWARE.  If software is modified to produce derivative works, such m
# odified software should be clearly marked, so as not to confuse it with the vers
# ion available from LANL.
# 
# Additionally, this program is free software; you can redistribute it and/or modi
# fy it under the terms of the GNU General Public License as published by the Free
#  Software Foundation; either version 3 of the License, or (at your option) any l
# ater version. Accordingly, this program is distributed in the hope that it will 
# be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHA
# NTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public Licens
# e for more details.
#
# [DESCRIPTIUON]
#
# This is a wrapper to run GOTTCHA profiling tool.
#
# 2014/04/20
#
# Po-E (Paul) Li
# B-11
# Los Alamos National Lab.
# po-e@lanl.gov
#

use Getopt::Long;
use FindBin qw($RealBin);
use strict;

# environment setup
my $ver        = "1.0b";
$ENV{PATH}     = "$RealBin:$RealBin/../ext/bin:$ENV{PATH}";
$ENV{PERL5LIB} = "$RealBin/../ext/lib/perl5:$ENV{PERL5LIB}";

$|=1;
my %opt;
my $res=GetOptions(\%opt,
    'input|i=s',
    'database|db|d=s',
    'dbLevel|l=s',
    'mode|m=s',
    'noPlasmidHit|n',
	'bwaOpt|b=s',
    'relAbu|r=s',
    'threads|t=s',
    'stDir|s=s',
    'outdir|o=s',
    'prefix|p=s',
    'fixL=i',
    'minQ=i', 
    'ascii=i',
    'minCov=f',
    'minMLHL=i',
    'cCov=f',
    'taxLvl=s',
    'minLen=i',
    'minHits=i',
    'dumpSam',
    'debug',
    'help|h|?') || &usage();

if ( $opt{help} ) { &usage(); }
if ( !defined $opt{input} && !-e $opt{stDir} ) { &usage("ERROR: No reads found."); }
if ( !defined $opt{database} ) { &usage("ERROR: Please specify a gottcha database."); }

my ($fn) = $opt{input} =~ /([^\/.]+)\.[^\/]+$/;
$fn ||= "gottcha_output";

my $time = time;
my $ct = &timeInterval($time);
print "[$ct] Starting GOTTCHA v$ver\n";

#-------GENERAL OPTIONS-------
my $INPUT      = $opt{input};
my $DB         = $opt{database};
my $DBLVL      = $opt{dbLevel};
my $THREADS    = $opt{threads};
my $noPHit     = $opt{noPlasmidHit};
my $PREFIX     = defined $opt{prefix}  ? $opt{prefix}  : $fn;
my $OUTDIR     = defined $opt{outdir}  ? $opt{outdir}  : ".";
my $TMPDIR     = "${PREFIX}_temp";
my $RELABU     = defined $opt{relAbu}  ? $opt{relAbu} : "LINEAR_DOC";
my $MODE       = defined $opt{mode}    ? $opt{mode} : "summary";
my $BWAMETHOD  = "mem";
my $BWA_OPT    = defined $opt{bwaOpt}  ? $opt{bwaOpt} : "-k 30 -T 0 -B 100 -O 100 -E 100";
#-------SPLITRIM OPTIONS-------
my $STDIR      = defined $opt{stDir}   ? $opt{stDir}   : "$OUTDIR/$TMPDIR/splitrim";
my $TRIM_FIXL  = defined $opt{fixL}    ? $opt{fixL}    : 30;
my $TRIM_MINQ  = defined $opt{minQ}    ? $opt{minQ}    : 20;
my $TRIM_ASCII = defined $opt{ascii}   ? $opt{ascii}   : 33;
#-------FILTER OPTIONS--------
#my $FIL_TAXL   = defined $opt{taxLvl}  ? $opt{taxLvl}  : "species" ;
my $FIL_MINC   = defined $opt{minCov}  ? $opt{minCov}  : 0.005;
my $FIL_MINM   = defined $opt{minMLHL} ? $opt{minMLHL} : 5;
my $FIL_CCOV   = defined $opt{cCov}    ? $opt{cCov}    : 0.006;
my $FIL_MINL   = defined $opt{minLen}  ? $opt{minLen}  : 100;
my $FIL_MINH   = defined $opt{minHits} ? $opt{minHits} : 10;
#-------DEBUG---------
my $DEBUG_MODE = $opt{debug};
my $LOGFILE    = "$OUTDIR/$PREFIX.gottcha.log";

my ($DBPATH) = $opt{database} =~ /^(.*)\/[^\/]+$/;
$DBPATH ||= ".";

my %db_level = (
        'superkingdom' => 10,
        'phylum'       => 20,
        'class'        => 30,
        'order'        => 40,
        'family'       => 50,
        'genus'        => 60,
        'species'      => 70,
        'strain'       => 80
        #'replicon'     => 90
);

# database level
$DBLVL = lc $DBLVL;
unless ( defined $db_level{$DBLVL} ){
    ($DBLVL) = $opt{database} =~ /\.(\w+)$/;
	$DBLVL = lc $DBLVL;
    if( defined $db_level{$DBLVL} ){
        print "[$ct] Auto set database level to ".uc($DBLVL).".\n";
    }
	else{
    	($DBLVL) = $opt{database} =~ /\.(\w+).\w+$/;
		$DBLVL = lc $DBLVL;
    	if( defined $db_level{$DBLVL} ){
        	print "[$ct] Auto set database level to ".uc($DBLVL).".\n";
		}
    	else{
        	die "ERROR: Please specify database level using --dbLevel option.\n";
    	}
	}
}

# number of threads
$ct = &timeInterval($time);
if( !defined $THREADS ){
	print "[$ct] Number of threads is not specified. Autodetecting CPU numbers...\n";
	my $autocpu = `cat /proc/cpuinfo | grep processor | wc -l`;
	chomp $autocpu;
	if( $autocpu > 0 ){
		print "[$ct] $autocpu processor(s) detected.\n";
		$THREADS = $autocpu;
	}
	else{
		print "[$ct] Failed to detect the number of processor(s). Use default value (2).\n";
		$THREADS = 2;
	}
}
print "[$ct] Number of threads: $THREADS\n";

#open (STDOUT, "| tee -a $LOGFILE");
$ct = &timeInterval($time);
`mkdir -p $OUTDIR/$TMPDIR`;
&executeCommand("rm $LOGFILE | touch $LOGFILE", "ERROR: Failed to create logfile: $LOGFILE."); 

# check running environment
$ct = &timeInterval($time);
print "[$ct] Checking running environment...\n";
&checkRunningEnv();
$ct = &timeInterval($time);
print "[$ct] Done. All required scripts and tools found.\n";

# run splitrim
my @fastqs = split /,/, $INPUT;
my @st_fastq_files;
my @st_stats_files;

# check if pre-splitrimmed FASTQ provided
if( defined $opt{stDir} ){
    print "[$ct] Pre-splitrimmed directory is specified to $STDIR. Skip split-trimming step.\n";
	foreach my $fastq ( @fastqs ){
		my ($p) = $fastq =~ /^([^\/]+)\.\w+$/;
    	my $file = &checkFileAbsence("$STDIR/${p}_splitrim.fastq", "$STDIR/${p}_splitrim.stats.txt");
    	if( $file ){
        	die "[$ct] ERROR: Can't find $file in $STDIR directory.\n";
    	}
		else{
			push @st_fastq_files, "$STDIR/${p}_splitrim.fastq";
			push @st_stats_files, "$STDIR/${p}_splitrim.stats.txt";
		}
	}
}
else{
	print "[$ct] Split-trimming with parameters fixL=$TRIM_FIXL, minQ=$TRIM_MINQ, ascii=$TRIM_ASCII.\n";
	foreach my $fastq ( @fastqs )
	{
		my ($p) = $fastq =~ /^([^\/]+)\.\w+$/;
		print "[$ct] Split-trimming: $fastq...\n";
		&executeCommand("splitrim                \\
		                   --inFile=$fastq       \\
		                   --fixL=$TRIM_FIXL     \\
		                   --recycle             \\
		                   --ascii=$TRIM_ASCII   \\
		                   --minQ=$TRIM_MINQ     \\
		                   --prefix=$p           \\
		                   --outPath=$STDIR"
		                , "Failed running splitrim. Please check $LOGFILE for detail.");
		push @st_fastq_files, "$STDIR/${p}_splitrim.fastq";
		push @st_stats_files, "$STDIR/${p}_splitrim.stats.txt";
		$ct = &timeInterval($time);
    	print "[$ct] Done splitrimming $fastq.\n";
	}
}

# merging
my $merged_stats;
foreach my $stats ( @st_stats_files ){
	open STATS, $stats or die "Can't open splitrim stats: $!\n";
	my @lines = <STATS>;
	my $tempidx=0;
	foreach my $line ( @lines[3..7] ){
		my @tmp = split(/\t/, $line);
		$merged_stats->{$tempidx}->{"$tmp[0]"}->{LINE} = $line;
		$merged_stats->{$tempidx}->{"$tmp[0]"}->{1} += $tmp[1];
		$merged_stats->{$tempidx}->{"$tmp[0]"}->{2} += $tmp[2];
		$tempidx++;
	}
	close STATS;
}

open STATS, ">$OUTDIR/$TMPDIR/${PREFIX}_splitrim.stats.txt" or die "Can't open splitrim stats: $!\n";
print STATS "\n\n\n\n";
foreach my $idx ( sort {$a<=>$b} keys %$merged_stats ){
	foreach my $header ( keys %{$merged_stats->{$idx}} ){
		if( $merged_stats->{$idx}->{$header}->{1} ){
			if( $header =~ /^Mean/ ){
				$merged_stats->{$idx}->{$header}->{1} /= scalar @fastqs;
				$merged_stats->{$idx}->{$header}->{2} /= scalar @fastqs;
			}

			printf STATS "%s\t%d\t%d\t(%.2f %%)\n",
				$header,
				$merged_stats->{$idx}->{$header}->{1},
				$merged_stats->{$idx}->{$header}->{2},
				$merged_stats->{$idx}->{$header}->{2}/$merged_stats->{$idx}->{$header}->{1}*100
			;
		}
		else{
			print STATS $merged_stats->{$idx}->{$header}->{LINE};
		}
	}
}
close STATS;

$ct = &timeInterval($time);
print "[$ct] Done merging splitrim stats.\n";

#print splitrim summary
open STATS, "$OUTDIR/$TMPDIR/${PREFIX}_splitrim.stats.txt" or die "Can't open splitrim stats: $!\n";
my @lines = <STATS>;
$lines[4] =~ s/=+/=============/;
foreach my $line ( @lines[3..8], ){
	$line =~ s/\cJ//;
	my @temp = map { &thousandsSep($_) } split(/\t/, $line);
	printf "%s  %20s  %20s  %s\n", @temp[0..2], defined $temp[3] ? $temp[3] : "";
}
close STATS;
print "\n";

# run bwa and profiling results
$ct = &timeInterval($time);
print "[$ct] Mapping split-trimmed reads to GOTTCHA database and profiling...\n";
my $sam_output = "$OUTDIR/$TMPDIR/$PREFIX.sam";
$sam_output = "$OUTDIR/$PREFIX.gottcha.sam" if $opt{dumpSam}; 
my $BWA_DEBUG = "| tee $sam_output" if $DEBUG_MODE || $opt{dumpSam};
my $extra_opts = "-noPlasmidHit" if $noPHit;
my $realTHREADS = $THREADS-1 < 1 ? 1 : $THREADS-2;
my $st_filelist = join " ", @st_fastq_files;

&executeCommand("cat $st_filelist | bwa $BWAMETHOD $BWA_OPT -t $realTHREADS $DB - $BWA_DEBUG \\
                 | profileGottcha.pl $extra_opts                              \\
                     --parsedDB=$DB.parsedGOTTCHA.dmp                         \\
                     --sam                                                    \\
                     --outdir=$OUTDIR/$TMPDIR                                 \\
                     --prefix=$PREFIX                                         \\
                     --trimStats=$OUTDIR/$TMPDIR/${PREFIX}_splitrim.stats.txt \\
                     --treeFile=$DBPATH/speciesTreeGI.dmp                     \\
                     --genomeVitals=$DBPATH/genomeVitals.dmp                  \\
					 --noFastqOut --method=1 > $OUTDIR/$TMPDIR/profileGottcha.log 2>&1"
                , "ERROR: Failed running BWA or profileGottcha.pl. Please check $LOGFILE for detail.");
&executeCommand("cat $OUTDIR/$TMPDIR/profileGottcha.log >> $LOGFILE");

#print bwa and profiling summary
my ($bwa_read,
	$bwa_read_raw,
	$bwa_processed,
	$bwa_mapped,
	$bwa_mapped_plasmid,
	$bwa_unmapped,
	$profile_taxa,
	$bwa_mapped_raw,
	$bwa_mapped_raw_plasmid) = (0,0,0,0,0,0,0,0,0);

open RUNNINGLOG, "$LOGFILE" or die "Can't open log file: $!\n";
while(<RUNNINGLOG>){
    ($bwa_mapped,$bwa_mapped_plasmid,$bwa_unmapped, $bwa_mapped_raw, $bwa_mapped_raw_plasmid) = ($1,$2,$3,$4,$5) if /Mapped split-trimmed reads: (\d+); Mapped split-trimmed reads to plasmids: (\d+); Unmapped split-trimmed reads: (\d+); Mapped raw reads: (\d+); Mapped raw reads to plasmids: (\d+)/;
    ($bwa_read_raw, $bwa_read) = ($1,$2) if /found (\d+) reads \(split-trimmed: (\d+)\)/;
    $bwa_processed = $1 if /Processed (\d+) reads/;
    $profile_taxa = $1 if /\.$DBLVL\.tsv\.ABU"\.\.\.done\. \((\d+) taxonomie/i;
}
close RUNNINGLOG;

$ct = &timeInterval($time);

print  "\n";
print  "                                       RAW         SPLIT-TRIMMED\n";
print  "                             =============         =============\n";
printf "# of Processed Reads: %20s  %20s\n",                &thousandsSep($bwa_read_raw), &thousandsSep($bwa_read);
printf "   # of Mapped Reads: %20s  %20s (genome)\n",       &thousandsSep($bwa_mapped_raw), &thousandsSep($bwa_mapped);
printf "   # of Mapped Reads: %20s  %20s (plasmid only)\n", &thousandsSep($bwa_mapped_raw_plasmid),  &thousandsSep($bwa_mapped_plasmid);
printf " # of Unmapped Reads: %20s  %20s\n",                &thousandsSep($bwa_read_raw-$bwa_mapped_raw), &thousandsSep($bwa_unmapped);
print  "\n";
print  "[$ct] Done profiling mapping results.\n";
print  "[$ct] The noPlasmidHit option is ON. $bwa_mapped_plasmid plasmid hits will be ignored.\n" if $noPHit;
print  "\n  $profile_taxa taxanomy(ies) found.\n\n";

if( $bwa_mapped == 0){
    $ct = &timeInterval($time);
    print "[$ct] No read mapped to $DBLVL-level signatures. Please try again with upper-level databases.\n";
    exit;
}

# run filterGottcha.pl
$ct = &timeInterval($time);
print "[$ct] Removing taxanomy(ies) that has insufficient coverage...\n";
&executeCommand("filterGottcha.pl                             \\
                   --strain=$DBPATH/variantStrainLookup.dmp   \\
                   --species=$DBPATH/variantSpeciesLookup.dmp \\
                   --genus=$DBPATH/genusLookupBySpecies.dmp   \\
                   --family=$DBPATH/familyLookupByGenus.dmp   \\
                   --order=$DBPATH/orderLookupByFamily.dmp    \\
                   --class=$DBPATH/classLookupByOrder.dmp     \\
                   --phylum=$DBPATH/phylumLookupByClass.dmp   \\
                   --taxlookup=$DBPATH/taxLookupBySpecies.dmp \\
                   --taxLevel=$DBLVL                          \\
                   --dir=$OUTDIR/$TMPDIR                      \\
                   --prefix=$PREFIX                           \\
                   --abuField=$RELABU                         \\
                   --method=1                                 \\
                   --minCov=$FIL_MINC                         \\
                   --minMLHL=$FIL_MINM                        \\
                   --cCov=$FIL_CCOV                           \\
                   --minLen=$FIL_MINL                         \\
                   --minHits=$FIL_MINH |tee $OUTDIR/$TMPDIR/$PREFIX.ABUX.out"
                , "[$ct] ERROR: Failed running filterGottcha.pl. Please check $LOGFILE for detail. ");

open ABUX, "$OUTDIR/$TMPDIR/$PREFIX.ABUX.out" or die "Can't open filterGottcha.pl output: $!\n";
my $flag=0;
my $filtered_tax_cnt=0;
while(<ABUX>){
    if(/^$DBLVL\t/i){
        $flag=1; next;
    }
    elsif( $flag && !/^===/ ){
	$filtered_tax_cnt++;
    }
    elsif( $flag && /^===/ ){
        last;
    }
}
close ABUX;

$ct = &timeInterval($time);
print "[$ct] Done filtering.\n";
print "\n  $filtered_tax_cnt taxanomy(ies) left.\n\n";

if( $filtered_tax_cnt == 0){
    $ct = &timeInterval($time);
    print "[$ct] No taxanomy found in $DBLVL level. Please try again with upper-level databases.\n";
    exit;
}

#mode
$ct = &timeInterval($time);
print "[$ct] Preparing result in $MODE mode...\n";

&executeCommand("convert_abu2list.pl --sig_lvl $DBLVL --abu_col $RELABU --sort_col $RELABU $OUTDIR/$TMPDIR/$PREFIX.ABUX.out > $OUTDIR/$PREFIX.gottcha.tsv");
$ct = &timeInterval($time);
print "[$ct] Done genereting a summary report ($OUTDIR/$PREFIX.gottcha.tsv).\n";

if( $MODE =~ /(all|full)/i ){
    &executeCommand("convert_abu2list.pl --abu_col $RELABU --sort_col $RELABU $OUTDIR/$TMPDIR/*.tsv.ABU > $OUTDIR/$PREFIX.gottcha_full.tsv");
    $ct = &timeInterval($time);
    print "[$ct] Done generating the report in full mode ($OUTDIR/$PREFIX.gottcha_full.tsv).\n";
}
if( $MODE =~ /(summary|full)/i ){
    &executeCommand("rm -rf $OUTDIR/$TMPDIR");
    $ct = &timeInterval($time);
    print "[$ct] Done cleaning temporary files.\n";
}
else{
    print "[$ct] All outputs stored in $OUTDIR/$TMPDIR directory.\n";
}

$ct = &timeInterval($time);
print "[$ct] Finished.\n";

#close STDOUT;
################################################################################
#
# SUBROUTINES
#
################################################################################

sub checkRunningEnv {
    my $ct = &timeInterval($time);
    &executeCommand("hash splitrim 2>/dev/null", "ERROR: splitrim not found. Please check your path or run INSTALL.sh.");
    &executeCommand("hash bwa 2>/dev/null", "ERROR: bwa not found. Please check your path or run INSTALL.sh.");
    &executeCommand("hash profileGottcha.pl 2>/dev/null", "ERROR: profileGottcha.pl not found. Please check your path or run INSTALL.sh.");
    &executeCommand("hash filterGottcha.pl 2>/dev/null", "ERROR: filterGottcha.pl not found. Please check your path or run INSTALL.sh.");
    &executeCommand("hash convert_abu2list.pl 2>/dev/null", "ERROR: convert_abu2list.pl not found. Please check your path or run INSTALL.sh.");

    my $file = &checkFileAbsence( "$DBPATH/genomeVitals.dmp", "$DBPATH/speciesTreeGI.dmp" );
   
    if( $file ){
        my $ct = &timeInterval($time);
        die "[$ct] ERROR: $file not found. Did you run mkGottchaTaxTree.pl or download one? Please check README for detail.";
    }
    
    $file = &checkFileAbsence( 
            "$DBPATH/phylumLookupByClass.dmp",
            "$DBPATH/classLookupByOrder.dmp",
            "$DBPATH/orderLookupByFamily.dmp",
            "$DBPATH/familyLookupByGenus.dmp",
            "$DBPATH/genusLookupBySpecies.dmp",
            "$DBPATH/variantSpeciesLookup.dmp",
            "$DBPATH/variantStrainLookup.dmp",
            "$DBPATH/taxLookupBySpecies.dmp"
    );
    
    if( $file ){
        my $ct = &timeInterval($time);
        die "[$ct] ERROR: $file not found. Did you run makeVariantTaxLookups.pl or download them? Please check README for detail.";
    }
 
    if( my $file = &checkFileAbsence(
            "$DB.sa",
            "$DB.bwt",
            "$DB.ann",
            "$DB.pac",
            "$DB.sa",
            ) ){
        my $ct = &timeInterval($time);
        print "[$ct] BWA index not found. Try to generate one ...\n";
        &executeCommand("bwa index $DB");
        $ct = &timeInterval($time);
        print "[$ct] Done generating BWA index.\n";
    }
    
    if( my $file = &checkFileAbsence( "$DB.parsedGOTTCHA.dmp" ) ){
        my $ct = &timeInterval($time);
        print "[$ct] $DB.parsedGOTTCHA.dmp does not exist. Try to generate one ...\n";
        &executeCommand( "profileGottcha.pl --db=$DB --make_dmp"
                           , "ERROR: Failed generating $DB.parsedGOTTCHA.dmp. Please check $LOGFILE for detail.");
        $ct = &timeInterval($time);
        print "[$ct] Done generating $DB.parsedGOTTCHA.dmp.\n";
    }
    return;
}

sub executeCommand {
    my ($command,$msg) = @_;
    $msg ||= "the command failed.\n";
    my $debug_cmd = "set -x; " if $DEBUG_MODE;
    my $exit_code = system("( $debug_cmd$command ) >> $LOGFILE 2>&1");
	my $ct = &timeInterval($time);
    unless ( $exit_code == 0 ){
        print "[$ct] $msg\n";
        exit 1;
    }
}

sub checkFileAbsence {
    my @file = @_;
    foreach my $file (@file){
        return "$file" if ! -r $file;
    }
    return 0;
}

sub timeInterval{
    my $strtime = shift;
    $strtime = time - $strtime;
    return sprintf "%02d:%02d:%02d", int($strtime / 3600), int(($strtime % 3600) / 60), int($strtime % 60);
}

sub thousandsSep {
	my $num = shift;
	$num =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1,/g;
	return $num;
}

sub usage {
my $msg = shift;
print <<__END__;
$msg

PROGRAM: GOTTCHA metagenomic taxonomic profiling tool

VERSION: $ver

DESCRIPTION:

Genomic Origin Through Taxonomic CHAllenge (GOTTCHA) is an 
annotation-independent and signature-based metagenomic taxonomic profiling tool
that has significantly smaller FDR than other profiling tools. This Perl script 
is a wrapper to run the GOTTCHA profiling tool with pre-computed signature 
databases. The procedure includes 3 major steps: split-trimming the input data, 
mapping reads to a GOTTCHA database using BWA, profiling/filtering the result.

USAGE: $0 [OPTIONS] --input <FASTQ> --database <DATABASE_PATH>

    --input|i    <STRING>  Input one or multiple FASTQ file(s). Use comma (,)
                           to separate multiple input files.
	                       
    --database|d <STRING>  The path of signature database. The database can be
                           in FASTA format or BWA index (5 files).

[OPTIONS]

  *** GENERAL OPTIONS ***

    --threads|t  <INT>     Number of threads [default: auto-detect]
    --dbLevel|l  <STRING>  Specify the taxonomic level of the input database 
                           (e.g. family, species, genus, strain, etc.). The
                           value will be auto-detected if the input database
                           ended with levels (e.g. GOTTCHA_db.species).
                           [default: none]
    --outdir|o   <STRING>  Output directory [default: ./]
    --prefix|p   <STRING>  Filename prefix of the output.
                           [default: <INPUT_FILENAME_PREFIX>]
    --relAbu|r   <STRING>  The field will be used to calculate relative
                           abundance. You can specify one of the following
                           fields: "LINEAR_LENGTH", "TOTAL_BP_MAPPED",
                           "HIT_COUNT", "LINEAR_DOC".
                           [default: LINEAR_DOC]
    --mode|m     <STRING>  You can specify one of the output mode:
                           "summary" : this mode will report a summary of
                                       profiling result to *.gottcha.tsv file.
                           "full"    : other than a summary, this mode will
                                       report unfiltered result to
                                       *.gottcha_full.tsv with more detail.
                           "all"     : other than two tables, this mode will 
                                       keep all output files that were 
                                       generated by each profiling step.
                           [default: summary]
    --noPlasmidHit|n       Ignore alignments that hit to plasmids
                           [default: null]
    --bwaOpt|b   <STRING>  BWA-MEM in this script is used to map input reads to 
                           GOTTCHA database. If you want to run it with your own
                           parameters, use this option to specify.
                           [default: "-k 30 -T 0 -B 100 -O 100 -E 100"]
    --stDir|s    <STRING>  Specify a directory contains pre-splitrimmed input
                           files. E.g. input file is "test.fastq", the script
                           will looking for "test_splitrim.fastq" and
                           "test_splitrim.stats.txt" in the specified directory.
    --dumpSam              Dump the mapping result in SAM format.

  *** OPTIONS FOR SPLIT-TRIMMING READS ***

    --minQ       <INT>     Minimum quality for a read to be considered valid
                           (0-41) [default: 20]
    --fixL       <INT>     Fixed length to which each trimmed read will be cut 
                           down to [default: 30]
    --ascii      <INT>     ASCII encoding of quality score (33 or 64) [default: 
                           33] 

  *** OPTIONS FOR FILTERING PROFILING RESULT ***

    --minCov     <FLOAT>   Minimum linear coverage to be considered valid in 
                           abundance calculation [default: 0.005]
    --minMLHL    <INT>     Minimum Mean-Linear-Hit-Length to be considered valid
                           in abundance calculation [default: 5]
    --cCov       <FLOAT>   Critical coverage below which --minMLHL will cause an
                           organism to fail [default: 0.006]
    --minLen     <INT>     Minimum unique length to be considered valid in
                           abundance calculation [default: 100]
    --minHits    <INT>     Minimum number of hits to be considered valid in 
                           abundance calculation [10]
    
    --help/h/?             display this help                   

__END__
exit 1;
}

