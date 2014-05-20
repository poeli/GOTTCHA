#!/usr/bin/env perl
#[COPYRIGHT]
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
# CHANGE LOG:
#   Version 0.9a:
#     05/12/2014 1.Add '--stDir' option for pre-splitrimmed input file
#
#     05/09/2014 1.This version used a pipeline version of
#                  profileGOTTCHA.pl.
#
#   Version 0.9:
#     04/30/2014 1.Provide bwaOpt option for user to use their own parameters to
#                  run BWA-MEM.
#                2.Use absolute path to run system calls
#                3.Provide relAbu option to choose column to calculate relative
#                  abundance.
#

use Getopt::Long;
use FindBin qw($RealBin);
use strict;

# environment setup
my $ver = "0.9a";
$ENV{PATH} = "$RealBin:$RealBin/../ext/bin:$ENV{PATH}";
$ENV{PERL5LIB} = "$RealBin/../ext/lib/perl5:$ENV{PERL5LIB}";

$|=1;
my %opt;
my $res=GetOptions(\%opt,
    'input|i=s',
    'database|d=s',
    'dbLevel|l=s',
    'mode|m=s',
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
    'debug',
    'help|h|?') || &usage();

if ( $opt{help} || !-e $opt{input} || ! defined $opt{database} ) { &usage(); }

my ($fn) = $opt{input} =~ /([^\/]+)\.[^\.]+$/;
$fn ||= "gottcha_output";

my $time = time;
my $ct = &timeInterval($time);
print "[$ct] Starting GOTTCHA v$ver\n";

#-------GENERAL OPTIONS-------
my $INPUT      = $opt{input};
my $DB         = $opt{database};
my $DBLVL      = $opt{dbLevel};
my $THREADS    = defined $opt{threads} ? $opt{threads} : 2;
my $PREFIX     = defined $opt{prefix}  ? $opt{prefix}  : $fn;
my $OUTDIR     = defined $opt{outdir}  ? $opt{outdir}  : ".";
my $TMPDIR     = "${PREFIX}_temp";
my $RELABU     = defined $opt{relAbu}  ? $opt{relAbu} : "LINEAR_DOC";
my $MODE       = defined $opt{mode}    ? $opt{mode} : "summary";
my $BWAMETHOD  = "mem";
my $BWA_OPT    = defined $opt{bwaOpt}  ? $opt{bwaOpt} : "-k 30 -T 0 -B 100 -O 100 -E 100";
#-------SPLITRIM OPTIONS-------
my $STDIR      = defined $opt{stDir}   ? $opt{stDir}   : "$TMPDIR/splitrim";
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
my $LOGFILE    = "$PREFIX.gottcha.log";

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
    if( defined $db_level{$DBLVL} ){
        print "[$ct] Auto set database level to ".uc($DBLVL).".\n";
    }
    else{
        die "ERROR: Please specify database level using --dbLevel option.\n";
    }
}

#open (STDOUT, "| tee -a $LOGFILE");

$ct = &timeInterval($time);
&executeCommand("rm $LOGFILE | touch $LOGFILE", "[$ct] ERROR: Failed to create logfile: $LOGFILE."); 

# check running environment
$ct = &timeInterval($time);
print "[$ct] Checking running environment...\n";
&checkRunningEnv();
$ct = &timeInterval($time);
print "[$ct] All required scripts and tools found.\n";

# run splitrim
$ct = &timeInterval($time);

# check if pre-splitrimmed FASTQ provided
if( defined $opt{stDir} ){
    print "[$ct] Pre-splitrimmed directory is specified to $STDIR. Skip split-trimming.\n";
    my $file = &checkFileAbsence("$STDIR/${PREFIX}_splitrim.fastq", "$STDIR/${PREFIX}_splitrim.stats.txt");
    if( $file ){
        die "[$ct] ERROR: Can't find $file in $STDIR directory.\n";
    }
}
else{
    print "[$ct] Split-trimming input reads (fixL=$TRIM_FIXL, minQ=$TRIM_MINQ, ascii=$TRIM_ASCII)\n";
    &executeCommand("splitrim              \\
                       --inFile=$INPUT     \\
                       --fixL=$TRIM_FIXL   \\
                       --recycle           \\
                       --ascii=$TRIM_ASCII \\
                       --minQ=$TRIM_MINQ   \\
                       --prefix=${PREFIX}  \\
                       --outPath=$STDIR"
                    , "[$ct] Failed running splitrim. Please check $LOGFILE for detail.");
    $ct = &timeInterval($time);
    print "[$ct] Done splitrim.\n";
}

# run bwa
$ct = &timeInterval($time);
print "[$ct] Mapping split-trimmed reads to GOTTCHA database and profiling...\n";
my $BWA_DEBUG = "| tee $OUTDIR/$TMPDIR/$PREFIX.sam" if $DEBUG_MODE;
&executeCommand("bwa $BWAMETHOD $BWA_OPT -t $THREADS $DB $STDIR/${PREFIX}_splitrim.fastq $BWA_DEBUG \\
                 | profileGottcha.pl                                          \\
                     --parsedDB=$DB.parsedGOTTCHA.dmp                         \\
                     --sam                                                    \\
                     --outdir=$OUTDIR/$TMPDIR                                 \\
                     --prefix=$PREFIX                                         \\
                     --trimStats=$STDIR/${PREFIX}_splitrim.stats.txt          \\
                     --treeFile=$DBPATH/speciesTreeGI.dmp                     \\
                     --genomeVitals=$DBPATH/genomeVitals.dmp                  \\
                     --noFastqOut --method=1"
                , "[$ct] ERROR: Failed running BWA or profileGottcha.pl. Please check $LOGFILE for detail.");
$ct = &timeInterval($time);
print "[$ct] Done result profiling.\n";

# run filterGottcha.pl
$ct = &timeInterval($time);
print "[$ct] Filtering profiling results...\n";
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
                   --minHits=$FIL_MINH > $OUTDIR/$TMPDIR/$PREFIX.ABUX.out"
                , "[$ct] ERROR: Failed running filterGottcha.pl. Please check $LOGFILE for detail. ");
$ct = &timeInterval($time);
print "[$ct] Done filtering.\n";

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
    &executeCommand("hash splitrim 2>/dev/null", "[$ct] ERROR: splitrim not found. Please check your path or run INSTALL.sh.");
    &executeCommand("hash bwa 2>/dev/null", "[$ct] ERROR: bwa not found. Please check your path or run INSTALL.sh.");
    &executeCommand("hash profileGottcha.pl 2>/dev/null", "[$ct] ERROR: profileGottcha.pl not found. Please check your path or run INSTALL.sh.");
    &executeCommand("hash filterGottcha.pl 2>/dev/null", "[$ct] ERROR: filterGottcha.pl not found. Please check your path or run INSTALL.sh.");
    &executeCommand("hash convert_abu2list.pl 2>/dev/null", "[$ct] ERROR: convert_abu2list.pl not found. Please check your path or run INSTALL.sh.");

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
        &executeCommand( "profileGottcha.pl --db=$DB --make_dmp --prefix=$DB"
                           , "[$ct] ERROR: Failed generating $DB.parsedGOTTCHA.dmp. Please check $LOGFILE for detail.");
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
    unless ( $exit_code == 0 ){
        print "$msg\n";
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

sub usage {
print <<__END__;
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

    --input|i    <STRING>  Input a single-ended FASTQ file.
    --database|d <STRING>  The path of signature database. The database can be
                           in FASTA format or BWA index (5 files).

[OPTIONS]

  *** GENERAL OPTIONS ***

    --threads|t  <INT>     Number of threads [default: 2]
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
    --bwaOpt|b   <STRING>  BWA-MEM in this script is used to map input reads to 
                           GOTTCHA database. If you want to run it with your own
                           parameters, use this option to specify.
                           [default: "-k 30 -T 0 -B 100 -O 100 -E 100"]
    --stDir|s    <STRING>  Specify a directory contains pre-splitrimmed input
                           files. E.g. input file is "test.fastq", the script
                           will looking for "test_splitrim.fastq" and
                           "test_splitrim.stats.txt" in the specified directory.
    --help/h/?             display this help                   

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

__END__
exit();
}

