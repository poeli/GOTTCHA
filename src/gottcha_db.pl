#!/usr/bin/env perl
#
#   gottcha_v1e13c.pl
#
# Created:   2011/05/09
# Updated:   2012/05/03
# Author:    Tracey Freitas, traceyf@lanl.gov
#
use strict;
use warnings;
no warnings 'portable';
use threads;
use threads::shared;
use Getopt::Long;
use Benchmark ':hireswallclock';
use Time::HiRes;
use File::Path;
use File::Basename;
use File::Copy;
use IO::Handle;
use YAML::XS;
use XML::Simple;
use Storable qw(store retrieve dclone);
use Statistics::Descriptive;
use Cwd qw(abs_path);
use Tie::IxHash;
use Algorithm::Combinatorics qw(variations_with_repetition);
use Devel::Size qw(size total_size);
use Sort::Key::Radix qw(rnkeysort ukeysort usort ssort);
use Sort::Key qw(keysort);
use Bit::Vector;
use feature 'switch';
use Fcntl qw(:flock SEEK_END);  # import LOCK_* and SEEK_END constants
#-------------------------------------------------------------------------------
# CONSTANTS (k-mer analysis)
my $OFFSET              = 1;
my $FIXEDWORDSIZE       = 0;
my $MINWORDSIZE         = 0;
my $MAXWORDSIZE         = 0;
my $STEPWORDSIZE        = 5;
my $DEFAULT_WORDSIZE    = 10;
my $DEFAULT_MINWORDSIZE = 2;
my $DEFAULT_MAXWORDSIZE = 25;
my $CRAZY_WORD_SIZE     = 1000; # formerly 200
my $LOADWORDS           = q{};
my $CLEANUPWORDS        = q{};
my $MINSTRINGPACKLEN    = 15; # if $k > than this number, use stringPack() on l-mers
my @nucs :shared = ("0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F");
#-------------------------------------------------------------------------------
# CONSTANTS (string matching)
my $PREFIX          = q{};
my $LOGFILENAME     = q{};
my $TAXIDXFILE      = q{};  # XML file specifying directory and subspecies names
my $LOADFASTA       = q{};
my $NOSTOREFASTA    = 0;
my $NTHREADS        = 1;    # No. of CPUs to use
my $PDIR            = ".";  # Default dir to where all partition files will be written
my $OUTDIR          = ".";  # Default dir to where all non-partition files will be written
my $DBLABEL         = "Bacteria"; 
my $UNIQUESBYREPLICONS = 0; # Export one 'uniques' file per replicon
my $UNIQUESBYORG       = 0; # Export one 'uniques' file per organism
my $UNIQUESSINGLE      = 0; # Export one 'uniques' file in total
my $COMPARESAMETAXLEVEL = 0;                            
my $TAXLEVEL            = "species"; # default taxonomic level comparison
my $SKIPNONSTDNUCS      = 0;   # Discard l-mers containing non-standard (non-A/T/G/C) nucleotides
my $MINUNIQUEFRAGLENGTH = 35;  # post-parsing filtering parameter of unique frags to report
my $MINUNIQUEFRAGLENRANGE = 35;
my $MIN_ORGS2COMPARE    = 2;   # After RestrictIn/Out, the min. no. of orgs to determine
                               # unique regions
my $NOSTATS             = 0;   # Flag. Produce stats (1) or not (0).
#my $OPTIMIZE            = q{}; # Optimize l-mer intersections for "mem" or "cpu"
my $STRPACKLEN          = 3;   # Length of each key in HREF $codon2char
my $PREFIXLENGTH        = 1;   # Default number of character prefixes for l-mer indexing
my $MINPARTITIONS       = 1;   # Default number of partitions to merge into larger partition
    # Set limits on $prefixLength so it doesn't get too large. $numPartitions tops out 
    # at 69 (64 codons + A/T/G/C/N). All combos are stored in memory, so this could crash 
    # the program if $prefixLength is too high.
my $MAXPREFIXLENGTH  = 8;
my $MAXPARTITIONS    = 69;     # 
my $MAXWRITETHREADS  = 20;
#-------------------------------------------------------------------------------
# CONSTANTS (statistics)
my $SKIPGC           = 0;
#-------------------------------------------------------------------------------
# GLOBALS (k-mer analysis)
my $offset           = $OFFSET;
my $wordSize         = q{};           # The string of word sizes input by user (comma-separated)
my $k                = 0;             # Word Size
my $bvLen            = 0;             # BitVector Length = 2*$k
my $baseNum          = 16;
my $fixedWordSize    = $FIXEDWORDSIZE;
my $minWordSize      = $MINWORDSIZE;
my $maxWordSize      = $MAXWORDSIZE;
my $stepWordSize     = $STEPWORDSIZE;
my $using_range      = 0;             # True with valid $minWordSize and $maxWordSize values
my $loadWords        = $LOADWORDS;    # Loads replicon binary lmer hashes from disk instead of
                                      #   calling lmerExtraction() sub
my $cleanupWords     = $CLEANUPWORDS; # Removes replicon binary lmer hashes from disk when done
my $prefix           = $PREFIX;
my $logfilename      = $LOGFILENAME;
my $wordSizeDetail        = q{};
my $uniqueFragLenDetail   = q{};
my @wordSizes             = ();
my @uniqueFragLenArray    = ();
my $logDetail             = q{};
my $LOGFILE               = q{};
my @reporting_filehandles = ();
#my %seen:shared;               # tracks l-mer processed replicons
my @intersectionJobs:shared;   # Each entry is a per-core, loadBalanced AREF of jobs
#my %rPartitionIntersectOpts :shared;# = &share({});
#my %lmerCoords:shared;         # KEY=org name; VAL = AREF of start positions
my %warning_log:shared;        # holds the warnings generated in each cycle and their frequency
#-------------------------------------------------------------------------------
# GLOBALS (partitioning)
my $partitionDir  = "partitions";
my $dmpExt        = "dmp";         # $partitionExt
my $prefixLength  = $PREFIXLENGTH; # number of prefixes to split into; user config
#my $numPartitions = $MINPARTITIONS;# prefixes are load-balanced across these no. of partitions
#my $char2posLookup= &share({});    # KEY: either ATGC(N) or the stringPacked() chars; VAL: counter
#my $pos2charLookup= &share({});
#my %prefix2partitionLookup :shared = ();
#my @charKeys                 : shared = ();
#my $maxCharKeys              : shared = 0;
#my %partitionHostCoords = ();
#my @partitionHostStartCoords : shared = ();
#my @partitionHostEndCoords   : shared = ();
my %bucketLookup:shared = ();       # KEY=>VAL:     partIdx=>numBucketsFilled (<= basenum**prefixLength)
my $MAX_BUCKETS = 0;            # 16^$prefixLength
#-------------------------------------------------------------------------------
# GLOBALS (benchmarking)
my $TEXT_SeqLoad                  = "LD";
my $TEXT_BitMasking               = "BM";    
my $TEXT_rcBitMasking             = "rcBM";
my $TEXT_Ranging                  = "R";
my $TEXT_rcRanging                = "rcR";
my $TEXT_Encoding                 = "E";
my $TEXT_rcEncoding               = "rcE";      
my $TEXT_kmerDecompositionFWD     = "kF";       #"DecompF";
my $TEXT_kmerDecompositionREV     = "kR";       #"DecompR";
my $TEXT_kmerDecompositionFWDwrap = "kFW";      #"DecompFW";
my $TEXT_kmerDecompositionREVwrap = "kRW";      #"DecompRW";
my $TEXT_Sort                     = "Sr";       # "Srt"
my $TEXT_Slice                    = "Sl";       # "Slice"
my $TEXT_Dereplication            = "Drp";      # "Derep"
my $TEXT_Prefixing                = "Prf";      # "Pref"
my $STORE_prefixBitMask           = "+PrfBM";   # "PrfBM"
my $STORE_kmerSuffixes            = "+Sf";      # "Suf"
my $STORE_kmerStarts              = "+St";      # "Start"
#my $STORE_TransformationVector    = "+Tv";
my $STORE_Level2BitMasks          = "+L2BM";
my $TEXT_TotalTime                = "tTime";
my $benchRez = 2;   # Benchmark resolution (after decimal point)
#-------------------------------------------------------------------------------
# GLOBALS (stored data structures)
my $seqBMfilename               = "seqBitMask.".$dmpExt;
my $kmers_eBV_hex_Filename      = "kmers_eBV_hex.".$dmpExt;
my $kmer_starts_Filename        = "kmer_starts.".$dmpExt;
my $uniqueSeqBMfilename         = "uniqueSeqBitMask.".$dmpExt;
#my $unionizedL2BMsFilename      = "uL2BMs.$dmpExt";    # Customized name in SUBTRACTION phase
#-------------------------------------------------------------------------------
# GLOBALS (string matching)
my %repliconTax = ();          # $replicon (Key) => { ORGNAME, REPLDESC, GENUS, SPECIES }
my %organisms   :shared;       # $organism (Key) => { $replicon => { SEQ    => $string,
                               #                                     LENGTH => $integer
                               #                                   },
                               #                              }
#my %contigHostPartitionLookup:shared = ();  # %hash = ( $repliconNick => { $partitionHostID => "" } );
my %orgLookupByContig:shared;
my %uniqueFragsByOrg:shared;   # shared Holds all the unique fragments, group by organism
my %fileNickLookup   :shared;  # KEY:replicon filename;           VALUE:replicon nickname (integer)
my %rFileNickLookup  :shared;  # KEY:replicon nickname (integer); VALUE:replicon filename
my %contigNickLookup :shared;  # KEY:contigID;                    VALUE:contig nickname (integer)
my %rContigNickLookup:shared;  # KEY:contig nickname (integer);   VALUE:contigID
my %contigID2filename:shared;  # KEY:contigID;                    VALUE:replicon filename
my %contigWordCount  :shared;  # KEY:contigID;                    VALUE:no. of l-mers extracted
my @processedOrgs  :shared;    # Processed orgs; taken from keys %uniqueFragsByOrg
my @fastaFilesOfUniques:shared;# Holds individual org's FASTA filenames of unique seqs
my @contigIDs      :shared;    # Holds full contig IDs/names inputted from XML file
my @contigFilenames:shared;    # Holds contig file basenames
my @orgNames       :shared;    # Holds full organism names inputted from XML file
my @processedContigs:shared;   # Holds contig nicknames after being l-mer processed
my @intersectionThreads = ();  # Holds the impending intersection jobIDs
my $numOrgs      = 0;
my $totalContigs = 0;
my %output_fh_hash :shared=(); # Holds the output filehandles for writing unique seqs to disk
#-------------------------------------------------------------------------------
my $totalProcessed :shared;    # scalar(@processedReplicons)
my $intersectCounter:shared=0; # Counter for completed lmer intersection jobs
my $taxIdxFile         = $TAXIDXFILE;
my $loadFasta          = $LOADFASTA;                                                    # not yet implemented
my $noStoreFasta       = $NOSTOREFASTA;
my $nThreads           = $NTHREADS;
my $pDir               = $PDIR;
my $outdir             = $OUTDIR;
my @valid_pDirs    : shared = ();   # Confirmed writeable directories for partition storage
my $num_pDirs      : shared = 0;
my %pDirNickLookup : shared = ();   # Use rNick to find its pDir code
my %pDirCodeLookup : shared = ();   # Use pDir code to find its home directory
my $seqDir             = q{}; # Populated in verifyInputOpts $outdir."/seqs";
my $fastadir           = "fasta";   # output subdir for FASTA files: $outdir/$fastadir
my $dbLabel            = $DBLABEL;
#my $skipNonStdNucs     = $SKIPNONSTDNUCS;
my $minUniqueFragLength= $MINUNIQUEFRAGLENGTH;
my $minUniqueFragLenRange = $MINUNIQUEFRAGLENRANGE;
my $wantErrorFreeSeqs  = 0;
my $noStats            = $NOSTATS;
#my $optimize           = $OPTIMIZE;
#my $optimization       = {
#    "mem" => 0,
#    "cpu" => 0,
#};
my $minStringPackLen   = $MINSTRINGPACKLEN;
my $stringPacking = q{};
my $updateDir     = "update";
my $storableDir   = "storable";
my $summaryFilename = q{};         # File of all user-specified-length FASTA sequences
#-------------------------------------------------------------------------------
# GLOBALS (output)
my $verbose     = 0;
my $vverbose    = 0;
my $help        = 0;
my $debug       = 0;
my $debug_kmerExtract = 1;          # kmer Extraction Debug
my $iDebug = 1;                     # Intersection Debug

my $exportWords = 0;
my $linebuffer      = 1000;
my $NUM_NUCLEOTIDES = 4;
my @day_of_week     = ("Sun","Mon","Tues","Wed","Thurs","Fri","Sat");
my $padLen = 9;     # pad zeroes up until these many places
#-------------------------------------------------------------------------------
# CONSTANTS (statistics)
my $skipGC           = $SKIPGC;
#-------------------------------------------------------------------------------
# GLOBALS (queues)
my $ripping_queue      = q{};
my $subtraction_queue  = q{};
my $ripCount         :shared = 0;
my $subtractionCount :shared = 0;
#-------------------------------------------------------------------------------
# GLOBALS (taxonomy)
my $taxTreeFile        = q{};         # Pre-computed Storable file of \%taxTree
my $taxLevel           = q{};         # Tax level to compare at (differences @ one level lower)
my $taxIdxLen          = 8;           # No. of characters used to index taxTree for matching
my $continueNoHostFound= 0;           # Continue program if no host match in TaxFile [default: STOP]
my $compareSameTaxLevel = $COMPARESAMETAXLEVEL;
my $originalCompareSameTaxLevel = 0;
my $doCompleteTax       = 0;          # user option; perform SPECIES intersection, then compute uniques for SPECIES -> PHYLUM
my %taxAbbr :shared;                  # maps full rank names to their abbreviations
tie %taxAbbr, "Tie::IxHash";
   %taxAbbr = (
#    species      => "S",
    genus        => "G",
    family       => "F",
    order        => "O",
    class        => "C",
    phylum       => "P",
    superkingdom => "SK",
);
my @ranks :shared;
   @ranks = keys %taxAbbr;

my %taxAbbrExt  :shared;
tie %taxAbbrExt, "Tie::IxHash";
    %taxAbbrExt = (
    subspecies   => "SS",
    species      => "S",
    genus        => "G",
    family       => "F",
    order        => "O",
    class        => "C",
    phylum       => "P",
    superkingdom => "SK",
);
my @ranksExt     :shared = keys %taxAbbrExt;
my @ranksExtAbbr :shared = values %taxAbbrExt;
my %revTaxAbbrExt:shared;
   @revTaxAbbrExt{@ranksExtAbbr} = @ranksExt;       # { SK => superkingdom, ... }

my $abbrTaxLevel:shared;

GetOptions(
    "taxIdxFile=s"              => \$taxIdxFile,
    "loadFasta=s"               => \$loadFasta,
    "noStoreFasta"              => \$noStoreFasta,
    "nThreads=i"                => \$nThreads,
    "wordSize:s"                => \$wordSize,         # mult. values comma-separated
    "minWordSize:i"             => \$minWordSize,     
    "maxWordSize:i"             => \$maxWordSize,
    "stepWordSize:i"            => \$stepWordSize,
    "loadWords=s"               => \$loadWords,
    "cleanupWords"              => \$cleanupWords,
    "compareSameTaxLevel"       => \$compareSameTaxLevel,
    "taxLevel=s"                => \$taxLevel,
    "taxTreeFile=s"             => \$taxTreeFile,
    "doCompleteTax"             => \$doCompleteTax,
    "taxIdxLength=i"            => \$taxIdxLen,             # no. of characters to index %taxTree for matching
    "continueNoHostFound"       => \$continueNoHostFound,   # proceed if genome cannot be host-matched
#    "skipNonStdNucs"            => \$skipNonStdNucs,
#    "minUniqueFragLength=i"     => \$minUniqueFragLength,
    "minUniqueFragLength=s"     => \$minUniqueFragLenRange, # OUTPUT...
    "returnErrorFreeSeqs"       => \$wantErrorFreeSeqs,
    "dbLabel=s"                 => \$dbLabel,
    "pdir=s"                    => \$pDir,
    "outdir=s"                  => \$outdir,
    "noStats"                   => \$noStats,
#    "optimize=s"                => \$optimize,
#    "minStringPackLength=i"     => \$minStringPackLen,
    "prefixLength=i"            => \$prefixLength,
#    "partitions=i"              => \$numPartitions,
    "skipGC"                    => \$skipGC,
    "verbose"                   => \$verbose,
    "vverbose"                  => \$vverbose,
    "help"                      => \$help,
    "debug"                     => \$debug,
    "exportWords"               => \$exportWords,
);

usage() if($help);

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
my $iter_time0 = new Benchmark;
#=-=-=-=-=-=-=-=-=-=-= START global script BENCHMARK =-=-=-=-=-=-=-=-=-=-=-=-=-#


#my $file = "/home/traceyf/uniques/benchmark24_packed/update/storable/rContigNickLookup.dmp";
#my %infilehash = %{ retrieve $file };
#print Dump(\%infilehash)."\n"; exit;

$verbose = 1 if($vverbose);

# Verify & Dipslay Input Options
processOptions();
my $iter_time0X = new Benchmark;    # Timer starting after removal of partition dir

# Completed Genomes: 1 FASTA filename = 1 replicon
# Draft Genomes:           one genome = 1 replicon
my $totalReplicons = 0;     

$originalCompareSameTaxLevel = $compareSameTaxLevel;    # Save for later

# Parses input from XML file
# Populates: %organisms         : main input data hash
#            @contigIDs         : Non-redundant list of all contig names/IDs
#            $totalContigs
#            $totalReplicons    : FASTA filename containing contig names/IDs
#            %orgLookupByContig : Hash lookup: use contigID to get OrgName
#            @orgNames          : Non-redundant list of all unique organisms
#            $numOrgs
populateGenomeSeqs($taxIdxFile) unless ($loadFasta);
displayWarnings();
#retrieveFasta() if($loadFasta);

stat_log("Loaded $totalContigs contigs [$totalReplicons replicons, $numOrgs organisms]");

# Attach full taxonomy lineage to organisms (Genus, Family, ... Superkingdom)
processTaxonomy();

#storeFasta() unless ($loadFasta || $noStoreFasta);
displayWarnings();  # called by sub storeFasta() and sub indexContigLmers()

#------------ NEW main() ------------------

# Prepare taxonomy analysis loop
my @ranks2analyze: shared = determineRanks2Analyze($doCompleteTax, $taxLevel);

# Populates globals: %fileNickLookup    # KEY:full repliconName; VALUE: replicon nickname
#                    %rFileNickLookup   # KEY:replicon nickname; KEY: full repliconName
#                    %contigNickLookup
#                    %rContigNickLookup
createRepliconNicknames();

# Populates globals: %partitionHostCoords       # KEY: host ID (int); VALUE: HREF to startCoord => endCoord
#                    @partitionHostStartCoords  # startCoords; array position is its host ID
#                    @partitionHostEndCoords    # endCoords;   array position is its host ID
#                    $char2posLookup            # KEY: ATGC(N) or stringPacked() chars; VAL: counter
#                    $pos2charLookup            # Inverse of $char2posLookup
#                    @charKeys
#                    $maxCharKeys               
#registerPartitions($char2codon, $minStringPackLen, $k, $prefixLength, $numPartitions, $skipNonStdNucs);

# Global variable "%lmerCoords" is populated in one of these subroutines.
# These subroutines create l-mer hashes with start position references for
# each replicon and stores them to disk under "$outdir/tmp".
if($loadWords) {
    # Check that all %wordsWithStarts binary hashes are on disk for $k
    # ****DISABLED****
    verifyWordsOnDisk($k);
    # ****DISABLED****
}
else {
    # Populates globals: @processedContigs (uses replicon nicknames, tho)
    #                    @processedContigs reset in here before (re)population
    #                    %contigID2filename
    #                    %lmerCoords
    extractKmersParallel();
    #prepForBowtie();
    
    # If a complete taxonomic analysis is desired, do not compare same tax levels
    # for ranks above SPECIES as this will just be repeating analyses
    if($doCompleteTax) {
        # Generate %lmerCoords at SPECIES level with $compareSameTaxLevel = 1,
        # BUT... only subtract based on $compareSameTaxLevel = 0;
        $compareSameTaxLevel = ($taxLevel =~ m/SPECIES/i) ? (1) : (0);
    }
    
    # The $taxLevel determines the organisms to be intersected. If this is "SPECIES"
    # then the intersections in %lmerCoords can be re-used at all higher tax levels.
    #lmerIntersectParallel();
    intersectKmersParallel(); #$kmerIntersectionOpts);

    #exportForBowtieComparison();
    #exit;

}

TAXLEVEL: foreach my $currTaxLevel (@ranks2analyze) {
    $taxLevel = $currTaxLevel;
    createDir($outdir."/".$fastadir."/".$taxLevel);

    # If $doCompleteTax and starting at SPECIES, then the finest granularity of intersections
    # were done (i.e. $compareSameTaxLevel = 1). Subtractions for all tax levels above SPECIES 
    # should be made with $compareSameTaxLevel = 0. The SPECIES level should have both done
    # if user originally supplied $compareSameTaxLevel = 1.
    #
    # If !$doCompleteTax, then leave $compareSameTaxLevel as is.
    $compareSameTaxLevel = 0 if($doCompleteTax);                # Needed for all non-SPECIES levels 
    $compareSameTaxLevel = 1 if($originalCompareSameTaxLevel);  # SPECIES-level not done yet if still == 1

    my $taxLevelText = $compareSameTaxLevel ? ("within") : ("among");
    stat_log("");
    stat_log("**********************************");
    stat_log("    TAX LEVEL: $taxLevel ($taxLevelText)");
    stat_log("**********************************");
    stat_log("");
    
    FRAGLENGTH: foreach my $currMinUniqueFragLen (@uniqueFragLenArray) {
        $minUniqueFragLength = $currMinUniqueFragLen; # override user's default
        %uniqueFragsByOrg    = ();                    # reset
        @processedOrgs       = ();                    # reset

        # Archive %lmerCoords here so they can be removed in sub subtractKmersParallel()
        # to save RAM; only needs to be done once if at species level.
        # UPDATE: 2012-10-09
        #archiveKmerIntersections() if(($taxLevel =~ m/SPECIES/i) && $originalCompareSameTaxLevel);
        
        # Global shared variable "%uniqueFragsByOrg" is populated in these subs
        # Organisms to be compared are determined by global $taxLevel
        # ***UPDATE***: Need to add filtering mechanism for appropriate $taxLevel-based
        #               subtraction, using %lmerCoords generated from the species-level intersections.
        #               Currently, $compareSameTaxLevel is set to 0 for all Tax Levels upon subtraction.
        #
        # This sub performs the same as storeUniqueFrags()
        subtractKmersParallel();

#        @processedOrgs       = keys %uniqueFragsByOrg;
#        @fastaFilesOfUniques = ();                   # reset

        # Write unique Frags to disk
        #storeUniqueFrags();     # uses @processedOrgs
        #storeUniqueLmers();     # uses @processedOrgs

        # DATE: 2012-02-19
        # Write stats of unique Frags to Disk
        storeUniqueFragStats() unless ($noStats);
        #exportUniqueCoords();

        archiveRunData() if(($taxLevel =~ m/SPECIES/i) && $originalCompareSameTaxLevel);

        resetOrgStats();
        
    } #FRAGLENGTH

    $originalCompareSameTaxLevel = 0;                           # Remove for next round

} #TAXLEVEL

#=-=-=-=-=-=-=-=-=-=-=-= END global script BENCHMARK =-=-=-=-=-=-=-=-=-=-=-=-=-#
STDOUT->autoflush(1);
my $iter_time1X = new Benchmark;
my $iter_time1  = new Benchmark;
stat_log("----------------------------------");
stat_log("See \"$logfilename\" for details of the run.");
stat_log("----------------------------------");
my $timeText0 = timestr(timediff($iter_time1, $iter_time0));
$timeText0 = $1 if($timeText0 =~ m/^\s*(\d+)\s+/);
$timeText0 = $1 if($timeText0 =~ m/(\d+\.?\d+)\s+/);
my $timeText1 = timestr(timediff($iter_time1X,$iter_time0X));
$timeText1 = $1 if($timeText1 =~ m/^\s*(\d+)\s+/);
$timeText1 = $1 if($timeText1 =~ m/(\d+\.?\d+)\s+/);
stat_log("Main output files");
stat_log("    LOG File                    : $outdir/");
stat_log("    STATISTICS Files            : $outdir/");
stat_log("    FASTA (Uniques/Shared) Files: $outdir/$fastadir/");
stat_log("    INCREMENTAL UPDATE Files    : $outdir/$updateDir/$storableDir/");
stat_log("Additional output files");
stat_log("    PARTITION DIR(s)            : [".join(",", @valid_pDirs)."]");
stat_log("----------------------------------");
stat_log("Total Script time:                               [$timeText0] wallsecs");
stat_log("Total Script time (after partition dir removal): [$timeText1] wallsecs");
stat_log("===================================================================");
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# ********** END MAIN() **********
################################################################################
# To do: 1. Update "die()" with "stat_log("...")"
#        2. format text within boundaries
#-------------------------------------------------------------------------------
sub checkWordSizes {

    my ($fixedWordSize, $minWordSize, $maxWordSize) = @_;

    if (!$fixedWordSize) {                                   
        if (!$minWordSize) {                           
            if (!$maxWordSize) {
                # Neither MIN nor MAX set; Set WORD_SIZE to default
                $fixedWordSize = $DEFAULT_WORDSIZE;
            }
            
            # MAX was set; 
            $minWordSize = ($DEFAULT_MINWORDSIZE < $maxWordSize) ? ($DEFAULT_MINWORDSIZE)    
                         :                                         ($maxWordSize);
    
            # MAX is less than MIN, so set it to MIN (equal to setting WORD_SIZE alone)
            $maxWordSize = ($maxWordSize < $minWordSize) ? ($minWordSize)    
                          :                                   ($maxWordSize);
                          
            if (($minWordSize > $CRAZY_WORD_SIZE) || ($maxWordSize > $CRAZY_WORD_SIZE)) {
                stat_log("YOU\'RE CRAZY!! Setting the word size greater than "
                        ."$CRAZY_WORD_SIZE will break your computer! Aborting...");
                die "\n";
            }
            if (($minWordSize < 1) || ($maxWordSize < 1)) {
                stat_log("Setting the word size less than 1 has no meaning.");
                die "\n";
            }
        }
        elsif (!$maxWordSize) {                        # MIN was set and MAX was not...
            $maxWordSize = ($DEFAULT_MAXWORDSIZE >= $minWordSize) 
                         ? ($DEFAULT_MAXWORDSIZE)   # ...set MAX to default if MIN <= DEFAULT_MAX
                         : ($minWordSize);          # ...else set MAX to current MIN
            if (($minWordSize > $CRAZY_WORD_SIZE) || ($maxWordSize > $CRAZY_WORD_SIZE)) {
                stat_log("YOU\'RE CRAZY!! Setting the word size greater than "
                        ."$CRAZY_WORD_SIZE will break your computer! Aborting...");
                die "\n";
            }
            if (($minWordSize < 1) || ($maxWordSize < 1)) {
                stat_log("Setting the word size less than 1 has no meaning.");
                die "\n";
            }
        }
        else {             # Both --min and --max SET; Check that --min <= --max
            if ($minWordSize > $maxWordSize) {
                stat_log("ERROR: \"min_word_size\" is GREATER than \"max_word_size!\" Aborting...");
                die "\n";
            }
            
            if (($minWordSize > $CRAZY_WORD_SIZE) || ($maxWordSize > $CRAZY_WORD_SIZE)) {
                stat_log("YOU\'RE CRAZY!! Setting the word size greater than "
                        ."$CRAZY_WORD_SIZE will break your computer! Aborting...");
                die "\n";
            }
                
            if (($minWordSize < 1) || ($maxWordSize < 1)) {
                stat_log("Setting the word size less than 1 has no meaning.");
                die "\n";
            }
        }
    }
    elsif ($fixedWordSize > $CRAZY_WORD_SIZE) {
        stat_log("YOU\'RE CRAZY!! Setting the word size greater than "
                ."$CRAZY_WORD_SIZE will break your computer! Aborting...");
        die "\n";
    }
    elsif ($fixedWordSize < 1) {
        stat_log("Setting the word size less than 1 has no meaning.");
        die "\n";
    }

    return ($fixedWordSize, $minWordSize, $maxWordSize);
}
################################################################################
#-------------------------------------------------------------------------------
sub verifyOutdir {

    # ----------------------------------------------------------------------------
    # Process the directories that will store the partition data. To minimize
    # the possibility of disk thrashing from many cores accessing the same disk
    # drive at the same time, multiple partition directories should be specified,
    # where each is controlled by a separate controller.
    # ----------------------------------------------------------------------------
    
    my @pDirs = split(/,/, $pDir);
    
    if(@pDirs) {
        # At least one non-empty pDir given
        foreach my $pDir (@pDirs) {
            if(-d $pDir) {
                die "*FATAL*: Cannot write to \"$pDir\"!\n" if(!-w $pDir);
                print "*WARNING*: Directory \"$pDir\" exists -- Data may be overwritten!\n";
            }
            #$pDir = abs_path($pDir);        # Resolve symbolic links
            createDir($pDir);
            push(@valid_pDirs, $pDir);
        }
    }
    else {
        # No pDir given: check if pwd is writeable or die
        my $cwd = abs_path();
        die "*FATAL*: Current working directory \"$cwd\" is not writeable!\n" if(!-w $cwd);
        
        # Add this directory to the list of valid pDirs, @pDirs
        print "*WARNING*: ALL PARTITIONS WILL BE WRITTEN TO THE CURRENT DIRECTORY \"$cwd\"!\n";
        push(@valid_pDirs, $cwd);
    }
    $num_pDirs = scalar(@valid_pDirs);
    $pDirCodeLookup{$_} = $valid_pDirs[$_] for (0..$#valid_pDirs);

    # ----------------------------------------------------------------------------
    # Process the $outdir
    # ----------------------------------------------------------------------------
    # Check if $outdir exists and if so, is it writeable?
    if(-d $outdir) {
        if(!-w $outdir) {
            print "Fatal: Unable to write to directory \'$outdir\'!\n";
            die "\n";
        }
    }
    # $outdir does not exist; possible to create?
    else {
        my $outdirText = ($outdir eq ".") ? ("Current working directory") : ("$outdir");
        print "$outdirText doesn't exist. Attempting to create...";
        createDir($outdir."/tmp");
        print"done\n";
  
    } #ELSE

    $seqDir = $outdir."/seqs";
    
    return;
}
################################################################################
#-------------------------------------------------------------------------------
sub verifyInputOpts {

    #------------------------------ LOGFILE ------------------------------------
    # Determine word sizes
    @wordSizes = split(/,/,$wordSize);
    if(scalar(@wordSizes) == 1) {
        $fixedWordSize  = $wordSizes[0];
        $wordSizeDetail = "l$fixedWordSize";
    }
    else {
        $wordSizeDetail = "lRANGE";
    }

    # Determine minimum unique fragment length sizes 
    @uniqueFragLenArray = split(/,/,$minUniqueFragLenRange);
    if(scalar(@uniqueFragLenArray) == 1) {
        $minUniqueFragLength = $uniqueFragLenArray[0];
        $uniqueFragLenDetail = "u$minUniqueFragLength";
    }
    else {
        $uniqueFragLenDetail = "uRANGE";
    }

    # partition, compression, and error info
    #my $partitionDetail   = "n$nThreads.p$prefixLength.m$numPartitions";
    my $partitionDetail   = "n$nThreads.p$prefixLength";
    my $compressionDetail = "s$minStringPackLen";
    #my $errorDetail       = ($skipNonStdNucs) ? ("N0") : ("N1");
    my $returnErrors      = ($wantErrorFreeSeqs) ? ("e0") : ("e1");

    # Logfile Initialization
    $logDetail    = $compressionDetail."."
                   .$wordSizeDetail."."
    #               .$errorDetail."."
                   .$uniqueFragLenDetail."."
                   .$returnErrors."."
                   .$partitionDetail;

    # Notify if compareSameTaxLevel used to create DB. 
    #   taxLevel0 = not compared at same level (e.g. species0)
    #   taxLevel1 = compared at same level     (e.g. order1)
    my $vs  = $compareSameTaxLevel ? "1" : "0";
    $prefix = "logfile.$dbLabel.$taxLevel$vs.$logDetail" if($dbLabel);

    $logfilename  = $logfilename ?  $logfilename
                  : $prefix      ?  $prefix.".log"
                  :                ($prefix = "logfile").".$logDetail.log";

    verifyOutdir();

    $logfilename  = "$outdir/$logfilename";
    clear_existing_output_file($logfilename);       # Delete any previous LOG file
                                                    # Start logging initial activity...
                                                   
    $LOGFILE = &init_log($logfilename);             # Initialize log file for writing
    @reporting_filehandles = (*STDOUT, $LOGFILE);   # Filehandles to report program's status

    # ----------
    my $lmerRangeText = join(":",@wordSizes);
    my $uniquesRangeText = join(":",@uniqueFragLenArray);
    stat_log("Using l-mer range: $lmerRangeText");
    stat_log("Using unique fragment range: $uniquesRangeText");

    # Threading with >1 CPUs?
    if($nThreads < $NTHREADS) {
        stat_log("Error: Invalid value for --nThreads [$nThreads]!");
        die "\n";
    }

    #------------------------------ GENOMES ------------------------------------

    if(!$taxIdxFile && !$loadFasta) {
        stat_log("Error: Please specify a taxIdx XML file indicating where all "
                 ."the FASTA files are located with the --taxIdxFile=<FILE> option.");
        #stat_log("This directory should be further subdivided by subdirectories "
        #         ."(named after the host organisms) where");
        #stat_log("the actual FASTA files are located.");
        die "\n";
    }

    if(!-e $taxIdxFile && !$loadFasta) {
        stat_log("Error: TaxIdx file \'$taxIdxFile\' does not exist!");
        die "\n";   
    }

    #------------------------------ WORDS --------------------------------------

    if($loadWords) {
        if(!-d $loadWords) {        # does directory exist?
            stat_log("Error: Directory \'$loadWords\' of stored words/l-mers does not exist.\n");
            die "\n";
        }
    }

    
    if($minUniqueFragLength < 1) {
        stat_log("Invalid \'$minUniqueFragLength\' value. Setting to minimum value of "
             ."$MINUNIQUEFRAGLENGTH-bp");
        $minUniqueFragLength = $MINUNIQUEFRAGLENGTH;
    }
    
    # Check valid l-mer values
    # If --word_size is set, then --min_word_size & --max_word_size are ignored.
    # If --min_word_size ONLY, then --max_word_size = DEFAULT or --min_word_size, whichever is greater.
    # If --max_word_size ONLY, then --min_word_size = DEFAULT or --max_word_size, whichever is less.
    # if --word_size is not set, but MIN & MAX set, *FATAL* if MIN > MAX

    ## Error checking: $fixedWordSize
    $fixedWordSize = $DEFAULT_WORDSIZE if (!$fixedWordSize && !$minWordSize && !$maxWordSize);
    ($fixedWordSize, $minWordSize, $maxWordSize) = checkWordSizes($fixedWordSize, $minWordSize, $maxWordSize);

    # FLAG for using range ($minWordSize and $maxWordSize) or using $fixedWordSize
    #my $using_range = ($fixedWordSize == 0) ? (1) : (0);
    $using_range = (($fixedWordSize == 0) && ($minWordSize != $maxWordSize)) ? (1) : (0);

    #--------------------------- OPTIMIZATION ----------------------------------

    # Specify optimization for l-mer intersections
#    if($optimize) {
#        $optimization->{cpu} = 1 if($optimize =~ m/^cpu/i);
#        $optimization->{mem} = 1 if($optimize =~ m/^mem/i);
#    }
#    else {
#        $optimization->{mem} = 1;          # ->{mem} is broken for some reason...
#        $optimization->{cpu} = 1;
#    }

    #--------------------------- TAXONOMY --------------------------------------

    # Loading a pre-computed taxonomy with $taxTreeFile takes prescedence over 
    # generating one from scratch with $namesTaxFile and $nodesTaxFile.

    # Are all present? TaxTreeFile takes prescedence...
    if(!$taxTreeFile) {
        stat_log("Error: Please specify species taxonomy file (usually \"speciesTreeGI.dmp\") with ");
        stat_log("       --taxTreeFile=<FILENAME>");
        die "\n"; 
    }
    elsif(!-e $taxTreeFile) {
        stat_log("FATAL: TAX TREE file \"$taxTreeFile\" does not exist!");
        die "\n";
    }
    
    # Taxonomic Options
    my $taxMatch = 0;

    foreach my $allowed (@ranksExt) {
        if($taxLevel =~ m/$allowed/i) {
            if($allowed eq "subspecies") {
                stat_log("Sorry, SUBSPECIES-level comparisons are not currently supported.");
                die "\n";
            }
            stat_log("Wanted Taxonomic Level ...FOUND: $allowed");
            $taxLevel = $allowed;
            $taxMatch = 1;
            last;
        }
    }
    die "Fatal: unrecognized Taxonomic Level: --taxLevel=$taxLevel\n"
#       ."Allowed Levels: {".join(",",@allowedTaxLevels)."}\n" if(!$taxMatch);
       ."Allowed Levels: {".join(",",@ranksExt)."}\n" if(!$taxMatch);

    # Add to processOptions() after TAXONOMY has been verified
    $doCompleteTax = ($taxLevel =~ m/SPECIES/i)     # only valid if SPECIES initially given;
                   ? ($doCompleteTax)               # return whatever user wanted for this var;
                   : (0);                           # not generating species-level intersections, so unable to create uniques
                                                    #     for all higher tax levels
    
    #--------------------------- COMPRESSION -----------------------------------

    # Any positive value is valid.
    # Note: stringPack() is disabled by default if($minStringPackLen > $wordSize).
    $k                = $wordSize;  # Assuming fixed wordSize (not array of values)
    $bvLen            = $k << 1;
    $minStringPackLen = ($minStringPackLen > 1)   ? ($minStringPackLen) : $MINSTRINGPACKLEN;
    $stringPacking    = ($k >= $minStringPackLen) ? ("ON")              : ("OFF");

    #---------------------------------------------------------------------------

    # Remove pre-existing partition directory, update dir, and FASTA dir
    # --> *must* remove FASTA dir because existing append feature will append to
    #     pre-existing FASTA files!
    removeDir("$_/$partitionDir") foreach (@valid_pDirs);
    createDir("$_/$partitionDir") foreach (@valid_pDirs);
    removeDir("$outdir/$updateDir");
    removeDir("$seqDir");
    createDir("$seqDir");
    if(-e $outdir."/".$fastadir) {
        print "*WARNING*: FASTA directory \"$outdir/$fastadir\" exists -- Exising files will be removed!\n";
    }
    removeDir($outdir."/".$fastadir);
    
    # UPDATE: 2012-10-09    Commented out archiveKmerIntersections() in main().
    print "*WARNING*: k-mer INTERSECTIONS will NOT be consolidated...\n";
    

    # ------- FIX THIS: TO BE INCLUDED IN sub displayInputOptions() ------ XXXXXXXX

    stat_log("Program: $0");
    stat_log("------------------------------------------------------------------");
    stat_log("Initial Parameters");
    stat_log("------------------------------------------------------------------");
    my     $fixedWordSize_text = ($fixedWordSize == 0)     ? ("n/a") : ($fixedWordSize);
    stat_log("    FIXED WORD SIZE = $fixedWordSize_text");
    my $minWordSize_text = ($minWordSize == 0) ? ("n/a") : ($minWordSize);
    stat_log("      MIN WORD SIZE = $minWordSize_text");
    my $maxWordSize_text = ($maxWordSize == 0) ? ("n/a") : ($maxWordSize);
    stat_log("      MAX WORD SIZE = $maxWordSize_text");

    $MAX_BUCKETS = scalar(@nucs)**$prefixLength;
    
}
################################################################################
#-------------------------------------------------------------------------------
sub displayInputOpts {
    my $stats_out_text = ($noStats) ? ("FALSE") : ("TRUE");
    my $cpuText = ($nThreads > 1) ? ("(PARALLEL)") : ("(SERIAL)");
    my $wordSizeText   = join(":",@wordSizes);
    my $uniqueSizeText = join(":",@uniqueFragLenArray);
    #my $errorBaseText  = ($skipNonStdNucs)      ? ("NO")        : ("YES");
    my $wantedBaseText = ($wantErrorFreeSeqs)   ? ("NO")        : ("YES");
#    my $optText        = $optimization->{mem}   ? "MEM"         : "CPU";
    my $comparingText  = ($compareSameTaxLevel) ? ("comparing") : ("*NOT* comparing");
    my $doCompleteText = ($doCompleteTax)       ? ("YES")       : ("NO");
    my $continueNoHostFoundText = ($continueNoHostFound)     
                                ? ("CONTINUE...") 
                                : ("STOP and report");
    
    stat_log("");
#    stat_log("Starting Parameters:");
    stat_log("    TAXIDXFILE                          = $taxIdxFile");
    stat_log("    OUTDIR                              = $outdir/");
    stat_log("    PARTITION DIRECTORY                 = $pDir/");
    stat_log("    SEQUENCE (BINARY) STORAGE DIRECTORY = $seqDir/");
    stat_log("    DB Label                            = $dbLabel");
    stat_log("    Using pre-computed WORDS            = $loadWords/") if($loadWords);
    stat_log(" TAXONOMY:");
    stat_log("    TAX TREE file                       = $taxTreeFile");
#    stat_log("    Output BINARY TAX TREE file        = $storableTaxFile");
#    stat_log("    Output ASCII  TAX TREE file        = $dumperTaxFile")  if($dumpASCIITaxFile);
    stat_log("    Comparison Tax Level                = ".uc($taxLevel)." ($comparingText within same ".uc($taxLevel).")");
    stat_log("    Creating uniques for ALL tax levels = $doCompleteText");
    stat_log("    No Match to Tax Host Action?        = $continueNoHostFoundText");
    stat_log(" OPTIMIZATION:");
#    stat_log("    Optimized for $optText");
    stat_log("    MAX_BUCKETS                         = $MAX_BUCKETS");
    stat_log("    [n] nThreads                        = $nThreads $cpuText");
    stat_log("    [l] WORD size                       = $wordSizeText");
    #stat_log("    [s] string pack length              = $minStringPackLen [stringPacking() is $stringPacking]");
    stat_log("    [p] prefix length (word indexing)   = $prefixLength");
    #stat_log("    [m] num partitions                  = $numPartitions");
    stat_log(" PROCESSING OPTIONS:");
    #stat_log("    [N] process WORDS w/error bases     = $errorBaseText");
    #stat_log("    [e] return fragments w/error bases  = $wantedBaseText");
    stat_log("    [u] return unique fragment lengths  = $uniqueSizeText");
    stat_log("    CALCULATE STATISTICS                = $stats_out_text");
    stat_log("    LOGFILE                             = $logfilename");
    stat_log("------------------------------------------------------------------");
}
################################################################################
#-------------------------------------------------------------------------------
sub processOptions {
    verifyInputOpts();      # does NOT verify  @wordSizes yet
    displayInputOpts();     # does NOT display @wordSizes yet
}
################################################################################
# In $orgsByFastaDir:
# ARGS: $taxIdxFile
#-------------------------------------------------------------------------------
sub populateGenomeSeqs {
    
    # Parse input XML file and populate %organisms
    loadOrganisms($_[0]);
    $totalContigs = scalar(@contigIDs);
    $numOrgs      = scalar(@orgNames);

	my $iter_time0 = new Benchmark;
    stat_log("Reading input data...");

    # Parse FASTA files and save them as Storable-formatted sequences (no header).
    # This will reduce memory usage and the sequences should be loaded quickly.
    my %contigIDlookup = ();
       @contigIDlookup{@contigIDs} = ();

    # (a) Open contig/replicon FASTA file
    # (b) Parse contig/replicon FASTA file
    # (c) Ensure $contigID matches parsed FASTA header 
    # (d) Write out binary version to tmp directory
    # (e) Store length, destination filename (ensure uniqueness)
    my $contigIdx = 0;
    ORG: foreach my $currOrg (keys %organisms) {
        FILE: foreach my $filename (keys %{ $organisms{$currOrg}->{REPL} }) {
            CONTIG: foreach my $contigID (keys %{ $organisms{$currOrg}->{REPL}->{$filename}->{CONTIG} }) {
                $contigIdx++;
                # Check contigID exists
                if(exists $contigIDlookup{$contigID}) {

            		my $iter_time0 = new Benchmark;
                    stat_log_("  ...[$contigIdx of $totalContigs] Parsing contig FASTA file \'$filename\' [$currOrg]...");
                    open my $FASTAFILE, '<', $filename || die "Fatal: populateGenomeSeqs(). Cannot open FASTA file "
                                                             ."\'$filename\' for read.\n";
                    # sub _parseFASTA() returns two references
                    my $root = $organisms{$currOrg}->{DIR}->{root};
                    my $dir = $organisms{$currOrg}->{DIR}->{dir};

                    #print "\n        Looking for ROOT=$root\tDIR=$dir\tFILENAME=$filename...\n";
                    
                    my ($seqRef, $contigHeaderRef) = _parseFASTA($root."/".$dir."/".$filename);

                    # Check contigID matches filename
                    if(${ $contigHeaderRef } ne $contigID) {
                        die "  **FATAL**: Mismatched FASTA file [$contigID] with contig header [".${ $contigHeaderRef }."]\n"
                           ."             Perhaps you neglected to include the full path to the data directories when generating your XML file?\n"
                           ."             Check the <root>:<name> field(s) in \"$taxIdxFile\" for this common error and ensure all sequence files\n"
                           ."             in their respective directories.\n";
                    }
                    
                    ## Destination directory
                    #my $destDir = "$outdir/tmp";
                    #if(!-d $destDir) {
                    #    stat_log_("$destDir doesn't exist. Attempting to create...");
                    #    my $mkdirCmd= `mkdir $destDir 2>&1`;
                    #    chomp $mkdirCmd;
                    #
                    #    if($mkdirCmd ne "") {
                    #        _stat_log("Fatal: Unable to create directory \"$destDir\"!");
                    #        die "\n";
                    #    }
                    #    _stat_log("done");
                    #} #DestDir
                    
                    # save Storable version of sequence: ENSURE UNIQUENESS OF FILENAME!!!
                    #my $seqStorableFilename = "$outdir/tmp/".basename($filename);
                    my $seqStorableFilename = $seqDir."/".basename($filename);              # ***POTENTIAL FOR CLOBBERING IF BASENAMES ARE THE SAME***
                    store clone($seqRef), $seqStorableFilename;

                    # record filename of Storable-formatted sequence (no header) and its length into %organisms
                    $organisms{$currOrg}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{SEQ}    = $seqStorableFilename;
                    $organisms{$currOrg}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{LENGTH} = length(${ $seqRef });

                    close $FASTAFILE;

                    #_stat_log("done. [".timestr(timediff($iter_time1, $iter_time0))."]");
                    
                    my $iter_time1 = new Benchmark;
                    my $statText = timestr(timediff($iter_time1,$iter_time0));
                    $statText = $1 if($statText =~ m/^\s*(\d+)\s+/);
                    $statText = $1 if($statText =~ m/(\d+\.?\d+)\s+/);
                    _stat_log("done in [$statText] wallsecs");

                }
                else {
                    die "Contig ID \"$contigID\" does not exist!!\n"
                }
            } #CONTIG
        } #FILE
    } #ORG

    # Sort @contigIDs in order of decreasing contig size
    my @sortedIdx = rnkeysort {   $organisms{ $orgLookupByContig{$contigIDs[$_]} }
                                  ->{REPL}->{ $contigID2filename{$contigIDs[$_]} }
                                ->{CONTIG}->{ $contigIDs[$_] }
                                ->{LENGTH}
                              } (0..$#contigIDs);
    @contigIDs = @contigIDs[@sortedIdx];
    @contigFilenames = @contigFilenames[@sortedIdx];
    
    #print "    There are ".(scalar(@contigIDs))." contigIDs\n";
    #print "    There are ".(scalar(@contigFilenames))." contig Filenames\n";
    #
    #print "    Contig\tFilename\n";
    #print "    ".$contigIDs[$_]."\t".$contigFilenames[$_]."\n" for (0..$#contigIDs);
    
    my $iter_timeY = new Benchmark;
    my $statText = timestr(timediff($iter_timeY,$iter_time0));
    $statText = $1 if($statText =~ m/(\d+\.?\d+)\s+/);
    $statText = $1 if($statText =~ m/^\s*(\d+)\s+/);
    stat_log("DONE parsing $totalContigs file(s) in [$statText] wallsecs.");
    stat_log("----------------------------------");
    return;
}
################################################################################
#-------------------------------------------------------------------------------
sub _parseFASTA {
    my $filename = shift;

    my $seq                  = q{}; # FASTA sequence
    my $geneID_already_found = 0;   # FLAG; TRUE when the first geneID is found
    my $header               = q{}; # FASTA header, from ">" up to first whitespace

    # Open genomic FASTA file for reading...
    open my $FASTAFILE, '<', $filename
        || die "Fatal: _parseFASTA(). Cannot open FASTA file for read: $!\n";
    
    STDOUT->autoflush(1);
    while (my $line = <$FASTAFILE>) {           # Read in one line...
        chomp $line;                               # Remove end-of-line marker
        if ($line =~ m/^[^(\s||\n)]*/) {           # Line not empty

            if ($line =~ /^>(\S+)\s+/) {              # Line is a geneID
                $geneID_already_found = 1;     # ... remember geneID found
                $header = $1;
            }
            elsif ($geneID_already_found) {        # is NOT a geneID, but one
                $line =~ m/^(\s*)(.*)$/;
                my $tentativeLine = $2;

                # Change all non-ATGC bases to N
                $tentativeLine =~ tr/ATGCN/N/c;
                $seq .= $tentativeLine;    # ...append to hash w/wspace
            }
        }
    }
    STDOUT->autoflush(0);

    close $FASTAFILE;

    return (\$seq, \$header);
}
################################################################################
# Populates global/shared variable:
#   %taxTree = ( $taxid1 => { S  => $species,
#                            SS => $subspecies,
#                            G  => $genus,
#                            F  => $family,
#                            O  => $order,
#                            C  => $class,
#                            P  => $phylum,
#                            SK => $superkingdom },
#                $taxid2 => { ... },
#                   ...
#              );
#
#   UPDATE: loads entire tax tree temporarily into memory, then matches organisms
#           to items in the tax tree, appending full taxonomic information to
#           %organisms. %taxTree is then dumped from memory.
#
#-------------------------------------------------------------------------------
sub processTaxonomy {

    #####################################################################
    # 1. Parsed *.xml file will contain the organisms and their TAX NAME, 
    #    TAX RANK, TAXID, and GI.   **A small number of orgs do NOT have a TAXID**
    # 2. Load "speciesTreeGI.dmp" to associate the {S} or {SS} name with 
    #    the rest of the tax tree [keep track of @unmatched]
    # 3. append the remaining tax info {G},{F},{O},{C},{P},{SK} to %organisms
    #    --see matchOrgs2TaxTree() for method of adding to %organisms
    #####################################################################
    # Load %taxTree = speciesTreeGI.dmp
    # Search %taxTree for each organisms' TAXID, which will be the "node" in the %taxTree
    # If match, add tax info to %organisms
    # If no match, add org to @unmatchedByTaxID
    # if(@unmatchedByTaxID), try old-school method to assign taxonomy
    # if(@unmatchedByName),  either ABORT, DELETE offending organisms, or create synthetic taxonomy



    # PROBLEM: TAXID in "speciesTreeGI.dmp" is the TAXID for the species, which will
    #          NOT correspond to those with strain-level information. A more complete
    #          process would store all strain-level TAXIDs with the species-level TAXID,
    #          but this will take more time to figure out, so right now, we just use
    #          a more brute force method.
    # Method #1:
    #          (1) create a temp hash of all the {S} and {SS} names as KEYS and their
    #              species TAXID as the VALUE;
    #          (2) search this temp tree with the orgName and, if found, use the VALUE
    #              (species-level TAXID) to get the final destination in "speciesTreeGI.dmp"
    # Method #2:
    #          (1) create a temp hash using the GIs in speciesTreeGI.dmp as KEYS and the
    #              species TAXID as the VALUE;
    #          (2) search this hash with one of the contig's GIs of the associated

    STDOUT->autoflush(1);

    my $fileDesc = "Species Taxonomy Tree (w/GIs)";
    stat_log("----------------------------------");
    my %taxTree = %{ loadFromDisk($taxTreeFile, $fileDesc) };
    
    ################################
    # Create temporary hash
    ################################
    stat_log_("Restructuring Taxonomic Tree...");
    my $timeX = new Benchmark;
    my %nameTaxidLookup = ();
       $nameTaxidLookup{SS} = ();
       $nameTaxidLookup{S}  = ();
    foreach my $node (keys %taxTree) {
        if(exists $taxTree{$node}->{SS}) {
            my @names = keys %{ $taxTree{$node}->{SS} };
            $nameTaxidLookup{SS}->{$_} = $node foreach(@names);
        }
        if(exists $taxTree{$node}->{S}) {
            my @names = keys %{ $taxTree{$node}->{S} };
            $nameTaxidLookup{S}->{$_} = $node foreach(@names);
        }
    }
    my $timeY = new Benchmark;
    my $restructText = timestr(timediff($timeY, $timeX));
    $restructText = $1 if($restructText =~ m/^\s*(\d+)\s+/);
    $restructText = $1 if($restructText =~ m/(\d+\.?\d+)\s+/);
    #_stat_log("done. [".timestr(timediff($timeY,$timeX))."]");
    _stat_log("done. [$restructText] wallsecs");
    
    ################################
    # Search temp hash with orgName
    ################################
    my $time1 = new Benchmark;
    stat_log_("Matching input organisms to $fileDesc ...");
    #my @unmatchedByTaxID = ();
    my @unmatchedByName  = ();
    
    ORG: foreach my $currOrg (keys %organisms) {
        #my $taxid = $organisms{$currOrg}->{TAXID};
        my $currOrgRank = $organisms{$currOrg}->{RANK};
        #print "RANK for $currOrg = $currOrgRank\n";

        #if(exists $taxTree{$taxid}) {
        if(exists $nameTaxidLookup{$currOrgRank}->{$currOrg}) {
            my $speciesTaxid = $nameTaxidLookup{$currOrgRank}->{$currOrg};
        
            # Add to $organisms{$currOrg}->{TAXTREE}
            #my $currOrgRank = $organisms{$currOrg}->{RANK};
            
            # Create reference lookup of valid taxonomies
            my %rankAbbrs = ();
               @rankAbbrs{@ranksExtAbbr} = ();  # { SS => (), S => (), ... };
            
            if($currOrgRank eq "SS") {
                $organisms{$currOrg}->{TAXTREE}->{SS} = $currOrg;
                delete $rankAbbrs{SS};
                # get sciname from "S"
                my $sciName = q{};
                my @names = keys %{ $taxTree{$speciesTaxid}->{S} };
                NAME: foreach my $name (@names) {
                    my $class = $taxTree{$speciesTaxid}->{S}->{$name};
                    if($class =~ m/scientific/i) {
                        $sciName = $name;
                        last NAME;
                    }
                } #NAME
                # Grab first name if no scientific name available
                $sciName = $names[0] if(!$sciName);
                $organisms{$currOrg}->{TAXTREE}->{S} = $sciName;
                delete $rankAbbrs{S};
            }
            elsif($currOrgRank eq "S") {
                $organisms{$currOrg}->{TAXTREE}->{S} = $currOrg;
                delete $rankAbbrs{S};
                delete $rankAbbrs{SS};  # No strain-level info, so just delete
            }
            else {
                # ... Organism name is neither a STRAIN/SUBSPECIES or SPECIES ...
                # ... This shouldn't happen, so we notify and abort ...
                stat_log("\n**FATAL**: Unsupported input taxonomic rank [$currOrgRank, $currOrg]");
                stat_log("\n           Only species \"S\" or subspecies/strain \"SS\" are supported");
                die      "           Abort.\n";
            }

            # Only non-SS and non-S ranks remain in %rankAbbrs and these are all SciNames: process them
            $organisms{$currOrg}->{TAXTREE}->{$_} = $taxTree{$speciesTaxid}->{$_} foreach (keys %rankAbbrs);

        } # exists
        else {
            push(@unmatchedByName, $currOrg);
        }

    } #ORG

    my $time2 = new Benchmark;
    my $timeText = timestr(timediff($time2, $time1));
    $timeText = $1 if($timeText =~ m/^\s*(\d+)\s+/);
    $timeText = $1 if($timeText =~ m/(\d+\.?\d+)\s+/);
    _stat_log("done. [$timeText] wallsecs");

    # Some organisms couldn't be TAXONOMY-matched with their SPECIES/STRAIN scientific name
    my $abortNoNameMatch   = 1;
    my $deleteNoNameMatch  = 0;
    my $createSyntheticTax = 0;
    if(@unmatchedByName) {

        stat_log("Unfortunately, the following ".(scalar(@unmatchedByName))." organism(s) could not be placed in the taxonomic tree by scientific name and will be removed:");
        stat_log("==============================================================");
        stat_log("    $_") foreach (@unmatchedByName);

        #if($abortNoNameMatch) {
        #    # ...
        #    die "not yet implemented abortNoNameMatch()...\n";
        #    # ...
        #} #abortNoNameMatch
        #elsif($deleteNoNameMatch) {

             wipeOrg(\@unmatchedByName);
            
        #} #deleteNoNameMatch
        #elsif($createSyntheticTax) {
        #    # ...
        #    die "not yet implemented createSyntheticTax()...\n";
        #    # ...
        #} #createSyntheticTax

    } #unmatchedByName
    else {
        stat_log("CONGRATULATIONS! All $numOrgs organisms were assigned "
                ."taxonomic identities at the species level or deeper.");
    }

#        stat_log("STOPPING.") unless ($continueNoHostFound);
#        die "\n"              unless ($continueNoHostFound);
#        stat_log("CONTINUING...");

    return;
}
################################################################################
# IF THERE IS ONLY ONE "root" NAME:
# ------------------------------------------------------------------------------
# taxIdx = ( linked => {
#              draft => {
#                  $kingdom => {
#                      root => { name => $name,
#                                org  => { $orgName1 => { acc      => $accession,
#                                                         date     => $date,
#                                                         dir      => $dirname,
#                                                         gi       => $gi,
#                                                         rank     => $rank,   [SS, S, ...]
#                                                         source   => $source, [plasmid, chr]
#                                                         stype    => $sType,  [GBK, GEN]
#                                                         taxid    => $taxid,
#                                                       },
#                                          $orgName2 => { ... },
#                                             ...
#                                        },
#                              },
#                              },
#                       },
#            unlinked => {  ... }
#          );
# ------------------------------------------------------------------------------
# IF THERE ARE MULTIPLE "root" NAMES:
# ------------------------------------------------------------------------------
# taxIdx = ( linked => {
#                draft => {
#                    $kingdom => {
#                        root => { 
#                            $rootDir => { org => { $orgName1 => { acc      => $accession,
#                                                                  date     => $date,
#                                                                  dir      => $dirname,
#                                                                  gi       => $gi,
#                                                                  rank     => $rank,   [SS, S, ...]
#                                                                  source   => $source, [plasmid, chr]
#                                                                  stype    => $sType,  [GBK, GEN]
#                                                                  taxid    => $taxid,
#                                                                },
#                                                   $orgName2 => { ... },
#                                                      ...
#                                        },
#                                 ...
#                                },
#                                },
#                         },
#                      },
#             unlinked => { ... },
#           );
#
# ******************************************************************************
# All organisms entered into the XML file, regardless of STATUS (completed or drafts)
# or KINGDOM, are subject to the process of "uniques" generation at the specified tax level.
# ******************************************************************************
#
# -loads %taxIdx
# -populates globally shared %organisms, but not its sequences
# ARGS: $taxIdxFile
################################################################################
sub loadOrganisms {

    STDOUT->autoflush(1);
    stat_log_("Importing XML file...");    
    my %taxIdx = %{ (XMLin($_[0])) };     # Local; goes out of scope at end of sub
    _stat_log("done.");

    standardizeTaxIdx(\%taxIdx);
    
    # XMLin parses the "root" data from the XML file slightly differently depending 
    # on the number of of root names in the XML file. If multiple root names, then 
    # these names become KEYS under {root}, and all organisms go under this $rootdir:
    #   {root}->{$rootdir}
    #   {root}->{$rootdir}->{org}->{$orgname}
    # If only ONE root, then this one name becomes the VALUE of {name} under {root}, and
    # the organisms go under {root} instead of {$rootdir}: 
    #   {root}->{name}->{$rootdir}
    #   {root}->{org}->{$orgname}

    # So just check for the fields "name" and "org" under the "root". If present, then
    # only 1 root. If not present, then we have multiple roots.

    #print Dump(clone(\%taxIdx))."\n";
    STDOUT->autoflush(1);

    stat_log_("Loading XML data into Perl hash...");
    #
    # All organisms entered into the XML file, regardless of STATUS (completed or drafts) or KINGDOM, 
    # will be subject to the process of "uniques" generation at the tax level desired.
    #
    STATUS: foreach my $status (keys %{ $taxIdx{linked} }) {
        #print "STATUS=$status\n";
        KINGDOM: foreach my $kingdom (keys %{ $taxIdx{linked}->{$status} }) {
            #print "KINGDOM=$kingdom\n";
            ROOTDIR: foreach my $rootDir (keys %{ $taxIdx{linked}->{$status}->{$kingdom}->{root} }) {
                #print "ROOTDIR=$rootDir\n";
                ORG: foreach (keys %{ $taxIdx{linked}->{$status}->{$kingdom}->{root}->{$rootDir}->{org} }) {
                    #print "ORG=$_\n";
                    $organisms{$_}            = &share({}) unless (exists $organisms{$_});
                    $organisms{$_}->{DIR}     = &share({}) unless (exists $organisms{$_}->{DIR});
                    $organisms{$_}->{TAXTREE} = &share({}) unless (exists $organisms{$_}->{TAXTREE});
                    $organisms{$_}->{STATS}   = &share({}) unless (exists $organisms{$_}->{STATS});
                    $organisms{$_}->{REPL}    = &share({}) unless (exists $organisms{$_}->{REPL});

                    $organisms{$_}->{STATS}->{GCMEAN}  = &share({}) unless (exists $organisms{$_}->{STATS}->{GCMEAN});
                    #$organisms{$_}->{STATS}->{GCSD}  = &share({}) unless (exists $organisms{$_}->{STATS}->{GCSD});
                    $organisms{$_}->{STATS}->{GCMIN}  = &share({}) unless (exists $organisms{$_}->{STATS}->{GCMIN});
                    $organisms{$_}->{STATS}->{GCMAX}  = &share({}) unless (exists $organisms{$_}->{STATS}->{GCMAX});
                    $organisms{$_}->{STATS}->{GCCOUNT}  = &share({}) unless (exists $organisms{$_}->{STATS}->{GCCOUNT});
                    $organisms{$_}->{STATS}->{GCTOT}  = &share({}) unless (exists $organisms{$_}->{STATS}->{GCTOT});
                    $organisms{$_}->{STATS}->{LMEAN}  = &share({}) unless (exists $organisms{$_}->{STATS}->{LMEAN});
                    #$organisms{$_}->{STATS}->{LSD}  = &share({}) unless (exists $organisms{$_}->{STATS}->{LSD});
                    $organisms{$_}->{STATS}->{LMIN}  = &share({}) unless (exists $organisms{$_}->{STATS}->{LMIN});
                    $organisms{$_}->{STATS}->{LMAX}  = &share({}) unless (exists $organisms{$_}->{STATS}->{LMAX});
                    $organisms{$_}->{STATS}->{LTOT}  = &share({}) unless (exists $organisms{$_}->{STATS}->{LTOT});
                    $organisms{$_}->{STATS}->{LCOUNT}  = &share({}) unless (exists $organisms{$_}->{STATS}->{LCOUNT});
                    $organisms{$_}->{STATS}->{LREMOVED}  = &share({}) unless (exists $organisms{$_}->{STATS}->{LREMOVED});
                    $organisms{$_}->{STATS}->{LPCTREMOVED}  = &share({}) unless (exists $organisms{$_}->{STATS}->{LPCTREMOVED});

                    $organisms{$_}->{STATS}->{uGCMEAN}  = &share({}) unless (exists $organisms{$_}->{STATS}->{uGCMEAN});
                    #$organisms{$_}->{STATS}->{uGCSD}  = &share({}) unless (exists $organisms{$_}->{STATS}->{uGCSD});
                    $organisms{$_}->{STATS}->{uGCMIN}  = &share({}) unless (exists $organisms{$_}->{STATS}->{uGCMIN});
                    $organisms{$_}->{STATS}->{uGCMAX}  = &share({}) unless (exists $organisms{$_}->{STATS}->{uGCMAX});
                    $organisms{$_}->{STATS}->{uGCCOUNT}  = &share({}) unless (exists $organisms{$_}->{STATS}->{uGCCOUNT});
                    $organisms{$_}->{STATS}->{uGCTOT}  = &share({}) unless (exists $organisms{$_}->{STATS}->{uGCTOT});
                    $organisms{$_}->{STATS}->{uLMEAN}  = &share({}) unless (exists $organisms{$_}->{STATS}->{uLMEAN});
                    #$organisms{$_}->{STATS}->{uLSD}  = &share({}) unless (exists $organisms{$_}->{STATS}->{uLSD});
                    $organisms{$_}->{STATS}->{uLMIN}  = &share({}) unless (exists $organisms{$_}->{STATS}->{uLMIN});
                    $organisms{$_}->{STATS}->{uLMAX}  = &share({}) unless (exists $organisms{$_}->{STATS}->{uLMAX});
                    $organisms{$_}->{STATS}->{uLTOT}  = &share({}) unless (exists $organisms{$_}->{STATS}->{uLTOT});
                    $organisms{$_}->{STATS}->{uLCOUNT}  = &share({}) unless (exists $organisms{$_}->{STATS}->{uLCOUNT});

                    push(@orgNames, $_);        # Non-redundant list of all organism (scientific) names (globally shared)

                    $organisms{$_}->{DIR}->{root} = $rootDir;
                    $organisms{$_}->{DIR}->{dir}  = $taxIdx{linked}->{$status}->{$kingdom}->{root}->{$rootDir}->{org}->{$_}->{dir};
                    $organisms{$_}->{RANK}        = $taxIdx{linked}->{$status}->{$kingdom}->{root}->{$rootDir}->{org}->{$_}->{rank};
                    $organisms{$_}->{TAXID}       = $taxIdx{linked}->{$status}->{$kingdom}->{root}->{$rootDir}->{org}->{$_}->{taxid};

                    my @filenames = keys %{ $taxIdx{linked}->{$status}->{$kingdom}->{root}->{$rootDir}->{org}->{$_}->{file} };        
                    FILENAME: foreach my $filename (@filenames) {
                        #print "FILENAME=$filename\n";
                        my $topology = $taxIdx{linked}->{$status}->{$kingdom}->{root}->{$rootDir}->{org}->{$_}->{file}->{$filename}->{topo};
                        my $replType = $taxIdx{linked}->{$status}->{$kingdom}->{root}->{$rootDir}->{org}->{$_}->{file}->{$filename}->{repltype};
                        my $nType    = $taxIdx{linked}->{$status}->{$kingdom}->{root}->{$rootDir}->{org}->{$_}->{file}->{$filename}->{ntype};
                        my $doRC     = $taxIdx{linked}->{$status}->{$kingdom}->{root}->{$rootDir}->{org}->{$_}->{file}->{$filename}->{rc};
                        $organisms{$_}->{REPL}->{$filename} = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename});
                        $organisms{$_}->{REPL}->{$filename}->{TOPO}     = $topology; # (L)inear or (C)ircular
                        $organisms{$_}->{REPL}->{$filename}->{REPLTYPE} = $replType; # CHR or PLASMID
                        $organisms{$_}->{REPL}->{$filename}->{NTYPE}    = $nType;    # Nucleotide type (DNA, RNA, ...)
                        $organisms{$_}->{REPL}->{$filename}->{RC}       = $doRC;     # perform RC? (0 or 1)
                        $organisms{$_}->{REPL}->{$filename}->{UNIQUE}   = 0;

                        # To add support for multiple contigs per filename (NOT YET IMPLEMENTED), replace "my $contigID =" with 
                        # "my @contigIDs = keys %{}". Currently, this sub assumes only one contig possible.
                        my @contigNames = keys %{ $taxIdx{linked}->{$status}->{$kingdom}->{root}->{$rootDir}->{org}->{$_}->{file}->{$filename}->{contig} };
                        CONTIG: foreach my $contigID (@contigNames) {
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}              = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG});
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID} = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID} );
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{PREFIX} = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{PREFIX} );
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{WORDS}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{WORDS} );
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{DESC}   = $taxIdx{linked}->{$status}->{$kingdom}->{root}->{$rootDir}->{org}->{$_}->{file}->{$filename}->{contig}->{$contigID}->{desc};
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{STATS}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{STATS});

                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{GCMEAN}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{GCMEAN});
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{GCSD}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{GCSD});
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{GCMIN}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{GCMIN});
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{GCMAX}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{GCMAX});
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{GCCOUNT}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{GCCOUNT});
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{GCTOT}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{GCTOT});
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{LMEAN}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{LMEAN});
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{LSD}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{LSD});
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{LMIN}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{LMIN});
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{LMAX}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{LMAX});
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{LTOT}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{LTOT});
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{LCOUNT}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{LCOUNT});
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{LREMOVED}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{LREMOVED});
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{LPCTREMOVED}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{LPCTREMOVED});
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{uGCMEAN}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{uGCMEAN});
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{uGCSD}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{uGCSD});
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{uGCMIN}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{uGCMIN});
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{uGCMAX}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{uGCMAX});
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{uGCCOUNT}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{uGCCOUNT});
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{uGCTOT}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{uGCTOT});
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{uLMEAN}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{uLMEAN});
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{uLSD}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{uLSD});
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{uLMIN}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{uLMIN});
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{uLMAX}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{uLMAX});
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{uLTOT}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{uLTOT});
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{uLCOUNT}  = &share({}) unless (exists $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{uLCOUNT});
                            
                            # Add GI for quick access to TAXONOMY info via GI lookup (assumes only ONE contig per filename)
                            $organisms{$_}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{GI} 
                                = $taxIdx{linked}->{$status}->{$kingdom}->{root}->{$rootDir}->{org}->{$_}->{file}->{$filename}->{contig}->{$contigID}->{gi};
                            
                            push(@contigIDs, $contigID);            # Full contigID, *NOT* nickname; Indices used in conjunction with @contigIDs
                            push(@contigFilenames, $filename);      # File basename (not full path); Indices used in conjunction with @contigIDs
                            $orgLookupByContig{$contigID} = $_;     # Full contigID, *NOT* nickname
                            #$orgLookupByFilename{$filename}= $_;    # Full organism scientific name
                            $contigID2filename{$contigID} = $filename;        # Globally shared lookup

                        } #CONTIG
                        
                        # Completed genomes: EACH filename constitutes a single replicon
                        $totalReplicons++              if($status eq "completed");               
                        $organisms{$_}->{REPL_COUNT}++ if($status eq "completed");

                    } #filename
                    
                    # Draft genomes: the sum of all filenames constitute a single replicon
                    $totalReplicons++              if($status eq "drafts");
                    $organisms{$_}->{REPL_COUNT}++ if($status eq "drafts");
                    
                } #ORG
            } #ROOTDIR
        } #kingdom
    } #status

    _stat_log("done.");
    return;
}
################################################################################
# ARGS: $orgName or \@orgNames
################################################################################
sub wipeOrg {

    # Add given org(2) to @orgs2wipe
    my @orgs2wipe = ();
    if(ref($_[0]) eq "ARRAY") {
        @orgs2wipe = @{ $_[0] };
    }
    else {
        push(@orgs2wipe, $_[0]);
    }

    foreach my $org (@orgs2wipe) {

        # Get orgName's contigIDs and filenames to wipe
        my @contigs2wipe = ();
        my @files2wipe   = keys %{ $organisms{$org}->{REPL} };
        push(@contigs2wipe, keys %{ $organisms{$org}->{REPL}->{$_}->{CONTIG} }) foreach (@files2wipe);

        # Reduce replicon count
        $totalReplicons -= $organisms{$org}->{REPL_COUNT};

        # Remove reference to org in %organisms
        delete $organisms{$org};    
        delete $orgLookupByContig{$_} foreach (@contigs2wipe);

        # Since the index in @contigIDs corresponds to the index in @contigFilenames,
        # we only need identify the index of one to delete from both arrays
        my %contigs = ();
        $contigs{ $contigIDs[$_] } = $_ foreach(0..$#contigIDs);

        # Match contigIDs to array indices
        my @matchedIdx = ();
        foreach my $contigID (@contigs2wipe) {
            push(@matchedIdx, $contigs{$contigID}) if(exists $contigs{$contigID});
        }
        
        my @unsharedContigIDs       = @{ clone(\@contigIDs) };
        my @unsharedContigFilenames = @{ clone(\@contigFilenames) };

        # Remove indices from arrays: since recurrent calls to splice() shorten
        # the array each time, we sort from greatest to smallest and work our way
        # from the end up to the beginning of the array
        foreach (sort {$b <=> $a} @matchedIdx) {
            splice(@unsharedContigIDs, $_, 1);
            splice(@unsharedContigFilenames, $_, 1);
        }
        
        @contigIDs = @unsharedContigIDs;
        @contigFilenames = @unsharedContigFilenames;
        
    } #ORG
    
    @orgNames     = keys %organisms;
    $numOrgs      = scalar(@orgNames);
    $totalContigs = scalar(@contigIDs); 
    
    return;
}
################################################################################
# GOAL: Standardize the %taxIdx hash by descending into the given hash and removing 
#       the {name} entries and use the actual $name as the hash KEY instead.
# Example: 
#   Change from: {root}->{name} = "/db/ftp"
#                {root}->{dir}  = $dir
#                {root}->{taxid}= $taxid
#            to: {root}->{"/db/ftp"}->{dir}  = $dir
#                {root}->{"/db/ftp"}->{taxid}= $taxid
# If only one entry, then {name} is present
# If multiple entries, then {name} is absent, but the key is the actual name needed
# ARGS: \%taxIdx
################################################################################
sub standardizeTaxIdx {

    #my %taxIdx = %{ $_[0] };

    foreach my $link (keys %{ $_[0] }) {
        #print "LINK=$link\n";
        foreach my $status (keys %{ $_[0]->{$link} }) {
            #print "STATUS=$status\n";
            foreach my $kingdom (keys %{ $_[0]->{$link}->{$status} }) {
                #print "KINGDOM=$kingdom\n";

                # Check whether {root} has multiple values or not
                my @roots = keys %{ $_[0]->{$link}->{$status}->{$kingdom}->{root} };
                my %rootCheck = ();
                   @rootCheck{@roots} = ();

                if(exists $rootCheck{name}) {
            
                    # Convert root's {name} to a KEY if {name} is present
                    my $rootName = $_[0]->{$link}->{$status}->{$kingdom}->{root}->{name};
                    $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$rootName} 
                        = &share({}) unless (exists $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$rootName});
                    $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$rootName}->{org} 
                        = &shared_clone(\%{ $_[0]->{$link}->{$status}->{$kingdom}->{root}->{org} });
                    delete $_[0]->{$link}->{$status}->{$kingdom}->{root}->{org};
                    delete $_[0]->{$link}->{$status}->{$kingdom}->{root}->{name};
                }
            
                # Process all roots since they are now in a consistent format
                ROOT: foreach my $root (keys %{ $_[0]->{$link}->{$status}->{$kingdom}->{root} }) {

                    #print "ROOT=$root\n";

                    # Check whether {org} has multiple values or not
                    my @orgs = keys %{ $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org} };
                    $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org} = &share({}) 
                        unless (exists $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org});
                    my %orgCheck = ();
                       @orgCheck{@orgs} = ();
                
                    # Convert single org's {name} into a KEY if {name} is present
                    if(exists $orgCheck{name}) {
                        my $orgName = $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{name};
                        $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$orgName} = &share({}) 
                            unless (exists $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$orgName});
                        $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$orgName}->{file} = &share({}) 
                            unless (exists $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$orgName}->{file});

                        $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$orgName}->{dir} 
                            = $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{dir};
                        $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$orgName}->{rank} 
                            = $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{rank};
                        $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$orgName}->{taxid} 
                            = $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{taxid};
                        $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$orgName}->{file} 
                            = &shared_clone(\%{ $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{file} });

                        delete $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{name};
                        delete $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{dir};
                        delete $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{rank};
                        delete $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{taxid};
                        delete $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{file};
                    }

                    # Process all orgs since they are now in a consistent format
                    ORG: foreach my $org (keys %{ $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org} }) {
                
                        #print "ORG=$org\n";
                
                        # Check wehther {file} has multiple values or not
                        my @files = keys %{ $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file} };
                        my %fileCheck = ();
                           @fileCheck{@files} = ();
                    
                        # Convert single file's {name} into a KEY if {name} is present
                        if(exists $fileCheck{name}) {
                            my $filename = $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{name};
                            $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{$filename} = &share({}) 
                                unless (exists $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{$filename});
                            $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{$filename}->{stype} 
                                = $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{stype};
                            $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{$filename}->{status} 
                                = $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{status};
                            $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{$filename}->{repltype} 
                                = $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{repltype};
                            $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{$filename}->{topo} 
                                = $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{topo};
                            $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{$filename}->{ntype} 
                                = $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{ntype};
                            $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{$filename}->{rc} 
                                = $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{rc};
                            $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{$filename}->{contig} 
                                = &shared_clone(\%{ $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{contig} });
                        
                            delete $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{stype};
                            delete $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{status};
                            delete $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{repltype};
                            delete $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{topo};
                            delete $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{ntype};
                            delete $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{rc};
                            delete $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{contig};
                            delete $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{name};
                        }
                
                        # Process all files since they are now in a consistent format
                        FILE: foreach my $filename (keys %{ $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file} }) {
                        
                            #print "FILENAME=$filename\n";
                        
                            my @contigs = keys %{ $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{$filename}->{contig} };
                            my %contigCheck = ();
                               @contigCheck{@contigs} = ();
                           
                            # Convert single contig's {name} into a KEY if {name} is present
                            if(exists $contigCheck{name}) {
                                my $contigID = $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{$filename}->{contig}->{name};
                                $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{$filename}->{contig}->{$contigID} = &share({}) 
                                    unless (exists $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{$filename}->{contig}->{$contigID});
                                $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{$filename}->{contig}->{$contigID}->{gi} 
                                    = $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{$filename}->{contig}->{gi};
                                $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{$filename}->{contig}->{$contigID}->{desc}
                                    = $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{$filename}->{contig}->{desc};
                                $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{$filename}->{contig}->{$contigID}->{acc}
                                    = $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{$filename}->{contig}->{acc};
                                $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{$filename}->{contig}->{$contigID}->{date}
                                    = $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{$filename}->{contig}->{date};

                                delete $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{$filename}->{contig}->{name};
                                delete $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{$filename}->{contig}->{gi};
                                delete $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{$filename}->{contig}->{desc};
                                delete $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{$filename}->{contig}->{acc};
                                delete $_[0]->{$link}->{$status}->{$kingdom}->{root}->{$root}->{org}->{$org}->{file}->{$filename}->{contig}->{date};
                            }

                        } #FILE

                    } #ORG
                
                } #ROOT

            } #KINGDOM
        } #STATUS
    } #LINK
   
    return;
    
}
################################################################################
# Function: Determine which taxonomic ranks to process in main() loop.
# ARGS: $doCompleteTax, $taxLevel
################################################################################
sub determineRanks2Analyze {

    my @toAnalyze = ();

    if($_[0]) {
        # SPECIES will be analyze TWICE if user supplies $doCompleteTax and $compareSameTaxLevel,
        # once for $compareSameTaxLevel = 1, and the other for = 0.
        push(@toAnalyze, "species") if($originalCompareSameTaxLevel && $doCompleteTax);
        push(@toAnalyze, "species", @ranks);    # (species, genus, family, order, class, phylum, superkingdom)
        pop(@toAnalyze);                        # (species, genus, family, order, class, phylum)
    }
    else {
        push(@toAnalyze, lc($_[1]));
    }

    #print "Will analyze: @toAnalyze"; <STDIN>;

    return @toAnalyze;
}
################################################################################
# Takes a replicon name and generates an integer index for it
#-------------------------------------------------------------------------------
sub createRepliconNicknames {

    # Increment $places and calculate $mod = $totalReplicons / (10**$n) until ($mod < 10)
    my $places = 0;
    for ($places = 1; $places < 10; $places++) {                   # Allow up into the billions
        my $quot = $totalReplicons/(10**$places);
        last if($quot < 10);
    }
    ++$places;

    my $pDirCode = 0;
    
    for(my $idx = 0; $idx < $totalContigs; $idx++) {
        my $filename = $contigFilenames[$idx];
        my $contigID = $contigIDs[$idx];
        #$counter = pad_zeroes(++$counter, $places);
        my $counter = pad_zeroes($idx, $places);
        $fileNickLookup{$filename} = $counter;
        $rFileNickLookup{$counter} = $filename;
        $contigNickLookup{$contigID} = $counter;
        $rContigNickLookup{$counter} = $contigID;
        
        # Assign ContigIDs to their pDir codes
        $pDirNickLookup{$counter} = $pDirCode;
        $pDirCode = ($pDirCode == ($num_pDirs-1))
                  ? (0)
                  : ($pDirCode+1);
                  
#        $contigHostPartitionLookup{$counter} = &share({});
    }
}
################################################################################
# Current implementation requires all replicon l-mer binary hashes to be found on
# disk, which could present a problem for those replicons who do not produce any
# valid l-mers for the given l-mer length (E.g. sequence is 150-bp long with an
# l-mer length of 151-bp)
#-------------------------------------------------------------------------------
sub verifyWordsOnDisk {
}
################################################################################
#    #####  #  #  #####  ####     #      ###  #####  #####    ###    ##   #
#    ##     #  #    #    #   #   # #    #       #      #     #   #   # #  #
#    ####    ##     #    ####   # # #  #        #      #    #     #  #  # #
#    ##     #  #    #    #  #   #   #   #       #      #     #   #   #   ##
#    #####  #  #    #    #   #  #   #    ###    #    #####    ###    #    #
################################################################################
################################################################################
#    #####  #  #  #####  ####     #      ###  #####  #####    ###    ##   #
#    ##     #  #    #    #   #   # #    #       #      #     #   #   # #  #
#    ####    ##     #    ####   # # #  #        #      #    #     #  #  # #
#    ##     #  #    #    #  #   #   #   #       #      #     #   #   #   ##
#    #####  #  #    #    #   #  #   #    ###    #    #####    ###    #    #
################################################################################
################################################################################
#    #####  #  #  #####  ####     #      ###  #####  #####    ###    ##   #
#    ##     #  #    #    #   #   # #    #       #      #     #   #   # #  #
#    ####    ##     #    ####   # # #  #        #      #    #     #  #  # #
#    ##     #  #    #    #  #   #   #   #       #      #     #   #   #   ##
#    #####  #  #    #    #   #  #   #    ###    #    #####    ###    #    #
################################################################################
# ******************************************************************************
# ENCODING/DECODING SUBS
# ******************************************************************************
################################################################################
#
#  DNA_string                 -----         stored/displayed    (LEFT->RIGHT)
#  bitVector                  -----         stored              (LEFT->RIGHT)
#
#  DNA_string->bitVector    dna2bv_2bit     stored              (LEFT->RIGHT)
#  bitVector->bin_string    ->to_Bin()      displayed           (RIGHT->LEFT)
#  bitVector->DNA_string    bv2dna_2bit     displayed           (LEFT->RIGHT)
#
#
################################################################################
# ARGS:    $bitMask, \$dna_seq, $doBenchmark(0/1)
# Returns: updates
#-------------------------------------------------------------------------------
sub dna2bm {

    my $iterX = new Benchmark if($_[2]);

    $_[0]->Bit_Off(
                   pos(${ $_[1] })-1
                  ) 
        while(${ $_[1] } =~ m/N/g);

    if($_[2]) {
        my $iterY = new Benchmark;
        print $_[0]->Size()."-bp vector bit-masked:  ".timestr(timediff($iterY,$iterX))."\n";
    }
    return;
}
################################################################################
# FUNCTION: This sub generates the reverse complement of $bit_vectorIN from $first
#           to $last1, storing the result in $bit_vectorOUT
#           ** READS FROM **
#
#          $_[0]                $_[1]
# ARGS: $bit_vector     $bit_vector->Size()-1
#-------------------------------------------------------------------------------
sub rc_bv {

    my $tmp_bv = $_[0]->Clone();

    # 1. Reverse
    $tmp_bv->Reverse($tmp_bv);
    
    # 2. Flip homogeneous bit pairs (00 => 11, 11 => 00)
    for (my $idx = 0; $idx < $_[1]; $idx+=2) {
        if($tmp_bv->bit_test($idx) == $tmp_bv->bit_test($idx+1)) {
            $tmp_bv->bit_flip($idx);
            $tmp_bv->bit_flip($idx+1);
        }
    }
    
    return $tmp_bv;
}
################################################################################
# FUNCTION: converts a DNA String into a Bit Vector object, a single
#           character (2 bits) at a time.
#           ** Bit vector must be initalized w/all cleared bits!! **
#
# ARGS: $bit_vector, \$dna_string, length($dna_string), 
#-------------------------------------------------------------------------------
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# OPTIMAL CODE FOR CONVERTING A DNA STRING INTO A BIT VECTOR
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
sub dna2bv_2bit7 {
    # Compile all the "C", "G", and "T" positions and perform bit flips all at once later
    my @toFlip = ();
    for my $idx (0..($_[2]-1)) {
        given (ord( substr(${ $_[1] },$idx,1) )) {                       # Grab first char and ord() it
            when (65)    { next;                                      }  # ord("A") = 65
            when (67)    { push(@toFlip, 2*$idx+1);                   }  # ord("C") = 67
            when (71)    { push(@toFlip, 2*$idx);                     }  # ord("G") = 71
            when (84)    { my $bitPos = 2*$idx;                          # ord("T") = 84
                           push(@toFlip, $bitPos, $bitPos+1);         }
            default      { next; #die "Invalid nucleotide!\n";        
                         }
        } #given
    } #for

    $_[0]->Index_List_Store(@toFlip);

    return;    
}
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
################################################################################
# FUNCTION: Converts a bit-vector into a string of DNA using 2-bit pairs
#           ** READS LEFT to RIGHT **
#
#           $_[0]         $_[1]         $_[2] = lastPositionInBitVector
# ARGS: $outputString, $bit_vector, $bit_vector->Size()-1
#-------------------------------------------------------------------------------
sub bv2dna_2bit {

    for (my $idx = 0; $idx < $_[2]; $idx+=2) {
        given ($_[1]->bit_test($idx).$_[1]->bit_test($idx+1)) {
            when ("00")     { $_[0] .= "A"; }
            when ("01")     { $_[0] .= "C"; }
            when ("10")     { $_[0] .= "G"; }
            when ("11")     { $_[0] .= "T"; }
        }
    }
    return;
}
################################################################################
# ARGS: $seqBM, $seqBM->Size()
#-------------------------------------------------------------------------------
sub rangeIt {
    my $start = 0;
    my %range = ();
    while (
           ($start < $_[1]) 
          &&
           (my ($min, $max) = $_[0]->Interval_Scan_inc($start))
          )
    {
        $start = $max + 2;
        #my $len = $max-$min+1;
        $range{$min}=$max;
    }
    #print "Given: ".$_[0]->to_Bin()."\n";
    #print "Range: \n";
    #print "    $_\t=> ".$range{$_}."\n" foreach(keys %range);
    #<STDIN>;
    
    return \%range;
}
################################################################################
# ARGS: $seqBM, $seqBM->Size(), $k
#-------------------------------------------------------------------------------
sub rangeIt_withK {
    my $start = 0;
    my %range = ();
    while (
           ($start < $_[1]) 
          &&
           (my ($min, $max) = $_[0]->Interval_Scan_inc($start))
          )
    {
        $start = $max + 2;
        my $len = $max-$min+1;
        $range{$min}=$len if($len >= $_[2]);   # Range is valid only if >= k
    }
    return \%range;
}
################################################################################
# FUNCTION: This sub expects the elements in @strings and @starts to coincide with
#           each other. That is, $strings[0] is related to $starts[0]. The first
#           string in @strings (i.e. $strings[0]) is added to @stringsDrp and its
#           start position added to @startsDrp as a single-element AREF at position
#           zero: $startsDrp[0] = [ $starts[0] ]. The next string (i.e. $strings[1]) 
#           is compared to the last unique string added to @stringsDrp, i.e. 
#           $stringsDrp[-1]. If it is a copy, then its start position is combined
#           with the start position of the last unique string, $stringsDrp[-1]. That
#           is, it is combined into the AREF: $startsDrp[0] = [ $starts[0], $starts[1] ].
#           We then look at the next string. If, however, the subsequent string is
#           unique, then we push this unique string to @stringsDrp and push its
#           start position as a single-element AREF to @startsDrp: 
#           $startsDrp[1] = [ $starts[1] ]. We procede to dereplicate the strings
#           in this manner.
#
#           Note: the starts within a AREF need to be sorted as well, so upon
#                 finding a new, unique string, we call usort() on the previous
#                 AREF of starts prior to pushing the unique string and its start
#                 to @stringsDrp and @startsDrp.
#           
# (In-place return):
#   1. @strings = array of unique strings
#   2. @starts  = array of AREFs, where the index to each AREF corresponds to the
#                 index to its corresponding string in @strings; AND, each AREF
#                 contains 1 or more integers representing the start position of
#                 said string in the original DNA sequence.
#      Note: the start positions are in the frame of the FWD strand. The subtraction
#            sub will deal with those strings that wrap around the end of the
#            contig.
#
# ARGS: \@strings, \@starts
################################################################################
sub dereplicate_inplace {
    my @stringsDrp = ();
    my @startsDrp  = ();
    
    $stringsDrp[0] = $_[0]->[0];            # Add first element of string array
    $startsDrp[0]  = [ $_[1]->[0] ];        # Add first element of starts array
    
    # First element is guaranteed a spot in the unique list. We check the rest...
    for (1..$#{ $_[0] }) {

        # New string!
        if($_[0]->[$_] ne $stringsDrp[-1]) {

            # 1. Sort previous AREF
            @{ $startsDrp[-1] } = usort @{ $startsDrp[-1] };                    # <------------------------------ Sort on @stringsDrp and then take slices for starts
            
            # 2. Add starts
            push(@stringsDrp, $_[0]->[$_]);
            push(@startsDrp, [ $_[1]->[$_] ]);

        }
        
        # Duplicate string!
        else {
            push(@{ $startsDrp[-1] }, $_[1]->[$_]);
        }
        
    }

    # Last AREF in @startsDrp hasn't had a chance to be sorted yet, so we sort now
    @{ $startsDrp[-1] } = usort @{ $startsDrp[-1] };
    
    # In-place substitution        
    @{ $_[0] } = @stringsDrp;
    @{ $_[1] } = @startsDrp;
    
    # Note: to return dereplicated arrays non-destructively, simply:
    #       return(\@stringsDrp, \@startsDrp);
    
    return;    
}
################################################################################
# FUNCTION: Partitions @array into buckets according to the prefixes in %findBucket.
#           %findBucket maps each $idxLen-character permutation to a unique position
#           in @indexed. Permutations that share the same prefix will have their
#           suffixes placed in an AREF at its specified location in @indexed.
#           The bucket size is then prepended to the list at position 0 to prevent
#           repeated renumeration.
# ARGS: \@array, \%findBucket
#-------------------------------------------------------------------------------
sub reindex_withCount {
    print "Indexing...";
    my @indexed = ();
    my $idxLen  = 0;
    
    # Get length of first key
    while(my($k,$v) = each %{ $_[1] }) {
        $idxLen = length($k);
        last;
    }

    # Split string into prefix and suffix, and index the suffix by its prefix
    #push(@{ $indexed[ $_[1]->{substr($_,0,$idxLen)} ] }, substr($_,$idxLen)) foreach (@{ $_[0] });

    foreach my $word (@{ $_[0] }) {
        my $prefix = substr($word,0,$idxLen);
        my $suffix = substr($word,$idxLen);
        push(@{ $indexed[ $_[1]->{ $prefix } ]}, $suffix);
    }
    
    
    # For each bucket, prepend the last position of its array given the addition 
    # of this new element. So if $indexed[$bucket] was originally 0..5285, the 
    # new array would be 0..5286, where $indexed[$bucket][0] = 5286, the index
    # of the last position in the array.
    foreach my $bucket (0..$#indexed) {
        next if(!$indexed[$bucket]);            # Skip empty buckets
        unshift(@{ $indexed[$bucket] }, scalar(@{ $indexed[$bucket] }) + 1);
    }
    
    @{ $_[0] } = @indexed;
    print "done.\n";
    return;
}
################################################################################
#   SUBROUTINES   SUBROUTINES   SUBROUTINES   SUBROUTINES   SUBROUTINES
################################################################################
# FUNCTION: Clips the prefix off the kmers, flips the corresponding bit in the
#           bitMask, and places all the kmer suffixes with the same prefix into
#           an AREF (bucket), whereby all suffixes remain in order. Their 
#           corresponding start positions are loaded into an AREF and stored in 
#           the same order as their suffixes. The number of suffixes present in
#           each bucket is then prepended (unshift) to the AREF.
#
#           An additional transformation vector (@trxVec) is created that 
#           maps the contents of each suffix bucket to its prefix position in the
#           bitMask. This transformation vector is valid for mapping both the 
#           buckets of suffixes and the buckets of start positions to the bitmask.
#
#                $_[0]         $_[1]           $_[2]       $_[3]
# ARGS:    \@kmers_eBV_hex, \@kmer_starts, $prefixLength, $seqIdx
#-------------------------------------------------------------------------------
sub prefixKmers {

    my @indexedKmers = ();
    my @indexedStarts = ();

    my $numBuckets = 0;
    
    foreach my $idx (0..$#{ $_[0] }) {
        my $prefix = substr($_[0]->[$idx], 0, $_[2]);       # Using the MSB ($_[2]*4) bits as prefix
        my $suffix = substr($_[0]->[$idx], $_[2]);
        my $bucket = hex($prefix);                          # Find prefix's home
        #print $_[0]->[$idx]." = $prefix + $suffix into bucket $bucket\n";
        push(@{ $indexedStarts[$bucket] }, $_[1]->[$idx]);  # Add its start to its home stack, too
        push(@{ $indexedKmers[$bucket]  }, $suffix);        # Add to its stack
    }

    # For each bucket, prepend the LAST POSITION of its array given the addition 
    # of this new element. So if $indexed[$bucket] was originally 0..5285, the 
    # new array would be 0..5286, where $indexed[$bucket][0] = 5286, the index
    # of the last position in the array.
    foreach my $bucket (0..$#indexedKmers) {
        next if(!$indexedKmers[$bucket]);           # Skip empty buckets
        $numBuckets++;                              # Track the no. of non-empty buckets
        
        # The no. of hex strings (not including the count) in current bucket. 
        # Also is used as the lastPos in the BitMask
        my $lastPos = scalar(@{ $indexedKmers[$bucket] });
        #print "Adding $lastPos to indexedKmers[$bucket]\n";
        unshift(@{ $indexedKmers[$bucket] }, $lastPos);
        #push(@{ $_[2] }, Bit::Vector->new($lastPos-1));
        #$_[2]->[$bucket] = Bit::Vector->new($lastPos);
    }

    #my $t1 = new Benchmark;
    #print "TIME: ".timestr(timediff($t1, $t0))."\n"; exit;

    #print "Press <ENTER> to START output..."; <STDIN>;
    #for (0..$#indexedKmers) {
    #    print $indexedKmers[$_]."," if($indexedKmers[$_]);
    #}
    #print "\nPress <ENTER> to continue (inside PREFIX)..."; <STDIN>; 



    # Return values in-place
    @{ $_[0] } = @indexedKmers;
    @{ $_[1] } = @indexedStarts;

    #$bucketLookup{ pad_zeroes($_[4],2) } = $numBuckets;
    $bucketLookup{ $_[3] } = $numBuckets;

    #print scalar(@indexedKmers)." kmers, ".scalar(@indexedStarts)." starts\tsize(IndexedKmers)=".total_size($_[0])."\tsize(IndexedStarts)=".total_size($_[1])."\n";

    ## ------------
    ## Display Data
    ## ------------
    #for (0..$#indexedKmers) {
    #    print "suffix=".join(",", @{ $indexedKmers[$_] })."\n";
    #    print "start =";
    #    foreach my $startAREF (@{ $indexedStarts[$_] }) {
    #        print join(",", @{ $startAREF })."\n       ";
    #    }
    #    print "\n";
    #}
    #
    
    return;

}
################################################################################
# Displaying binary strings of bitMasks (BM) or encoded BitVectors (eBV) should
# be read from right to left, since the LSB and MSB are written in the flavor
# of American numbers, where the largest digit is on the LEFT and the least most
# digit is to the RIGHT.
#
# So if I have a sequence such as "AAAnAAAAA"   =   5'-AAAnAAAAA-3'
#   Bitmask version (memory):      111011111
#   Bitmask version (displayed):   111110111
#
# Reverse Complement would be     "TTTTTnTTT"   =   5'-TTTTTnTTT-3'
# rcBitmask version (memory):      111110111
# rcBitmask version (displayed):   111011111
#
#    012345678
# 5'--AAANAAAAA--3'     seq
# 3'--TTTNTTTTT--5'
#     876543210
# MSB-111011111-LSB
# LSB-111011111-MSB
#
#
# 5'--AAANAAAAA--3'     seq
# 5'--TTTTTNTTT--3'     rcSeq
# MSB-111011111-LSB
# MSB-111110111-LSB
#
#
# Another Example:
# Given the DNA sequence "ACGTnATCGTn"  =   5'-ACGTnATCGTn-3'
#   BitMask (mem):        11110111110
#   BitMask (disp):       01111101111
# Reverse Complement is  "nACGATnACGT"  =   5'-nACGATnACGT-3'
# rcBitMask (mem):        01111101111
# rcBitMask (disp):       11110111110
#
# ------------------------------------------------------------------------------
# FUNCTION:
# ARGS: 0: \@arefOfIndices
#       1: \%kmerExtractionOpts = ( SEQFILES => $seqFiles,
#                                   PREFLEN  => $prefixLength,
#                                   K        => $k,
#                                   BVLEN    => $bvLen,
#                                   THREADS  => $numThreads,
#                                 )
# ------------------------------------------------------------------------------
# ARGS: \@arefOfIndices of @contigIDs, $totalJobs (across all nodes)
#-------------------------------------------------------------------------------
sub indexContigKmers {

    STDOUT->autoflush(1);
    #my @contigs   = @contigIDs[@{ $_[0] }];       # slice off the contigIDs we were given
    #my @filenames = @contigFilenames[@{ $_[0] }]; # slice off the contigID filenames, too
    #my $totalJobs = $_[1];

    foreach my $contigIdx (@{ $_[0] }) {
        #print "Index $contigIdx of ".join(",", @{ $_[0] })."\n";
        
        my $indexTime0 = new Benchmark;
        
        # Index for @contigIDs is the same for @contigFilenames
        my $currContig   = $contigIDs[$contigIdx];
        my $currFilename = $contigFilenames[$contigIdx]; # Full basename as loaded in from XML file
        #my $seqDir       = "$outdir/seqs";               # formerly $tmpdir = "$outdir/tmp";
        my $currOrg      = $orgLookupByContig{$currContig};
        my $rNick        = $fileNickLookup{$currFilename};  # Same as $contigNickLookup{$currContig}
        my $pDir         = $pDirCodeLookup{ $pDirNickLookup{$rNick} };
        my $repliconDir  = $pDir."/".$partitionDir."/".$rNick;
        
        # Topology is "C"ircular or "L"inear. Circular topologies get wrap-around processing
        my $doWrap       = ($organisms{$currOrg}->{REPL}->{$currFilename}->{TOPO} eq "C")
                         ? (1)
                         : (0);
        #my $gc_text   = q{};
 
        #my $repliconDir  = $seqIdx;   # Will be global, user-determined
        #my $storeDir = $outdir."/".$partitionDir."/".$repliconDir;
        #createDir($storeDir);
        #$partDirs{$seqIdx} = $storeDir;
        
            my $loadTime0 = new Benchmark;
        my $seq         = ${ retrieve $organisms{$currOrg}->{REPL}->{$currFilename}->{CONTIG}->{$currContig}->{SEQ} };
            my $loadTime1 = new Benchmark;
            my $loadTimeText = sciExpand($1, $benchRez) if(timestr(timediff($loadTime1,$loadTime0)) =~ m/^\s*(\S+)\s+/);
            my $benchmarkString .= "[$TEXT_SeqLoad=$loadTimeText";

        my $seqLen      = length($seq);

        #---------------------------------------------------------------------------
        # No shreds possible if sequence length < l-mer length (i.e. l), so return;
        #---------------------------------------------------------------------------
        if ($seqLen < $k) {
            my $warning_text = "Warning: length of sequence ($seqLen) is < WORD_SIZE ($k)";
            {
                #lock(%warning_log);
                $warning_log{$warning_text}++;
            }
            return;
        }

        createDir($repliconDir);    # Valid replicon, so create its output directory

        my @kmers_eBV_hex = ();     # hex strings of bitVectors
        my @kmer_starts   = ();     # start positions w.r.t. FWD strand in bitMask (DNA string)
    
        # ------------------------------------------------------------------------------
        # Generate BitMask (BM) for contig sequence
        # ------------------------------------------------------------------------------
        my $seqBM = Bit::Vector->new($seqLen);        # Match exact length of DNA
        $seqBM->Fill();                               # Fill with 1's
            my $bmTime0 = new Benchmark;
        dna2bm($seqBM, \$seq, 0);                     # Mark N's with 0's
            my $bmTime1 = new Benchmark;
            #print "".timestr(timediff($bmTime1,$bmTime0))."\n";
            my $bmTimeText = sciExpand($1, $benchRez) if(timestr(timediff($bmTime1,$bmTime0)) =~ m/^\s*(\S+)\s+/); 
            #print "[BitMasking=$bmTimeText]\n";
            #$kmerExtractionTime += $bmTimeText;
            $benchmarkString .= "|$TEXT_BitMasking=$bmTimeText";
        # Faster to call scalar than repeatedly calling ->Size().
        # Note: $seqBM->Size() is the same as $rcSeqBM->Size() so we just use this one.
        my $seqBM_Size = $seqBM->Size();
            # SEQUENCE BITMASK
            my $store8 = new Benchmark;
        store $seqBM, $repliconDir."/".$seqBMfilename;
            my $store9 = new Benchmark;
            my $level2store = sciExpand($1, $benchRez) if(timestr(timediff($store9, $store8)) =~ m/^\s*(\S+)\s+/);
            #$kmerExtractionTime += $level2store;
            $benchmarkString .= "|$STORE_Level2BitMasks=$level2store";
        #store \$seq, $seqDir."/".$seqFilename;

    
        # ------------------------------------------------------------------------------
        # Generate the Reverse Complement BitMask (rcSeqBM): used for shredding 
        # ------------------------------------------------------------------------------
            my $revTime0 = new Benchmark;
        my $rcSeqBM = $seqBM->Clone();
        $rcSeqBM->Reverse($rcSeqBM);                        # only need to reverse a bitmask
            my $revTime1 = new Benchmark;
            my $revTimeText = sciExpand($1, $benchRez) if(timestr(timediff($revTime1,$revTime0)) =~ m/^\s*(\S+)\s+/);
            #print "[rcBitMasking=$revTimeText]\n";
            #$kmerExtractionTime += $revTimeText;
            $benchmarkString .= "|$TEXT_rcBitMasking=$revTimeText";
        # ------------------------------------------------------------------------------
        # Determine valid (non-N) range in BitMask (seqBM and rcSeqBM)
        # ------------------------------------------------------------------------------
            my $rangeTime0 = new Benchmark;
        my %shredRange = %{ rangeIt_withK($seqBM, $seqBM_Size, $k) };
            my $rangeTime1 = new Benchmark;
            my $rangeTimeText = sciExpand($1, $benchRez) if(timestr(timediff($rangeTime1,$rangeTime0)) =~ m/^\s*(\S+)\s+/);
            #print "[Ranging=$rangeTimeText]\n";
            #$kmerExtractionTime += $rangeTimeText;
            $benchmarkString .= "|$TEXT_Ranging=$rangeTimeText";
        #my %rcShredRange = %{ rangeIt_withK($rcSeqBM, $seqBM_Size, $k) };     # Direct, via iterating over entire BitMask
            my $rcRangeTime0 = new Benchmark;
        my %rcShredRange = ();
        while(my($start,$len) = each %shredRange) {         # Indirect, via calculation from existing range
            $rcShredRange{$seqLen-$len-$start} = $len;      # Note: Given start1=>len1 from %shredRange,
        }                                                   #       (start2,len2) from %rcShredRange are
                                                            #       start2 = seqLen - (len1 + start1)
                                                            #       len2   = len1
            my $rcRangeTime1 = new Benchmark;
            my $rcRangeTimeText = sciExpand($1, $benchRez) if(timestr(timediff($rcRangeTime1,$rcRangeTime0)) =~ m/^\s*(\S+)\s+/);
            #print "[rcRanging=$rcRangeTimeText]\n";
            #$kmerExtractionTime += $rcRangeTimeText;
            $benchmarkString .= "|$TEXT_rcRanging=$rcRangeTimeText";

        # ------------------------------------------------------------------------------
        # Generate Encoded BitVector (seq_eBV) for contig sequence
        # ------------------------------------------------------------------------------
            my $encodeTime0 = new Benchmark;
        my $seq_eBV = Bit::Vector->new(2*$seqLen);
        dna2bv_2bit7($seq_eBV, \$seq, $seqLen);       # Generate eBV from seq
            my $encodeTime1 = new Benchmark;
            my $encodeTimeText = sciExpand($1, $benchRez) if(timestr(timediff($encodeTime1,$encodeTime0)) =~ m/^\s*(\S+)\s+/);
            #print "[Encoding=$encodeTimeText]\n";
            #$kmerExtractionTime += $encodeTimeText;
            $benchmarkString .= "|$TEXT_Encoding=$encodeTimeText";

        # Faster if we put in a scalar instead of repeatedly calling ->Size()
        my $seq_eBV_Size = $seq_eBV->Size();

        # Don't need the actual sequence anymore...
        undef $seq;

        # ------------------------------------------------------------------------------
        # Generate Reverse Complement of the Encoded BitVector (rcSeq_eBV)
        # ------------------------------------------------------------------------------
            my $rcTime0 = new Benchmark;
        my $rcSeq_eBV = rc_bv($seq_eBV, $seq_eBV_Size-1);
            my $rcTime1 = new Benchmark;
            my $rcTimeText = sciExpand($1, $benchRez) if(timestr(timediff($rcTime1, $rcTime0)) =~ m/^\s*(\S+)\s+/);
            #print "[rcEncoding=$rcTimeText]\n";
            #$kmerExtractionTime += $rcTimeText;
            $benchmarkString .= "|$TEXT_rcEncoding=$rcTimeText";

        # ------------------------------------------------------------------------------
        # Shred FWD strand of $seq_eBV according to %shredRange
        # ------------------------------------------------------------------------------
        # Given $start and $stop of a valid shred range, we shred it inline here
            my $fwdShred0 = new Benchmark;
        while(my($start,$len) = each %shredRange) {                     # bitMask coords
            my $gap = $start + $len - $k;                               # bitMask coords 
            for(my $offset = $start; $offset <= $gap; ++$offset) {      # bitMask coords
                # -----------------
                # Forward Direction
                # -----------------
                my $fwd = Bit::Vector->new($bvLen);
                $fwd->Interval_Copy($seq_eBV, 0, 2*$offset, $bvLen);
                push(@kmers_eBV_hex, $fwd->to_Hex());
                push(@kmer_starts, $offset);                                     # bitMask coords
            }
        }
            my $fwdShred1 = new Benchmark;
            my $fwdShredTimeText = sciExpand($1, $benchRez) if(timestr(timediff($fwdShred1, $fwdShred0)) =~ m/^\s*(\S+)\s+/);
            #print "[kmerDecompositionFWD=$fwdShredTimeText]\n";
            #$kmerExtractionTime += $fwdShredTimeText;
            $benchmarkString .= "|$TEXT_kmerDecompositionFWD=$fwdShredTimeText";

        # ------------------------------------------------------------------------------
        # Shred REV strand of $rcSeq_eBV according to %shredRange
        # ------------------------------------------------------------------------------
        # Note: Currently, I cannot combine with FWD shred's while loop because $start/$len 
        #       in %rcShredRange is different from %shredRange. There is linear relationship between 
        #       the two so it is possible to determine that linear map and combine the two, but I
        #       haven't done that yet.
        #
        # Note: Mapping function is the following:
        #       For any ($start1, $stop1) pair in %shredRange, the corresponding values in %rcShredRange are
        #       ($start2, $stop2)
        #       where   $start2 = $lastPos - $stop1
        #               $stop2  = $lastPos - $start1

            my $revShred0 = new Benchmark;
        my $rcConstOffset = $seqLen-$k;
        while(my($start,$len) = each %rcShredRange) {
            my $gap = $start + $len - $k;
            for(my $offset = $start; $offset <= $gap; ++$offset) {
                # -----------------
                # Reverse Direction
                # -----------------
                my $rev = Bit::Vector->new($bvLen);
                $rev->Interval_Copy($rcSeq_eBV, 0, 2*$offset, $bvLen);
                push(@kmers_eBV_hex, $rev->to_Hex());
                push(@kmer_starts, $rcConstOffset-$offset);  # Simplify: $const = $seqLen - $k
            }
        }
            my $revShred1 = new Benchmark;
            my $revShredTimeText = sciExpand($1, $benchRez) if(timestr(timediff($revShred1, $revShred0)) =~ m/^\s*(\S+)\s+/);
            #print "[kmerDecompositionREV=$revShredTimeText]\n";
            #$kmerExtractionTime += $revShredTimeText;
            $benchmarkString .= "|$TEXT_kmerDecompositionREV=$revShredTimeText";

        # ------------------------------------------------------------------------------
        # Handle Wrap-Arounds
        # ------------------------------------------------------------------------------
        if($doWrap) {

            my $wrapConst  = $k-1;
            my $wrapLength = $wrapConst << 1;   # same as 2*($k-1) = 2*$k-2

                my $fwdWrapShred0 = new Benchmark;
            # --------------------------------------------------------------------------
            # Generate FWD Wrap-Around and shred
            # --------------------------------------------------------------------------
            # -------------------
            # Wrap-Around BitMask
            # -------------------
            my $leadBM = Bit::Vector->new($wrapConst);
            my $tailBM = Bit::Vector->new($wrapConst);
            $tailBM->Interval_Copy($seqBM, 0,                      0, $wrapConst);
            $leadBM->Interval_Copy($seqBM, 0, $seqBM_Size-$wrapConst, $wrapConst);
            my $wrapBM = Bit::Vector->Concat_List($tailBM, $leadBM);
            # ---------------------
            # Wrap-Around BitVector
            # ---------------------
            my $lead_eBV = Bit::Vector->new($wrapLength);   # same as 2*($k-1)
            my $tail_eBV = Bit::Vector->new($wrapLength);   # same as 2*($k-1)
            $tail_eBV->Interval_Copy($seq_eBV, 0,                            0, $wrapLength);
            $lead_eBV->Interval_Copy($seq_eBV, 0, $seq_eBV_Size-$wrapLength, $wrapLength);
            my $wrap_eBV = Bit::Vector->Concat_List($tail_eBV, $lead_eBV);
            # -----------------
            # Shred Wrap-Around
            # -----------------
            # Pos (in FWD frame) will be at least:  length($seqBM)-$k+1
            #                    will never exceed: length($seqBM)-1
            # Note: length(wrapAround) = 2(k-1)
            #       gap = length(wrapAround) - k
            #           = 2k-2-k 
            #           = k-2 (which is a const, so we can pull out of the while loop)
            my %wrapShredRange  = %{ rangeIt_withK($wrapBM, $wrapLength, $k) };
            my $gap             = $k - 2;
            my $wrapConstOffset = $seqLen-$wrapConst;
    
            while(my($start,$len) = each %wrapShredRange) {
                for(my $offset = $start; $offset <= $gap; ++$offset) {
                    my $fwd = Bit::Vector->new($bvLen);
                    $fwd->Interval_Copy($wrap_eBV, 0, 2*$offset, $bvLen);
                    #my $tmp = $fwd->Clone(); $tmp->Reverse($tmp);
                    push(@kmers_eBV_hex, $fwd->to_Hex());
                    push(@kmer_starts, $wrapConstOffset+$offset);   # Note: seqLen-k+1 = const, so
                                                                    # push(@kmer_starts, const+$offset);
                }
            }
                my $fwdWrapShred1 = new Benchmark;
                my $fwdWrapShredTimeText = sciExpand($1, $benchRez) if(timestr(timediff($fwdWrapShred1, $fwdWrapShred0)) =~ m/^\s*(\S+)\s+/);
                #print "[kmerDecompositionFWDwrap=$fwdWrapShredTimeText]\n";
                #$kmerExtractionTime += $fwdWrapShredTimeText;
                $benchmarkString .= "|$TEXT_kmerDecompositionFWDwrap=$fwdWrapShredTimeText";
        
                my $revWrapShred0 = new Benchmark;
            # --------------------------------------------------------------------------
            # Generate REV Wrap-Around and shred
            # --------------------------------------------------------------------------
            # -------------------
            # Wrap-Around BitMask
            # -------------------
            my $rcLeadBM = Bit::Vector->new($wrapConst);
            my $rcTailBM = Bit::Vector->new($wrapConst);
            $rcTailBM->Interval_Copy($rcSeqBM, 0,                      0, $wrapConst);
            $rcLeadBM->Interval_Copy($rcSeqBM, 0, $seqBM_Size-$wrapConst, $wrapConst);  # seqBM->Size() = rcSeqBM->Size()
            my $rcWrapBM = Bit::Vector->Concat_List($rcTailBM, $rcLeadBM);    
            # ---------------------
            # Wrap-Around BitVector
            # ---------------------
            my $rcLead_eBV = Bit::Vector->new($wrapLength);  # same as 2*($k-1)=2*$k-2
            my $rcTail_eBV = Bit::Vector->new($wrapLength);  # same as 2*($k-1)=2*$k-2
            $rcTail_eBV->Interval_Copy($rcSeq_eBV, 0,                            0, $wrapLength);
            $rcLead_eBV->Interval_Copy($rcSeq_eBV, 0, $seq_eBV->Size()-$wrapLength, $wrapLength);   # seq_eBV->Size() = rcSeq_eBV->Size()
            my $rcWrap_eBV = Bit::Vector->Concat_List($rcTail_eBV, $rcLead_eBV);
            # -----------------
            # Shred Wrap-Around
            # -----------------
            # Note: $gap is the same for the $rcWrap_eBV as it is for $wrap_eBV
            # Note: $gap isn't modified, so we use the same value here
            my %rcWrapShredRange = %{ rangeIt_withK($rcWrapBM, $wrapLength, $k) };
            # Note: $rcWrapConstOffset = $seqLen-$k+1+$gap = $seqLen-$k+1+($k-2) = $seqLen-1;
            my $rcWrapConstOffset = $seqLen-1;
            while(my($start,$len) = each %rcWrapShredRange) {
                for(my $offset = $start; $offset <= $gap; ++$offset) {
                    my $rev = Bit::Vector->new($bvLen);
                    $rev->Interval_Copy($rcWrap_eBV, 0, 2*$offset, $bvLen);
                    push(@kmers_eBV_hex, $rev->to_Hex());
                    push(@kmer_starts, $rcWrapConstOffset-$offset);   # Note: $const1 = $seqLen-$k+1
                                                                      # Note: $const2 = $gap = $k-2
                                                                      # Note: $const  = $const1 + $const2
                }
            }
                my $revWrapShred1 = new Benchmark;
                my $revWrapShredTimeText = sciExpand($1, $benchRez) if(timestr(timediff($revWrapShred1, $revWrapShred0)) =~ m/^\s*(\S+)\s+/);
                #print "[kmerDecompositionREVwrap=$revWrapShredTimeText]\n";
                #$kmerExtractionTime += $revWrapShredTimeText;
                $benchmarkString .= "|$TEXT_kmerDecompositionREVwrap=$revWrapShredTimeText";
    
        }   # DOWRAP

        ################################################################################
        # At this point, we have the kmers and their start positions from both FWD and
        # REV strands, along with the (optional) wrap-arounds
        ################################################################################
        #print "kmers:\n";
        ##print $_->to_Bin()."\n" foreach (@kmers);
        #print $_."\n" foreach (@kmers_eBV_hex);
        #print "kmer_starts:\n";
        #print $_."\n" foreach (@kmer_starts);

        #-------------------------------------------------------------------------------
        # SORT (kmers and kmer_starts)
        #-------------------------------------------------------------------------------
        # -> need sorted indices because both @kmers_eBV_hex and @kmer_starts must correlate;
        # -> ukeysort is the fastest method I could find. First we convert the hex strings into
        #    their integer equivalents. Then we sort these values as unsigned ints (ukeysort),
        #    thus generating a list of sorted indices.
        #
        # FAST ...
        #@kmers_eBV_hex_sortIdx = skeysort { $kmers_eBV_hex[$_] } (0..$#kmers_eBV_hex);
    
        # FASTER ...
        #@kmers_eBV_hex_sortIdx = skeysort { pack "H*", $kmers_eBV_hex[$_] } (0..$#kmers_eBV_hex);
    
        # FASTEST ...
        # ...code block so @sortedIndices will go out of scope when done...
        {
            my $sortTime0     = new Benchmark;
        #my @sortedIndices = ukeysort { hex $kmers_eBV_hex[$_] } (0..$#kmers_eBV_hex);
        my @sortedIndices = keysort { $kmers_eBV_hex[$_] } (0..$#kmers_eBV_hex);
        #die "Index count mismatch!\n" if(scalar(@sortedIndices1) != scalar(@sortedIndices2));
        #for(0..$#sortedIndices1) {
        #    print $sortedIndices1[$_]." != ".$sortedIndices2[$_]."\n" if($sortedIndices1[$_] != $sortedIndices2[$_]);
        #}

            my $sortTime1     = new Benchmark;
            my $sortTimeText  = sciExpand($1, $benchRez) if(timestr(timediff($sortTime1,$sortTime0)) =~ m/^\s*(\S+)\s+/);
            #print "[Sort=$sortTimeText]\n";
            #$kmerExtractionTime += $sortTimeText;
            $benchmarkString .= "|$TEXT_Sort=$sortTimeText";

        # Use slice() to create sorted (i) kmers and sorted (ii) starts
            my $sortSliceTime0 = new Benchmark;
        @kmers_eBV_hex     = @kmers_eBV_hex[@sortedIndices];    # in-place sort
        @kmer_starts       = @kmer_starts[@sortedIndices];      # in-place sort
            my $sortSliceTime1 = new Benchmark;
            my $sortSliceTimeText = sciExpand($1, $benchRez) if(timestr(timediff($sortSliceTime1,$sortSliceTime0)) =~ m/^\s*(\S+)\s+/);
            #print "[Slice=$sortSliceTimeText]\n";
            #$kmerExtractionTime += $sortSliceTimeText;
            $benchmarkString .= "|$TEXT_Slice=$sortSliceTimeText";
        #print $kmers_eBV_hex[$_]." => ".$kmer_starts[$_]."\n" for(0..$#kmers_eBV_hex);
        }

        #-------------------------------------------------------------------------------
        # DEREPLICATE IN-PLACE (kmers and kmer_starts)
        #-------------------------------------------------------------------------------
        # We leave the value at pos 0 untouched and compare everything else in the
        # sorted array to it.
            my $derepTime0 = new Benchmark;
        dereplicate_inplace(\@kmers_eBV_hex, \@kmer_starts);
            my $derepTime1 = new Benchmark;
            
            my $numKmers = scalar(@kmers_eBV_hex);
            
            #print "Generated ".$numKmers." unique kmers\n" if($debug_kmerExtract);
            my $derepTimeText = sciExpand($1, $benchRez) if(timestr(timediff($derepTime1,$derepTime0)) =~ m/^\s*(\S+)\s+/);
            #print "[Dereplication=$derepTimeText]\n";
            #$kmerExtractionTime += $derepTimeText;
            $benchmarkString .= "|$TEXT_Dereplication=$derepTimeText";
        #print $kmers_eBV_hex[$_]." => ".join(",",@{ $kmer_starts[$_] })."\n" for(0..$#kmers_eBV_hex);
    
        # These values are needed only for the statistics sub
        $organisms{$currOrg}->{REPL}->{$currFilename}->{CONTIG}->{$currContig}->{WORDS}->{COUNT} = $numKmers;
        $organisms{$currOrg}->{REPL}->{$currFilename}->{CONTIG}->{$currContig}->{WORDS}->{GCMEAN}= 0;
        $organisms{$currOrg}->{REPL}->{$currFilename}->{CONTIG}->{$currContig}->{WORDS}->{GCSD}  = 0;
        $organisms{$currOrg}->{REPL}->{$currFilename}->{CONTIG}->{$currContig}->{WORDS}->{GCMIN} = 0;
        $organisms{$currOrg}->{REPL}->{$currFilename}->{CONTIG}->{$currContig}->{WORDS}->{GCMAX} = 0;    
        #-------------------------------------------------------------------------------
        # PREFIX & SUFFIX (kmers and kmer_starts)
        #-------------------------------------------------------------------------------
        # -----------------------------------
        # Structures needed for Intersection:
        # -----------------------------------
        # @kmers_eBV_hex:   Kmer suffixes partitioned into buckets based on common prefixes [1..n]
        # @kmer_starts:     Corresponding start positions for the kmer suffix buckets above [1..n]
        # @kmersLevel2BitMasks: BitMasks for the kmers in each bucket                       [1..n]

            my $prefixTime0 = new Benchmark;
        #my @kmersLevel2BitMasks = ();
        prefixKmers(\@kmers_eBV_hex, \@kmer_starts, $prefixLength, $rNick);
        #prefixKmers(\@kmers_eBV_hex, \@kmer_starts, $prefixLength, $contigIdx);
        #prefixKmers2(\@kmers_eBV_hex, \@kmer_starts, \@kmersLevel2BitMasks, $prefixLength, $contigIdx);
            my $prefixTime1 = new Benchmark;
            my $prefixTimeText = sciExpand($1, $benchRez) if(timestr(timediff($prefixTime1, $prefixTime0)) =~ m/^\s*(\S+)\s+/);
            #print "[Prefixing=$prefixTimeText]\n";
            #$kmerExtractionTime += $prefixTimeText;
            $benchmarkString .= "|$TEXT_Prefixing=$prefixTimeText";

        #print "Press <ENTER> to see prefix BitMask..."; <STDIN>;
        #print $prefixBitMask->to_Bin()."\n";
    
        #print "Press <ENTER> to see how many buckets there are..."; <STDIN>;
        my $totalBuckets = $bucketLookup{$rNick};
        #my $totalBuckets = $bucketLookup{$contigIdx};
        #for (0..$#kmers_eBV_hex) {
        #    $totalBuckets++ if($kmers_eBV_hex[$_]);
        #}



        #print $bucketLookup{$contigIdx}." buckets (".sprintf("%.$res", $bucketLookup{$contigIdx}/($baseNum**$prefixLength)*100)." \% theor.)\n" if($debug_kmerExtract);
    
        #print "Press <ENTER> to see what the buckets look like..."; <STDIN>;
        #for (0..$#kmers_eBV_hex) {
        #    print "    Bucket $_: ".$kmers_eBV_hex[$_]->[0]." elements\n";
        #}
    
        #print "Level2 BitMasks:\n";
        #print "  ".$_->to_Bin()."\n" foreach (@kmersLevel2BitMasks);
        #print "kmers_suffixes:\n";
        #print "    ".join(",", @{ $kmers_eBV_hex[$_] })."\n" foreach (0..$#kmers_eBV_hex);
        #print "kmer_starts:    @kmer_starts\n";

        #my $repliconDir  = pad_zeroes($contigIdx,2);   # Will be global, user-determined
        
            # KMER SUFFIXES
            my $store2 = new Benchmark;
        store \@kmers_eBV_hex, $repliconDir."/".$kmers_eBV_hex_Filename;
            my $store3 = new Benchmark;
            my $kmersStore = sciExpand($1, $benchRez) if(timestr(timediff($store3, $store2)) =~ m/^\s*(\S+)\s+/);
            #$kmerExtractionTime += $kmersStore;
            $benchmarkString .= "|$STORE_kmerSuffixes=$kmersStore";
            undef @kmers_eBV_hex;
        
            # KMER STARTS
            my $store4 = new Benchmark;
        store \@kmer_starts, $repliconDir."/".$kmer_starts_Filename;
            my $store5 = new Benchmark;
            my $kmerStartsStore = sciExpand($1, $benchRez) if(timestr(timediff($store5, $store4)) =~ m/^\s*(\S+)\s+/);
            #$kmerExtractionTime += $kmerStartsStore;
            $benchmarkString .= "|$STORE_kmerStarts=$kmerStartsStore";
            undef @kmer_starts;
        
        # We want the contigIDs in the same sorted order (decreasing replicon size)
        # so we put their nicknames in the same positions as in @contigIDs. Empty
        # array elements corresponding to contigIDs (filtered out by length < k)
        # will be removed in main thread.
        #lock(@processedContigs);
        $processedContigs[$contigIdx] = $rNick;     # Add nickname to same position

        my $indexTime1 = new Benchmark;
        my $indexTimeText = timestr(timediff($indexTime1, $indexTime0));
        if($indexTimeText =~ m/^\s*(\S+)\s+/) {
            $indexTimeText = sciExpand($1, $benchRez);
        }
        elsif($indexTimeText =~ m/(\d+\.?\d+)\s+/) {
            $indexTimeText = $1;
        }
        
        $benchmarkString .= "|$TEXT_TotalTime=$indexTimeText]";

        my $res = "2f";
        #print $totalBuckets." buckets (".sprintf("%.$res", $totalBuckets/($baseNum**$prefixLength)*100)." \% theor.)\n" if($debug_kmerExtract);
        #my $pctBuckets = sprintf("%.$res", $totalBuckets/($baseNum**$prefixLength)*100);
        my $pctBuckets = sprintf("%.$res", $totalBuckets/($MAX_BUCKETS)*100);
        
        # Thread ### (of ###): Streptococcus pyogenes strain XYZ, gi|######|ref|XXXXXXX (256 buckets:100%)
        my $currCount    = 0;
        {
            lock($ripCount);
            $ripCount++;
            $currCount = $ripCount;
        }
        $benchmarkString = "Thread $currCount (of ".$_[1]."): \"$numKmers\" $k-mers in [$indexTimeText] wallsecs\n"
                          ."                          "
                          .$benchmarkString."\n"
                          ."                          "
                          ."[$totalBuckets/$MAX_BUCKETS buck:$pctBuckets%%] //$currOrg////$currContig//\n";
#                          .$benchmarkString."\n";
        #stat_log("Thread $currCount (of ".$_[1]."):  $ripOutput");


        
        #print $benchmarkString."\n";
        stat_log("$benchmarkString");

        #print "Index Time = ".timestr(timediff($indexTime1, $indexTime0))."\n";
        
    } #contigIdx

    return;
}
################################################################################
# ARGS: none (all globals)
#-------------------------------------------------------------------------------
sub extractKmersParallel {
    
    # Create temp dir for Storable versions of FASTA files
    my $tmpdir = "$outdir/tmp";
    if(!-d $tmpdir) {                           # tmpdir doesn't already exist
        my $mkdirCmd= `mkdir $tmpdir 2>&1`;     # make tmpdir
        chomp $mkdirCmd;
        if($mkdirCmd ne "") {
            stat_log("Fatal: Unable to create directory \"$tmpdir\"!");
            die "\n";
        }
    }
    
    STDOUT->autoflush(1);
    stat_log("----------------------------------");
    stat_log("Decomposing contigs into words of length $k...");
    @processedContigs = ();           # Global; Need to reset between subsequent iterations of $k
    
    # Init benchmark
    my $iter_time0 = new Benchmark;

    my $totalContigIDs  = scalar(@contigIDs);
    #my $loadBalancedContigs_aref = loadBalanceLinear(\@contigIDs, $nThreads);
    my $loadBalancedIndices = loadBalanceLinearPaired([0..$#contigIDs], $nThreads);
    
    #print "Total of $totalContigIDs\n";
    #foreach my $aref (@{ $loadBalancedIndices}) {
    #    print join(",", @{ $aref })."\n";
    #}
    
    my @jobs = ();
    foreach (@{ $loadBalancedIndices }) {
        my $tid = threads->new(\&indexContigKmers, $_, $totalContigIDs);
        push(@jobs, $tid);
    }
    
    $_->join() foreach (@jobs);

    # Sort order in @processedContigs is maintained.
    # Scan @processedContigs for empty indices (those that failed extraction); remove
    # ...Scan must happen from end-bounds to start-bounds
    # **** Splice not implemented for shared arrays ****
    #for (reverse 0..$#processedContigs) {
    #    splice(@processedContigs, $_, 1) if(!defined $processedContigs[$_]);
    #}
    # **************************************************

    my @processedContigsUpdated : shared = ();
    foreach (0..$#processedContigs) {
        push(@processedContigsUpdated, $processedContigs[$_]) if(defined $processedContigs[$_]);
    }
    @processedContigs = @processedContigsUpdated;
    
    # End benchmark
    my $iter_time1 = new Benchmark;
    my $timeText = timestr(timediff($iter_time1, $iter_time0));
    $timeText = $1 if($timeText =~ m/^\s*(\d+)\s+/);
    $timeText = $1 if($timeText =~ m/(\d+\.?\d+)\s+/);
    
    stat_log("TOTAL EXTRACTION TIME: [$timeText] wallsecs");

    displayWarnings();
    
}
##############################################################################################
#                                                                                                                
#    #####  #   #  #####  #####  ####    ####  #####    ###  #####  #####    ###    ##   #                                                                                              
#      #    ##  #    #    #      #   #  #      #       #       #      #     #   #   # #  #                                      
#      #    # # #    #    ####   ####    ###   ####   #        #      #    #     #  #  # #                                                 
#      #    #  ##    #    #      #  #       #  #       #       #      #     #   #   #   ##                                      
#    #####  #   #    #    #####  #   #  ####   #####    ###    #    #####    ###    #    #                                                       
#
##############################################################################################
##############################################################################################
#                                                                                                                
#    #####  #   #  #####  #####  ####    ####  #####    ###  #####  #####    ###    ##   #                                                                                              
#      #    ##  #    #    #      #   #  #      #       #       #      #     #   #   # #  #                                      
#      #    # # #    #    ####   ####    ###   ####   #        #      #    #     #  #  # #                                                 
#      #    #  ##    #    #      #  #       #  #       #       #      #     #   #   #   ##                                      
#    #####  #   #    #    #####  #   #  ####   #####    ###    #    #####    ###    #    #                                                       
#
##############################################################################################
##############################################################################################
#                                                                                                                
#    #####  #   #  #####  #####  ####    ####  #####    ###  #####  #####    ###    ##   #                                                                                              
#      #    ##  #    #    #      #   #  #      #       #       #      #     #   #   # #  #                                      
#      #    # # #    #    ####   ####    ###   ####   #        #      #    #     #  #  # #                                                 
#      #    #  ##    #    #      #  #       #  #       #       #      #     #   #   #   ##                                      
#    #####  #   #    #    #####  #   #  ####   #####    ###    #    #####    ###    #    #                                                       
#
##############################################################################################
################################################################################
# ARGS: $rOrg, $qOrg, $taxLevel
# Returns TRUE if tax level is the same
# IMPORTANT: NO TAXONOMIC RESOLUTION PROVIDED BELOW THE SPECIES LEVEL.
#            ALL FURTHER DIVISIONS (SUBSPECIES, STRAINS, ...) ARE COMBINED
#            INTO THE "subspecies" CATEGORY. SERIOUS UPGRADES ARE REQUIRED
#            TO PROVIDE FURTHER RESOLUTION:
#            (i)   program used to generate "bySpecies_full4.dmp"
#            (ii)  %taxAbbr
#            (iii) this sub
#-------------------------------------------------------------------------------
sub isSameTaxLevel {

    # Each organism has their complete taxonomy tree associated with their name,
    # including SK, P, C, O, F, G, S and SS. To match these, just probe with an 
    # "eq" at the appropriate tax level. Tax info also exists for the subspecies
    # (SS) rank, however this rank contains the merged tax info of divisions 
    # deeper than subspecies alone. Comparisons should be made up to the species
    # (S) level only.

    my $wantedRank = $taxAbbrExt{$_[2]};

    my $rRank = $organisms{$_[0]}->{TAXTREE}->{$wantedRank};
    my $qRank = $organisms{$_[1]}->{TAXTREE}->{$wantedRank};

    if(!defined $rRank) {
        stat_log("**WARNING**: RANK ".$_[2]." not defined for \"".$_[0]."\"!");
        return 0;
    }
    if(!defined $qRank) {
        stat_log("**WARNING**: RANK ".$_[2]." not defined for \"".$_[1]."\"!");
        return 0;
    }

    my $result = ($rRank eq "Unassigned") ? (0)
               : ($qRank eq "Unassigned") ? (0)
               : ($rRank eq $qRank)       ? (1)
               :                            (0);

    return $result;
}
################################################################################
# FUNCTION: Intersect sets of kmers from the prescribed organisms.
# ARGS: 0: \@arefOfIndices
# ------------------------------------------------------------------------------
sub intersectKmers {

    my %localResults = ();
    
    # Thus, if($bucketSwap), 
    #
    # "Natural" order is REF > QRY, therefore the "natural" output order is 
    #     $largestArray   = REF
    #     $smallestArray  = QRY
    # and
    #     $largestBucket  = REF
    #     $smallestBucket = QRY
    #
    # if($bucketSwap): 
    #                   if($arraySwap): orientation maintained  => $largestBucket writes to REF
    #                                                              $smallestBucket writes to QRY
    #                                                              -----------------------------
    #                            else : orientation reversed    => $smallestBucket writes to REF
    #                                                              $largestBucket writes to QRY
    #           else :
    #                   if($arraySwap): orientation reversed    => $smallestBucket writes to REF
    #                                                              $largestBucket writes to QRY
    #                                                              -----------------------------
    #                            else : orientation maintained  => $largestBucket writes to REF
    #                                                              $smallestBucket writes to QRY
    # 
    #   $arraySwap | $bucketSwap | $largestBucket writes to: | $smallestBucket writes to:
    #  ------------|-------------|---------------------------|----------------------------
    #       0             0   (XOR=0)       REF                         QRY         $destination = ($arraySwap->Xor($arraySwap, $bucketSwap))
    #       0             1   (XOR=1)       QRY                         REF                      ? (\@qryOutput)
    #       1             0   (XOR=1)       QRY                         REF                      : (\@refOutput);
    #       1             1   (XOR=0)       REF                         QRY
    
    #my @partDirs = keys %partDirs;
    
    #print "AREF=".join(",", @{ $_[0] })."\n";
    my $threadTime0 = new Benchmark;
    my $refOutputText = q{};        # Text displayed about ref intersections just before returning to main thread
    # Cycle through the slice of @processedContigs containing rContigNick's
    REF: for my $rIdx (@{ $_[0] }) {
        #my $rDir = $outdir."/".$partitionDir."/".pad_zeroes($rIdx,2);
        my $rNick     = $processedContigs[$rIdx];
        #my $rDir      = $outdir."/".$partitionDir."/".$rNick;
        my $pDir      = $pDirCodeLookup{ $pDirNickLookup{$rNick} };
        my $rDir      = $pDir."/".$partitionDir."/".$rNick;
        my $rContigID = $rContigNickLookup{ $rNick };
        my $rOrg      = $orgLookupByContig{ $rContigID };
        my $rFilename = $rFileNickLookup{$rNick};
        my $numKmers  = $organisms{$rOrg}->{REPL}->{$rFilename}->{CONTIG}->{$rContigID}->{WORDS}->{COUNT};
        
        my $refTime0 = new Benchmark;
#        print "Loading REF (rIdx=$rIdx, nick=$rNick) hex strings...";
        my $rkmers0   = new Benchmark;
        my @rKmers    = @{ retrieve $rDir."/".$kmers_eBV_hex_Filename };
        my $rkmers1   = new Benchmark;
        my $rLoadTime = timestr(timediff($rkmers1, $rkmers0));
           $rLoadTime = sciExpand($1, $benchRez) if($rLoadTime =~ m/^\s*(\S+)\s+/);
           $rLoadTime = $1 if($rLoadTime =~ m/^(\d+\.?\d+)\s+/);
           
           
#        print "done [".sciExpand($1, $benchRez)."]\n" if(timestr(timediff($rkmers1, $rkmers0)) =~ m/^\s*(\S+)\s+/);

        # Loop through all of ref's buckets, creating bitMasks for each, if suffixes exist in that bucket
#        print "Generating BitMasks for REF ($rNick)...";
        my @rL2BM = ();
        my $rBM0  = new Benchmark;
        for(0..$#rKmers) {
            #print "***Creating REF bitmask [$_] of length ".$rKmers[$_]->[0]."***\n";
            $rL2BM[$_] = Bit::Vector->new($rKmers[$_]->[0]) if($rKmers[$_]);
        }
        my $rBM1  = new Benchmark;
#        print "done [".sciExpand($1, $benchRez)."]\n" if(timestr(timediff($rBM1, $rBM0)) =~ m/^\s*(\S+)\s+/);

        #my $rBM0      = new Benchmark;
        #my @rL2BM     = @{ retrieve $rDir."/kmersLevel2BitMasks.dmp" };
        #my $rBM1      = new Benchmark;
        #print "done [".sciExpand($1, $benchRez)."]\n" if(timestr(timediff($rBM1, $rBM0)) =~ m/^\s*(\S+)\s+/);

        # $organisms{$org}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{BASENUM} = 16;
        # $organisms{$org}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{PREFLEN} = 4;
        # $organisms{$org}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{BUCKETS} = $M;
        #my $maxBuckets = scalar(@rKmers);       # or (BASENUM**PREFLEN)
        my $rBuckets   = $bucketLookup{$rNick};  # ->{BUCKETS};  No. of buckets filled

        #print "Checking if kmers are sorted...\n";
        #for my $bucket (0..$#rKmers) {
        #    next if(!$rKmers[$bucket]);
        #    for (2..$#{ $rKmers[$bucket] }) {
        #        if($rKmers[$bucket]->[$_] lt $rKmers[$bucket]->[$_-1]) {
        #            print "ORDERING ERROR!!: \$rKmers[$bucket]->[$_] is NOT > \$rKmers[$bucket]->[".($_-1)."]\t".$rKmers[$bucket]->[$_]." is NOT > ".$rKmers[$bucket]->[$_-1]."\n";
        #        }
        #    }
        #}

        my $totalQryTime = 0;
        my $totalQryLoaded = 0;
        
        # $rIdx alone forces self-vs-self intersections ($rIdx+1 to skip self)
        #QRY: for my $qIdx ($rIdx+1..$#contigIDs) {
        QRY: for my $qIdx ($rIdx+1..$#processedContigs) {
            
            my $qNick     = $processedContigs[$qIdx];
            my $qContigID = $rContigNickLookup{$qNick};
            my $qOrg      = $orgLookupByContig{$qContigID};

            # Eliminate comparisons of contigs from the same organism unless we explicitly want to compare
            next QRY if($qOrg eq $rOrg);
            next QRY if(!$compareSameTaxLevel && isSameTaxLevel($rOrg,$qOrg,$taxLevel));
            my $pDir      = $pDirCodeLookup{ $pDirNickLookup{$qNick} };
            my $qDir      = $pDir."/".$partitionDir."/".$qNick;
#            print "    Loading QRY (qIdx=$qIdx, nick=$qNick) hex strings...";
            my $qkmers0   = new Benchmark;
            my @qKmers    = @{ retrieve $qDir."/".$kmers_eBV_hex_Filename };
            my $qkmers1   = new Benchmark;
#            print "done [".sciExpand($1, $benchRez)."]\n" if(timestr(timediff($qkmers1, $qkmers0)) =~ m/^\s*(\S+)\s+/);
            my $qLoadTime = timestr(timediff($qkmers1, $qkmers0));
               $qLoadTime = sciExpand($1, $benchRez) if($qLoadTime =~ m/^\s*(\S+)\s+/);
               $qLoadTime = $1 if($qLoadTime =~ m/^(\d+\.?\d+)\s+/);
            $totalQryTime += $qLoadTime;
            ++$totalQryLoaded;
               
               
            # Loop through all of ref's buckets, creating bitMasks for each, if suffixes exist in that bucket
            my @qL2BM = ();
#            print "Generating BitMasks for QRY ($qNick)...";
            my $qBM0      = new Benchmark;
            for(0..$#qKmers) {
                #print "***Creating QRY bitmask [$_] of length ".$qKmers[$_]->[0]."***\n";
                $qL2BM[$_] = Bit::Vector->new($qKmers[$_]->[0]) if($qKmers[$_]);                
            }
            my $qBM1      = new Benchmark;
#            print "done [".sciExpand($1, $benchRez)."]\n" if(timestr(timediff($qBM1, $qBM0)) =~ m/^\s*(\S+)\s+/);

            #my $qBM0      = new Benchmark;
            #my @qL2BM     = @{ retrieve $qDir."/kmersLevel2BitMasks.dmp"  };
            #my $qBM1      = new Benchmark;
            #print "done [".sciExpand($1, $benchRez)."]\tProcessing QRY..." if(timestr(timediff($qBM1, $qBM0)) =~ m/^\s*(\S+)\s+/);

            my $qTime0 = new Benchmark;
            
            my $qBuckets = $bucketLookup{$qNick};    # No. of buckets filled
            #my $t0 = new Benchmark;

            # Given two kmer arrays (@rKmers, @qKmers), we select the smaller of the 
            # two (assume @rKmers) and iterate over its existing AREFs as long as one 
            # exists in the larger kmer array (assume @qKmers) as well. If one exists, 
            # we Bit_On the Level1 BitMask of the REF (rLevel1BitMask) and QRY
            # (qLevel1BitMask), and immediately perform a Level2 intersection. REF
            # Level2 intersections are then recorded in their respective Level2 BitMasks
            # in @rKmersLevel2BitMasks, and those in QRY are recorded in their Level2
            # BitMasks in @qKmersLevel2BitMasks2.
            #
            # For INCREMENTAL UPDATE purposes, we require all pairwise intersection
            # data to be recorded. Therefore, the Level1 BitMask for both REF and QRY
            # are recorded:
            #   $rDir/          @rKmers
            #      |__$qDir/    $rLevel1BitMask         // specific to r_vs_q
            #                   @rKmersLevel2BitMasks   // specific to r_vs_q
            #   $qDir/          @qKmers
            #      |__$rDir/    $qLevel1BitMask         // specific to q_vs_r
            #                   @qKmersLevel2BitMasks   // specific to q_vs_r
            #-------------------------------------------------------------------
            # Given 2 arrays, @N and @M, we let: size(N) > size(M)
            #-------------------------------------------------------------------
            #print "============ Starting bucket-based intersection =============\n" if($iDebug);

            # We are given 2 arrays of data, and 2 arrays of output to fill. If
            # the 2 arrays of data are swapped (for efficient computation), we increment
            # our arraySWAP counter. If the buckets are swapped, we bit-flip the bucketSWAP
            # counter. Therefore, if the 2 arrays of data are swapped:

            # Determine which array has the most Level1 buckets
            my $largestArray;
            my $smallestArray;
            my $arraySwap  = Bit::Vector->new(1);       # 1 if @rKmers swapped for @qKmers
            my $bucketSwap = Bit::Vector->new(1);       # 1 if the array buckets have been swapped
            my $netSwap    = Bit::Vector->new(1);       # the result, Xor($arraySwap,$bucketSwap)
            
            if($rBuckets >= $qBuckets) {
                #print "##REF ($rBuckets) has more buckets than QRY ($qBuckets) => NO SWAP\n" if($iDebug);
                $largestArray  = \@rKmers; 
                $smallestArray = \@qKmers;
                $arraySwap->Bit_Off(0);
            }
            else {
                #print "##QRY ($qBuckets) has more buckets than REF ($rBuckets) => SWAP\n" if($iDebug);
                $largestArray  = \@qKmers; 
                $smallestArray = \@rKmers;
                $arraySwap->Bit_On(0);
            }

            # Loop through each bucket and perform Level2 intersections if necessary
            my $bucket = 0;
            
            # or:
            #       set $maxBuckets = smaller of($rBuckets, $qBuckets)
            #           while($bucket != $maxBuckets) ...
            #       and only increase $bucket when exists($smallestArray->[$bucket])
            #
            
            #for (0..$#qL2BM) {
            #    print ">>>>>>>>>>>>>>>>>$_: ".$qL2BM[$_]->Size()." bits\n" if($qL2BM[$_]);
            #}
            
            while ($bucket != $MAX_BUCKETS) {

                # --------------------------------------------------------------
                # 1. We skip empty buckets by iterating through the smallest 
                #    array for efficiency
                # --------------------------------------------------------------
                if(!$smallestArray->[$bucket] || !$largestArray->[$bucket]) {
                    ++$bucket;
                    #print "##Skipping empty bucket...\n" if($iDebug);
                    next;
                }

                # --------------------------------------------------------
                # 2. Check for *any* overlap among corresponding buckets; 
                #    Next bucket if no overlap at all
                #    No overlap if {minM} > {maxN} or {minN} > {maxM}
                # --------------------------------------------------------
                if(
                   #___1stStringInM____   ____lastStringInN_____
                   ($smallestArray->[$bucket][1] gt $largestArray->[$bucket][-1]) 
                   #___1stStringInN____   ____lastStringInM_____
                || ($largestArray->[$bucket][1]  gt $smallestArray->[$bucket][-1])
                  ) {
                    ++$bucket;
                    #print "##Skipping non-overlapping buckets...\n" if($iDebug);
                    next;
                }

                # ------------------------------------------
                # 3. Find out which Level2 array is largest
                # ------------------------------------------
                my $largestBucket  = q{};       # AREF
                my $smallestBucket = q{};       # AREF
                #print "##lastPos(array):  LargestArray->[bucket$bucket][0] = ".$largestArray->[$bucket][0]."\n" if($iDebug);
                #print "##lastPos(array): SmallestArray->[bucket$bucket][0] = ".$smallestArray->[$bucket][0]."\n" if($iDebug);
                if($largestArray->[$bucket][0] > $smallestArray->[$bucket][0]) { 
                    $largestBucket = \@{ $largestArray->[$bucket] }; $smallestBucket = \@{ $smallestArray->[$bucket] };
                    $bucketSwap->Bit_Off(0);
                }
                else { 
                    $largestBucket = \@{ $smallestArray->[$bucket] }; $smallestBucket = \@{ $largestArray->[$bucket] };
                    $bucketSwap->Bit_On(0);
                }
                #print "##********************************************\n" if($iDebug);
                #print "## Largest array's bucket=$bucket has ".(scalar(@{ $largestArray->[$bucket] })-1)." elements\n" if($iDebug);
                #print "##Smallest array's bucket=$bucket has ".(scalar(@{ $smallestArray->[$bucket] })-1)." elements\n" if($iDebug);

                my $destLargestBucket  = q{};               # Either \@qryOutput or \@refOutput
                my $destSmallestBucket = q{};               # Either \@refOutput or \@qryOutput
                
                $netSwap->Xor($arraySwap, $bucketSwap);
                if($netSwap->bit_test(0)) {      # largest writes to QRY
                    #if($iDebug) {
                    #    if($bucketSwap) { print "##Bucket SWAPPED\n"; }
                    #    else            { print "##Bucket *not* SWAPPED!\n"; }
                    #}
                    $destLargestBucket  = \$qL2BM[$bucket];    # REF output
                    $destSmallestBucket = \$rL2BM[$bucket];    # QRY output
                    #print "## destLargestBucket: binary(\"".${ $destLargestBucket }->to_Bin()."\"): ".join(",", @{ $qKmers[$bucket] })."\n" if($iDebug);
                    #print "##destSmallestBucket: binary(\"".${ $destSmallestBucket }->to_Bin()."\"): ".join(",", @{ $rKmers[$bucket] })."\n" if($iDebug);
                }
                else {                                              # largest writes to REF
                    #if($iDebug) {
                    #    if($bucketSwap) { print "##Bucket SWAPPED\n"; }
                    #    else            { print "##Bucket *not* SWAPPED!\n"; }
                    #}
                    $destLargestBucket  = \$rL2BM[$bucket];    # QRY output
                    $destSmallestBucket = \$qL2BM[$bucket];    # REF output
                    #print "## destLargestBucket = binary(\"".${ $destLargestBucket }->to_Bin()."\"): ".join(",", @{ $rKmers[$bucket] })."\n" if($iDebug);
                    #print "##destSmallestBucket = binary(\"".${ $destSmallestBucket }->to_Bin()."\"): ".join(",", @{ $qKmers[$bucket] })."\n" if($iDebug);
                }

                #print "##Press <ENTER>..."; <STDIN>;

                # ----------------------------------------
                # 4. Initialize limits for LARGEST bucket
                # ----------------------------------------
                #    LARGEST
                my $firstN = 1;                         # we take pos1 because pos0 is the array size
                my $lastN  = $largestBucket->[0]+1;     # position is 1+ lastPos in bitMask

                #print " arraySWAP=".$arraySwap->bit_test(0)."\n";
                #print "bucketSWAP=".$bucketSwap->bit_test(0)."\n";

                # -----------------------------------------
                # 5. Initialize limits for SMALLEST bucket
                # -----------------------------------------
                #    SMALLEST
                #TRYING TO DETERMINE LIMITS OF M BY FINDING FIRSTPOS(N) AND LASTPOS(N) IN M[]:
                #-firstM: if FIRSTPOS(N) > LASTPOS(M), then dont proceed further => no overlap
                #         otherwise, return position (in M) of next closest value (in M) to FIRSTPOS(N).
                #-lastM:  if LASTPOS(N) < FIRSTPOS(M), then dont proceed further => no overlap
                #         otherwise, return position (in M) of next closest value (in M) to LASTPOS(N).
                # Find limits of overlap with smallest array using {imin} and {imax} of the LARGEST array
                my $firstM = 1;
                my $lastM  = 1;
                if( ($largestBucket->[$firstN] gt $smallestBucket->[-1]) || ($largestBucket->[-1] lt $smallestBucket->[1]) ) { ++$bucket; next }
                else {
                    $firstM = binary_search_overlap($largestBucket->[1],  $smallestBucket, 1, $smallestBucket->[0]-1);
                    # Sub binary_search_overlap() will never return more than the position of the last element in the array.
                    # For the zipper intersection algorithm, we need lastM to be 1 + the last elemetn in the array.
                    $lastM  = binary_search_overlap($largestBucket->[-1], $smallestBucket, 1, $smallestBucket->[0]-1)+1;
                }

                #print "##FirstM = $firstM\n" if($iDebug);
                #print "##LastM  = $lastM\n"  if($iDebug);
                #print "##FirstN = $firstN\n" if($iDebug);
                #print "##LastN  = $lastN\n"  if($iDebug);
                ##print "Press <ENTER>..."; <STDIN>;

                # firstN is the position in the STRING  =>  position in BitMask is firstN-1
                # firstM is the position in the STRING  =>  position in BitMask is firstM-1
                #print "##bucket$bucket: SmallestBucket: ".join(",", @{ $smallestBucket })."\n" if($iDebug);
                #print "##bucket$bucket:  LargestBucket: ".join(",", @{ $largestBucket })."\n"  if($iDebug);

                # Note: firstM/lastM/firstN/lastN correspond to the positions in the ARRAYS, which is
                #       one-off from the corresponding position in the bitMask, 
                #           i.e.    BMpos = ARRAYpos-1;
                my @tmpResultsM = ();
                my @tmpResultsN = ();
                while ($firstM != $lastM && $firstN != $lastN) {
                    if   ($smallestBucket->[$firstM] lt $largestBucket->[$firstN])  { ++$firstM; }
                    elsif($largestBucket->[$firstN]  lt $smallestBucket->[$firstM]) { ++$firstN; }
                    else {
                    
                        #print "## lastPos(destLargestBucket: $bucket)=".(${ $destLargestBucket }->Size()-1)."=?$lastN\tAttempting to Bit_On(".($firstN-1)."): ".$smallestBucket->[$firstM]." = ".$largestBucket->[$firstN]."\n" if($iDebug);
                        #print "##lastPos(destSmallestBucket: $bucket)=".(${ $destSmallestBucket }->Size()-1)."=?$lastM\tAttempting to Bit_On(".($firstM-1).")\n" if($iDebug);
                    
                        #$destLargestBucket->Bit_On($firstN++);
                        #$destSmallestBucket->Bit_On($firstM++);
                        ${ $destLargestBucket }->Bit_On($firstN-1); ++$firstN;
                        ${ $destSmallestBucket }->Bit_On($firstM-1); ++$firstM;
                    }
                    #print "##next element...\n" if($iDebug);
                }
                #print "##Exiting...\n" if($iDebug);
                #print "## destLargestBucket = [".${ $destLargestBucket }->Size()."]binary(\"".${ $destLargestBucket }->to_Bin()."\")\n" if($iDebug);
                #print "                       (".join(",", @tmpResultsN).")\n" if($iDebug);
                #print "##destSmallestBucket = [".${ $destSmallestBucket }->Size()."]binary(\"".${ $destSmallestBucket }->to_Bin()."\")\n" if($iDebug);
                #print "                       (".join(",", @tmpResultsM).")\n" if($iDebug);
                ++$bucket;
            }
            #-------------------------------------------------------------------

            #my $t1 = new Benchmark;
            #print "And() + Ranging TIME: ".timestr(timediff($t1, $t0)).": ".scalar(@rIndices)." = ".scalar(@qIndices)." elements in PrefixBitMask\n";

            # -------------------
            # ACCUMULATE RESULTS
            # -------------------
            # @qL2BM contains the results of the REF_vs_QRY intersections
            # @rL2BM contains the results of the REF_vs_QRY intersections
            $localResults{$rNick}->{$qNick} = \@rL2BM;
            $localResults{$qNick}->{$rNick} = \@qL2BM;

            #print "ON bitMask positions: ".join(",",@rIndices)."\n";            
            #print "qIndices: ".join(",",@qIndices)."\n";
            #<STDIN>;
            
            my $qTime1 = new Benchmark;
#            print "done [".sciExpand($1, $benchRez)."]\n" if(timestr(timediff($qTime1, $qTime0)) =~ m/^\s*(\S+)\s+/);

        } #QRY

        my $aveQryLoadTime = ($totalQryLoaded > 0) ? sprintf("%.2f", $totalQryTime/$totalQryLoaded) : (0);
       
        my $refTime1 = new Benchmark;
        my $refTimeText = timestr(timediff($refTime1, $refTime0));
           $refTimeText = sciExpand($1, $benchRez) if($refTimeText =~ m/^\s*(\S+)\s+/);
           $refTimeText = $1 if($refTimeText =~ m/^(\d+\.?\d+)\s+/);

        $refOutputText .= "//$rContigID////$rOrg// $numKmers k-mers [rLD=$rLoadTime|qLD=$aveQryLoadTime|QRY=$totalQryLoaded] [$refTimeText] wallsecs\n                          ";
        #print "\tQRY TIME: ".sciExpand($1, $benchRez)."\n" if(timestr(timediff($refTime1, $refTime0)) =~ m/^\s*(\S+)\s+/);
        

    } #REF

    # We cannot thread-share bitVectors, so we dump groups of them to disk and, 
    # when done, we load them ALL back into the main thread, and proceed to the
    # SUBTRACTION phase.
    #
    # -----
    # DUMP
    # -----
    my $store0 = new Benchmark;
    while(my($rNick, $href) = each %localResults) {
        #store $href->{$_}, $outdir."/".$partitionDir."/".$rNick."/L2BM_$_.dmp" foreach (keys %{ $href });
        store $href->{$_}, $pDirCodeLookup{ $pDirNickLookup{$rNick} }."/".$partitionDir."/".$rNick."/L2BM_$_.dmp" foreach (keys %{ $href });
    }
    #foreach my $rNick (%localResults) {
    #    store $localResults{$rNick}->{$_}, $outdir."/".$partitionDir."/".$rNick."/L2BM_$_.dmp" foreach (keys %{ $localResults{$rNick} });
    #}
    my $store1 = new Benchmark;
    my $storeTimeText = timestr(timediff($store1, $store0));
       $storeTimeText = sciExpand($1, $benchRez) if($storeTimeText =~ m/^\s*(\S+)\s+/);
       $storeTimeText = $1 if($storeTimeText =~ m/^(\d+\.?\d+)\s+/);
       
    my $threadTime1 = new Benchmark;
    my $threadTimeText = timestr(timediff($threadTime1, $threadTime0));
       $threadTimeText = sciExpand($1, $benchRez) if($threadTimeText =~ m/^\s*(\S+)\s+/);
       $threadTimeText = $1 if($threadTimeText =~ m/^(\d+\.?\d+)\s+/);

    my $currCounter = 0;
    {
        lock($intersectCounter);
        $intersectCounter++;
        $currCounter = $intersectCounter;
    }
       
    stat_log("Thread $currCounter (of ".$_[1]."): INTERSECTION BitMasks stored in [$storeTimeText] wallsecs; [$threadTimeText] wallsecs\n"
            ."                          ".$refOutputText."");
#            ."Stored INTERSECTION BitMasks to disk in [".$storeTimeText."] wallsecs");

    
#    {
#        lock(%intersectionResults);
#        foreach my $rIdx (keys %localResults) {
#            $intersectionResults{$rIdx} = $
#        }
#    }
    # ...or...
    #
    # foreach my $idx (keys %localResults) {
    #     store $href->{$_}, $outdir."/".$partitionDir."/".$idx."/L2BM_$_.dmp" foreach (keys %{ $localResults{$idx} });
    # }
    #
    
    return;

}
################################################################################
# ARGS: $kmerIntersectionOpts
# ------------------------------------------------------------------------------
sub intersectKmersParallel {

    stat_log("----------------------------------");
    stat_log("L-mer Intersection");

    my $time0 = new Benchmark;                     # Sub() Benchmark

    # %bucketLookup = ($contigNick => $numBuckets) is populated in sub prefixKmers
    # so numContigs also equal scalar(@processedContigs)
    my $numContigs = scalar(@processedContigs);
    #my $loadBalancedContigs_aref = loadBalanceLinear(\@contigIDs, $nThreads);
    
    my $loadBalancedIndices = loadBalanceLinearPaired([0..$#processedContigs], $nThreads);
    my $numJobs = scalar(@{ $loadBalancedIndices });
    #print "LoadBalIndices = ".join(",", @{ $loadBalancedIndices })."\n";
    
    my @jobs = ();
    foreach my $arefOfIndices (@{ $loadBalancedIndices }) {
        #print "arefOfIndices=".join(",", @{ $arefOfIndices })."\n";
        #my $jobID = threads->new(\&intersectKmers, $arefOfIndices, $_[0]);
        my $jobID = threads->new(\&intersectKmers, $arefOfIndices, $numJobs);
        push(@jobs, $jobID);
    }

    $_->join() foreach(@jobs);

    my $time2 = new Benchmark;
    my $intersectionTimeText = timestr(timediff($time2,$time0));
       $intersectionTimeText = $1 if($intersectionTimeText =~ m/^\s*(\d+)\s+/);
       $intersectionTimeText = $1 if($intersectionTimeText =~ m/(\d+\.?\d+)\s+/);
    stat_log("TOTAL INTERSECTION TIME: [$intersectionTimeText] wallsecs");
    displayWarnings();
    
    return;
}
########################################################################################
#                                                                                                                
#     ####  #   #  ####   #####  ####     #      ###  #####  #####    ###    ##   #                                                                                              
#    #      #   #  #   #    #    #   #   # #    #       #      #     #   #   # #  #                                      
#     ###   #   #  ####     #    ####   #####  #        #      #    #     #  #  # #                                                 
#        #  #   #  #   #    #    #  #   #   #   #       #      #     #   #   #   ##                                      
#    ####    ###   ####     #    #   #  #   #    ###    #    #####    ###    #    #                                                       
#
########################################################################################
########################################################################################
#                                                                                                                
#     ####  #   #  ####   #####  ####     #      ###  #####  #####    ###    ##   #                                                                                              
#    #      #   #  #   #    #    #   #   # #    #       #      #     #   #   # #  #                                      
#     ###   #   #  ####     #    ####   #####  #        #      #    #     #  #  # #                                                 
#        #  #   #  #   #    #    #  #   #   #   #       #      #     #   #   #   ##                                      
#    ####    ###   ####     #    #   #  #   #    ###    #    #####    ###    #    #                                                       
#
########################################################################################
########################################################################################
#                                                                                                                
#     ####  #   #  ####   #####  ####     #      ###  #####  #####    ###    ##   #                                                                                              
#    #      #   #  #   #    #    #   #   # #    #       #      #     #   #   # #  #                                      
#     ###   #   #  ####     #    ####   #####  #        #      #    #     #  #  # #                                                 
#        #  #   #  #   #    #    #  #   #   #   #       #      #     #   #   #   ##                                      
#    ####    ###   ####     #    #   #  #   #    ###    #    #####    ###    #    #                                                       
#
########################################################################################
################################################################################
# ARGS: \@arefOfIndices
# ------------------------------------------------------------------------------
sub subtractKmers {
    
    # Load all of REF's @L2BM's (or just the orgs that require intersection)
    my $threadTime0 = new Benchmark;

    my $prefix            = "GOTTCHA";
    my $vs                = $compareSameTaxLevel ? "1" : "0";
    my $taxLevelDetail    = "$taxLevel$vs";
    my $contigText        = "c$totalContigs";
    my $lmerText          = "k$k";
    my $uniqueText        = "u$minUniqueFragLength";
    my $partitionDetail   = "n$nThreads.p$prefixLength";
    my $dbLabelDetail     = $contigText."."
                           .$lmerText."."
                           .$uniqueText."."
                           .$partitionDetail;
                           
    my $refDbLabelDetail  = $contigText."."
                           .$lmerText."."
                           ."u1."
                           .$partitionDetail;

    my $referenceFilename = $contigText."."
                           .$lmerText."."
                           ."u1."                       # <--- for profiling later
                           .$partitionDetail;

    my $ext             = "fna";      # FASTA file extension 
    
    my %localStats = ();

    # Loop through this thread's indices
    # Loop through set of indicies for @processedContigs
    REF: for my $rIdx (@{ $_[0] }) {
        my $rNick = $processedContigs[$rIdx];
        my $pDir = $pDirCodeLookup{ $pDirNickLookup{$rNick} };
        my $rDir = $pDir."/".$partitionDir."/".$rNick;
        my $rContigID = $rContigNickLookup{$rNick};
        my $rFilename = $rFileNickLookup{$rNick};
        my $rOrg      = $orgLookupByContig{$rContigID};
        my $orgDir    = $organisms{$rOrg}->{DIR}->{dir};
        my $totalUniqueLen = 0;

        ## DEBUG
        #print "******** Processed Contigs ($rIdx): ".join(",", @processedContigs)." ********\n";
        
        $uniqueFragsByOrg{$rOrg} = &share({}) unless (exists $uniqueFragsByOrg{$rOrg});
        $uniqueFragsByOrg{$rOrg}->{$rContigID} = &share({}) unless (exists $uniqueFragsByOrg{$rOrg}->{$rContigID});

        my $outFilename1 = "$outdir/$fastadir/$taxLevel/$prefix.$taxLevelDetail.$dbLabel.$dbLabelDetail.$orgDir.$ext";      # uX
        my $outFilename2 = "$outdir/$fastadir/$taxLevel/$prefix.$taxLevelDetail.$dbLabel.$refDbLabelDetail.$orgDir.$ext";   # u1

        # Needed for statistics
        my %uniqueFrags : shared = ();
        
        my @uL2BMs = ();        # Stores the unions of each bitmask into each bucket

        # BEGIN_LOCAL_BLOCK
        {
            my @L2BMs = ();

            QRY: foreach my $qIdx (0..$#processedContigs) {
                ## DEBUG
                #print "----------------> query $qIdx (CSTL: $compareSameTaxLevel)\n";
                next if($qIdx == $rIdx);                    # Skip identical contigs
                my $qNick     = $processedContigs[$qIdx];
                my $qContigID = $rContigNickLookup{$qNick};     #************************* FIX: 2012-10-15 ******************
                my $qOrg      = $orgLookupByContig{$qContigID}; #************************* FIX: 2012-10-15 ******************

                next QRY if(!$compareSameTaxLevel && isSameTaxLevel($rOrg, $qOrg, $taxLevel));
                #if(!$compareSameTaxLevel && isSameTaxLevel($rOrg, $qOrg, $taxLevel)) {
                #    print "$rOrg ($rContigID) is SKIPPING $qOrg ($qContigID) @ $taxLevel...\n";
                #    next QRY;
                #}
                #else {
                #    print "$rOrg ($rContigID) is SUBTRACTING $qOrg ($qContigID) @ $taxLevel...\n";
                #}

               
                my $filename = $rDir."/L2BM_$qNick.dmp";    # Load REF_vs_QRY file from *REF* directory
                next if(!-e $filename);
                #if(!-e $filename) {                             #************************* FIX: 2012-10-15 ******************
                #    print "**** WARNING ****: $filename is NOT found...!\n";
                #    next QRY;
                #}
                my $l2bm_aref = retrieve $filename;
                push(@L2BMs, $l2bm_aref);                   # add an AREF of BitMasks                
            }

            # Unionize the bitmasks (for identical buckets) across all @L2BMs
            #
            # NOTE: If the first bucket is empty, we will have problems! We need to protect against this...
            #       Either create an empty bitmask of the correct length and compare $_+1 against an exists() function
            #       or some other means...
            #
            foreach my $bucketIdx (0..$#{ $L2BMs[0]}) {
                my $uBM = q{};      # Unionized BitMask across all equivalent buckets
                my $idx = 0;        # Location of first extant bucket
                
                # Find first extant bucket among all ordered arrays
                for(0..$#L2BMs) {
                    if($L2BMs[$_]->[$bucketIdx]) {
                        $uBM = $L2BMs[$_]->[$bucketIdx]->Clone();
                        $idx = $_;
                        last;
                    }
                }
                
                # If no bucket has a representative BitMask, then neither will the UNION
                next unless $uBM;
                
                for($idx+1..$#L2BMs) {
                    next if(!$L2BMs[$_]->[$bucketIdx]);
                    $uBM->Union($uBM, $L2BMs[$_]->[$bucketIdx]);
                }                
                $uL2BMs[$bucketIdx] = $uBM;     # Store
            } #bucketIdx
        } #END_LOCAL_BLOCK

        # Load $seqBM
        my $seqBM = retrieve $rDir."/".$seqBMfilename;

        # Load @kmer_starts
        my @kmer_starts = @{ retrieve $rDir."/".$kmer_starts_Filename };

        # Range the ON bits foreach bucket
        foreach my $bucketIdx (0..$#uL2BMs) {
        
            next if(!$uL2BMs[$bucketIdx]);
        
            # Returns hash of (min,max) values, where min = first ON bit, and
            # max = first OFF bit - 1.
            # %range specifies which indices in @kmer_starts to use for turning
            # OFF bits in the $seqBM.
            my %range = %{ rangeIt($uL2BMs[$bucketIdx], $uL2BMs[$bucketIdx]->Size()) };
            next if(!%range);
           
            # Correlate each ON bit in the BitMask with its corresponding start position in @start_kmers.
            # Use Interval_Empty(min,max) on $seqBM for each start position to turn
            # corresponding bits OFF, from (startPos0 + $k) to (startPosN + $k)
            #my $lastPos      = $uL2BMs[$bucketIdx]->Size()-1;
            my $wrapConst    = $k - 1;
            my $seqBMlastPos = $seqBM->Size()-1;
            my $boundary     = $seqBMlastPos - $wrapConst;
            #my $extension    = $k - ($seqBMlastPos - $start + 1) = $k - 1 - $seqBMlastPos + $start
            #                 = $wrapConst - $seqBMlastPos + $start
            #                 = -(-$wrapConst + $seqBMlastPos) + $start
            #                 = -($boundary) + $start

            foreach my $min (usort keys %range) {
                my $max = $range{$min};
            }
            
            foreach my $min (usort keys %range) {
                my $max = $range{$min};
                for($min..$max) {
                    foreach my $startPos (@{ ${ $kmer_starts[$bucketIdx] }[$_] }) {
                        if($startPos <= $boundary) {
                            $seqBM->Interval_Empty($startPos, $startPos + $wrapConst);        # <<<<<<<<<<<<<<<<<<< CHECK THIS
                        }
                        else {
                            $seqBM->Interval_Empty($startPos, $seqBMlastPos);
                            $seqBM->Interval_Empty(0, $startPos - $boundary - 1);  
                        }
                    }
                } #for
                
            } #foreach

        } #bucketIdx

        # Save resulting $seqBM
        store $seqBM, $rDir."/".$uniqueSeqBMfilename;

        # Load $seq
        my $seqStorableFilename = $organisms{$rOrg}->{REPL}->{$rFilename}->{CONTIG}->{$rContigID}->{SEQ};
        my $seq = ${ retrieve $seqStorableFilename };
        my $desc = $organisms{$rOrg}->{REPL}->{$rFilename}->{CONTIG}->{$rContigID}->{DESC};
    
        # Use $seqBM to extract unique strings from $seq of a minimum size
        my %uniqueRange = %{ rangeIt($seqBM, $seqBM->Size()) };
               
        # Pre-existing FASTA files will be appended to. Necessary so that each replicon can write to 
        # the appropriate file, but between runs of this program, this can cause a problem if the 
        # entire directory is not cleared first.
        open my $fh1, '>>', $outFilename1 || die "Cannot open filehandle for \"$outFilename1\" - $!\n";
        open my $fh2, '>>', $outFilename2 || die "Cannot open filehandle for \"$outFilename2\" - $!\n";
        fileLock($fh1);
        fileLock($fh2);
        while(my($start, $stop) = each %uniqueRange) {
            my $fragLen = $stop-$start+1;
            my $frag    = substr($seq, $start, $fragLen);
            my $header  = $rContigNickLookup{$rNick} #$contigIDs[$rIdx]
                        .pad_zeroes($start+1,9)
                        ."|".pad_zeroes($start+$fragLen,9)      # (start) + fraglen - 1
                        ."|".pad_zeroes($fragLen,9);

            # Print entire fragment, regardless of length
            print $fh2 ">".$header."| ".$desc."\n$frag\n";
            $totalUniqueLen += $fragLen;

            next unless ($fragLen >= $minUniqueFragLength);
            # Print u=minUniqueFragLength fragments
            print $fh1 ">".$header."| ".$desc."\n$frag\n";

            # UPDATE LOCAL STATS
            $localStats{$rOrg}->{$rFilename} += $totalUniqueLen;

            $uniqueFrags{$header} = $frag;
        }
        fileUnlock($fh1);
        fileUnlock($fh2);
        close($fh1);
        close($fh2);
        
        # ---------------------------------------------------------------------------------------------
        # We cannot close these filehandles within the REF loop, because init_file_write() initializes
        # an overwrite of the file. If we change init_file_write() to append, then we should be okay.
        # ---------------------------------------------------------------------------------------------
        #close $output_fh_hash{$orgDir}->{FH_u1};
        #close $output_fh_hash{$orgDir}->{FH_uX};

        # Update this globally shared hash to use for statistics
        {
            lock(%uniqueFragsByOrg);
            $uniqueFragsByOrg{$rOrg}->{$rContigID} = \%uniqueFrags;
        }
        
        # Store the @uL2BMs for future incremental update
        # UPDATE: 2012-10-09
        #store \@uL2BMs, $rDir."/".$unionizedL2BMsFilename;
        my $destDir = "$outdir/$updateDir/$storableDir/$rNick";
        createDir($destDir);
        store \@uL2BMs, $destDir."/"."uL2BMs_$taxLevel$compareSameTaxLevel\_$rNick.$dmpExt";
        
    } #REF

    # UPDATE %organisms WITH {UNIQUE}, the total unique length
    {
        lock(%organisms);
        while(my ($org, undef) = each %localStats) {
            while( my ($replFilename, undef) = each %{ $localStats{$org} } ) {
                $organisms{$org}->{REPL}->{$replFilename}->{UNIQUE} = $localStats{$org}->{$replFilename};
            }
        }        
        
    }
        
    my $currCounter = 0;
    {
        lock($subtractionCount);
        $subtractionCount++;
        $currCounter = $subtractionCount;
    }

    my $threadTime1 = new Benchmark;
    my $threadTimeText = timestr(timediff($threadTime1, $threadTime0));
       $threadTimeText = sciExpand($1, $benchRez) if($threadTimeText =~ m/^\s*(\S+)\s+/);
       $threadTimeText = $1 if($threadTimeText =~ m/^(\d+\.?\d+)\s+/);
       
    stat_log("Thread $currCounter (of ".$_[1]."): [$threadTimeText] wallsecs");
    
    return;
    
}
################################################################################
# ARGS: none
# ------------------------------------------------------------------------------
sub subtractKmersParallel {
    
    my $subtractTime0 = new Benchmark;
    #my $numContigs = scalar(keys %bucketLookup);
    my $numContigs = scalar(@processedContigs);

    stat_log("----------------------------------");
    stat_log("Processing unique regions for all $totalContigs contigs with $nThreads threads "
            ."(k=$k, u=$minUniqueFragLength)...");

    for (0..$#processedContigs) {
        my $rOrg      = $orgLookupByContig{ $rContigNickLookup{ $processedContigs[$_] } };
        my $orgDir    = $organisms{$rOrg}->{DIR}->{dir};
        $output_fh_hash{$orgDir} = &share({});
        #$output_fh_hash{$orgDir}->{FH_uX}       = &share({});
        #$output_fh_hash{$orgDir}->{FILENAME_uX} = &share({});
        #$output_fh_hash{$orgDir}->{FH_u1}       = &share({});
        #$output_fh_hash{$orgDir}->{FILENAME_u1} = &share({});
    }

    # Balance across the total no. of contigIDs (before ripping)
    my $loadBalancedIndices = loadBalanceLinearPaired([0..$#processedContigs], $nThreads);
    my $numJobs = scalar(@{ $loadBalancedIndices});
    
    my @jobs = ();
    foreach my $arefOfIndices (@{ $loadBalancedIndices }) {
        my $jobID = threads->new(\&subtractKmers, $arefOfIndices, $numJobs);
        push(@jobs, $jobID);
    }

    $_->join() foreach(@jobs);

    # End benchmark
    my $subtractTime1 = new Benchmark;
    my $intersectionTimeText = timestr(timediff($subtractTime1, $subtractTime0));
       $intersectionTimeText = $1 if($intersectionTimeText =~ m/^\s*(\d+)\s+/);
       $intersectionTimeText = $1 if($intersectionTimeText =~ m/(\d+\.?\d+)\s+/);
    stat_log("TOTAL SUBTRACTION TIME: [$intersectionTimeText] wallsecs");

    return;
}
################################################################################
################################################################################
# OTHER SUBS
################################################################################
################################################################################
################################################################################
#   1. A single file indicating stats for each replicon's unique fragments (count, 
#      length, stdev, histogram data)
#      + for the entire organism
#      + for the entire Genus
#      + for the entire dataset
#   2. A single histogram data file detailing the stats for each replicon's sequence 
#      (histogram data: seq1 => length)
#
        # ---- SUMMARY of STATS --------------------------------------------------------
        # ---> l-mer length = #
        # ---> total number of organisms represented (mean # of replicons per org +/- stdev)
        # ---> number of replicons representing (# of fragments, mean length in bp +/- stdev)
        # ---> total # of fragments (mean length in bp +/- stdev)
        # ---> total length in bp
        #
        # ---- DETAILS of STATS --------------------------------------------------------
        # ---> l-mer length = #
        # ---> Organism: Staphylococcus aureus (total # of replicons)
        #      "ORIG_LEN" "RAW_UNI_LEN" "FILTERED_UNI_LEN" "Raw%ofTot" "Filt%ofTot" "%Filtered:Raw"
        # ---> ----"CP001844": 18759 9348 8311 49.83 44.30 88.91
#
# DB Features:
#   1. Breakdown of unique sequences per replicon 
#   2. Breakdown of unique sequences per organism
#   3. Breakdown of unique sequences per taxon (species, genus, ...)
#
#   ARGS: \%uniqueFragsByOrg
#-------------------------------------------------------------------------------
sub storeUniqueFragStats {

    STDOUT->autoflush(1);

    stat_log("----------------------------------");
    stat_log("Database statistics");
    my $time0 = new Benchmark;

    # ------------------------------------------
    # BEGIN: THREADED PORTION
    # ------------------------------------------
    # Process core stats across multiple threads
    # Using global: @orgNames, but check exists $uniqueFragsByOrg{$org} in case 
    #               an org didn't make it out of the XML file...
    my $numOrgs = scalar(@orgNames);
    my $loadBalIndices = loadBalanceStatic($numOrgs, $nThreads);
    my @jobs = ();
    foreach my $arefOfIndices (@{ $loadBalIndices }) {
        my $tid = threads->new(\&_storeUniqueFragStats, $arefOfIndices, $numOrgs);
        push(@jobs, $tid);
    }
    
    # process all stats and updates %organisms just before returning
    my $statIdx   = 0;
    my $totalJobs = scalar(@jobs);
    foreach my $tid (@jobs) {
        my $returnText = $tid->join();
        $statIdx++;
        stat_log("Thread $statIdx (of $totalJobs):  $returnText");
    }
    # ------------------------------------------
    # END: THREADED PORTION
    # ------------------------------------------

    # End benchmark
    my $time1 = new Benchmark;
    my $dbText= timestr(timediff($time1, $time0));
    $dbText = $1 if($dbText =~ m/^\s*(\d+)\s+/);
    $dbText = $1 if($dbText =~ m/(\d+\.?\d+)\s+/);
    stat_log("TOTAL STATISTICS CALCULATION TIME: [$dbText] wallsecs");

    # Print to file
    printStatsToFile();

    return;
    
}
################################################################################
# Function: calculates stats of select organisms, updating a local hash until the
#           end of the sub, at which point it will update the global %organisms
#           just before returning.
# ARGS: \@arefOfIndices to @orgNames
#-------------------------------------------------------------------------------
sub _storeUniqueFragStats {

    my %localStats = ();
    STDOUT->autoflush(1);
    my $iter0 = new Benchmark;
        
    my @orgs = @orgNames[@{ $_[0] }];       # Get subset of organism names provided
    ORG: foreach my $org (@orgs) {
    
        my ($orgLenCount,  $orgLenTot,  $orgLenMin,  $orgLenMax) = (0,0,0,0,0);
        my ($orgGCCount,   $orgGCTot,   $orgGCMin,   $orgGCMax)  = (0,0,0,0,0);
        my ($orguLenCount, $orguLenTot, $orguLenMin, $orguLenMax)= (0,0,0,0,0);
        my ($orguGCCount,  $orguGCTot,  $orguGCMin,  $orguGCMax) = (0,0,0,0,0);
        my ($contigCount,$totGenomeLen) = (0,0);

        #$uniqueFragsByOrg{$org}->{$contigID}->{$header} = $fragment;

        next ORG unless (exists $uniqueFragsByOrg{$org});
        foreach my $contigID (keys %{ $uniqueFragsByOrg{$org} } ) {
                                 #                 _________contigNick________
            my $cFilename = $rFileNickLookup{ $contigNickLookup{$contigID} };
            $contigCount++;          # to help track MIN/MAX values
            $totGenomeLen += $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{LENGTH};
            
            my %uany_len  = ();
            my %uany_gc   = ();
            my %uuser_len = ();
            my %uuser_gc  = ();

            my @headers = keys %{ $uniqueFragsByOrg{$org}->{$contigID} };
            foreach (@headers) {
                
                # contigID:   gi|XXXXX|xxx|XXXXXXX|
                # header:     start|stop|length|
                # longHeader: contigID|start|stop|length| description
                my $longHeader = $contigID.$_." ".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{DESC};
                my $seq = $uniqueFragsByOrg{$org}->{$contigID}->{$_};
               
                # *** CURRENTLY ALL REPORTED STATS INCLUDE SEQUENCES WITH ERROR BASES ***
                # Collect stats on ALL unique fragments
                my $seqLen = length($seq);
                $uany_len{$seqLen}++;
                $uany_gc{ pctGC($seq) }++;
                # Collect stats on unique fragments, where length >= u
                if($seqLen >= $minUniqueFragLength) {
                    $uuser_len{$seqLen}++;
                    $uuser_gc{ pctGC($seq) }++;
                }
            } # foreach(@headers)

        $localStats{$org}->{STATS}->{LGENOME}  = $totGenomeLen ? $totGenomeLen : 0;
        $localStats{$org}->{STATS}->{LREMOVED} = sprintf("%.2f", ($totGenomeLen - $orgLenTot));             # <---- here's the problem: $orgLenTot is never updated, and always zero
        my $orgLenPctRemoved = ($totGenomeLen == 0) ? (0) : ((1 - $orgLenTot/$totGenomeLen)*100);
        $localStats{$org}->{STATS}->{LPCTREMOVED} = sprintf("%.2f", $orgLenPctRemoved);

            # PROCESS UANY_GC
            {
                my @unfilteredGCs = ();
                while (my ($gc, $freq) = each %uany_gc) {
                    push(@unfilteredGCs, $gc) for(1..$freq);
                }
                my ($meanGC, $stdevGC, $minGC, $maxGC, $sumGC) = doStats(\@unfilteredGCs);
                $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCMEAN} = $meanGC  ? sprintf("%.2f", $meanGC)  : 0;
                $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCSD}   = $stdevGC ? sprintf("%.2f", $stdevGC) : 0;
                $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCMIN}  = $minGC   ? sprintf("%.2f", $minGC)   : 0;
                $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCMAX}  = $maxGC   ? sprintf("%.2f", $maxGC)   : 0;
                $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCCOUNT}= (@unfilteredGCs) ? (scalar(@unfilteredGCs)) : (0);
                $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCTOT}  = $sumGC   ? sprintf("%.2f", $sumGC)   : 0;
            }

            # PROCESS UANY_LEN
            {
                my @unfilteredLengths = ();
                while (my ($len, $freq) = each %uany_len) {
                    push(@unfilteredLengths, $len) for(1..$freq);
                }
                my ($meanLen, $stdevLen, $minLen, $maxLen, $sumLen) = doStats(\@unfilteredLengths);
                $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LMEAN}  = $meanLen  ? sprintf("%.2f", $meanLen) : 0;
                $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LSD}    = $stdevLen ? sprintf("%.2f", $stdevLen): 0;
                $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LMIN}   = $minLen   ? $minLen                   : 0;
                $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LMAX}   = $maxLen   ? $maxLen                   : 0;
                $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LTOT}   = $sumLen   ? $sumLen                   : 0;
                $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LCOUNT} = (@unfilteredLengths) ? (scalar(@unfilteredLengths)) : (0);    # redundant; same value as ->{LENcount}
            }
            
            # PROCESS DELETIONS
            {
                my $contigLen  = $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{LENGTH}
                               ? $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{LENGTH}
                               : 0;
                my $lenRemoved = $contigLen - $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LTOT};
                my $pctRemoved = ($contigLen == 0) ? (0) : (($lenRemoved/$contigLen)*100);
                $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LREMOVED}    = $lenRemoved ? sprintf("%.2f", $lenRemoved) : 0;
                $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LPCTREMOVED} = $pctRemoved ? sprintf("%.2f", $pctRemoved) : 0;
            }

            # PROCESS UUSER_GC
            {
                my @filteredGCs = ();
                while (my ($gc, $freq) = each %uuser_gc) {
                    push(@filteredGCs, $gc) for(1..$freq);
                }
                my $uTotalGC  = 0;
                   $uTotalGC += $_ foreach (@filteredGCs);
                my ($umeanGC, $ustdevGC, $uminGC, $umaxGC, $usumGC) = doStats(\@filteredGCs);
                $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCMEAN} = $umeanGC  ? sprintf("%.2f", $umeanGC) : 0;
                $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCSD}   = $ustdevGC ? sprintf("%.2f", $ustdevGC): 0;
                $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCMIN}  = $uminGC   ? sprintf("%.2f", $uminGC)  : 0;
                $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCMAX}  = $umaxGC   ? sprintf("%.2f", $umaxGC)  : 0;
                $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCCOUNT}= (@filteredGCs) ? scalar(@filteredGCs) : 0;
                $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCTOT}  = $usumGC   ? sprintf("%.2f", $usumGC)  : 0;
            }

            # PROCESS UUSER_LEN        
            {
                my @filteredLengths = ();
                while (my ($len, $freq) = each %uuser_len) {
                    push(@filteredLengths, $len) for(1..$freq);
                }
                my ($umeanLen, $ustdevLen, $uminLen, $umaxLen, $usumLen) = doStats(\@filteredLengths);
                $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLMEAN} = $umeanLen  ? sprintf("%.2f", $umeanLen)  : 0;
                $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLSD}   = $ustdevLen ? sprintf("%.2f", $ustdevLen) : 0;
                $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLMIN}  = $uminLen   ? $uminLen : 0;
                $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLMAX}  = $umaxLen   ? $umaxLen : 0;
                $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLTOT}  = $usumLen   ? $usumLen : 0;
                $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLCOUNT}= (@filteredLengths) ? scalar(@filteredLengths) : 0;  # redundant; same value as ->{uGCcount}
                
            }
            ###################################################
            # Incremental addition of replicon data to organism
            ###################################################
            if($contigCount == 1) {
                $orgLenMin = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LMIN};
                $orgLenMax = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LMAX};
                $orgGCMin  = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCMIN};
                $orgGCMax  = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCMAX};
                # ----
                $orguLenMin = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLMIN};
                $orguLenMax = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLMAX};
                $orguGCMin  = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCMIN};
                $orguGCMax  = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCMAX};
            }
            # ---- all fragments
            # (LENGTH)
            $orgLenCount += $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LCOUNT};
            $orgLenTot   += $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LTOT};
            $orgLenMin    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LMIN}
                if($localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LMIN} < $orgLenMin);
            $orgLenMax    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LMAX}
                if($localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LMAX} > $orgLenMax);
            # (%GC)
            $orgGCCount += $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCCOUNT};
            $orgGCTot   += $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCTOT};
            $orgGCMin    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCMIN}
                if($localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCMIN} < $orgGCMin);
            $orgGCMax    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCMAX}
                if($localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCMAX} > $orgGCMax);
            # ---- user-selected fragments
            # (LENGTH)
            $orguLenCount += $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLCOUNT};
            $orguLenTot   += $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLTOT};      # <----------------------- HERE!
            $orguLenMin    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLMIN}
                if($localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLMIN} < $orguLenMin);
            $orguLenMax    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLMAX}
                if($localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLMAX} > $orguLenMax);
            # (%GC)
            $orguGCCount += $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCCOUNT};
            $orguGCTot   += $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCTOT};
            $orguGCMin    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCMIN}
                if($localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCMIN} < $orguGCMin);
            $orguGCMax    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCMAX}
                if($localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCMAX} > $orguGCMax);
            # ----

        } # foreach(@replicons)

        ###################################################
        # Update statistics for ORGANISM
        # ..apply seqLength ranges instead of exact no., since seqs of length 1 or 2 mean nothing
        ###################################################
        $localStats{$org}->{STATS}->{LGENOME}    = $totGenomeLen;
        $localStats{$org}->{STATS}->{LREMOVED}   = sprintf("%.2f", ($totGenomeLen - $orgLenTot));
        my $orgLenPctRemoved = ($totGenomeLen == 0) ? (0) : ((1 - $orgLenTot/$totGenomeLen)*100);
        $localStats{$org}->{STATS}->{LPCTREMOVED}= sprintf("%.2f", $orgLenPctRemoved);

        # All fragments
        # (LENGTH)
        $localStats{$org}->{STATS}->{LCOUNT}= $orgLenCount;
        $localStats{$org}->{STATS}->{LTOT}  = $orgLenTot;
        my $orgLenMean = ($orgLenCount == 0) ? (0) : ($orgLenTot/$orgLenCount);
        $localStats{$org}->{STATS}->{LMEAN} = sprintf("%.2f", $orgLenMean);
        $localStats{$org}->{STATS}->{LMIN}  = $orgLenMin;
        $localStats{$org}->{STATS}->{LMAX}  = $orgLenMax;
        # All fragments
        # (LENGTH)
        $localStats{$org}->{STATS}->{GCCOUNT}= $orgGCCount;
        $localStats{$org}->{STATS}->{GCTOT}  = $orgGCTot;
        my $orgGCMean = ($orgGCCount == 0) ? (0) : ($orgGCTot/$orgGCCount);
        $localStats{$org}->{STATS}->{GCMEAN} = sprintf("%.2f", $orgGCMean);
        $localStats{$org}->{STATS}->{GCMIN}  = $orgGCMin;
        $localStats{$org}->{STATS}->{GCMAX}  = $orgGCMax;
        # User-selected fragments
        # (LENGTH)
        $localStats{$org}->{STATS}->{uLCOUNT}= $orguLenCount;
        $localStats{$org}->{STATS}->{uLTOT}  = $orguLenTot;
        my $orguLenMean = ($orguLenCount == 0) ? (0) : ($orguLenTot/$orguLenCount);
        $localStats{$org}->{STATS}->{uLMEAN} = sprintf("%.2f", $orguLenMean);
        $localStats{$org}->{STATS}->{uLMIN}  = $orguLenMin;
        $localStats{$org}->{STATS}->{uLMAX}  = $orguLenMax;
        # User-selected fragments
        # (LENGTH)
        $localStats{$org}->{STATS}->{uGCCOUNT}= $orguGCCount;
        $localStats{$org}->{STATS}->{uGCTOT}  = $orguGCTot;
        my $orguGCMean = ($orguGCCount == 0) ? (0) : ($orguGCTot/$orguGCCount);
        $localStats{$org}->{STATS}->{uGCMEAN} = sprintf("%.2f", $orguGCMean);
        $localStats{$org}->{STATS}->{uGCMIN}  = $orguGCMin;
        $localStats{$org}->{STATS}->{uGCMAX}  = $orguGCMax;

    } # foreach(@orgs)

    my $iter1 = new Benchmark;
    my $statText = timestr(timediff($iter1,$iter0));
       $statText = $1 if($statText =~ m/^\s*(\d+)\s+/);
       $statText = $1 if($statText =~ m/(\d+\.?\d+)\s+/);
    my $returnText = "calculated in [$statText] wallsecs... ";

    # Update GLOBAL %organisms before returning
    {
        #lock %organisms;
        foreach my $org (keys %localStats) {
            $organisms{$org}->{STATS}->{uGCMEAN}    = $localStats{$org}->{STATS}->{uGCMEAN};
            $organisms{$org}->{STATS}->{uGCMIN}     = $localStats{$org}->{STATS}->{uGCMIN};
            $organisms{$org}->{STATS}->{uGCMAX}     = $localStats{$org}->{STATS}->{uGCMAX};
            $organisms{$org}->{STATS}->{uGCCOUNT}   = $localStats{$org}->{STATS}->{uGCCOUNT};
            $organisms{$org}->{STATS}->{uGCTOT}     = $localStats{$org}->{STATS}->{uGCTOT};
            $organisms{$org}->{STATS}->{uLMEAN}     = $localStats{$org}->{STATS}->{uLMEAN};
            $organisms{$org}->{STATS}->{uLMIN}      = $localStats{$org}->{STATS}->{uLMIN};
            $organisms{$org}->{STATS}->{uLMAX}      = $localStats{$org}->{STATS}->{uLMAX};
            $organisms{$org}->{STATS}->{uLCOUNT}    = $localStats{$org}->{STATS}->{uLCOUNT};
            $organisms{$org}->{STATS}->{uLTOT}      = $localStats{$org}->{STATS}->{uLTOT};
            
            $organisms{$org}->{STATS}->{GCMEAN}     = $localStats{$org}->{STATS}->{GCMEAN};
            $organisms{$org}->{STATS}->{GCMIN}      = $localStats{$org}->{STATS}->{GCMIN};
            $organisms{$org}->{STATS}->{GCMAX}      = $localStats{$org}->{STATS}->{GCMAX};
            $organisms{$org}->{STATS}->{GCCOUNT}    = $localStats{$org}->{STATS}->{GCCOUNT};
            $organisms{$org}->{STATS}->{GCTOT}      = $localStats{$org}->{STATS}->{GCTOT};
            $organisms{$org}->{STATS}->{LMEAN}      = $localStats{$org}->{STATS}->{LMEAN};
            $organisms{$org}->{STATS}->{LMIN}       = $localStats{$org}->{STATS}->{LMIN};
            $organisms{$org}->{STATS}->{LMAX}       = $localStats{$org}->{STATS}->{LMAX};
            $organisms{$org}->{STATS}->{LCOUNT}     = $localStats{$org}->{STATS}->{LCOUNT};
            $organisms{$org}->{STATS}->{LTOT}       = $localStats{$org}->{STATS}->{LTOT};

            $organisms{$org}->{STATS}->{LGENOME}    = $localStats{$org}->{STATS}->{LGENOME};
            $organisms{$org}->{STATS}->{LREMOVED}   = $localStats{$org}->{STATS}->{LREMOVED};
            $organisms{$org}->{STATS}->{LPCTREMOVED}= $localStats{$org}->{STATS}->{LPCTREMOVED};
            foreach my $cFilename (keys %{ $localStats{$org}->{REPL} }) {
                foreach my $contigID (keys %{ $localStats{$org}->{REPL}->{$cFilename}->{CONTIG} }) {
                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCMEAN}
                    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCMEAN};
                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCSD}  
                    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCSD};
                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCMIN}
                    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCMIN};
                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCMAX}
                    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCMAX};
                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCCOUNT}
                    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCCOUNT};
                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCTOT}
                    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCTOT};
                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LMEAN}
                    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LMEAN};
                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LSD}
                    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LSD};
                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LMIN}
                    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LMIN};
                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LMAX}
                    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LMAX};
                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LTOT}
                    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LTOT};
                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LCOUNT}
                    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LCOUNT};
                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LREMOVED}
                    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LREMOVED};
                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LPCTREMOVED}
                    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LPCTREMOVED};


                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCMEAN}
                    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCMEAN};
                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCSD}
                    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCSD};
                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCMIN}
                    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCMIN};
                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCMAX}
                    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCMAX};
                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCCOUNT}
                    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCCOUNT};
                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCTOT}
                    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCTOT};
                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLMEAN}
                    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLMEAN};
                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLSD}
                    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLSD};
                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLMIN}
                    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLMIN};
                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLMAX}
                    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLMAX};
                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLTOT}
                    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLTOT};
                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLCOUNT}
                    = $localStats{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLCOUNT};
                } #CONTIGID
            } #CFILENAME
        } #ORG
    } #LOCK
    # --------------------------------------


    my $iter2 = new Benchmark;
    my $statText2 = timestr(timediff($iter2,$iter1));
       $statText2 = $1 if($statText2 =~ m/^\s*(\d+)\s+/);
       $statText2 = $1 if($statText2 =~ m/(\d+\.?\d+)\s+/);
    $returnText .= "updated in [$statText2] wallsecs";
    
    return $returnText;
 
}
################################################################################
#  -----------------------------
#  DESCRIPTION OF THE STATISTICS
#  -----------------------------
#    Length_Replicon = total number of bases in the replicon 
#    numUniqueWords  = total number of words(l-mers) found in the replicon
#    Mean%GC_UniqueWords = the mean %GC of all the unique words found in the replicon
#    STDEV = standard deviation of the Mean%GC
#    COV = coefficient of variation = STDEV/MEAN*100
#    MIN = smallest %GC value found
#    MAX = largest %GC value found
#    oWordCoverage = Overlapping Word Coverage = (word length)*(numUniqueWords)
#    LengthTotalFrags = sum total lenth (in bp) of all unique fragments (unfiltered)
#    n/oWordCoverage = non-overlapping word coverage = the length (in bp) of all non-overlapping parts of the unique words = (Length_Replicon) - (LengthTotalFrags) [unfiltered]
#    %n/oWordCoverage = the % of Length_Replicon which is the non-overlapping word coverage
#    %n/o2oWordCoverage = the % non-overlapping to overlapping coverage = ratio of the non-overlapping to the overlapping unique word lengths of the replicon
#    uLengthTotalFrags = sum total length (in bp) of all unique fragments (filtered for min.length)
#    u_n/oWordCoverage = the length (in bp) of all non-overlapping parts of the unique words = (Length_Replicon) - (uLengthTotalFrags) [filtered for min.length]
#    u_%n/oWordCoverage = the % of Length_Replicon which is the user's (length-filtered) non-overlapping word coverage
#    u_%n/o2oWordCoverage = ratio of the non-overlapping (filtered for min.length) to the overlapping unique word lengths of the replicon
#    MeanLength_UnfiltFrags = mean length of all unique fragments (unfiltered)
#    numUnfiltFrags = the count of the unique fragments (unfiltered)
#    Mean%GC_UnfiltFrags = the mean %GC of all unique fragments (unfiltered)
#    MeanLength_FiltFrags = the mean length of all unique fragments (filtered for min.length)
#    Mean%GC_FiltFrags = the mean %GC of all unique fragments (filtered for min.length)
################################################################################
sub printStatsToFile {

    #    %summaryTaxTree = (
    #       $rank    =>  {  $name1  => { COUNT => count(),
    #                                    PREV  => { RANKNAME => { $name1 = count(),
    #                                                             $name2 = count(),
    #                                                               ...
    #                                                           },
    #                                               LENGTH   => length(RANKNAME()),
    #                                             },
    #                                  },
    #                       $name2  => { ... },
    #                         ...
    #                    },
    #        ...
    #    );
    #
    # see "test_displayTax5.pl" for additional info.

    my %summaryTaxTree = ();

    ## Take relevant parameters from logfile                                     EDIT HERE
    #my $statsCore = $logfilename;
    #my $statsCore = $summaryFilename;
    #if($statsCore =~ m/(GOTTCHA)\.(.+)\.ALL.fna$/) {
    #    $statsCore = ".".$2;
    #}
    #else {
    #    $statsCore = "";
    #}

    
    my $prefix            = "GOTTCHA";
    my $vs                = $compareSameTaxLevel ? "1" : "0";
    my $taxLevelDetail    = "$taxLevel$vs";
    my $contigText        = "c$totalContigs";
    my $compressionDetail = "s$minStringPackLen";
    my $lmerText          = "k$k";
    my $uniqueText        = "u$minUniqueFragLength";
    my $partitionDetail   = "n$nThreads.p$prefixLength";
    my $dbLabelDetail     = $contigText."."
                           .$compressionDetail."."
                           .$lmerText."."
                           .$uniqueText."."
                           .$partitionDetail;
    
    my $statsCore = "$taxLevelDetail.$dbLabel.$dbLabelDetail";
    
    my $statsByContigFilename = $outdir."/statistics_byContig.$statsCore.csv";
    my $statsByOrgFilename    = $outdir."/statistics_byOrganism.$statsCore.csv";
    my $dbStatsFilename       = $outdir."/statistics_GOTTCHA.$statsCore.csv";

    open my $STATS, '>', $statsByContigFilename || die "Cannot open statistics file \"$statsByContigFilename\" for write.\n";
    stat_log("----------------------------------");
    stat_log_("Statistics file (by contig) created:       \"$statsByContigFilename\"...");
    my $iterX = new Benchmark;
    print $STATS "SUPERKINGDOM\tPHYLUM\tCLASS\tORDER\tFAMILY\tGENUS\tSPECIES\tORGANISM\tCONTIG\t"
                ."Length_Contig(bp)\t"
                ."TotalNucsRemoved(bp)\t"
                ."\%ContigRemoved\t"
                # Words
                ."numUniqueWords\t"
                ."Mean\%GC_UniqueWords\t"
                ."+/-\t"
                ."STDEV\t"
                ."COV\t"
                ."MIN\t"
                ."MAX\t"

                ."oWordCoverage\t"
                ."LengthTotalFrags\t"
                ."n/oWordCoverage\t"
                ."\%n/oWordCoverage\t"
                ."\%n/o2oWordCoverage\t"

                ."uLengthTotalFrags\t"
                ."u_n/oWordCoverage\t"
                ."\%u_n/oWordCoverage\t"
                ."\%u_n/o2oWordCoverage\t"
                # Unique fragments, all (unfiltered)
                ."MeanLength_UnfiltFrags\t"
                ."+/-\t"
                ."STDEV\t"
                ."numUnfiltFrags\t"
                ."MIN\t"
                ."MAX\t"
                ."Mean\%GC_UnfiltFrags\t"
                ."+/-\t"
                ."STDEV\t"
                ."COV\t"
                ."MIN\t"
                ."MAX\t"
                # Unique fragments, user-selected (filtered)
                ."MeanLength_FiltFrags\t"
                ."+/-\t"
                ."STDEV\t"
                ."numFiltFrags\t"
                ."MIN\t"
                ."MAX\t"
                ."Mean\%GC_FiltFrags\t"
                ."+/-\t"
                ."STDEV\t"
                ."COV\t"
                ."MIN\t"
                ."MAX\n";
                
    ORG: foreach my $org (@orgNames) {
        next ORG unless (exists $organisms{$org});
        
        my $s  = $organisms{$org}->{TAXTREE}->{S}; 
#        my $ss = $organisms{$org}->{TAXTREE}->{SS};
        my $g  = $organisms{$org}->{TAXTREE}->{G};
        my $f  = $organisms{$org}->{TAXTREE}->{F};
        my $o  = $organisms{$org}->{TAXTREE}->{O};
        my $c  = $organisms{$org}->{TAXTREE}->{C};
        my $p  = $organisms{$org}->{TAXTREE}->{P};
        my $sk = $organisms{$org}->{TAXTREE}->{SK};
#        print $STATS $org."|$s|$g|$f|$o|$c|$p|$sk\n";

        foreach my $cFilename (keys %{ $organisms{$org}->{REPL} }) {
            
            foreach my $contigID (keys %{ $organisms{$org}->{REPL}->{$cFilename}->{CONTIG} }) {

            # non-overlapping word coverage = (contig_length) - (total_unique_frag_length)
            my $NOWC  = $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{LENGTH} - $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LTOT};   # <---- HERE! #<---- ERROR:  LTOT must NOT be > contig length
            my $uNOWC = $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{LENGTH} - $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLTOT};  # <---- HERE! #<---- ERROR: uLTOT must not be > contig length

            # overlapping word coverage = (word_size) * (num_words)
            my  $OWC  = ($k)*($organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{WORDS}->{COUNT});

            # percent Non-Overlapping Word Coverage = NOWC / contig_length * 100
            my $pctNOWC  = ($organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{LENGTH} == 0)
                         ? (0)
                         : sprintf("%.2f",  $NOWC/$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{LENGTH}*100);
            my $pctuNOWC = ($organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{LENGTH} == 0)
                         ? (0)
                         : sprintf("%.2f", $uNOWC/$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{LENGTH}*100);
            
            # percent non-overlapping to overlapping word coverage = NOWC / OWC * 100
            my $pctNOWC2OWC  = ($OWC == 0) 
                             ? (0) 
                             : sprintf("%.2f",  $NOWC/$OWC*100);
            my $pctuNOWC2OWC = ($OWC == 0) 
                             ? (0)
                             : sprintf("%.2f", $uNOWC/$OWC*100);

            # coeff. of variance (all fragments %GC)
            my $COV = ($organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCMEAN} == 0)
                    ? (0)
                    : (sprintf("%.2f", $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCSD}/
                                       $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCMEAN}*100));
            # coeff. of variance (user's %GC)
            my $uCOV = ($organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCMEAN} == 0)
                     ? (0)
                     : (sprintf("%.2f", $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCSD}/
                                        $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCMEAN}*100));
            # coeff. of variance (word's %GC)
            my $wCOV = ($organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{WORDS}->{GCMEAN} == 0)
                     ? (0)
                     : (sprintf("%.2f", $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{WORDS}->{GCSD}/
                                        $organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{WORDS}->{GCMEAN}*100));

            print $STATS "$sk\t$p\t$c\t$o\t$f\t$g\t$s\t";
            print $STATS $org."\t".$contigID
                        # l-mers
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{LENGTH}        # REPLICON_LENGTH
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LREMOVED}
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LPCTREMOVED}

                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{WORDS}->{COUNT}# NUM_UNIQUE_WORDS
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{WORDS}->{GCMEAN}."\t+/-"
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{WORDS}->{GCSD}
                        ."\t".$wCOV
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{WORDS}->{GCMIN}
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{WORDS}->{GCMAX}
                        ."\t".$OWC                                               #       OVERLAPPING WORD COVERAGE
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LTOT} # ALL_FRAGS_TOTAL_LENGTH                            #<------------------
                        ."\t".$NOWC                                              #   NON-OVERLAPPING WORD COVERAGE
                        ."\t".$pctNOWC                                           # % NON-OVERLAPPING WORD COVERAGE
                        ."\t".$pctNOWC2OWC                                       # % NON-OVERLAPPING TO OVERLAPPING WORD COVERAGE
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLTOT}# USR_FRAGS_TOTAL_LENGTH                            #<------------------
                        ."\t".$uNOWC                                             #   user's NON-OVERLAPPING WORD COV.
                        ."\t".$pctuNOWC                                          # % user's NON-OVERLAPPING WORD COV.
                        ."\t".$pctuNOWC2OWC                                      # % user's NON-OVERLAPPING to OVERLAPPING WORD COV.
                        # all unique fragments (Length)
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LMEAN}."\t+/-"
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LSD}
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LCOUNT}
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LMIN}
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{LMAX}
                        # all unique fragments (%GC)
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCMEAN}."\t+/-"
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCSD}
                        ."\t".$COV
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCMIN}
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{GCMAX}
                        # user's unique fragments (Length)
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLMEAN}."\t+/-"
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLSD}
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLCOUNT}
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLMIN}
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uLMAX}
                        # user's unique fragments (%GC)
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCMEAN}."\t+/-"
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCSD}
                        ."\t".$uCOV
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCMIN}
                        ."\t".$organisms{$org}->{REPL}->{$cFilename}->{CONTIG}->{$contigID}->{STATS}->{uGCMAX}
                        
                        ."\n";
                        
            } #CONTIGID
        } #CFILENAME
    } #ORG
    print $STATS "";
    close $STATS;
    my $iterY = new Benchmark;
    my $timetext = timestr(timediff($iterY, $iterX));
       $timetext = $1 if($timetext =~ m/^\s*(\d+)\s+/);
       $timetext = $1 if($timetext =~ m/(\d+\.?\d+)\s+/);
    _stat_log("done in [$timetext] wallsecs.");


    ################################################
    # Save Organism-based Statistics
    ################################################
    open my $STATS_BYORG, '>', $statsByOrgFilename;
    $iterX = new Benchmark;
    stat_log_("Statistics file (by organism) created:     \"$statsByOrgFilename\"...");
    print $STATS_BYORG "SUPERKINGDOM\tPHYLUM\tCLASS\tORDER\tFAMILY\tGENUS\tSPECIES\tORGANISM\t";
    print $STATS_BYORG "TotalContigLength(bp)\t";
    print $STATS_BYORG "TotalNucsRemoved(bp)\t";
    print $STATS_BYORG "\%GenomeRemoved\t";
    print $STATS_BYORG "TotalFragLength,Unfilt(bp)\t";
    print $STATS_BYORG "TotalFrags,Unfilt\t";
    print $STATS_BYORG "MeanLength,Unfilt(bp)\t";
    print $STATS_BYORG "MinLength,Unfilt(bp)\t";
    print $STATS_BYORG "MaxLength,Unfilt(bp)\t";
    print $STATS_BYORG "Mean\%GC,Unfilt\t";
    print $STATS_BYORG "Min\%GC,Unfilt\t";
    print $STATS_BYORG "Max\%GC,Unfilt\t";
    print $STATS_BYORG "TotalFragLength,Filter(bp)\t";
    print $STATS_BYORG "TotalFrags,Filter\t";
    print $STATS_BYORG "MeanLength,Filter(bp)\t";
    print $STATS_BYORG "MinLength,Filter(bp)\t";
    print $STATS_BYORG "MaxLength,Filter(bp)\t";
    print $STATS_BYORG "Mean\%GC,Filter\t";
    print $STATS_BYORG "Min\%GC,Filter\t";
    print $STATS_BYORG "Max\%GC,Filter\n";

    ORG: foreach my $org (@orgNames) {
        next ORG if(!exists $organisms{$org});
        my $s  = $organisms{$org}->{TAXTREE}->{S}; 
        my $g  = $organisms{$org}->{TAXTREE}->{G};
        my $f  = $organisms{$org}->{TAXTREE}->{F};
        my $o  = $organisms{$org}->{TAXTREE}->{O};
        my $c  = $organisms{$org}->{TAXTREE}->{C};
        my $p  = $organisms{$org}->{TAXTREE}->{P};
        my $sk = $organisms{$org}->{TAXTREE}->{SK};

        print $STATS_BYORG "$sk\t$p\t$c\t$o\t$f\t$g\t$s\t$org\t";

        print $STATS_BYORG $organisms{$org}->{STATS}->{LGENOME}."\t";
        print $STATS_BYORG $organisms{$org}->{STATS}->{LREMOVED}."\t";
        print $STATS_BYORG $organisms{$org}->{STATS}->{LPCTREMOVED}."\t";

        print $STATS_BYORG $organisms{$org}->{STATS}->{LTOT}."\t";
        print $STATS_BYORG $organisms{$org}->{STATS}->{LCOUNT}."\t";
        print $STATS_BYORG $organisms{$org}->{STATS}->{LMEAN}."\t";
        print $STATS_BYORG $organisms{$org}->{STATS}->{LMIN}."\t";
        print $STATS_BYORG $organisms{$org}->{STATS}->{LMAX}."\t";
        # All fragments
        # (LENGTH)
        print $STATS_BYORG $organisms{$org}->{STATS}->{GCMEAN}."\t";
        print $STATS_BYORG $organisms{$org}->{STATS}->{GCMIN}."\t";
        print $STATS_BYORG $organisms{$org}->{STATS}->{GCMAX}."\t";
        # User-selected fragments
        # (LENGTH)
        print $STATS_BYORG $organisms{$org}->{STATS}->{uLTOT}."\t";
        print $STATS_BYORG $organisms{$org}->{STATS}->{uLCOUNT}."\t";
        print $STATS_BYORG $organisms{$org}->{STATS}->{uLMEAN}."\t";
        print $STATS_BYORG $organisms{$org}->{STATS}->{uLMIN}."\t";
        print $STATS_BYORG $organisms{$org}->{STATS}->{uLMAX}."\t";
        # User-selected fragments
        # (LENGTH)
        print $STATS_BYORG $organisms{$org}->{STATS}->{uGCMEAN}."\t";
        print $STATS_BYORG $organisms{$org}->{STATS}->{uGCMIN}."\t";
        print $STATS_BYORG $organisms{$org}->{STATS}->{uGCMAX}."\n";
    } #ORGS
    
    close $STATS_BYORG;
    $iterY = new Benchmark;
       $timetext = timestr(timediff($iterY, $iterX));
       $timetext = $1 if($timetext =~ m/^\s*(\d+)\s+/);
       $timetext = $1 if($timetext =~ m/(\d+\.?\d+)\s+/);
    _stat_log("done in [$timetext] wallsecs.");


    ################################################
    # Save Global DB Statistics
    ################################################
    my @ranksExt = values(%taxAbbrExt);

    ORG: foreach my $org (sort @orgNames) {
        next ORG if(!exists $organisms{$org});
        RANK: foreach my $rankIdx (0..(scalar(@ranksExt)-1)) {
            my $rank = $ranksExt[$rankIdx];
            if(!$organisms{$org}->{TAXTREE}->{$rank}) {
                $organisms{$org}->{TAXTREE}->{$rank} = "Unassigned";
                #print "going to next rank...\n";
                next RANK;
            }
            $organisms{$org}->{TAXTREE}->{$rank} =~ s/^\s+//g;
            $organisms{$org}->{TAXTREE}->{$rank} =~ s/\s+$//g;
            my $name = $organisms{$org}->{TAXTREE}->{$rank};
            next RANK if($name eq "Unassigned");                                # <--- 2012-03-19. Fixes error: "Invalid value for shared scalar..." below

            # Total no. of unique SUBTAXA (i.e. total count of ->{PREV}
            $summaryTaxTree{$rank}->{$name}->{COUNT}++;     
            
            # Add rank NAME to {PREV}->{RANK}
            $summaryTaxTree{$rank}->{$name}
                                  ->{PREV}
                                  ->{RANKNAME}
                                  ->{ $organisms{$org}->{TAXTREE}->{$ranksExt[$rankIdx-1]} }++ 
                if($rankIdx > 0);
            
            # Add $organisms{$org}->{STATS}->{LGENOME} if at SPECIES level
            if($rank eq "SS") {
                $summaryTaxTree{$rank}->{$name}->{PREV}->{LENGTH}               = $organisms{$name}->{STATS}->{LGENOME};
                $summaryTaxTree{$rank}->{$name}->{PREV}->{STATS}->{GCCOUNT}     = $organisms{$name}->{STATS}->{GCCOUNT}; # No. of %GC calcs done (i.e. no. of unfilt. frags)
                $summaryTaxTree{$rank}->{$name}->{PREV}->{STATS}->{GCMAX}       = $organisms{$name}->{STATS}->{GCMAX}; # Max %GC of unfilt. frags
                $summaryTaxTree{$rank}->{$name}->{PREV}->{STATS}->{GCMEAN}      = $organisms{$name}->{STATS}->{GCMEAN}; # Mean %GC of unfilt. frags
                $summaryTaxTree{$rank}->{$name}->{PREV}->{STATS}->{GCMIN}       = $organisms{$name}->{STATS}->{GCMIN}; # Min %GC of unfilt. frags
                $summaryTaxTree{$rank}->{$name}->{PREV}->{STATS}->{GCTOT}       = $organisms{$name}->{STATS}->{GCTOT}; # Sum of all %GC values for all unfilt. frags
                $summaryTaxTree{$rank}->{$name}->{PREV}->{STATS}->{LCOUNT}      = $organisms{$name}->{STATS}->{LCOUNT}; # Same as GCCOUNT
                $summaryTaxTree{$rank}->{$name}->{PREV}->{STATS}->{LMAX}        = $organisms{$name}->{STATS}->{LMAX}; # Max Length out of all unfilt. frags
                $summaryTaxTree{$rank}->{$name}->{PREV}->{STATS}->{LMEAN}       = $organisms{$name}->{STATS}->{LMEAN}; # Mean Length out of all unfilt. frags
                $summaryTaxTree{$rank}->{$name}->{PREV}->{STATS}->{LMIN}        = $organisms{$name}->{STATS}->{LMIN}; # Min Length out of all unfilt. frags
                $summaryTaxTree{$rank}->{$name}->{PREV}->{STATS}->{LPCTREMOVED} = $organisms{$name}->{STATS}->{LPCTREMOVED}; # % of the Genome removed
                $summaryTaxTree{$rank}->{$name}->{PREV}->{STATS}->{LREMOVED}    = $organisms{$name}->{STATS}->{LREMOVED}; # Total bp removed from all unfilt. frags
                $summaryTaxTree{$rank}->{$name}->{PREV}->{STATS}->{LTOT}        = $organisms{$name}->{STATS}->{LTOT}; # Total bp of all unfilt. frags
                $summaryTaxTree{$rank}->{$name}->{PREV}->{STATS}->{uGCCOUNT}    = $organisms{$name}->{STATS}->{uGCCOUNT};
                $summaryTaxTree{$rank}->{$name}->{PREV}->{STATS}->{uGCMAX}      = $organisms{$name}->{STATS}->{uGCMAX};
                $summaryTaxTree{$rank}->{$name}->{PREV}->{STATS}->{uGCMEAN}     = $organisms{$name}->{STATS}->{uGCMEAN};
                $summaryTaxTree{$rank}->{$name}->{PREV}->{STATS}->{uGCMIN}      = $organisms{$name}->{STATS}->{uGCMIN};
                $summaryTaxTree{$rank}->{$name}->{PREV}->{STATS}->{uGCTOT}      = $organisms{$name}->{STATS}->{uGCTOT};
                $summaryTaxTree{$rank}->{$name}->{PREV}->{STATS}->{uLCOUNT}     = $organisms{$name}->{STATS}->{uLCOUNT};
                $summaryTaxTree{$rank}->{$name}->{PREV}->{STATS}->{uLMAX}       = $organisms{$name}->{STATS}->{uLMAX};
                $summaryTaxTree{$rank}->{$name}->{PREV}->{STATS}->{uLMEAN}      = $organisms{$name}->{STATS}->{uLMEAN};
                $summaryTaxTree{$rank}->{$name}->{PREV}->{STATS}->{uLMIN}       = $organisms{$name}->{STATS}->{uLMIN};
                $summaryTaxTree{$rank}->{$name}->{PREV}->{STATS}->{uLTOT}       = $organisms{$name}->{STATS}->{uLTOT};
                #print "Adding ".$organisms{$name}->{STATS}->{LGENOME}."-bp from $name to $name(RANK=$rank)\n";
            }
        } #RANK
    } #ORG
    
    #------------------------------------
    # Loop through all ranks, skipping SS
    # foreach rankName found (S->G->...->SK), get all ->{PREV}->{RANKNAME} names ($prevRankName)
    # foreach $prevRankName, grab their lengths from $prevRank
    # **Must process these in order from S->G->F->O->C->P->SK

    foreach my $rankIdx (0..(scalar(@ranksExt)-1)) {
        my $rank = $ranksExt[$rankIdx];
        next if($rank eq "SS");     # Already processed these
        
        foreach my $taxon (keys %{ $summaryTaxTree{$rank} }) {
            my @prevTaxa = keys %{ $summaryTaxTree{$rank}->{$taxon}->{PREV}->{RANKNAME} };
            foreach my $prevTaxonIdx (0..(scalar(@prevTaxa)-1)) {
                my $prevTaxon = $prevTaxa[$prevTaxonIdx];
                my $prevRank   = $ranksExt[$rankIdx-1];
                my $prevLength = $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{LENGTH}
                               ? $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{LENGTH}
                               : 0;
                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{LENGTH} += $prevLength;

                # Total length of unfilt. fragments                
                my $prevGCcount    = $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{GCCOUNT}
                                   ? $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{GCCOUNT}
                                   : 0;
                my $prevGCtot      = $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{GCTOT}
                                   ? $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{GCTOT}
                                   : 0;
                my $prevLenCount   = $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{LCOUNT}
                                   ? $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{LCOUNT}
                                   : 0;
                my $prevLenRemoved = $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{LREMOVED}
                                   ? $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{LREMOVED}
                                   : 0;
                my $prevLenTot     = $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{LTOT}
                                   ? $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{LTOT}
                                   : 0;
                my $prevuGCcount   = $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{uGCCOUNT}
                                   ? $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{uGCCOUNT}
                                   : 0;
                my $prevuGCtot     = $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{uGCTOT}
                                   ? $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{uGCTOT}
                                   : 0;
                my $prevuLenCount  = $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{uLCOUNT}
                                   ? $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{uLCOUNT}
                                   : 0;
                my $prevuLenTot    = $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{uLTOT}
                                   ? $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{uLTOT}
                                   : 0;
                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{GCCOUNT} += $prevGCcount;
                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{GCTOT}   += $prevGCtot;
                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{LCOUNT}  += $prevLenCount;
                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{LREMOVED}+= $prevLenRemoved;
                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{LTOT}    += $prevLenTot;
                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{uGCCOUNT}+= $prevuGCcount;
                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{uGCTOT}  += $prevuGCtot;
                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{uLCOUNT} += $prevuLenCount;
                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{uLTOT}   += $prevuLenTot;

                # For min, max, averages, and %'s, set "=" instead of "+="
                my $prevGCmax  = $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{GCMAX}
                               ? $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{GCMAX}
                               : 0;
                my $prevGCmin  = $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{GCMIN}
                               ? $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{GCMIN}
                               : 0;
 #               my $prevGCmean    = ($summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{GCCOUNT} > 0)
 #                                 ? ($summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{GCTOT}/
 #                                    $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{GCCOUNT})
 #                                 : 0;
                my $prevLenMax  = $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{LMAX}
                                ? $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{LMAX}
                                : 0;
                my $prevLenMin  = $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{LMIN}
                                ? $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{LMIN}
                                : 0;
                my $prevLenMean = $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{LMEAN}
                                ? $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{LMEAN}
                                : 0;
                my $prevLenPctRemoved = $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{LPCTREMOVED}
                                      ? $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{LPCTREMOVED}
                                      : 0;
                my $prevuGCmax  = $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{uGCMAX}
                                ? $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{uGCMAX}
                                : 0;
                my $prevuGCmin  = $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{uGCMIN}
                                ? $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{uGCMIN}
                                : 0;
                                
                $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{uGCCOUNT} = 0 
                    if(!defined $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{uGCCOUNT});
                my $prevuGCmean = ($summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{uGCCOUNT} > 0)
                                ? ($summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{uGCTOT}/
                                   $summaryTaxTree{$prevRank}->{$prevTaxon}->{PREV}->{STATS}->{uGCCOUNT})
                                :  0;
               
#                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{GCMAX} = $prevGCmax if($prevTaxonIdx == 0);
#                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{GCMAX}      = ($prevGCmax > $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{GCMAX})
#                                                                                ? ($prevGCmax)
#                                                                                : ($summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{GCMAX});
#
#                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{GCMIN} = $prevGCmin if($prevTaxonIdx == 0);
#                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{GCMIN}      = ($prevGCmin < $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{GCMIN})
#                                                                                ? ($prevGCmin)
#                                                                                : ($summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{GCMIN});
#
#                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{GCMEAN} = $prevGCmean;
#                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{GCMEAN} = 0 unless ($summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{GCMEAN}); 
#                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{GCMEAN}     = ($summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{GCMEAN})
#                                                                                ? ($prevGCmean)
#                                                                                : ($summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{GCMEAN});

                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{LMAX} = $prevLenMax if($prevTaxonIdx == 0);
                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{LMAX} = ($prevLenMax > $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{LMAX})       # <----- HERE!
                                                                          ? ($prevLenMax)
                                                                          : ($summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{LMAX});
                                                                                
                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{LMIN} = $prevLenMin if($prevTaxonIdx == 0);
                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{LMIN} = ($prevLenMin < $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{LMIN})       # <----- HERE!
                                                                          ? ($prevLenMin)
                                                                          : ($summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{LMIN});
                
#                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{LMEAN}      = $prevLenMean;

                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{LPCTREMOVED}= $prevLenPctRemoved;               # <---------- not being utilized!!

#                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{uGCMAX} = $prevuGCmax if($prevTaxonIdx == 0);
#                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{uGCMAX}     = ($prevuGCmax > $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{uGCMAX}) 
#                                                                                ? ($prevuGCmax)
#                                                                                : ($summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{uGCMAX});
#                                                                                
#                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{uGCMIN} = $prevuGCmin if($prevTaxonIdx == 0);
#                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{uGCMIN}     = ($prevuGCmin < $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{uGCMIN})
#                                                                                ? ($prevuGCmin)
#                                                                                : ($summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{uGCMIN});
#                                                                                
 #               $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{uGCMEAN}    = $prevuGCmean;

                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{uLMAX} = $prevLenMax if($prevTaxonIdx == 0);
#                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{uLMAX} = ($prevTaxonIdx == 0) 
#                                                                           ? $prevLenMax 
#                                                                           : 0;
                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{uLMAX} = ($prevLenMax > $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{uLMAX})       # <----- HERE!
                                                                           ? ($prevLenMax)
                                                                           : ($summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{uLMAX});
                                                                                
                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{uLMIN} = $prevLenMin if($prevTaxonIdx == 0);
                $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{uLMIN} = ($prevLenMin < $summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{uLMIN})       # <----- HERE!
                                                                           ? ($prevLenMin)
                                                                           : ($summaryTaxTree{$rank}->{$taxon}->{PREV}->{STATS}->{uLMIN});
                #----------------

            } #PREVTAXON
        } #TAXON
    } #RANKIDX


    #------------------------------------
    # Write to file
    $iterX = new Benchmark;
    stat_log_("Statistics file on entire GOTTCHA DB created: \"$dbStatsFilename\"...");
    open my $DBSTATS, '>', $dbStatsFilename || die "Cannot open $dbStatsFilename for write!\n";

    # Print headers to $DBSTATS
    print $DBSTATS "RANK\t#_UNIQUE_TAXA\tTAXON_ORDER_IN_RANK\tCONSTITUENT_TAXA\t"
                  ."SIZE_OF_TAXA(bp)\t"
#                  ."MEAN\%GC,TAXA\t"                                            # <--- not yet computed
                  ."TOTAL_BP_REMOVED,TAXA(bp)\t\%TAXA_REMOVED\t"
#                  ."MEAN\%GC_REMOVED,TAXA\t"                                    # <--- not yet computed
                  ."#_SUBTAXA\tORDINAL_POS_IN_SUBTAXA\tSUBTAXA\tSIZE_OF_SUBTAXA(bp)\t"
#                  ."MEAN\%GC,SUBTAXA\t"                                         # <--- not yet computed
                  ."\%SIZE_OF_TAXA\tTOTAL_BP_REMOVED,SUBTAXA(bp)\t\%SUBTAXA_REMOVED\t"
#                  ."MEAN\%GC_REMOVED,SUBTAXA\t"                                 # <--- not yet computed
                  ."TOTAL_BP_REMAINING(bp)\t\%SUBTAXA_REMAINING\t"
#                  ."MEAN\%GC_REMAINING,SUBTAXA\t"                               # <--- not yet computed
                  ."MIN_FRAG_LEN,UNFILT.(bp)\tMAX_FRAG_LEN,UNFILT.(bp)\t"
#                  ."MEAN\%GC,UNFILT.\t"                                         # <--- not yet computed
                  ."MIN_FRAG_LEN,FILT.(bp)\tMAX_FRAG_LEN,FILT.(bp)\t"
#                  ."MEAN\%GC,FILT.\t"                                           # <--- not yet computed
                  ."\n";
    
    
    RANK: foreach my $rankIdx (1..(scalar(@ranksExt)-1)) {
        my $rank = $ranksExt[$rankIdx];
    
        # Total no. of (REDUNDANTLY) represented taxons in current rank
        my @uniqueNames = sort keys %{ $summaryTaxTree{$rank} };  # Staphylococcus, or Streptococcus
        my $numUniqueNames = scalar(@uniqueNames);                # 2
        
        # Total no. of (UNIQUELY) represented taxons in current rank
        my @copies = @{ $summaryTaxTree{$rank} }{@uniqueNames};
        my $numNonuniqueNames  = 0;
           $numNonuniqueNames += $_ foreach (@copies);            # 3+1 = 4

        # All organism's represented in current rank           
        my @orgsInRank = ();
        my %uniqueRankNames = ();
           @uniqueRankNames{@uniqueNames} = ();
        foreach my $org (@orgNames) {
            next if(!exists $organisms{$org});
            # Save organism's name if it's rank name is in @uniqueNames
            push(@orgsInRank, $org) if(exists $uniqueRankNames{ $organisms{$org}->{TAXTREE}->{ $rank } });
        } 

        # Display the unique taxon NAMES within the current rank
        UNIQUENAME: foreach my $uniqueNameIdx (0..(scalar(@uniqueNames)-1)) {
            my $uniqueName = $uniqueNames[$uniqueNameIdx];
            
            # If exists(), grab the $summaryTaxTree{$rank}->{$uniqueName} organisms and the associated data
            my $currLength         = $summaryTaxTree{$rank}->{$uniqueName}->{PREV}->{LENGTH};
            my $currMeanGC         = ($summaryTaxTree{$rank}->{$uniqueName}->{PREV}->{STATS}->{GCCOUNT} > 0)
                                   ? (sprintf("%.2f", ($summaryTaxTree{$rank}->{$uniqueName}->{PREV}->{STATS}->{GCTOT}/
                                                       $summaryTaxTree{$rank}->{$uniqueName}->{PREV}->{STATS}->{GCCOUNT})))
                                   : (0);
            my @prevUniqueNames    = keys %{ $summaryTaxTree{$rank}->{$uniqueName}->{PREV}->{RANKNAME} };
            my $numPrevUniqueNames = scalar(@prevUniqueNames);
            my $currLenMIN         = $summaryTaxTree{$rank}->{$uniqueName}->{PREV}->{STATS}->{LMIN}; 
            my $currLenMAX         = $summaryTaxTree{$rank}->{$uniqueName}->{PREV}->{STATS}->{LMAX};
            my $currLenRemoved     = $summaryTaxTree{$rank}->{$uniqueName}->{PREV}->{STATS}->{LREMOVED};
#            my $currLenPctRemoved  = sprintf("%.2f", $summaryTaxTree{$rank}->{$uniqueName}->{PREV}->{STATS}->{LPCTREMOVED});
            my $currLenPctRemoved  = ($currLength == 0) 
                                   ? (0)
                                   : sprintf("%.2f", $currLenRemoved/$currLength*100);
            my $currLentot         = sprintf("%.2f", $summaryTaxTree{$rank}->{$uniqueName}->{PREV}->{STATS}->{LTOT});
            if ($uniqueNameIdx == 0) {
                print $DBSTATS $revTaxAbbrExt{$rank}."\t"                                                  # RANK
                     .$numUniqueNames."\t"                                                                 # TAXON_NAME
                     .($uniqueNameIdx+1)."\t"                                                              # TAXON_ORDER_IN_RANK
                     ."";
            }
            else {
                print $DBSTATS $revTaxAbbrExt{$rank}."\t"                                                  # RANK
                     ."\t"                                                                                 # (spacer)
                     .($uniqueNameIdx+1)."\t"                                                              # TAXON_ORDER_IN_RANK
                     ."";                                                                                  
            }
            foreach my $prevUniqueNameIdx (0..($numPrevUniqueNames-1)) {
                my $prevUniqueName = $prevUniqueNames[$prevUniqueNameIdx];

                # Get length of previous taxa
                my $prevRank       = $ranksExt[$rankIdx-1];
                my $prevLength     = $summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{LENGTH}
                                   ? $summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{LENGTH}
                                   : 0;
                my $prevLenRemoved = $summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{STATS}->{LREMOVED}
                                   ? $summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{STATS}->{LREMOVED}
                                   : 0;
#                my $prevLenPctRemoved = sprintf("%.2f", $summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{STATS}->{LPCTREMOVED});    #<---- not updating properly so computing on the fly...
               my $prevLenPctRemoved = ($prevLength == 0) 
                                     ? (0)
                                     : sprintf("%.2f", $prevLenRemoved/$prevLength*100);
                
                # Get unique fragment lengths 
                $summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{STATS}->{GCCOUNT} = 0 
                    if (!defined $summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{STATS}->{GCCOUNT});
                my $unfiltMeanGC = ($summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{STATS}->{GCCOUNT} > 0)
                                 ? (sprintf("%.2f", $summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{STATS}->{GCTOT}/
                                                    $summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{STATS}->{GCCOUNT}))
                                 : (0);
                my $subtaxaPctMakeupOfTaxa = ($currLength == 0) 
                                           ? (0)
                                           : sprintf("%.2f", $prevLength/$currLength*100);
                my $pctSubtaxaRemaining    = ($prevLength == 0) 
                                           ? (0)
                                           : sprintf("%.2f", $summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{STATS}->{LTOT}/$prevLength*100);
                
                # Display data
                print $DBSTATS $revTaxAbbrExt{$rank}."\t\t\t" if($prevUniqueNameIdx > 0);                         
                if($prevUniqueNameIdx == 0) {   
                    print $DBSTATS $uniqueName."\t"                                                        # CONSTITUENT_TAXA
                         .$currLength."\t"                                                                 # SIZE_OF_TAXA
#                         ."\t"                                                                             # MEAN_%GC,TAXA (currently UNTRACKED)
                         .$currLenRemoved."\t"                                                             # TOTAL_BP_REMOVED,TAXA
                         .$currLenPctRemoved."\t"                                                          # %TAXA_REMOVED
#                         ."\t"                                                                             # MEAN%GC_REMOVED,TAXA
                         .$numPrevUniqueNames."\t"                                                         # #_SUBTAXA
                         ."";
                }
               print $DBSTATS "\t\t\t\t\t" if($prevUniqueNameIdx > 0);
                print $DBSTATS "".($prevUniqueNameIdx+1)."\t"                                               # SUBTAXON'S ORDINAL POSITION IN SUBTAXA
                     .$prevUniqueName."\t"                                                                  # SUBTAXON NAME
                     .$prevLength."\t"                                                                      # SIZE_OF_SUBTAXA
#                     ."\t"                                                                                  # MEAN%GC,SUBTAXA (currently UNTRACKED)
                     .$subtaxaPctMakeupOfTaxa."\t"                                                          # %SIZE_OF_TAXA
                     ."";
                     
                $prevLenRemoved    = 0 if(!defined $prevLenRemoved);
                $prevLenPctRemoved = 0 if(!defined $prevLenPctRemoved);
                $pctSubtaxaRemaining = 0 if(!defined $pctSubtaxaRemaining);
                $summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{STATS}->{LTOT} = 0
                    if(!defined $summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{STATS}->{LTOT});
                print $DBSTATS $prevLenRemoved."\t"                                                         # TOTAL_BP_REMOVED,SUBTAXA
                     .$prevLenPctRemoved."\t"                                                               # %SUBTAXA_REMOVED
#                     ."\t"                                                                                  # MEAN%GC_REMOVED,SUBTAXA
                     .$summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{STATS}->{LTOT}."\t"           # TOTAL_BP_REMAINING
                     .$pctSubtaxaRemaining."\t"                                                             # %SUBTAXA_REMAINING
#                     ."\t"                                                                                  # MEAN%GC_REMAINING,SUBTAXA
                     ."";
                
                # Display fragment data (use hash to store column headers)
                $summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{STATS}->{LMIN} = 0
                    if(!defined $summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{STATS}->{LMIN});
                $summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{STATS}->{LMIN} = 0
                    if(!defined $summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{STATS}->{LMAX});
                $summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{STATS}->{LMIN} = 0
                    if(!defined $summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{STATS}->{uLMIN});
                $summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{STATS}->{LMIN} = 0
                    if(!defined $summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{STATS}->{uLMAX});
                print $DBSTATS
                      $summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{STATS}->{LMIN}."\t"           # MIN_FRAG_LEN,UNFILT.(bp)
                     .$summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{STATS}->{LMAX}."\t"           # MAX_FRAG_LEN,UNFILT.(bp)
#                     .$unfiltMeanGC."\t"                                                                    # MEAN%GC,UNFILT.
#                print $DBSTATS $summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{STATS}->{uLCOUNT}."\t";
#                print $DBSTATS $summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{STATS}->{uLTOT}."\t";  
                     .$summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{STATS}->{uLMIN}."\t"          # MIN_FRAG_LEN,FILT.(bp)
                     .$summaryTaxTree{$prevRank}->{$prevUniqueName}->{PREV}->{STATS}->{uLMAX}."\t"          # MAX_FRAG_LEN,FILT.(bp)
#                     ."\t"                                                                                  # MEAN%GC,FILT.
                     ."\n";
            } #FOREACH
        } #UNIQUENAME
        print $DBSTATS "\n";
    } #RANK

    close ($DBSTATS);
    $iterY = new Benchmark;
    $timetext = timestr(timediff($iterY, $iterX));
    $timetext = $1 if($timetext =~ m/^\s*(\d+)\s+/);
    $timetext = $1 if($timetext =~ m/(\d+\.?\d+)\s+/);
    _stat_log("done in [$timetext] wallsecs.");

    return;
}
################################################################################
# ARGS: \@data
# Returns: $mean, $stdev, $min, $max
# If \@data is emtpy, returns 0 for all
#-------------------------------------------------------------------------------
sub doStats {

    my $count = scalar(@{ $_[0] });
    return (0,0,0,0,0) if($count == 0);

    my ($mean, $stdev, $min, $max, $sum);
    $stdev = 0;
    $mean  = 0;
    $sum   = 0;
    
    $max = $_[0]->[0];
    $min = $_[0]->[0];
    foreach (@{ $_[0] }) {
        $max = $_ if ($_ > $max);       # COMPLETE
        $min = $_ if ($_ < $min);       # COMPLETE
        $sum += $_;
    }
    $mean = $sum;
    $mean /= $count;                    # COMPLETE
    
    $stdev += ($_ - $mean)**2 foreach (@{ $_[0] });
    if($count == 1) {
        $stdev = 0;
    }
    else {
        $stdev  = sqrt($stdev/($count-1));
    }
    
    return ($mean, $stdev, $min, $max, $sum);
}
################################################################################
# Function: This sub is called between iterations of TAXLEVEL in main() to reset
#           all statistics.
################################################################################
sub resetOrgStats {

    my @level1fields = ("uGCMEAN","uGCMIN","uGCMAX","uGCCOUNT","uGCTOT","uLMEAN","uLMIN","uLMAX","uLCOUNT","uLTOT"
                       , "GCMEAN", "GCMIN", "GCMAX", "GCCOUNT", "GCTOT", "LMEAN", "LMIN", "LMAX", "LCOUNT", "LTOT"
                       ,"LGENOME", "LREMOVED", "LPCTREMOVED");

    my @level2fields = ("uGCMEAN","uGCMIN","uGCMAX","uGCCOUNT","uGCTOT","uLMEAN","uLMIN","uLMAX","uLCOUNT","uLTOT"
                       , "GCMEAN", "GCMIN", "GCMAX", "GCCOUNT", "GCTOT", "LMEAN", "LMIN", "LMAX", "LCOUNT", "LTOT"
                       ,"LREMOVED", "LPCTREMOVED", "GCSD","LSD","uGCSD","uLSD");

    foreach my $org (keys %organisms) {
        $organisms{$org}->{STATS}->{$_} = 0 foreach (@level1fields);
        foreach my $filename (keys %{ $organisms{$org}->{REPL} }) {
            foreach my $contigID (keys %{ $organisms{$org}->{REPL}->{$filename}->{CONTIG} }) {
                $organisms{$org}->{REPL}->{$filename}->{CONTIG}->{$contigID}->{STATS}->{$_} = 0 foreach (@level2fields);
            } #CONTIGID
        } #FILENAME
    } #ORG
    
    return;

}
################################################################################
#-------------------------------------------------------------------------------
sub archiveKmerIntersections {
    stat_log("----------------------------------");
    stat_log_("Archiving k-mer Intersections for incremental update function...");

    my $time0 = new Benchmark;
    
    my $datadir    = $outdir."/".$updateDir;
    my $dirPrefix  = $datadir."/".$storableDir;

    # This subroutine is called before the subtraction process is called, so
    # the only datastructures we copy from the partition directory to the
    # update directory are the following: 
    
    foreach my $rNick (@processedContigs) {

        my $destDir = $dirPrefix."/".$rNick;
        createDir($destDir);
        
        my $origDir = $pDirCodeLookup{ $pDirNickLookup{$rNick} }."/".$partitionDir."/".$rNick;
        
    # 1. "L2BM_###.dmp"
        foreach my $L2BMfile (glob "$origDir/L2BM_*.dmp") { 
            copy($L2BMfile, $destDir) || stat_log("Copy of L2BM_* {$rNick} failed: $!");
        }
    
    # 2. "kmers_eBV_hex.dmp"
        copy($origDir."/".$kmers_eBV_hex_Filename, $destDir) || stat_log("Copy of kmers_eBV_hex {$rNick} failed: $!");
    
    # 3. "kmer_starts.dmp"
        copy($origDir."/".$kmer_starts_Filename, $destDir) || stat_log("Copy of kmer_starts {$rNick} failed: $!");
    
    # 4. "seqBitMask.dmp"
        copy($origDir."/".$seqBMfilename, $destDir) || stat_log("Copy of (pre-subtraction) seqBM {$rNick} failed: $!");
    
    }

    my $time1 = new Benchmark;

    #my $timeTxt = q{};
    my $timediff = timestr(timediff($time1,$time0));
    if($timediff =~ m/^\s*(\S+)\s+/) {
        $benchRez = 5;
        $timediff = sciExpand($1, $benchRez);
    }
    elsif ($timediff =~ m/(\d+\.?\d+)\s+/) {
        $timediff = $1;
    }
    #my $loadTimeText = sciExpand($1, $benchRez) if(timestr(timediff($time1,$time0)) =~ m/^\s*(\S+)\s+/);
    _stat_log("done in [$timediff] wallsecs");

#    
#    dump2disk(\%lmerCoords, $dirPrefix."/lmerCoords.dmp", "lmerCoords", "bin");
#
    return;
}
################################################################################
# Function: Stores all datastructures, constants, and user options that are 
#           essential for incrementally updating the uniques database.
#-------------------------------------------------------------------------------
sub archiveRunData {

    stat_log("----------------------------------");
    stat_log("Archiving run data and user parameters for Incremental Update...");

    my $time0 = new Benchmark;

    my $dataDir    = $outdir."/".$updateDir;
    my $dirPrefix  = $dataDir."/".$storableDir;
    createDir($dirPrefix);
    my $configfile = $dataDir."/update.config";

    foreach my $rNick (@processedContigs) {
        my $origDir = $pDirCodeLookup{ $pDirNickLookup{$rNick} }."/".$partitionDir."/".$rNick;
        my $destDir = $dirPrefix."/".$rNick;
        
        # 5. "uniqueSeqBitMask.dmp"
        copy($origDir."/".$uniqueSeqBMfilename, $destDir) || stat_log("Copy of (post-subtraction) seqBM {$rNick} failed: $!");
        
        # UPDATE: 2012-10-09    Already stored in SUBTRACTION PHASE
        # 6. "uL2BMs.dmp"
        #copy($origDir."/".$unionizedL2BMsFilename, $destDir) || stat_log("Copy of (post-subtraction) Unionized L2BMs {$rNick} failed: $!");
    }
    
    # Storable Datastructures
    dump2disk(\%orgLookupByContig, $dirPrefix."/orgLookupByContig.dmp", "orgLookupByContig", "bin");
    dump2disk(\%fileNickLookup, $dirPrefix."/fileNickLookup.dmp", "fileNickLookup", "bin");
    dump2disk(\%rFileNickLookup, $dirPrefix."/rFileNickLookup.dmp", "rFileNickLookup", "bin");
    dump2disk(\%contigNickLookup, $dirPrefix."/contigNickLookup.dmp", "contigNickLookup", "bin");
    dump2disk(\%rContigNickLookup, $dirPrefix."/rContigNickLookup.dmp", "rContigNickLookup", "bin");
    dump2disk(\%contigID2filename, $dirPrefix."/contigID2filename.dmp", "contigID2filename", "bin");
    dump2disk(\@contigIDs, $dirPrefix."/contigIDs.dmp", "contigIDs", "bin");
    dump2disk(\@contigFilenames, $dirPrefix."/contigFilenames.dmp", "contigFilenames", "bin");
    dump2disk(\@processedContigs, $dirPrefix."/processedContigs.dmp", "processedContigs", "bin");
    dump2disk(\@processedOrgs, $dirPrefix."/processedOrgs.dmp", "processedOrgs", "bin");
#    dump2disk(\@fastaFilesOfUniques, $dirPrefix."/fastaFilesOfUniques.dmp", "fastaFilesOfUniques", "bin");
    dump2disk(\@orgNames, $dirPrefix."/orgNames.dmp", "orgNames", "bin");
    #dump2disk(\@mergedIndicesRefs, $dirPrefix."/mergedIndiciesRefs.dmp", "mergedIndiciesRefs", "bin");
    #dump2disk(\@mergedPartitionRefs, $dirPrefix."/mergedPartitionRefs.dmp", "mergedPartitionRefs", "bin");
    dump2disk(\%taxAbbr, $dirPrefix."/taxAbbr.dmp", "taxAbbr", "bin");
    dump2disk(\%taxAbbrExt, $dirPrefix."/taxAbbrExt.dmp", "taxAbbrExt", "bin");
    dump2disk(\%pDirCodeLookup, $dirPrefix."/pDirCodeLookup.dmp", "pDirCodeLookup", "bin");
    dump2disk(\%pDirNickLookup, $dirPrefix."/pDirNickLookup.dmp", "pDirNickLookup", "bin");
    dump2disk(\@valid_pDirs, $dirPrefix."/valid_pDirs.dmp", "valid_pDirs", "bin");
#    dump2disk($codon2char, $dirPrefix."/codon2char.dmp", "codon2char", "bin");
#    dump2disk($char2codon, $dirPrefix."/char2codon.dmp", "char2codon", "bin");
#    dump2disk(\%partitionHostCoords, $dirPrefix."/partitionHostCoords.dmp", "partitionHostCoords", "bin");
#    dump2disk(\@partitionHostStartCoords, $dirPrefix."/partitionHostStartCoords.dmp", "partitionHostStartCoords", "bin");
#    dump2disk(\@partitionHostEndCoords,   $dirPrefix."/partitionHostEndCoords.dmp", "partitionHostEndCoords", "bin");
#    dump2disk(\%prefix2partitionLookup, $dirPrefix."/prefix2partitionLookup.dmp", "prefix2partitionLookup", "bin");
#    dump2disk(\%contigHostPartitionLookup, $dirPrefix."/contigHostPartitionLookup.dmp","contigHostPartitionLookup","bin");
#    dump2disk(\@charKeys, $dirPrefix."/charKeys.dmp","charKeys","bin");
#    dump2disk($char2posLookup, $dirPrefix."/char2posLookup.dmp", "char2posLookup", "bin");
#    dump2disk($pos2charLookup, $dirPrefix."/pos2charLookup.dmp", "pos2charLookup", "bin");
    dump2disk(\%organisms, $dirPrefix."/organisms.dmp", "organisms", "bin");
    dump2disk(\%uniqueFragsByOrg, $dirPrefix."/uniqueFragsByOrg.dmp", "uniqueFragsByOrg", "bin");
    #dump2disk(\@kmers_eBV_hex, $dirPrefix."/kmers_eBV_hex", "kmers_eBV_hex", "bin");   

    open my $CONFIG, '>', $configfile;
    # CONSTANTS and USER OPTIONS
    print $CONFIG "[CONSTANTS / USER OPTIONS]\n";
    print $CONFIG "totalReplicons=$totalReplicons\n"
                 ."totalContigs=$totalContigs\n"
                 ."numOrgs=$numOrgs\n"
#                 ."pDir=$pDir\n"
                 ."outdir=$outdir\n"
                 ."storableDir=$storableDir\n"
                 ."seqDir=$seqDir\n"
                 ."updateDir=$updateDir\n"
                 ."dataDir=$dataDir\n"
                 ."dirPrefix=$dirPrefix\n"
                 ."dmpExt=$dmpExt\n"
                 ."dblabel=$dbLabel\n"
#                 ."skipNonStdNucs=$skipNonStdNucs\n"
                 ."minUniqueFragLength=$minUniqueFragLength\n"
#                 ."wantErrorFreeSeqs=$wantErrorFreeSeqs\n"
#                 ."optimize=$optimize\n"
#                 ."minStringPackLen=$minStringPackLen\n"
#                 ."maxCharKeys=$maxCharKeys\n"
                 ."partitionDir=$partitionDir\n"
#                 ."partitionExt=$partitionExt\n"
                 ."prefixLength=$prefixLength\n"
#                 ."partitions=$numPartitions\n"
                 ."padLed=$padLen\n"
                 ."taxLevel=$taxLevel\n"
                 ."compareSameTaxLevel=$compareSameTaxLevel\n"
                 ."currWordSize=$k\n";
    print $CONFIG "\n";
    print $CONFIG "[DATASTRUCTURE LOCATIONS]\n";
    print $CONFIG "ext=dmp\n"
                 ."dir=$dirPrefix\n";
    print $CONFIG "[DATASTRUCTURES]\n";
    print $CONFIG "orgLookupByContig\n"
                 ."fileNickLookup\n"
                 ."rFileNickLookup\n"
                 ."contigNickLookup\n"
                 ."rContigNickLookup\n"
                 ."contigID2filename\n"
                 ."contigIDs\n"
                 ."contigFilenames\n"
                 ."processedContigs\n"
                 ."processedOrgs\n"
#                 ."fastaFilesOfUniques\n"
                 ."orgNames\n"
#                 ."partitionHostCoords\n"
#                 ."partitionHostStartCoords\n"
#                 ."partitionHostEndCoords\n"
#                 ."contigHostPartitionLookup\n"
#                 ."prefix2partitionLookup\n"
#                 ."charKeys\n"
#                 ."char2posLookup\n"
#                 ."pos2charLookup\n"
#                 ."mergedIndicesRefs\n"
#                 ."mergedPartitionRefs\n"
                 ."taxAbbr\n"
                 ."taxAbbrExt\n"
                 ."pDirCodeLookup\n"
                 ."pDirNickLookup\n"
                 ."valid_pDirs\n"
#                 ."codon2char\n"
#                 ."char2codon\n"
                 ."organisms\n"
#                 ."lmerCoords\n"
                 ."uniqueFragsByOrg\n"
                 ."seqBM\n"
                 ."kmers_eBV_hex\n"
                 ."kmer_starts\n"
                 ."L2BM_\n"
#                 .$seqBMfilename."\n"
#                 .$kmers_eBV_hex_Filename."\n"
#                 .$kmer_starts_Filename."\n"
#                 .$uniqueSeqBMfilename."\n"
                 ;

    close $CONFIG;

    my $time1 = new Benchmark;
    my $timeText = timestr(timediff($time1,$time0));
    $timeText = $1 if($timeText =~ m/^\s*(\d+)\s+/);
    $timeText = $1 if($timeText =~ m/(\d+\.?\d+)\s+/);
    stat_log("Done in [$timeText] wallsecs.");    
    
    return;
}
################################################################################
# ARGS: $string
#-------------------------------------------------------------------------------
sub pctGC {

    my $seq = shift;

    # Delete all G's and C's, and determine %GC
    my $len = length($seq);            # Get length of expanded string before deleting
    
    return 0 if($len == 0);    
    
    my $count = ($seq =~ tr/ATN//c);
    return ($count/$len*100);

#    return (1-length($seq)/$len)*100;    # 3/100-sec slower than the above bitch

}
################################################################################
# ARGS: $filename, $description
################################################################################
sub loadFromDisk {

    die "Attempted to load non-existent file \"".$_[0]."\" from disk. Abort\n" if(!-e $_[0]);

    stat_log_("Loading ".$_[1]." from disk [".$_[0]."]...");
    my $iter0 = new Benchmark;
    my $href = retrieve $_[0];
    my $iter1 = new Benchmark;
    my $timeText = timestr(timediff($iter1, $iter0));
    $timeText = $1 if($timeText =~ m/^\s*(\d+)\s+/);
    $timeText = $1 if($timeText =~ m/(\d+\.?\d+)\s+/);
    _stat_log("done. [$timeText] wallsecs");
    
    return $href;
}
################################################################################
#        $_[0]    $_[1]        $_[2]     $_[3]
# ARGS: \%hash, $filename, $description, $mode      $mode = "bin" or "ascii"
################################################################################
sub dump2disk {

    STDOUT->autoflush(1);

    my $time0 = new Benchmark;

    if($_[3] =~ m/^bin/i) {
        #print "-> Storing datastructure ".$_[2]." to disk in BINARY format as \"".$_[1]."\"...";
        stat_log_("Storing datastructure ".$_[2]." to disk (BINARY) as \"".$_[1]."\"...");
        store clone($_[0]), $_[1];
    }
    elsif($_[3] =~ m/^ascii/i) {
        #print "-> Storing datastructure ".$_[2]." to disk in ASCII (text) format as \"".$_[1]."\"...";
        stat_log_("Storing datastructure ".$_[2]." to disk (ASCII) as \"".$_[1]."\"...");
        open my $OUTFILE, '>', $_[1];
        print $OUTFILE Dump clone($_[0]);
        close $OUTFILE;
    }
    else {
        #print "-> Attempting to store datastructure \"".$_[2]."\" to \"".$_[1]."\" in an unknown format: skipping...\n";
        stat_log("**Warning**: Attempting to store datastructure \"".$_[2]."\" to \"".$_[1]."\" in an unknown format: skipping...");
        return;
    }

    my $time1 = new Benchmark;
    #my $timeTxt = timestr(timediff($time1,$time0));
    #   $timeTxt = $1 if($timeTxt =~ m/^\s*(\d+)\s+/);
    
    my $timeTxt = q{};
    my $timediff = timestr(timediff($time1,$time0));
    if($timediff =~ m/^\s*(\S+)\s+/) {
        $benchRez = 5;
        $timeTxt = sciExpand($1, $benchRez);
    }
    elsif ($timediff =~ m/(\d+\.?\d+)\s+/) {
        $timeTxt = $1;
    }
    #my $loadTimeText = sciExpand($1, $benchRez) if(timestr(timediff($time1,$time0)) =~ m/^\s*(\S+)\s+/);
    _stat_log("done in [$timeTxt] wallsecs");
    return;
}
################################################################################
# ARGS: $filename
################################################################################
sub checkFile {
    if($_[0]) {
        die "Fatal: File \"".$_[0]."\" does not exist. Abort.\n" if(!-e $_[0]);
        return 1;                                                       # Given && Exists
    }
    else {
        return -1;                                                      # !Given
    }
}    
################################################################################
# Used to overcome Storable's limitations on storing shared hashes
# ------------------------------------------------------------------------------
# Usage:
#   store clone( \%hash ), "$0.bin";
#   my $retrieved :shared = shared_clone( retrieve "$0.bin" );
#
#-------------------------------------------------------------------------------
sub clone {
    my $ref = shift;
    ref( $ref ) or return $ref;
    ref( $ref ) eq 'SCALAR' 
        and return \(''.$$ref);
    ref( $ref ) eq 'HASH'   
        and return { map+( $_, clone( $ref->{$_} )), keys  %{ $ref } };
    ref( $ref ) eq 'ARRAY'  
        and return [ map       clone( $ref->[$_] ) , 0 .. $#{ $ref } ];
    die "\n";
}
################################################################################
# usage: lock($fh)
#-------------------------------------------------------------------------------
sub fileLock {
    my ($fh) = @_;
    flock($fh, LOCK_EX) or die "ERROR OBTAINING FILE LOCK! - $!\n";
    
    # seek to end in case some thread appended while waiting...
    seek($fh, 0, SEEK_END) or die "CANNOT SEEK TO END OF FILE! - $!\n";
}
################################################################################
# usage: unlock($fh)
#-------------------------------------------------------------------------------
sub fileUnlock {
        my ($fh) = @_;
        flock($fh, LOCK_UN) or die "ERROR RELEASING FILE LOCK! - $!\n";
}
################################################################################
# COMPLETE
#-------------------------------------------------------------------------------
sub init_file_write {
  my $fname = shift;
  open my $fh, '>', $fname || die "Fatal: cannot open \'$fname\' for writing.\n";
  return \*$fh;
}
################################################################################
# COMPLETE
#-------------------------------------------------------------------------------
sub init_file_append {
  my $fname = shift;
  open my $fh, '>>', $fname || die "Fatal: cannot open \'$fname\' for writing.\n";
  return \*$fh;
}
################################################################################
#-------------------------------------------------------------------------------
sub _init_file_read {
  my $fname = shift;
  open my $fh, '<', $fname || die "Fatal: cannot open \'$fname\' for reading.\n";
  return \*$fh;
}
################################################################################
# Display contents of global variable %warning_log, then resets it
#-------------------------------------------------------------------------------
sub displayWarnings {
    
    my $idx = 0;
    
    foreach my $warning_text (keys %warning_log) {
        $idx++;
        my $count = $warning_log{$warning_text};
        my $count_text = ($count == 1) ? ("")
                       :                 ("($count times)");
        my $notice     = ($count == 1) ? ("    ***A WARNING WAS") 
                       :                 ("    ***SEVERAL WARNINGS WERE");
        stat_log("$notice GENERATED $count_text***") if($idx == 1);
        stat_log("    $warning_text");
    }
    
    %warning_log = ();
    
    return;
}
################################################################################
# ARGS: $number, $sig
#-------------------------------------------------------------------------------
sub sciExpand {
    my $n = $_[0];
    return sprintf("%.${_[1]}f", $n) unless $n =~ /^(.*)e([-+]?)(.*)$/;
    my ($num, $sign, $exp) = ($1, $2, $3);
    #my $sig = $sign eq '-' ? ".".($exp-1+length($num)) : '';
    #return sprintf("%${sig}f", $n);
    return sprintf("%.${_[1]}f",$n);

    #my $n = $_[0];
    #return $n unless $n =~ /^(.*)e([-+]?)(.*)$/;
    #my ($num, $sign, $exp) = ($1, $2, $3);
    #my $sig = $sign eq '-' ? ".".($exp-1+length($num)) : '';
    #return sprintf("%${sig}f", $n);

}
################################################################################
#-------------------------------------------------------------------------------
sub removeDir {
    stat_log("Removing directory $_[0]...");
    eval { rmtree($_[0]) };
    if($@) {
        stat_log("Could not remove directory $_[0]: $@") if($@);
        return 0;
    }
    return;
}
################################################################################
#-------------------------------------------------------------------------------
sub createDir {
    eval { mkpath($_[0]) };
    die "Couldn't create ".$_[0].": $@\n" if($@);
    return;
}
################################################################################
#-------------------------------------------------------------------------------
sub pad_zeroes {
  my ($num, $opt_len) = @_;
  die "Cannot pad zeroes for less that the tens place!\n" if($opt_len < 2);
#  my $optimal_length = 3;   # equivlanet to the hundreds place
  $num =~ s/^(\d+)$/("0"x($opt_len - length$1)).$1/e;
  return $num;
}
################################################################################
# FUNCTION: Performs a binary search for ELEM in a **sorted ARRAY** (increasing), 
#           between positions START and END in ARRAY. There are 7 return statements. 
# CONDITIONS FOR RETURN:
# RETURNS:  1. END, if ELEM out-of-bounds(upper limit)      [req's sorted array]
#           2. START, if ELEM out-of-bounds(lower limit)    [req's sorted array]
#           3. LAST_POS ($p), if ELEM found in ARRAY
#
#               if ELEM < ARRAY[LAST_POS]:
#           4.          0 if LAST_POS is out-of-bounds(lower limit)
#           5.   LAST_POS if LAST_POS within bounds         [lower bounds]
#
#               if ELEM > ARRAY[LAST_POS]:
#           6.        END if LAST_POS >= END
#           7. LAST_POS+1 if LAST_POS within bounds         [upper bounds]
#         $_[0]    $_[1]   $_[2]    $_[3]
# ARGS: $element, \@array, $start, $end(i.e. $top)
# Note: $top = scalar(@array) - 1, which is the last position in the array
#-------------------------------------------------------------------------------
sub binary_search_overlap {
   my $bot = $_[2];     # = 0;                      to start from scratch;
   my $top = $_[3];     # = scalar(@{ $_[1] })-1;   to end at last element in array

    #print "Looking for ".$_[0]." in list from positions $bot (val=".$_[1]->[$bot].") to $top (val=".$_[1]->[$top].").\n";

    # ELEM is above ARRAY[END] (for sorted arrays)
    # So return highest position in array: $end
    # >>>> SHOULD I JUST RETURN A NEG NUMBER TO INDICATE OUT OF BOUNDS, SO DO NOT SEARCH FURTHER?? <<<<
    # If looking for MINIMUM value and ELEM is out-of-bounds (upper), then yes, don't proceed further.
    return $_[3] if($_[0] gt $_[1]->[ $_[3] ]); 
                                                   
    # ELEM is below ARRAY[START] (for sorted arrays);
    # So return lowest position in array: $start
    # >>>> SHOULD I JUST RETURN A NEG NUMBER TO INDICATE OUT OF BOUNDS, SO DO NOT SEARCH FURTHER?? <<<<
    # If looking for MAXIMUM value and ELEM is out-of-bounds (lower), then yes, don't proceed further.
    return $_[2] if($_[0] lt $_[1]->[ $_[2] ]);

    my $p = 0;
    while($top >= $bot)
    {

        # From: http:// googleresearch.blogspot.com/2006/06/extra-extra-read-all-
        #       about-it-nearly.html#!/2006/06/extra-extra-read-all-about-it-nearly.html
        # Most all binary search algorithms implement:
        #
        #       int mid = (low + high) / 2;
        #
        # This binary search fails if the sum of low (START) and high (END) is 
        # greater than the maximum positive int value (2^31 - 1). The sum overflows
        # to a negative value, and the value stays negative when divided by two.
        # In C this causes an array index out of bounds with unpredictable results.
        # This bug can manifest itself for arrays whose length (in elements) is 
        # 2^30 or greater (roughly a billion elements).
        #
        # Solutions:
        #
        # Perl:     int mid = low + ((high - low) / 2);
        #  or       int mid = (low + high) >> 1;
        # C/C++:    mid = ((unsigned int)low + (unsigned int)high)) >> 1;       
        #
        # Antoine Trux, Principal Member of Engineering Staff at Nokia Research
        # Center Finland for pointing out that the original proposed fix for C and 
        # C++, was not guaranteed to work by the relevant C99 standard 
        # (INTERNATIONAL STANDARD - ISO/IEC - 9899 - Second edition - 1999-12-01, 
        # 3.4.3.3), which says that if you add two signed quantities and get an 
        # overflow, the result is undefined. So the fix requires unsigned integers.

        $p = $top + $bot >> 1;               # add up and divide by 2
        
        #print "p=$p [VAL=".$_[1]->[$p]."]\tKEY=$_[0]\n";
        
        if ($_[1]->[$p] lt $_[0]) {              # ARRAY[p] < ELEM
            $bot = ++$p; 
        }
        elsif ($_[1]->[$p] gt $_[0]) {           # ARRAY[p] > ELEM
            $top = --$p; 
        }
        else {                                  # ARRAY[p] == ELEM
            return $p; 
        }
    }
 
    # Do one last comparison
    if($_[0] lt $_[1]->[$p]) {          # ELEM < nearest VAL @p, 
        return 0 if($p <= 0);           # ... return position 0 if p < 0
        return $p;                      # ... else return p
    }
    else {                              # ELEM > nearest VAL @p
        return $_[3] if($p >= $_[3]);   # ... return max position if p > max position
        return $p+1;                    # ... else return the next position above p
    }
}
################################################################################
# Given:    (1) Total number of datapoints
#           (2) Total number of CPUs available to process the datapoints
# Function: Given the size of the array (datapoints), the sub will determine the 
#           appropriate number of array indices to process per node, while
#           minimizing the difference in load between nodes. These indices will
#           be grouped into individual arrays -- where each array of indices will
#           be processed per CPU node -- and these arrays will be returned as a
#           single AREF.
# Returns:  AREF containing all the load-balanced arrays of indices.
# ARGS:     $totalPoints, $nThreads
#-------------------------------------------------------------------------------
sub loadBalanceStatic {
   
    my @jobs = ();
   
    my $n = int($_[0]/$_[1]); # no. of times $nThreads are fully loaded 
    my $r = $_[0] % $_[1];    # leftover datapoints

    my $unbalDatapoints = $n;
    my $balDatapoints   = $n+1;
    my $balNodes        = $r;
    my $unbalNodes      = $_[1] - $r;

    if($r == 0) {
        $balDatapoints   = $n;
        $balNodes        = $_[1]; 
        $unbalNodes      = 0;
        $unbalDatapoints = 0;
    }
    $unbalNodes = 0 if($unbalDatapoints == 0);

    # Returns the balanced indices of the input array
    #
    # // BALANCED
    for (my $idx = 0; $idx < $balNodes; $idx++) {
        my $startIdx = $idx*$balDatapoints;
        my $stopIdx  = $balDatapoints*($idx+1)-1;
        push(@jobs, [$startIdx..$stopIdx]);   # Each entry is a per-core, loadBalanced AREF of jobs
    }

    # // UNBALANCED
    my $balStop = $balNodes*$balDatapoints;

    for (my $idx = 0; $idx < $unbalNodes; $idx++) {
        my $startIdx = $idx*$unbalDatapoints + $balStop;
        my $stopIdx  = $unbalDatapoints*($idx+1)-1 + $balStop;
        push(@jobs, [$startIdx..$stopIdx]);
    }

    return \@jobs;
}
################################################################################
# PURPOSE:  We call this function in order to:
#           (1) minimize the repetitive loading of the largest of arrays; and
#           (2) optimally distribute the workload;
#
#           By loading the largest arrays
# FUNCTION: Given an array of indices where the indices correspond to an array
#           whose contents are sorted in order of decreasing size or complexity,
#           this function will load-balance these indices across the specified
#           number of threads so as to minimize the difference in total complexity
#           assigned to each thread. The data at the last index will not be
#           compared to any other datasets, so it is considered to be a no-op,
#           and is excluded from the load-balancing process.
#
#           Example: Eight contigs sorted in decreasing order of size from 0..7.
#                    Each contig idx must be compared with all other contigs of
#                    indices greater than itself.
#                                                                       Times Loading Data Element
#                 PAIRWISE_INTERACTIONS           INDEX     NO.JOBS     0   1   2   3   4   5   6   7
#           ----------------------------------    -----     -------    --- --- --- --- --- --- --- ---
#            0:1, 0:2, 0:3, 0:4, 0:5, 0:6, 0:7      0          7        1   1   1   1   1   1   1   1
#                 1:2, 1:3, 1:4, 1:5, 1:6, 1:7      1          6            1   1   1   1   1   1   1
#                      2:3, 2:4, 2:5, 2:6, 2:7      2          5                1   1   1   1   1   1
#                           3:4, 3:5, 3:6, 3:7      3          4                    1   1   1   1   1
#                                4:5, 4:6, 4:7      4          3                        1   1   1   1
#                                     5:6, 5:7      5          2                            1   1   1
#                                          6:7      6          1                                1   1
#                                                                      === === === === === === === ===
#                                                                       1   2   3   4   5   6   7   7
#            Thread0     Thread1     Thread2     Thread3
#            idx(jobs)   idx(jobs)   idx(jobs)   idx(jobs) 
#            -------     -------     -------     -------
#            0 (7)       1 (6)       2 (5)       3 (4)
#            6 (1)       5 (2)       4 (3)
#     Total
#     Jobs:    7           7           7           7
#
#
# RETURNS:  array of AREFs of indices
# ARGS:     \@array_of_indices, $nThreads
# ------------------------------------------------------------------------------
sub loadBalanceLinearPaired {
    my @loadBalanced = ();
    my $relativePos  = 0;
    my $flip         = -1;
    
    my $idx = 0;
    
    # We use scalar()-1 because the data at the last index will not be compared
    # to any other dataset, and is essentially a no-op, so it is excluded.
    #while($idx < scalar(@{ $_[0] })-1) {
    while($idx < scalar(@{ $_[0] })) {
        if($idx % $_[1] == 0) {
            $relativePos = 0;
            $flip = $flip*(-1);
        }
        if($flip == 1) {                 # $flip = 1
            push(@{ $loadBalanced[$relativePos] }, $idx);
        }
        else {                      # $flip = -1
            push(@{ $loadBalanced[$_[1] - $relativePos -1] }, $idx);
        }
        ++$relativePos;
        ++$idx;
    }
    return \@loadBalanced;
}
################################################################################
# Given:    STRING. LOG filename to open.
# Function: Initializes the log file and returns the filehandle
# Returns:  FILEHANDLE_REF. The open LOG filehandle
#------------------------------------------------------------------------------#
sub init_log() {                                                                                    # <----------------- misspelling. should be called like:    init_log($logfilename)
    my ($logfilename) = @_;

    open(my $LOGFILE, ">$logfilename") || die "Cannot write to log file \"$logfilename\": $!\n";
    return \*$LOGFILE;
}
################################################################################
# Given:    FILEHANDLE, STRING.
# Function: Print's STRING with a timestamp to a filehandle
# Returns:  Nothing.
#------------------------------------------------------------------------------#
sub print_timestamp() {
    my ($filehandle, $text_to_print) = @_;
    my ($sec, $min, $hour, $day, $month, $year, $wday) = &get_time();

    printf {$filehandle} ("$day_of_week[$wday] %02d/%02d/%04d at %02d:%02d:%02d",
            "    -->$text_to_print",$month+1, $day, $year+1900, $hour, $min, $sec);
}
################################################################################
# Given:    FH_REF, STRING. The string $phrase is written to the open filehandle,
#           FH_REF. FH_REF must have already been open, say, by sub init_log().
# Function: Print's the STRING to the FH_REF.
# Returns:  Nothing.
#------------------------------------------------------------------------------#
# REWRITE: sub can be re-written in terms of sub _begin_phase() and sub get_time()
sub stat_log {
#    my ($LOGFILE, $phrase) = @_;
    my ($phrase) = @_;
    my ($START_SEC, $START_MIN, $START_HOUR, $START_DAY,
        $START_MONTH, $START_YEAR, $START_WDAY          ) = &get_time();

    # Print to *STDOUT and $LOGFILE
    foreach my $filehandle (@reporting_filehandles) {
#        printf {$filehandle} ("$day_of_week[$START_WDAY] %02d/%02d/%04d (".
#            "%02d:%02d:%02d) -->$phrase\n",$START_MONTH+1, $START_DAY,
        printf {$filehandle} ("$day_of_week[$START_WDAY] %02d/%02d/%04d(".
            "%02d:%02d:%02d) $phrase\n",$START_MONTH+1, $START_DAY,
            $START_YEAR+1900, $START_HOUR, $START_MIN, $START_SEC);
    }
}
################################################################################
#  stat_log () will print a line: prefixed by date/time, suffixed by "\n"
#  stat_log_() will print a line: prefixed by date/time
# _stat_log () will print a line:                        suffixed by "\n"
#
# Usage: 
#   stat_log_("print ");
#   _stat_log("this\n");
#   stat_log("stuff\n");
#
# should print:
#   Wed 08/10/2011 (16:41:48) -->print this
#   Wed 08/10/2011 (16:41:48) -->stuff
#
#------------------------------------------------------------------------------#
sub stat_log_ {
#    my ($LOGFILE, $phrase) = @_;
    my ($phrase) = @_;
    my ($START_SEC, $START_MIN, $START_HOUR, $START_DAY,
        $START_MONTH, $START_YEAR, $START_WDAY          ) = &get_time();

    # Print to *STDOUT and $LOGFILE
    foreach my $filehandle (@reporting_filehandles) {
#        printf {$filehandle} ("$day_of_week[$START_WDAY] %02d/%02d/%04d (".
#            "%02d:%02d:%02d) -->$phrase",$START_MONTH+1, $START_DAY,
        printf {$filehandle} ("$day_of_week[$START_WDAY] %02d/%02d/%04d(".
            "%02d:%02d:%02d) $phrase",$START_MONTH+1, $START_DAY,
            $START_YEAR+1900, $START_HOUR, $START_MIN, $START_SEC);
    }
}
################################################################################
#------------------------------------------------------------------------------#
sub _stat_log_ {
#    my ($LOGFILE, $phrase) = @_;
    my ($phrase) = @_;
    my ($START_SEC, $START_MIN, $START_HOUR, $START_DAY,
        $START_MONTH, $START_YEAR, $START_WDAY          ) = &get_time();

    # Print to *STDOUT and $LOGFILE
    foreach my $filehandle (@reporting_filehandles) {
        printf {$filehandle} ("$phrase");
    }
}
################################################################################
#------------------------------------------------------------------------------#
sub _stat_log {
#    my ($LOGFILE, $phrase) = @_;
    my ($phrase) = @_;
    my ($START_SEC, $START_MIN, $START_HOUR, $START_DAY,
        $START_MONTH, $START_YEAR, $START_WDAY          ) = &get_time();

    # Print to *STDOUT and $LOGFILE
    foreach my $filehandle (@reporting_filehandles) {
        printf {$filehandle} ("$phrase\n");
    }
}
################################################################################
# Given:    LOG filehandle
# Function: Close log filehandle
# Returns:  Nothing.
#------------------------------------------------------------------------------#
sub close_log() {
    my ($LOGFILE) = @_;

    close ($LOGFILE) || die "Cannot close log file \"$logfilename\"\n";

}
################################################################################
# Given:    STRING. The name of a file to delete.
# Function: Delete an existing file from a previous run to prevent problems with
#           file write/update on the same filename.
# Returns:  Nothing.
#------------------------------------------------------------------------------#
sub clear_existing_output_file {

    my ($filename_to_delete) = @_;
    if (-e $filename_to_delete) {
        unlink ($filename_to_delete) or die "Cannot delete existing file ",
            "$filename_to_delete!\n";
    }
}
################################################################################
# Given:    Nothing.
# Function: Get time and date.
# Returns:  INTEGER x 7. INTEGERs represent the time and date;
#------------------------------------------------------------------------------#
sub get_time() {
    # Get current date & time
    my ($SEC, $MIN, $HOUR, $DAY, $MONTH, $YEAR, $WDAY) = (localtime)[0,1,2,3,4,5,6];
    # Return time and date
    return ($SEC, $MIN, $HOUR, $DAY, $MONTH, $YEAR, $WDAY);
}
################################################################################
# Given n things to choose from, choose r of them; Order matters; Repetition OK
# P(n,r) = n^r
sub permute_repOK {
    my ($choices, $replicates) = @_;
    return ($choices**$replicates);
}
################################################################################
# Given n things to choose from, choose r of them; Order matters; NO Repetition
# P(n,r) = n! / (n-r)!
sub permute_NOrep() {
    my ($choices, $replicates) = @_;
    my $numerator   = &fac($choices);
    my $denominator = &fac($choices - $replicates);
    return ($numerator/$denominator);
}
################################################################################
# Given n things to choose from, choose r of them; Order doesn't matter; Repetition OK
# C(n,r) = [r + (n-1)]! / [r!(n-1)!]
sub combine_repOK() {
    my ($choices, $replicates) = @_;
    my $numerator = &fac($replicates + $choices - 1);
    my $denominator = &fac($replicates) * &fac($choices - 1);
    return ($numerator/$denominator);
}
################################################################################
# Given n things to choose from, choose r of them; Order doesn't matter; NO Repetition
# C(n,r) = n! / [r!(n-r)!]
sub combine_NOrep() {
    my ($choices, $replicates) = @_;
    my $numerator = &fac($choices);
    my $denominator = &fac($replicates) * &fac($choices - $replicates);
    return ($numerator/$denominator);
}
################################################################################
# Returns the factorial of a number
# Note: 1! = fac(1) = 1; 
#       0! = fac(0) = 1; 
#      -4! = (-1)4*3*2*1 = -24 *NOT Implemented YET*
#     3.2! = *illegal*
sub fac() {
    my $num = shift;
#    my $sign = ($num < 0) ? (-1) : (1);
#    my $num = abs($num) if ($sign == -1);

    die "Illegal Operation!\n"
       ."Cannot take factorial of non-integer number \'$num\'\n" if (($num - int($num)) != 0);
    die "Illegal Operation!\n"
       ."Cannot take factorial of negative number \'$num\'\n" if ($num < 0);
    return ($num) if ($num < 2);
    return $num * &fac($num-1);
}
################################################################################
#-------------------------------------------------------------------------------
################################################################################
#-------------------------------------------------------------------------------
sub usage {
    my $help = <<ENDHELP;

    Usage: $0   REQUIRED    [OPTIONS]

    Objective: Given a set of genome sequences, generate a sequence database
               of DNA fragments that are unique amongst the genome sequences
               analyzed. Defaults are in square brackets [].

    REQUIRED
        --taxIdxFile=               Taxonomy Index file (in XML format), specifying all the input
                                    organisms and their sequence files.
        --fastaext=                 Extension found on FASTA sequence files.
        
        *******DEPRECATED*******
        --loadFasta=                For faster processing times, this overrides the 
                                    --orgsByFastaDir option and loads two binary hashes 
                                    from disk that contain all relevant sequence and
                                    taxonomy data.
        --noStoreFasta              Disables storing of binary hashes of sequence and
                                    taxonomy info to disk [false]
        *******DEPRECATED*******


        (discrete values)
        --wordSize=X,Y,Z            Explicitly specifies single or multiple l-mer lengths, each 
                                    separated by a comma (no spaces). 
                                    Disables both --minWordSize and --maxWordSize
                                    NOTE: wordSize >= minUniqueFragLength

        (range of values)
        --minWordSize=              MIN l-mer length; used with --maxWordSize and --stepWordSize [NOT YET IMPLEMENTED]
        --maxWordSize=              MAX l-mer length; used with --minWordSize and --stepWordSize [NOT YET IMPLEMENTED]
        --stepWordSize=             Used with --minWordSize and --maxWordSize, this value
                                    specifies the gradations between the two. If MIN=5, MAX=10,
                                    STEP=5, then words of length 5 and 10 are created. If
                                    STEP=3, then words of length 5, 8, and 10 are created. 
                                    [NOT YET IMPLEMENTED]
        --loadWords                 Directory containing the pre-computed words/l-mers
                                    for each replicon, encoded by their nicknames and value
                                    of l. E.g. "00001.100.dmp" (replicon 00001, l=100)

    [OPTIONS]
        --nThreads=                 Number of threads to use [1]

        INPUT RESTRICTIONS
       --------------------------
        --minUniqueFragLength=      Minimum length of unique gene fragments to report. This 
                                    only filters the reported nucmer output and is different 
                                    from the "nucmer -l" option. [20]
        --returnErrorFreeSeqs       Return only those unique gene fragments that do not contain
                                    any errors.
        TAXONOMY
       --------------------------
        --taxTreeFile=              Perl Storable file (denovoTaxTree.dmp) containing the taxonomic organization
                                    of all organisms in --orgsByFastaDir, indexed by the Species name. 
                                    To generate this file:
                                    1. ncbiTax2RankLevelTable.pl -name names.dmp -node nodes.dmp > outfile
                                    2. perl parseTax2b.pl outfile
        --taxIdxLength=             No. of characters to index the \%taxTree for organism identification.
                                    Higher numbers yields higher performance up to 8.
        --taxLevel=                 Taxanomic level to restrict comparison to. Allowed values:
                                    "species","genus","family","order","class","phylum", or "superkingdom"
        --compareSameTaxLevel       Compare genomes from the same taxnomic level (Species,
                                    Genus, etc.) when determining unique regions. [no]
        --doCompleteTax             Generate uniques for all taxonomic levels, starting from SPECIES and up to
                                    PHYLUM. Requires --taxLevel=species for this to be enabled. [no]
        --continueNoHostFound       Continue processing even if an input genome has no organism match found
                                    in the Taxonomy Tree file specified by --taxTreeFile [halt]

        MEMORY/CPU OPTIMIZATION
       --------------------------
        --optimize=                 Optimize run for either memory ("mem") or for CPU ("cpu").  
                                    If "mem", all data is retained in memory for shorter run time.
                                    If the run fails, specify "cpu" to swap RAM for disk. [cpu] 
                                    See --minStringPackLength, --prefixLength, and --numPartitions
                                    for additional configuration options.
        --minStringPackLength=      L-mers above this length (in bp) will be compressed by
                                    approximately 1/3 using the stringPack() function. [15]
                                    Note: the stringPack() function is disabled by default 
                                    if(minStringPackLen > wordSize).
        --prefixLength=             Index each replicon's l-mers by its first N characters
                                    (where N = prefixLength), thus creating N files (partitions)
                                    for each replicon. 
                                    If stringPack() is used, then (currently) 64^N partitions
                                    are created (at most)
                                        Ex. For N=1, 64 partitions are created, at most. 
                                            For N=2, 4906 partitions. 
                                    If stringPack() is *not* used (see --minStringPackLength), 
                                    then at most 4^N partitions are created. 
                                        Ex. For N=1, 4 partitions are created, at most. 
                                            For N=2, 16 partitions. 
                                    Specifying partitions > 1 swaps RAM for disk space. It
                                    effectively lowers the RAM requirement for the l-mer comparisons
                                    but increases CPU time because of disk access. [1] 
                                    See --numPartitions for additional options.
        --partitions=               Merges the individual partitions into these many "merged" partitions. 
                                    The merging is effectively load-balanced, where differences between 
                                    each partition are minimized. If stringPacking is enabled, there 
                                    are 69^N partitions, where N is specified by --partitions=<INTEGER>.
                                    If --skipNonStdNucs is enabled, there are only 64^N partitions -- the
                                    extra 5 come from: A/T/G/C/N. If stringPacking is *NOT* enabled,
                                    there are 5^N partitions. If --skipNonStdNucs is enabled, then there
                                    are only 4^N.
                                    If --partitions=1, then all partitions will be merged into a single
                                    file, effectively reducing the amount of disk access by using more
                                    RAM and thus decreasing run time. [1]
        STATISTICS
       --------------------------
        --skipGC                    Disables \%GC calculations on l-mers [0]

        OUTPUT
       --------------------------
        --pDir=                     Directory to where all partitions are written. Separate multiple
                                    directories with commas (no spaces). If giving multiple directories,
                                    they should be on separate drives with separate controllers. [./]
        --outdir=                   Directory to where all other results are written. [./]
        --dbLabel=                  Label to add to the output DB file. [Bacteria]
        --noStats                   Do not generate statistics file [generate file]
        --verbose                   Output verbosity [no]
        --vverbose                  Increased output verbosity [no]
        --help                      Prints this message
        --debug                     Debugging condition


ENDHELP
print "$help\n";
exit;
}
