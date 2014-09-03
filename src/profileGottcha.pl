#!/usr/bin/env perl
#
# /*                            
#     Program:      profileDTAR_v1.1.pl 
#     Author:       Tracey Freitas, Po-E Li
#     Created:      2011-09-19
#     Email:        traceyf@lanl.gov
#     Last Updated: 2013-01-15 (v2.10)
#                   2013-04-10 (v2.14)
#                   2013-04-11 (v2.14b)
#                   2013-04-11 (v2.14b)
#                   2014-05-08 (v2.15)
#
#     2015-05-08:   (1) SAM file can be inputted through STDIN.
#                   (2) Merging format checking and parsing to the same loop for SAM input
#
#     2014-04-20:   Default parsedDTAR file changes filename to *.parsedGOTTCHA.dmp
#
#     2013-03-14:   SAM file parsing taking alignment code (field #2) into account
#
#     2013-04-02:   (1) Changed all report file extensions from *.csv to *.tsv (I hate you, Matt)
#                   (2) Additional report file was generated specifically for REPLICONs (*.replicon.tsv)
#                   (3) STRAIN report (*.strain.csv) actually contains data on the entire strain, rather 
#                   than just individual replicons (all replicons)
#                   (4) Updated each report (*.tsv files) at given RANK:
#                   -> added "# of Subranks" (i.e. # SPECIES under the given GENUS)
#                   -> added "# of Genome Projects" (i.e. # of STRAINS under the given GENUS)
#                   -> added "FULL_REFDB_LENGTH" (the full length of the GI's supporting the entry)
#                   -> removed internally scaled parameters
#                   -> moved up "HIT_COUNT" and "TOTAL_BP_MAPPED" to a slot just after "LINEAR_COV"
#                   -> renamed "FOLD_COV_UNIQUE_SAMPLE" to "LINEAR_DOC", where DOC = depth-of-coverage
#                   -> renamed "FOLD_COV_UNIQUE_REFDB"  to "UREF_DOC", where UREF = unique reference
#                   -> renamed "REF_CMAX" to "UREF_CMAX"
#                   -> added all the above parameters for the BEST SUBRANK, as determined by its LINEAR_COV
#                   (5) Undef'd as many of the large data structures as early as possible to reduce max RAM usage
#
#     2013-04-03:   (1) added contig stats (count, mean, stdev, min, max)
#                   (2) added histogram string
#                   (3) require that each giFrag entry has the start, stop, and length fields,
#                       as in "gi|XXXXXXX|xx|XXXXXXXX|#########|#########|#########|"
#                       This holds for the contents of GOTTCHA's *.fna and *.parsedDTAR.dmp files, the SAM files, etc.
#
#     2013-04-10    (1) Added new output file types
#                   -> <PREFIX>.replicon.reads.coords.csv
#                   -> <PREFIX>.replicon.contig.HistoByEntry.csv
#                   -> <PREFIX>.replicon.contig.HistoByGI.csv
#                   -> <PREFIX>.strain.contig.coords.tsv
#                   -> <PREFIX>.MAPPED.trimmed.fastq        [reads that mapped to >=1 GOTTCHA entry]
#                   -> <PREFIX>.UNMAPPED.trimmed.fastq      [reads that did NOT map to >=1 GOTTCHA entry]
#                   (2) Implemented ability to capture full-read hit counts (FULL_HIT_COUNT)
#
#     2013-04-16    (1) Added N50, N60, N70, N80, N90, N95 stats
#
#     2013-06-09    Modified output file types:
#                   -> <PREFIX>.MAPPED.trimmed.fastq has been expanded to include the organism's taxonomic name
#                      at the level in which GOTTCHA is being run at. For example:
#                      Streptococcus pyogenes NZ131 will be:
#                      <PREFIX>.Firmicutes.MAPPED.trimmed.fastq                 [Phylum-GOTTCHA]
#                      <PREFIX>.Bacilli.MAPPED.trimmed.fastq                    [Class-GOTTCHA]
#                      <PREFIX>.Lactobacillales.MAPPED.trimmed.fastq            [Order-GOTTCHA]
#                      <PREFIX>.Streptococcaceae.MAPPED.trimmed.fastq           [Family-GOTTCHA]
#                      <PREFIX>.Streptococcus.MAPPED.trimmed.fastq              [Genus-GOTTCHA]
#                      <PREFIX>.Streptococcus_pyogenes.MAPPED.trimmed.fastq     [Species-GOTTCHA]
#                      <PREFIX>.Streptococcus_pyogenes_NZ131.trimmed.fastq      [Strain-GOTTCHA]
#
#
#
#
# $1"|"$2"|"$3"|"$4"|"$5"|"$6"|"$7"|"$8"|"$9"|"$10"|"$11"|"$12"|"$13"|"$14"|"$15"|"$16"|"$17"|"$18"|"$19"|"$20"|"$21"|"$22"|"$23"|"$24"|"$25"|"$26"|"$27"|"$28"|"$29"|"$30"|"$31"|"$32"|"$33"|"$34"|"$35"|"$36"|"$37"|"$38
# 
# /*
#        
#
# =============================================================================
#
#
# Given a database of unique fragments, the "%non-overlapping coverage" metric 
# is used to quantify just how much of these unique fragments are present in a
# given sample. 
#
# This program calculates the %non-overlapping coverage of reference sequences, 
# organized by GI/strain, given a reference FASTA DB file and a mapping results 
# file (nucmer, bwa, or blast tabular format). The mapping file may be either 
# reads or contigs whose lengths are greater than or equal to the minimum unique
# fragment length reported in the uniques DB used (usually, this is 30-35 bp).
#
#
#
# //////
# BUGS:
# //////
#   1. In sub parseSAM(): # Omit 3'-overhanging reads
#      |--> this obviates the use of local alignment! Comment out and profile results
#           to show that this line really rejects local alignments with 3'-overhanging reads
#
#
# ======
# TO DO:
# ======
#
#
# A. use updated trimmed START:STOP positions in FASTQ entries to track: 
#   (1) no. of full-length reads, and
#   (2) the average fraction of reads mapped (see "TO DO" section in splitrim_v1j3d.d)
#
# B. export FASTQ headers (trimmed and/or full-reads... haven't decided yet) that have 
#    a mapping entry in the given SAM file into a "MAPPED.fastqHeaders.txt" file, and
#    organized/indexed by GI.
#    ... need a separate script that will parse the "MAPPED.fastqHeaders.txt" file along 
#        with the "speciesTreeGI.dmp" and "genomeVitals.dmp" files to organize the 
#        FASTQ headers by STRAIN, SPECIES, GENUS, FAMILY, ORDER, CLASS, and PHYLUM.
#
# C. 
#
# -Clean up the *.coordDetails.tsv file for easy parsing (perhaps make a Storable
#  version for fast loading?)
#
#
#
# -Added "Observed Abundance" per replicon. 
#
#                                         TOTAL_BP_MAPPED
#                            -----------------------------------------
#                            (UNIQUE_REFDB_LENGTH / FULL_REFDB_LENGTH)
#   Observed Abundance = -------------------------------------------------
#                                        TOTAL_INPUT_BASES
#
# 
#                             (TOTAL_BP_MAPPED * FULL_REFDB_LENGTH) 
#                      = -----------------------------------------------
#                           (UNIQUE_REFDB_LENGTH * TOTAL_INPUT_BASES)
#
#  where FULL_REFDB_LENGTH is obtained from genomeVitals file: $genomeVitals_HREF->{GI}->{$gi}->{SIZE}
#
# -Enable one thread per SAM file (only if headers are trimmed)
# -Enable multiple threads per SAM file 
#
# ---------
# Caveats: 
# ---------
# 1. Mapping results to a FAMILY-level uniques DB are *NOT* the same as the mapping 
#    results to a SPECIES-level uniques DB that has been rolled up to the FAMILY level. 
# 2. HIT_COUNT is only meaningful when mapping reads to the unqiues DB. It is 
#    relatively meaningless if mapping contigs.
#
#
# REQUIREMENTS:
#   --genomeVitals=<genomeVitals.dmp>      Created with "getGenomeVitals_v1.2.pl".
#                                          Parses *.gbk files and NCBI's "gi_taxid_nucl.dmp"
#                                          for vital data: GI, TAXID, ACC, ORG, REPL, SIZE, DATE.
#
#   --treeFile=<speciesTreeGI.dmp>         Created in "dtar_editor_v1.21b.pl".
#                                          Parses NCBI's "names.dmp" and "nodes.dmp" files to create
#                                          a hash (indexed by TAXID) representing each organism's
#                                          "complete" taxonomic tree.
#
#   --db=<DB.fna>                          MultiFASTA "uniques" database of fragments whose entries
#                                          follow the format: gi|######|XXX|XXXXXX|start|stop|
#
#   --map=<MAPFILE>                        Either a SAM- or NUCMER-formatted mapping file.
#
# ======
# To Do: 
# ======
#   1. Modify Output:
#      STRAIN | GI | NON-OVERLAPPING_MAPPED_LENGTH | UNIQUE_DB_LENGTH | %COVERAGE | 
#      REPLICON_ORIGINAL_LENGTH | %REDUCTION | No._of_UNIQUE_FRAGS ("|"-separated list of integers)
#
#   2. Incorporate Depth of Coverage (mean fold coverage).
#
#   3. add user option to restrict qualifying HITS to completely mapped contigs/reads
#      and/or completely mapped unique fragments.
#
#   4. Implement --blast input format options
#
#   5. Multithreaded benefits: prefetch the (i)  $dbFile  (sub parseMultiFASTA) 
#                                       and (ii) $mapFile (sub parseXXXXXX)
#
#   6. User specifies the rank level that the uniques DB was created at, so that
#      inappropriate ROLLUPS are avoided.
#      Example: a FAMILY-level uniques DB will only map reads to the FAMILY level
#               so STRAIN-, SPECIES-, and GENUS-level rollups should be omitted.
#
#   7. Add --parsedMAP=<parsedMAP.dmp> to quickly load a parsed MAP file (for subsequent analyses)
#      -> must also add in a storeMap() function
#
use strict;
use warnings;
use Getopt::Long;
use Storable;
use YAML;
use Benchmark;
use File::Path;
use IO::Handle;
#use feature 'switch';

my $MINPCTID=98;

my $dbFile     = q{};
my $parsedDTAR = q{};
my $make_dmp   = 0;
my $mapFile = q{-};
my $nucmer  = 0;
my $blast   = 0;
my $sam     = 0;
my $minPctID= $MINPCTID;
my $genomeVitalsFile = q{};
my $treeFile        = q{};
my $trimStats       = q{};
my @trimStatsFiles  = ();
my $inputReads      = 0;          # Total no. of input reads (from trimming)
my $inputBases      = 0;          # Total no. of input bases (from trimming)
my $detailFile      = q{};
my $noVitalsFile    = 0;
my $outdir          = ".";
my $prefix          = q{};
my $noMappedFastq   = 0;
my $noUnmappedFastq = 0;
my $noFastqOut      = 0;
#---------------------------
my $findTopHit       = 0;
my $method           = 0;
my %methods          = (1 => q{}, 2 => q{}, 3 => q{});
my @wantedMethods    = ();      # populated in sub checkInput()
my $fieldSep         = "\t";
my $minCoverage      = 0.005; #0.038 / 0.005 / 0.001
my $minHits          = 5;
my $minLength        = 100;
my $minMLHL          = 10;
my $criticalCoverage = 0.006; # 0.5;
my %abundanceOptions = ();
#---------------------------
my $help    = 0;
my $verbose = 0;
my $debug   = 0;
my $ddebug  = 0;
my $lineBuffer = 1000000;
my $countdown  = 10;            # Time delay when map format believed to be inconsistent w/selection
my @rankAbbr   = ("SS","S","G","F","O","C","P","SK");
my @validRanks = ("strain","species","genus","family","order","class","phylum");
my %validRanks = ();
   @validRanks{@validRanks} = ();
my $rankAbbr = {strain  => "SS",
                species => "S",
                genus   => "G",
                family  => "F",
                order   => "O",
                class   => "C",
                phylum  => "P",
                superkingdom => "SK"};
my %revRankAbbr = ();
while(my ($k, $v) = each %{ $rankAbbr }) {
    $revRankAbbr{$v} = $k;
}
my @hitTreeRanks = ("GI","SS","S","G","F","O","C","P","SK");

GetOptions(
    #=========================================
    # Input Options: Database
    #=========================================
    "db=s"              => \$dbFile,
    "make_dmp"          => \$make_dmp,
    "parsedDB=s"        => \$parsedDTAR,
    #-----------------------------------------
    # Input Options: Map & Reference Files
    #-----------------------------------------
    "map=s"             => \$mapFile,
    "genomeVitals=s"    => \$genomeVitalsFile,  # GI-linked data (genomeVitals.dmp)
    "treeFile=s"        => \$treeFile,          # contains tax tree with GI
    "trimStats=s"       => \$trimStats,         # Comma-separated *.stats.txt files containing "No. of READS" and "No. of BASES"
    #-----------------------------------------
    # Mapping file type:
    #-----------------------------------------
    "nucmer"            => \$nucmer,            # format: using nucmer COORDS file
    "blast"             => \$blast,             #         using BLAST tabular format
    "sam"               => \$sam,               #         using BWA format
    #-----------------------------------------
    # Selection criteria:
    #-----------------------------------------
    "minID"             => \$minPctID,          # Minimum %ID to filter out
    #=========================================
    # Output Options
    #=========================================
    "outdir=s"          => \$outdir,
    "prefix=s"          => \$prefix,
    "noUnmappedFastq"   => \$noUnmappedFastq,   # if TRUE, suppresses writing of umapped reads to disk
    "noMappedFastq"     => \$noMappedFastq,     # if TRUE, suppresses writing of mapped reads to disk
    "noFastqOut"        => \$noFastqOut,        # if TRUE, suppresses writing of both mapped and unmapped reads to disk
    #=========================================
    # Post-Processing: Top Hit
    #=========================================
    "topHit"            => \$findTopHit,            # if TRUE, returns data on the ONE most probable organism in sample
    #=========================================
    # Post-Processing: Abundance 
    #=========================================
    "method=s"          => \$method,           # 1 or 2 or 3 or 1,2 or 1,3 or 2,3 or 1,2,3
    "field=s"           => \$fieldSep,         # field separator [default: tab, "\t"]
    #-----------------------------------------
    # Abundance Filter 1: Coverage-based
    #-----------------------------------------
    # Criteria: PASS if((best_LINEAR_COV > minCov) 
    #                && (best_HIT_COUNT  > minHits) 
    #                && (best_LINEAR_LEN > minLength))
    "minCov=f"          => \$minCoverage,      # FILTER1a: minimum best_LINEAR_COVERAGE
    "minHits=i"         => \$minHits,          # FILTER1b: minimum best_HIT_COUNT
    "minLen=i"          => \$minLength,        # FILTER1c: minimum best_LINEAR_LENGTH
    #-----------------------------------------
    # Abundance Filter 2: Stacking-based (static)
    #-----------------------------------------
    # Criteria: PASS if((MLHL < minMLHL) 
    #                && (best_LINEAR_COV > cCov))
    "minMLHL=f"         => \$minMLHL,          # FILTER2a: mean Linear Hit Length
    "cCov=f"            => \$criticalCoverage, # FILTER2b: critical coverage
    #-----------------------------------------
    # Abundance Filter 3: Stacking-based (dynamic)
    #-----------------------------------------
    # Plot (best_LINEAR_COV) vs. (best_HIT_COUNT) vs. (best_MLHL)
    # ...
    #=========================================
    # Miscellaneous Options
    #=========================================
    "help"              => \$help,
    "verbose"           => \$verbose,
    "debug"             => \$debug,
    "ddebug"            => \$ddebug,
);

my $iter_timeX = new Benchmark;

# ----------------------------------------
# Error checks
# ----------------------------------------
checkInput();

my %giHisto = ();           # KEY: $gi      VAL: $contigLength => $count
my %giFragHisto = ();       # KEY: $giFrag  VAL: $contigLength => $count

# ----------------------------------------------------------------------------------
# 1. Calculate total unique length, per GI, from the uniques DB
# ----------------------------------------------------------------------------------
# KEY: gi (from the uniques DB)
# VAL: total length of GI in uniques DB
my $uniqueGIlengths = $dbFile ? calculateUniqueLengthsByGI($dbFile)
                    :           retrieveParsedDB($parsedDTAR);

# ----------------------------------------------------------------------------------
# 2. Parse $mapFile
# ----------------------------------------------------------------------------------
my $gi2mapFrags = parseMapFile($mapFile);

# ----------------------------------------------------------------------------------
# 3. Merge hits for each entry in the uniques DB
# ----------------------------------------------------------------------------------
#   Holds START|STOP positions all uniquely mapped "unique" fragments
my $gi2mergedFrags = mergeOverlappingHits2($gi2mapFrags);    # Contains only gi's from the map file

# ----------------------------------------------------------------------------------
# 4. Associate GI with Strain name using "genomeVitals.dmp"
# ----------------------------------------------------------------------------------
# KEY1: GI->{ORG,REPL,TAXID,SIZE,DATE,ACC} (from all bacterial genbank files; created 2011-09-27)
# KEY2: ORG->{$org}->{REPL}->{$repl}->{ACC,GI,DATE,SIZE}
# KEY3: REPL->{ORG,ACC,GI,DATE,SIZE}
# KEY4: ACC->{$acc}->{$gi}
# KEY5: TAXID->{$taxid}->{$gi}
my $genomeVitals = q{};
if(!$noVitalsFile) {
    my $vitalsIter0 = new Benchmark;
    print "->Retrieving Genome Vitals from disk [".$genomeVitalsFile."]...";
    $genomeVitals = retrieve $genomeVitalsFile if(!$noVitalsFile);  
    my $vitalsIter1 = new Benchmark;
    print "done. ".timestr(timediff($vitalsIter1,$vitalsIter0))."\n";
}

# ----------------------------------------------------------------------------------
# 5. Count up lengths per GI
# ----------------------------------------------------------------------------------
# KEY: gi (from map file)
# VAL: total non-overlapping length mapped to GI
my $giLengths = sumMappedLengths($gi2mergedFrags, $genomeVitals, $gi2mapFrags, \%giFragHisto, \%giHisto);
                                                            
# ----------------------------------------------------------------------------------
# 6. Calculate %non-overlapped coverage
# ----------------------------------------------------------------------------------
#                                                        LINEAR_LENGTHs   UNIQUE_REFDB_LINEAR_LENGTHs
my $nonOverlappedCoverage = calculateNonOverlappedCoverage($giLengths,     $uniqueGIlengths,           $genomeVitals);

# Free the RAM-sucking variables!
undef $giLengths;
undef $uniqueGIlengths;

my $opts = {
             DBFILE              => $dbFile,
             PARSEDDTAR          => $parsedDTAR,
             #UGILENGTHS_HREF     => $uniqueGIlengths,
             MAPFILE             => $mapFile,
             GI2MAPFRAGS_HREF    => $gi2mapFrags,
             GI2MERGEDFRAGS_HREF => $gi2mergedFrags,
             #GILENGTHS_HREF      => $giLengths,
             NOVITALS            => $noVitalsFile,
             GVITALSFILE         => $genomeVitalsFile,
             GVITALS_HREF        => $genomeVitals,
             N_O_COV_HREF        => $nonOverlappedCoverage,
             TREEFILE            => $treeFile,
             INPUTREADS          => $inputReads,
             INPUTBASES          => $inputBases,
           };

# ----------------------------------------------------------------------------------
# 7. Attach tax tree data to %nonOverlappedCoverage
# ----------------------------------------------------------------------------------
#    **Only do this for a STRAIN-level Uniques DB**
#computeRepliconData($nonOverlappedCoverage, $gi2mergedFrags);
computeRepliconData($opts);

# ----------------------------------------------------------------------------------
# 8. Attach tax tree data to %nonOverlappedCoverage and calculate taxonomic rollups.
# ----------------------------------------------------------------------------------
#computeTaxRollups($nonOverlapedCoverage, $treeFile, $genomeVitals, $gi2mergedFrags);
computeTaxRollups($opts);

# ----------------------------------------------------------------------------------
# 9. Calculate Relative Abundances
# ----------------------------------------------------------------------------------
calcAbundance(\%abundanceOptions);

# ----------------------------------------------------------------------------------
# 10. Find Top Hit
# ----------------------------------------------------------------------------------
findTopHit() if($findTopHit);

my $iter_timeY = new Benchmark;

print "------------------\n";
print "TOTAL SCRIPT TIME:\t$mapFile\t".timestr(timediff($iter_timeY, $iter_timeX))."\n";
print "------------------\n";

################################################################################
##############################  SUBROUTINES ####################################
################################################################################
################################################################################
# ARGS: $string, \%results
sub getIdxOfHeader {
    
    my $wantedIdx   = 0;
    my $found       = 0;
    
    foreach my $headerIdx (keys %{ $_[1] }) {
        if($_[1]->{$headerIdx}->{NAME} eq $_[0]) {
            $wantedIdx = $headerIdx;
            $found = 1;
            last;
        }
    }
    die "*FATAL*: \"LINEAR_DOC\" could not be found!\n" if(!$found);

    return $wantedIdx;
}
################################################################################
# ARGS: \%abundanceOptions
#    $abundanceOptions{TOPHIT}     = $findTopHit;        # bool
#    $abundanceOptions{ABU_METHOD} = $method;            # Comma-delimited string
#    $abundanceOptions{FIELD_SEP}  = $fieldSep;          # "\t"
#    $abundanceOptions{MIN_COV}    = $minCoverage;       # float
#    $abundanceOptions{MIN_HITS}   = $minHits;           # int
#    $abundanceOptions{MIN_LEN}    = $minLength;         # ulong
#    $abundanceOptions{MIN_MLHL}   = $minMLHL;           # float
#    $abundanceOptions{CCOV}       = $criticalCoverage;  # float
################################################################################
sub calcAbundance {

    RANK: foreach my $cRank (@validRanks) {         # process in order
        my $gottchaTable = $prefix 
                         ? ($outdir."/".$prefix.".".$cRank.".tsv") 
                         : ($outdir."/".$cRank.".tsv");
        if(!-e $gottchaTable) {
            print "*WARNING*: GOTTCHA table \"".$gottchaTable."\" is not found! Skipping...";
            next RANK;
        }

        #print Dump(\%abundanceOptions)."\n";

        # 1. Parse GOTTCHA table
        my %gottchaResults = %{ parseTable($gottchaTable, $fieldSep) };              # what if the field separator differs

        # 2. Select Abundance Method
        foreach (@wantedMethods) {
            abundanceMethod1(\%gottchaResults, $_[0]) if $_==1;
            #calcAbundance1(\%fullResults) if(%fullResults);
            abundanceMethod2(\%gottchaResults, $_[0]) if $_==2;
            #calcAbundance2(\%fullResults) if(%fullResults); 
            abundanceMethod3(\%gottchaResults, $_[0]) if $_==3;
            #calcAbundance3(\%fullResults) if(%fullResults); 
        } #foreach(@wantedMethods)

        exportUpdatedTable(\%gottchaResults, $gottchaTable.".ABU");

    } #RANK
}
################################################################################
# ARGS: \%results, \%abundanceOptions
#    $abundanceOptions{TOPHIT}     = $findTopHit;        # bool
#    $abundanceOptions{ABU_METHOD} = $method;            # Comma-delimited string
#    $abundanceOptions{FIELD_SEP}  = $fieldSep;          # "\t"
#    $abundanceOptions{MIN_COV}    = $minCoverage;       # float
#    $abundanceOptions{MIN_HITS}   = $minHits;           # int
#    $abundanceOptions{MIN_LEN}    = $minLength;         # ulong
#    $abundanceOptions{MIN_MLHL}   = $minMLHL;           # float
#    $abundanceOptions{CCOV}       = $criticalCoverage;  # float
################################################################################
sub abundanceMethod1 {

    my %abuOpts = %{ $_[1] };

    # Get last number in hash; increment as new column for RelAbundance1
    my $nextIdx      = (sort {$a <=> $b} keys %{ $_[0] })[-1]+1;
    my $totalEntries = (sort {$a <=> $b} keys %{ $_[0]->{0}->{VAL} })[-1]+1;
    
    #print "NextIdx      = $nextIdx\n";
    #print "TotalEntries = $totalEntries\n";
    
    $_[0]->{$nextIdx}->{NAME} = "REL_ABUNDANCE";
    
    # Find the LINEAR_DOC column header
    my $linearDOCidx = getIdxOfHeader("LINEAR_DOC", $_[0]);

    # Add FILTER here
    #     LINEAR_LENGTH   > 100
    #     HIT_COUNT       > 10
    #     best_LINEAR_COV > 0.005
    # ...
    
    my @qualIdx = ();          # Qualifying indices
    my $relAbundanceSum = 0;
    for (0..$totalEntries-1) {
        $_[0]->{$nextIdx}->{VAL}->{$_} = 0;                 # Initialize

        # FILTER
        if (
                     ($_[0]->{23}->{VAL}->{$_} < $abuOpts{MIN_COV})     # best_LINEAR_COV
                  || ($_[0]->{24}->{VAL}->{$_} < $abuOpts{MIN_HITS})    # best_HIT_COUNT
                  || ($_[0]->{20}->{VAL}->{$_} < $abuOpts{MIN_LEN})     # best_LINEAR_LENGTH
                ){
        	$_[0]->{$nextIdx}->{VAL}->{$_} = "filtered";
                next;
        }

        # A low MEAN_LINEAR_HIT_LENGTH with a low coverage suggests background or unknown bleed-through
        if(
                ($_[0]->{33}->{VAL}->{$_} < $abuOpts{MIN_MLHL})         # best_MEAN_LINEAR_HIT_LENGTH
             && ($_[0]->{23}->{VAL}->{$_} < $abuOpts{CCOV})             # best_LINEAR_COV
               ){
		$_[0]->{$nextIdx}->{VAL}->{$_} = "filtered";
                next;
	}

        push(@qualIdx, $_);
        $relAbundanceSum += $_[0]->{$linearDOCidx}->{VAL}->{$_};
    }

    #for (0..$totalEntries-1) {
    #    my $relAbundance               = sprintf("%.15f", $_[0]->{$linearDOCidx}->{VAL}->{$_} / $relAbundanceSum);
    #    $_[0]->{$nextIdx}->{VAL}->{$_} = $relAbundance;
    #    $maxRelAbundance = $relAbundance if($relAbundance > $maxRelAbundance);
    #}

    # Loop through each entry, creating its new abundance value
    my $maxRelAbundance = 0;               # Set some arbitrary max
    foreach (@qualIdx) {
        my $relAbundance = sprintf("%.15f", $_[0]->{$linearDOCidx}->{VAL}->{$_} / $relAbundanceSum);
        $_[0]->{$nextIdx}->{VAL}->{$_} = $relAbundance;
        $maxRelAbundance = $relAbundance if($relAbundance > $maxRelAbundance);
    }

    # Add new entry to results hash (must increment if any new fields to be added)
    $nextIdx++;
    $_[0]->{$nextIdx}->{NAME} = "SCALED_REL_ABUNDANCE";
    # Scale each relative Abundance value to the smallest value
    for (0..$totalEntries-1) {
	if( $_[0]->{$nextIdx-1}->{VAL}->{$_} eq 'filtered' ){
        	$_[0]->{$nextIdx}->{VAL}->{$_} = "filtered";
		next;
	}
        $_[0]->{$nextIdx}->{VAL}->{$_} = "0";
        $_[0]->{$nextIdx}->{VAL}->{$_} = sprintf("%.15f", $_[0]->{$nextIdx-1}->{VAL}->{$_} / $maxRelAbundance) if $maxRelAbundance;
    }
    
}
################################################################################
# ARGS: \%results
################################################################################
sub abundanceMethod2 {
}
################################################################################
# ARGS: \%results
################################################################################
sub abundanceMethod3 {
}
################################################################################
# ARGS: \%results, $filename
################################################################################
sub exportUpdatedTable {

    print "Saving updated table to \"".$_[1]."\"...";
    open my $OUTFILE, '>', $_[1];
    
    my @sortedHeaderIndices = sort {$a <=> $b} keys %{ $_[0] };
    print $OUTFILE join("\t", map { $_[0]->{$_}->{NAME} } @sortedHeaderIndices)."\n";
    my @sortedEntryIndices = sort { $a <=> $b } keys %{ $_[0]->{0}->{VAL} };
    my $taxa_cnt=0;

    foreach my $entryIdx (@sortedEntryIndices) {
        foreach my $sortedHeaderIdx (0..$#sortedHeaderIndices) {
            if(exists $_[0]->{$sortedHeaderIdx}->{VAL}->{$entryIdx}) {
                print $OUTFILE $_[0]->{$sortedHeaderIdx}->{VAL}->{$entryIdx}."\t";
            }
            else {
                print $OUTFILE "\t";
            }
        }
        print $OUTFILE "\n";
	$taxa_cnt++;
    }
        
    close $OUTFILE;
    
    print "done. ($taxa_cnt taxonomies)\n";
    
}
################################################################################
# ARGS: $filename, $fieldSep
################################################################################
sub parseTable {

    my %results = ();

    my @requiredHeaders = ("LINEAR_COV","LINEAR_DOC","best_LINEAR_COV", "best_LINEAR_DOC");
    my $numFields = 0;
    
    # Need to check for one of the following in each table:
    my @rankNames = ("STRAIN-LEVEL_REPLICON","STRAIN","SPECIES","GENUS","FAMILY","ORDER","CLASS","PHYLUM");

    open my $INFILE, '<', $_[0];
    
    # Process just the header
    while(my $line = <$INFILE>) {
        chomp $line;
        next if($line eq "");
        my @colHeaders = split($_[1], $line);
        my $numFields = scalar(@colHeaders);
        die "Invalid table format or field separator \"".$_[1]."\"" if(scalar(@colHeaders) <= 1);
        
        # Check for required field headers
        my %colHeaders = ();
           @colHeaders{@colHeaders} = ();
        my @missingHeaders = ();
        foreach my $currHeader (@requiredHeaders) {
            push(@missingHeaders, $currHeader) if(!exists $colHeaders{$currHeader});
        }
        if(@missingHeaders) {
            print "*FATAL*: The following column headers are required in your table, but are missing:\n";
            print "     ".$_."\n" foreach (@missingHeaders);
            die "\n";
        }
        
        # Populate results with column headers
        $results{$_}->{NAME} = $colHeaders[$_] for (0..$#colHeaders);
        
        last;
    }

    # Process the rest of the table
    my $entryNum = 0;
    LINE: while(my $line = <$INFILE>) {
        chomp $line;
        next if($line eq "");
        my @colHeaders = split($_[1], $line);
        #print join("|",@colHeaders)."\n";
        die "Invalid table format or field separator \"".$_[1]."\"" if(scalar(@colHeaders) < $numFields);
        
        $results{$_}->{VAL}->{$entryNum} = $colHeaders[$_] for (0..$#colHeaders);
        $entryNum++;
    }
    close $INFILE;

    #die "*FATAL*: No valid entries detected in table!*\n" if($entryNum < 1);
    #print Dump(\%results)."\n"; exit;

    return \%results;

}
################################################################################
# ARGS:
################################################################################
sub findTopHit {
    RANK: foreach my $cRank (@validRanks) {         # process in order
        my $outfile = $prefix ? ($outdir."/".$prefix.".".$revRankAbbr{$_[2]}.".tsv") : ($outdir."/".$revRankAbbr{$_[2]}.".tsv");
    }
}
################################################################################
################################################################################
################################################################################
# ARGS: \%opts
################################################################################
sub computeTaxRollups {

    if($_[0]->{TREEFILE} && -e $_[0]->{TREEFILE} && !$_[0]->{NOVITALS}) {

        print "===== Extended Taxonomic Rank-level Analysis =====\n";
        attachTaxTree($_[0]->{N_O_COV_HREF}, $_[0]->{TREEFILE}, $_[0]->{GVITALS_HREF});
        RANK: foreach my $cRank (@validRanks) {         # process in order
            # -----------------------------------------------------------
            # Removed option to specify desired tax-rollups: 05-02-2013
            # -----------------------------------------------------------
            #next RANK if(!exists $wanted{$cRank});      # only want user-specified ranks
            computeData($_[0]->{N_O_COV_HREF}, $_[0]->{GI2MERGEDFRAGS_HREF}, $rankAbbr->{$cRank}, $_[0]->{GI2MAPFRAGS_HREF}, $_[0]->{GVITALS_HREF});
        }    
    }
    else {
        print "->No Taxonomy tree datafile provided. Exiting...\n";
        exit;
    }
    
    return;
}
################################################################################
#
# ------------------
# 2013-03-28 UPDATE
# ------------------
# Including STRAIN-level data within this sub, where a strain's unique DB length
# is comprised of it's CHR + PLASMID lengths.
# 
#
#
#
# %nonOverlappedCoverage =
# KEY: gi (from map file)
# VAL: {$gi}->{LINLEN}:   total non-overlapping (LINEAR) length
#      {$gi}->{USIZE}:    total unique reference length of GI in uniques DB
#      {$gi}->{COV}:      {LINLEN}/{USIZE}*100 = non-overlapping %coverage
#      {$gi}->{REPLICON}: strain's replicon name associated with the $gi
#      {$gi}->{STRAIN}:   strain's name associated with the $gi
#   
################################################################################
#              $_[0]                   $_[1]       $_[2]        $_[3]         $_[4]
# ARGS: \%nonOverlappedCoverage, \%gi2mergedFrags, $cRank, \%gi2mapFrags, \%genomeVitals
################################################################################
sub computeData {

    # Initialize filename based on rank-level analysis desired
    my $outfile = $prefix ? ($outdir."/".$prefix.".".$revRankAbbr{$_[2]}.".tsv") : ($outdir."/".$revRankAbbr{$_[2]}.".tsv");
    open my $OUTFILE, '>', $outfile;

    my $iter0 = new Benchmark;
    print "->Rolling up results for rank ".uc($revRankAbbr{$_[2]})." [$outfile]...";

    # Headers
    # ==========================================================================
    # column 01:
    #
    # RANKNAME      (REPLICON)  = replicon name (source + plasmid/chr)
    #               (STRAIN)    = strain name
    #               (SPECIES)   = species name
    #               (GENUS)     = genus name
    #                ...
    # ==========================================================================
    # column 02:
    #
    # NUM_SUBRANKS              = no. of distinct subranks for the current rank
    #                             (E.g.  the no. of SPECIES under the current GENUS)
    # ==========================================================================
    # column 03:
    #
    # GPROJ_ENTRIES             = no. of genome projects (i.e. STRAINS) under this RANK NAME
    # ==========================================================================
    # column 04:
    #
    # LINEAR_LENGTH             = N/O_LENGTH 
    #                           = non-overlapping length 
    #                           = no. of non-overlapping bases covering the unique DB
    # ==========================================================================
    # column 05:
    #
    # UNIQUE_DB_LENGTH          = no. of unique bases for this organism
    # ==========================================================================
    # column 06:
    #
    # FULL_REFDB_LENGTH         = no. of bases in full reference
    # ==========================================================================
    # column 07:
    #
    # LINEAR_COV                = LINEAR_LENGTH / UNIQUE_DB_LENGTH
    # ==========================================================================
    # column 08:
    #
    # HIT_COUNT                 = no. of hits recruited to genome
    # ==========================================================================
    # column 09:
    #
    # FULL_HIT_COUNT            = no. of full-length read hits recruited to genome
    # ==========================================================================
    # column 10:
    #
    # TOTAL_BP_MAPPED           = sum total of all hit lengths recruited to genome
    #                           = hit1.length + hit2.length + ... hitX.length
    # ==========================================================================
    # column 11:
    #
    #                           [formerly FOLD_COV_UNIQUE_SAMPLE]
    # LINEAR_DOC                = linear depth-of-coverage
    #                           = fold coverage of sample's LINEAR_LENGTH 
    #                           = TOTAL_BP_MAPPED / LINEAR_LENGTH
    # ==========================================================================
    # column 12:
    #
    #                           [formerly FOLD_COV_UNIQUE_REFDB]
    # UREF_DOC                  = unique reference's depth-of-coverage
    #                           = fold coverage of reference's UNIQUE_DB_LENGTH
    #                           = TOTAL_BP_MAPPED / UNIQUE_DB_LENGTH
    # ==========================================================================
    # column 13:
    #
    # UREF_CMAX                 = MAX COVERAGE OF REFDB POSSIBLE, GIVEN SAMPLE INPUT BASES
    #                           = Cmax = L0/l0 
    #                           = TOTAL_INPUT_BASES / UNIQUE_DB_LENGTH
    # --------------------------------------------------------------------------
    # <<<< REMOVED >>>>
    # LINEAR_COV_SCALING_FACTOR = UNIQUE_DB_LENGTH / MAX(UNIQUE_DB_LENGTH)
    # --------------------------------------------------------------------------
    # <<<< REMOVED >>>>
    # SCALED_LINEAR_COV         = LINEAR_LENGTH / UNIQUE_DB_LENGTH * LINEAR_COV_SCALING_FACTOR
    #                           = LINEAR_COV * LINEAR_COV_SCALING_FACTOR
    # --------------------------------------------------------------------------
    # <<<< REMOVED >>>>
    # HIT_SCALING_FACTOR        = MIN(UNIQUE_DB_LENGTH) / UNIQUE_DB_LENGTH
    # --------------------------------------------------------------------------
    # <<<< REMOVED >>>>
    # SCALED_HIT_COUNT          = HIT_COUNT * [ MIN(UNIQUE_DB_LENGTH) / UNIQUE_DB_LENGTH ]
    #                           = HIT_COUNT * HIT_SCALING_FACTOR
    # ==========================================================================
    # column 14:
    # ---------
    # FRAC_HITS_POSSIBLE        = HIT_COUNT / TOTAL_INPUT_READS
    # ==========================================================================
    # column 15:
    #
    # FRAC_BASES_POSSIBLE       = TOTAL_BP_MAPPED / TOTAL_INPUT_BASES
    # ==========================================================================
    # column 16:
    #
    # MEAN_HIT_LENGTH           = TOTAL_BP_MAPPED / HIT_COUNT
    # ==========================================================================
    # column 17:
    #
    # MEAN_LINEAR_HIT_LENGTH    = LINEAR_LENGTH / HIT_COUNT
    # ==========================================================================

    #///////////////////////////////////////////////////////////////////////////
    # **** BEST SUBRANK **** BEST SUBRANK **** BEST SUBRANK **** BEST SUBRANK ****
    #///////////////////////////////////////////////////////////////////////////

    # ==========================================================================
    # column 18:
    #
    # best_SUBRANK              = name of the best subrank (determined by the highest LINEAR_COV)
    # ==========================================================================
    # column 19:
    #
    # best_NUM_SUBRANKS         = no. of subranks supporting current "SUBRANK"
    #
    #                      {SS} = no. of GI entries supporting this strain
    #                      {S}  = no. of strains supporting this species
    #                      {G}  = no. of species supporting this genus
    #                      {F}  = no. of genera supporting this family
    #                      {O}  = no. of families supporting this order 
    #                      {C}  = no. of orders supporting this class
    #                      {P}  = no. of classes supporting this phylum
    # ==========================================================================
    # column 20:
    #
    # best_GPROJ_ENTRIES        = no. of genome projects (i.e. STRAINS) under this best_SUBRANK
    #
    #                      {SS} = no. of genome projects supporting this strain = 1
    #                      {S}  = no. of genome projects supporting this species
    #                      {G}  = no. of genome projects supporting this genus
    #                      {F}  = no. of genome projects supporting this family
    #                      {O}  = no. of genome projects supporting this order
    #                      {C}  = no. of genome projects supporting this class
    #                      {P}  = no. of genome projects supporting this phylum
    # ==========================================================================
    # column 21:
    #
    # best_LINEAR_LENGTH
    # ==========================================================================
    # column 22:
    #
    # best_UNIQUE_DB_LENGTH
    # ==========================================================================
    # column 23:
    #
    # best_FULL_REFDB_LENGTH
    # ==========================================================================
    # column 24:
    #
    # best_LINEAR_COV
    # ==========================================================================
    # column 25:
    #
    # best_HIT_COUNT
    # ==========================================================================
    # column 26:
    #
    # best_FULL_HIT_COUNT
    # ==========================================================================
    # column 27:
    #
    # best_TOTAL_BP_MAPPED
    # ==========================================================================
    # column 28:
    #
    # best_LINEAR_DOC (a.k.a. Abundance)
    # ==========================================================================
    # column 29:
    #
    # best_UREF_DOC
    # ==========================================================================
    # column 30:
    #
    # best_UREF_CMAX
    # ==========================================================================
    # column 31:
    #
    # best_FRAC_HITS_POSSIBLE
    # ==========================================================================
    # column 32:
    #
    # best_FRAC_BASES_POSSIBLE
    # ==========================================================================
    # column 33:
    #
    # best_MEAN_HIT_LENGTH
    # ==========================================================================
    # column 34:
    #
    # best_MEAN_LINEAR_HIT_LENGTH
    # ==========================================================================
    # column 35:
    #
    # CONTIG_COUNT              = No. of contiguous fragments
    #                             (after mapping & generating non-overlapping fragments)
    # ==========================================================================
    # column 36:
    #
    # CONTIG_MEAN_LEN           = Mean length of contigs (bp)
    # ==========================================================================
    # column 37:
    #
    # CONTIG_STDEV_LEN          = Standard deviation of contig lengths (bp)
    # ==========================================================================
    # column 38:
    #
    # CONTIG_MINLEN             = Length of smallest contig
    # ==========================================================================
    # column 39:
    #
    # CONTIG_MAXLEN             = Length of largest contig
    # ==========================================================================
    # column 40:
    #
    # CONTIG_HISTOGRAM(LEN:FREQ)
    #                           = Contig Length Histogram
    #                             (in the format contigLength:frequency)
    # ==========================================================================
    # <DISABLED>
    # column 41:
    # 
    # CONTIG_N50:N60:N70:N80:N90:N95
    #                           =
    # ==========================================================================
    print $OUTFILE uc($revRankAbbr{$_[2]})."\t"       # rank name
                  #."GI_ENTRIES\t".                  
        ."NUM_SUBRANKS\t".                  "GPROJ_ENTRIES\t".              "LINEAR_LENGTH\t"
        ."UNIQUE_DB_LENGTH\t".              "FULL_REFDB_LENGTH\t".          "LINEAR_COV\t"
        ."HIT_COUNT\t".                     "FULL_HIT_COUNT\t".             "TOTAL_BP_MAPPED\t".            "LINEAR_DOC\t"
        ."UREF_DOC\t".                      "UREF_CMAX\t".                  "FRAC_HITS_POSSIBLE\t"
        ."FRAC_BASES_POSSIBLE\t".           "MEAN_HIT_LENGTH\t".            "MEAN_LINEAR_HIT_LENGTH\t"
        #------------------------------------------------------------------------------------------
        ."best_SUBRANK\t"                  
        ."best_NUM_SUBRANKS\t".             "best_GPROJ_ENTRIES\t".         "best_LINEAR_LENGTH\t"
        ."best_UNIQUE_DB_LENGTH\t".         "best_FULL_REFDB_LENGTH\t".     "best_LINEAR_COV\t"
        ."best_HIT_COUNT\t".                "best_FULL_HIT_COUNT\t".        "best_TOTAL_BP_MAPPED\t".       "best_LINEAR_DOC\t"
        ."best_UREF_DOC\t".                 "best_UREF_CMAX\t".             "best_FRAC_HITS_POSSIBLE\t"
        ."best_FRAC_BASES_POSSIBLE\t".      "best_MEAN_HIT_LENGTH\t".       "best_MEAN_LINEAR_HIT_LENGTH\t"
        #------------------------------------------------------------------------------------------
        ."CONTIG_COUNT\t".                  "CONTIG_MEAN_LEN\t".            "CONTIG_STDEV_LEN\t"
        ."CONTIG_MINLEN\t".                 "CONTIG_MAXLEN\t".              "CONTIG_HISTOGRAM(LEN:FREQ)\t"
        ."CONTIG_N50:N60:N70:N80:N90:N95\n";

    # Generate %hitTree (\%nonOverlappedCoverage, \%gi2mergedFrags, \%gi2mapFrags, \%genomeVitals)
    my $hitTree = genTree($_[0], $_[1], $_[3], $_[4]);

    # Write out strain contig coordinates
    writeCoords($hitTree, "SS", "strain");

    # Calculate coverage, hit_count, etc. using the GIs in $hitTree{$gi} and $
    my @rankNames = keys %{ $hitTree->{$_[2]} };       # Depending on tax rank, contains all STRAIN names,
                                                       # SPECIES names, GENUS names, FAMILY names, etc.

    RANKNAME: foreach my $rankName (sort @rankNames) {
        # The number of subranks supporting the current rank
        my $numSubranks            = scalar(keys %{ $hitTree->{$_[2]}->{$rankName}->{CHILDREN}->{NAMES} });
        my $entries                = ($_[2] eq "SS")                              # total no. of Genome Project entries in rank (a count of strains)
                                   ? scalar(keys %{ $hitTree->{$_[2]}->{$rankName}->{CHILDREN}->{NAMES} })
                                   : $hitTree->{$_[2]}->{$rankName}->{GPROJ};

        my $linearLength           = $hitTree->{$_[2]}->{$rankName}->{LINLEN};       # total non-overlapping LINEAR length
        my $uniqueDBlength         = $hitTree->{$_[2]}->{$rankName}->{USIZE};        # total unique DB length
        
        if(!exists $hitTree->{$_[2]}->{$rankName}->{COV}) {
            print "NOT EXIST!\n";
            print "    rankName = \"$rankName\"\n";
            print Dump($hitTree->{$_[2]}->{$rankName})."\n";
            <STDIN>;
        }
                
        my $linearCoverage         = sprintf("%.15f", $hitTree->{$_[2]}->{$rankName}->{COV});       # total linear coverage
        my $hitCount               = $hitTree->{$_[2]}->{$rankName}->{HITCOUNT};     # total hit count
        my $fullReadHitCount       = 0; #$hitTree->{$_[2]}->{$rankName}->{FULLREADHITS}; # total hit count (in terms of full-length reads)
        my $totalBpMapped          = $hitTree->{$_[2]}->{$rankName}->{MAPPED};       # total bases mapped (sum of all hits)
        my $fullRefDBlength        = $hitTree->{$_[2]}->{$rankName}->{REFSIZE};      # total bases in the COMPLETE reference (not just the unique bases)
        my $linearDOC              = sprintf("%.15f", $totalBpMapped / $linearLength);
        my $uRefDOC                = sprintf("%.15f", $totalBpMapped / $uniqueDBlength);
        # Max Coverage of RefDB Possible in Sample = Cmax = L0/l0 = (Total Possible Bases) / (unique refDB Bases)
        my $uRefCmax               = sprintf("%.15f", $inputBases/$uniqueDBlength);
        # Fraction of Total Hits Possible = uN = Ni/N0 = (HIT_COUNT, genome i) / (Total Potential Hits)
        my $fracHitsPossible       = ($inputReads == 0) ? (0.0) : sprintf("%.15f", $hitCount / $inputReads);
        # Fraction of Total Bases Possible = uL = Li/L0 = (TOTAL_BP_MAPPED, genome i) / (Total Potential Bases)
        my $fracBasesPossible      = ($inputBases == 0) ? (0.0) : sprintf("%.15f", $totalBpMapped / $inputBases);
        my $mean_hit_length        = sprintf("%.2f", $totalBpMapped / $hitCount);
        my $mean_linear_hit_length = sprintf("%.2f", $linearLength  / $hitCount);

        # ----------------------------------------------------------------------
        # Calculate values for the BEST Subrank
        # ----------------------------------------------------------------------
        # Get previous tax rank abbr
        my $prevTreeRank = q{};
        RANKIDX: foreach my $rankIdx (1..$#hitTreeRanks) {      # Search from "SS" on up
            if($hitTreeRanks[$rankIdx] eq $_[2]) {
                $prevTreeRank = $hitTreeRanks[$rankIdx-1];
                last RANKIDX;
            }
        }
        my $bestName                = $hitTree->{$_[2]}->{$rankName}->{CHILDREN}->{BEST};
        # The number of subranks supporting the current rank
        my $bestNumSubranks         =  scalar( keys %{ $hitTree->{$prevTreeRank}->{$bestName}->{CHILDREN}->{NAMES} } );
        my $bestNumEntries          = (($prevTreeRank eq "SS") || ($prevTreeRank eq "GI"))
                                    ? (scalar( keys %{ $hitTree->{$prevTreeRank}->{$bestName}->{CHILDREN}->{NAMES} } ))
                                    : $hitTree->{$prevTreeRank}->{$bestName}->{GPROJ};
        my $bestLinearLength        = $hitTree->{$prevTreeRank}->{$bestName}->{LINLEN};
        my $bestUniqueDBlength      = $hitTree->{$prevTreeRank}->{$bestName}->{USIZE};
        my $bestFullRefDBlength     = $hitTree->{$prevTreeRank}->{$bestName}->{REFSIZE};
        my $bestLinearCoverage      = sprintf("%.15f", $hitTree->{$prevTreeRank}->{$bestName}->{COV});
        my $bestHitCount            = $hitTree->{$prevTreeRank}->{$bestName}->{HITCOUNT};
        my $bestFullReadHitCount    = 0; #$hitTree->{$prevTreeRank}->{$bestName}->{FULLREADHITS};
        my $bestTotalBpMapped       = $hitTree->{$prevTreeRank}->{$bestName}->{MAPPED};
        my $bestLinearDOC           = sprintf("%.15f", $bestTotalBpMapped / $bestLinearLength);
        my $bestURefDOC             = sprintf("%.15f", $bestTotalBpMapped / $bestUniqueDBlength);
        my $bestURefCmax            = sprintf("%.15f", $inputBases / $bestUniqueDBlength);
        my $bestFracHitsPossible    = ($inputReads == 0) ? (0.0) : sprintf("%.15f", $bestHitCount / $inputReads);
        my $bestFracBasesPossible   = ($inputBases == 0) ? (0.0) : sprintf("%.15f", $bestTotalBpMapped / $inputBases);
        my $bestMeanHitLength       = sprintf("%.15f", $bestTotalBpMapped / $bestHitCount);
        my $bestMeanLinearHitLength = sprintf("%.15f", $bestLinearLength  / $bestHitCount);

        print $OUTFILE $rankName."\t"                               # Rank Name
                      .$numSubranks."\t"                            # No. of Subranks (i.e. the no. of SPECIES under the current GENUS)
                      .$entries."\t"                                # No. of GPROJ entries
                      .$linearLength."\t"
                      .$uniqueDBlength."\t"
                      .$fullRefDBlength."\t"                        # ADDED: 2013-04-02
                      .$linearCoverage."\t"
                      #.$fullLengthHitCount."\t"                     #
                      .$hitCount."\t"
                      .$fullReadHitCount."\t"
                      .$totalBpMapped."\t"
                      .$linearDOC."\t"
                      .$uRefDOC."\t"
                      .$uRefCmax."\t"
                      .$fracHitsPossible."\t"
                      .$fracBasesPossible."\t"
                      .$mean_hit_length."\t"
                      .$mean_linear_hit_length."\t"
                      #----------------------------
                      .$bestName."\t"
                      .$bestNumSubranks."\t"
                      .$bestNumEntries."\t"
                      .$bestLinearLength."\t"
                      .$bestUniqueDBlength."\t"
                      .$bestFullRefDBlength."\t"
                      .$bestLinearCoverage."\t"
                      .$bestHitCount."\t"
                      .$bestFullReadHitCount."\t"
                      .$bestTotalBpMapped."\t"
                      .$bestLinearDOC."\t"
                      .$bestURefDOC."\t"
                      .$bestURefCmax."\t"
                      .$bestFracHitsPossible."\t"
                      .$bestFracBasesPossible."\t"
                      .$bestMeanHitLength."\t"
                      .$bestMeanLinearHitLength."\t";

        # Sort descending
        my @contigLengths = sort {$b <=> $a} keys %{ $hitTree->{$_[2]}->{$rankName}->{HISTO} };

        # ----------------------------------
        # Calculate contig length statistics
        # ----------------------------------
        my ($totalCount, $currFreq, $mean, $stdev, $min, $max, $sum) = (0,0,0,0,0,0,0);
        $min = $contigLengths[-1];       # req's the array to be decr. sorted
        $max = $contigLengths[0];        # req's the array to be decr. sorted
        foreach (@contigLengths) {
            $currFreq    = $hitTree->{$_[2]}->{$rankName}->{HISTO}->{$_};
            $totalCount += $currFreq;
            $sum        += ($_)*($currFreq);
        }
        
        # ----------------------------------
        # Calculate N50 stats
        # ----------------------------------
        my ($n50, $n60, $n70, $n80, $n90, $n95) = (0,0,0,0,0,0);     
        my $n50cutoff = (0.50)*$sum;
        my $n60cutoff = (0.60)*$sum;
        my $n70cutoff = (0.70)*$sum;
        my $n80cutoff = (0.80)*$sum;
        my $n90cutoff = (0.90)*$sum;
        my $n95cutoff = (0.95)*$sum;
        my $tmpSum    = 0;
        foreach (@contigLengths) {   
            $currFreq = $hitTree->{$_[2]}->{$rankName}->{HISTO}->{$_};
            $tmpSum += ($_)*($currFreq);
            $n50 = $_ if(!$n50 && $tmpSum >= $n50cutoff);
            $n60 = $_ if(!$n60 && $tmpSum >= $n60cutoff);
            $n70 = $_ if(!$n70 && $tmpSum >= $n70cutoff);
            $n80 = $_ if(!$n80 && $tmpSum >= $n80cutoff);
            $n90 = $_ if(!$n90 && $tmpSum >= $n90cutoff);
            $n95 = $_ if(!$n95 && $tmpSum >= $n95cutoff);
        }

        $mean = $sum;
        $mean /= $totalCount;
        
        if($totalCount == 1) {
            $stdev = 0;
        }
        else {
            #$stdev += ($_ - $mean)**2 foreach (@contigLengths);
            foreach my $cLen (@contigLengths) {
                $stdev += ($cLen - $mean)**2 for (1..$hitTree->{$_[2]}->{$rankName}->{HISTO}->{$cLen});
            }
            $stdev = sqrt($stdev/($totalCount-1));
        }
        $mean  = sprintf("%.1f", $mean);
        $stdev = sprintf("%.1f", $stdev);
        
        print $OUTFILE $totalCount."\t"
                      .$mean."\t"
                      .$stdev."\t"
                      .$min."\t"
                      .$max."\t";

        # ----------------------------------
        # Print out the histogram string:    LENGTH1:FREQ1,LENGTH2:FREQ2,...
        # ----------------------------------
        #my @contigLengths = sort {$a <=> $b} keys %{ $hitTree->{$_[2]}->{$rankName}->{HISTO} };
        print $OUTFILE $contigLengths[0].":";
        for (1..$#contigLengths) {
            print $OUTFILE join(",",$hitTree->{$_[2]}->{$rankName}->{HISTO}->{$contigLengths[$_-1]},$contigLengths[$_]).":";
        }
        print $OUTFILE $hitTree->{$_[2]}->{$rankName}->{HISTO}->{$contigLengths[$#contigLengths]}."\t";

        # ----------------------------------
        # Print out the N50:N60:N70:N80:N90:N95 string
        # ----------------------------------
        print $OUTFILE $n50.":".$n60.":".$n70.":".$n80.":".$n90.":".$n95."\n";

    } #rankNames
    
    close $OUTFILE;
    my $iter1 = new Benchmark;
    print "done. ".timestr(timediff($iter1,$iter0))."\n";

    return;
}
################################################################################
# ARGS: \%opts
################################################################################
sub computeRepliconData {

    # Initialize filename based on rank-level analysis desired
    my $outfile = $prefix ? ($outdir."/".$prefix.".replicon.tsv") : ($outdir."/"."replicon.tsv");
    open my $OUTFILE, '>', $outfile;

    print "===== REPLICON-level Analysis =====\n";
    my $iter0 = new Benchmark;
    print "->Writing REPLICON-level results to disk [$outfile]...";

    # Headers
    # --------------------------------------------------------------------------
    # STRAIN-LEVEL_REPLICON     = Name of the replicon
    # --------------------------------------------------------------------------
    # GI
    # --------------------------------------------------------------------------
    # GPROJ_ENTRIES
    # --------------------------------------------------------------------------
    # LINEAR_LENGTH             = N/O_LENGTH 
    #                           = non-overlapping length 
    #                           = no. of non-overlapping bases covering the unique DB
    # --------------------------------------------------------------------------
    # UNIQUE_DB_LENGTH          = no. of unique bases for this organism
    # --------------------------------------------------------------------------
    # FULL_REFDB_LENGTH         = no. of bases of full-length replicon (not just uniques)
    # --------------------------------------------------------------------------
    # LINEAR_COV                = LINEAR_LENGTH / UNIQUE_DB_LENGTH
    # --------------------------------------------------------------------------
    # HIT_COUNT                 = no. of hits recruited to genome
    # --------------------------------------------------------------------------
    # FULL_HIT_COUNT            = no. of full-length read hits recruited to genome
    # --------------------------------------------------------------------------
    # TOTAL_BP_MAPPED           = sum total of all hit lengths recruited to genome
    #                           = hit1.length + hit2.length + ... hitX.length
    # --------------------------------------------------------------------------
    # "ABUNDANCE"
    # LINEAR_DOC                = depth of coverage of LINEAR_LENGTH 
    # (FOLD_COV_UNIQUE_SAMPLE)  = TOTAL_BP_MAPPED / LINEAR_LENGTH
    # --------------------------------------------------------------------------
    # UREF_DOC                  = TOTAL_BP_MAPPED / UNIQUE_DB_LENGTH
    # (FOLD_COV_UNIQUE_REFDB)
    # --------------------------------------------------------------------------
    # "UNIQUE REFERENCE MAX COVERAGE IN SAMPLE"
    # UREF_CMAX                 = TOTAL_INPUT_BASES / UNIQUE_DB_LENGTH
    # (REF_CMAX)
    # --------------------------------------------------------------------------
    # FRAC_HITS_POSSIBLE        = HIT_COUNT / TOTAL_INPUT_READS
    # --------------------------------------------------------------------------
    # FRAC_BASES_POSSIBLE       = TOTAL_BP_MAPPED / TOTAL_INPUT_BASES
    # --------------------------------------------------------------------------
    # MEAN_HIT_LENGTH           = TOTAL_BP_MAPPED / HIT_COUNT
    # --------------------------------------------------------------------------
    # MEAN_LINEAR_HIT_LENGTH    = LINEAR_LENGTH / HIT_COUNT
    # --------------------------------------------------------------------------
    # CONTIG_COUNT
    # --------------------------------------------------------------------------
    # CONTIG_MEAN_LEN
    # --------------------------------------------------------------------------
    # CONTIG_STDEV_LEN
    # --------------------------------------------------------------------------
    # CONTIG_MINLEN
    # --------------------------------------------------------------------------
    # CONTIG_MAXLEN
    # --------------------------------------------------------------------------
    # CONTIG_HISTO(LEN:FREQ)
    # --------------------------------------------------------------------------
    # <DISABLED>
    # CONTIG_N50:N60:N70:N80:N90:N95
    # --------------------------------------------------------------------------
    #
    #
    print $OUTFILE "STRAIN-LEVEL_REPLICON\t" if(!$noVitalsFile);
    print $OUTFILE "GI\t".                     "GPROJ_ENTRIES\t".             "LINEAR_LENGTH\t"
                  ."UNIQUE_DB_LENGTH\t".       "FULL_REFDB_LENGTH\t".         "LINEAR_COV\t"
                  ."HIT_COUNT\t".              "FULL_HIT_COUNT\t".            "TOTAL_BP_MAPPED\t".           "LINEAR_DOC\t"
                  ."UREF_DOC\t".               "UREF_CMAX\t".                 "FRAC_HITS_POSSIBLE\t"
                  ."FRAC_BASES_POSSIBLE\t".    "MEAN_HIT_LENGTH\t".           "MEAN_LINEAR_HIT_LENGTH\t"
                  #------------------------------------------------------------------------------------------
                  ."CONTIG_COUNT\t".           "CONTIG_MEAN_LEN\t".           "CONTIG_STDEV_LEN\t"
                  ."CONTIG_MINLEN\t".          "CONTIG_MAXLEN\t".             "CONTIG_HISTO(LEN:FREQ)\t"
                  ."CONTIG_N50:N60:N70:N80:N90:N95\n";

    # Find the MAX and MIN unique reference lengths for proper scaling in next step
    my @gis = sort { $_[0]->{N_O_COV_HREF}->{$b}->{COV} 
                 <=> $_[0]->{N_O_COV_HREF}->{$a}->{COV} 
                   } keys %{ $_[0]->{N_O_COV_HREF} };
    my $minRefLen = $_[0]->{N_O_COV_HREF}->{$gis[0]}->{USIZE};  # init with first value
    my $maxRefLen = $_[0]->{N_O_COV_HREF}->{$gis[0]}->{USIZE};  # init with first value
    for (1..$#gis) {
        $minRefLen = $_[0]->{N_O_COV_HREF}->{$gis[$_]}->{USIZE} if($_[0]->{N_O_COV_HREF}->{$gis[$_]}->{USIZE} < $minRefLen);
        $maxRefLen = $_[0]->{N_O_COV_HREF}->{$gis[$_]}->{USIZE} if($_[0]->{N_O_COV_HREF}->{$gis[$_]}->{USIZE} > $maxRefLen);
    }

    # Cycle through GIs in %nonOverlappedCoverage
    foreach my $gi (@gis) {
    
        my $replicon = q{};
           $replicon = $_[0]->{N_O_COV_HREF}->{$gi}->{REPLICON} 
               if(!$noVitalsFile && exists $_[0]->{N_O_COV_HREF}->{$gi});

        my $fullReadHitCount = 0; #$_[0]->{GI2MAPFRAGS_HREF}->{$gi}->{FULLREADHITS};
        my $totalBpMapped    = $_[0]->{GI2MAPFRAGS_HREF}->{$gi}->{TOTALBPMAPPED};
        my $uniqueDBlength   = $_[0]->{N_O_COV_HREF}->{$gi}->{USIZE};

        my $covScalingFactor     = $uniqueDBlength / $maxRefLen;
        my $scaledLinearCoverage = $_[0]->{N_O_COV_HREF}->{$gi}->{LINLEN} / $uniqueDBlength * $covScalingFactor;

        #my $hitScalingFactor     = $minRefLen / $uniqueDBlength;
        #my $scaledHitCount       = $_[0]->{GI2MERGEDFRAGS_HREF}->{$gi}->{HITCOUNT} * $hitScalingFactor;
        #   $scaledHitCount       = ($scaledHitCount > 1) 
        #                         ? $scaledHitCount 
        #                         : 1;


        my $linearDOC = sprintf("%.4f", $totalBpMapped / $_[0]->{N_O_COV_HREF}->{$gi}->{LINLEN});
        my $uniqueDOC = sprintf("%.4f", $totalBpMapped / $uniqueDBlength);
        my $mean_hit_length = 
            sprintf("%.4f", $totalBpMapped / $_[0]->{GI2MERGEDFRAGS_HREF}->{$gi}->{HITCOUNT});
        my $mean_linear_hit_length =
            sprintf("%.4f", $_[0]->{N_O_COV_HREF}->{$gi}->{LINLEN} / $_[0]->{GI2MERGEDFRAGS_HREF}->{$gi}->{HITCOUNT});

        # Fraction of Total Hits Possible = uN = Ni/N0 = (HIT_COUNT, genome i) / (Total Potential Hits)
        my $fracHitsPossible  = ($inputReads == 0)
                              ? (0)
                              : $_[0]->{GI2MERGEDFRAGS_HREF}->{$gi}->{HITCOUNT} / $inputReads;

        # Fraction of Total Bases Possible = uL = Li/L0 = (TOTAL_BP_MAPPED, genome i) / (Total Potential Bases)
        my $fracBasesPossible = ($inputBases == 0)
                              ? (0)
                              : $totalBpMapped / $inputBases;
        
        # Max Coverage of RefDB Possible in Sample = Cmax = L0/l0 = (Total Possible Bases) / (refDB Bases)
        my $uRefCmax =sprintf("%.4f", $inputBases/$uniqueDBlength);

        # Observed Abundance = (TOTAL_BP_MAPPED * FULL_REFDB_LENGTH) / (UNIQUE_DB_LENGTH * TOTAL_POSSIBLE BASES)
        #my $obsAbundance1 = ($inputBases == 0)
        #                  ? (0)
        #                  : ($totalBpMapped * $_[0]->{GVITALS_HREF}->{GI}->{$gi}->{SIZE}) / ($uniqueDBlength * $inputBases);      

        #my $obsAbundance2 = ($inputBases == 0)
        #                  ? (0)
        #                  : ($fracBasesPossible * $uniqueDBlength / $_[0]->{GVITALS_HREF}->{GI}->{$gi}->{SIZE});


        print $OUTFILE $replicon."\t" if(!$noVitalsFile);
        print $OUTFILE $gi."\t"
                      ."1\t"                                                # supposed to be the "No. of Genome Projects" but since its a replicon, it's always one
                      .$_[0]->{N_O_COV_HREF}->{$gi}->{LINLEN}."\t"          # non-overlapping length (i.e. Linear Length)
                      .$uniqueDBlength."\t"                                 # unique DB length
                      .$_[0]->{GVITALS_HREF}->{GI}->{$gi}->{SIZE}."\t"      # full reference length of replicon, not just unique length
                      .$_[0]->{N_O_COV_HREF}->{$gi}->{COV}."\t"             # linear coverage
                      .$_[0]->{GI2MERGEDFRAGS_HREF}->{$gi}->{HITCOUNT}."\t" # no. of counted hits to reference
                      .$fullReadHitCount."\t"                               # no. of full read hits to reference
                      .$totalBpMapped."\t"                                  # total BP mapped
                      .$linearDOC."\t"                                      # linear depth-of-coverage
                      .$uniqueDOC."\t"                                      # unique genome depth-of-coverage
                      .$uRefCmax."\t"                                       # max coverage of unique DB given sample input
                      .$fracHitsPossible."\t"                               # fraction of all possible hits
                      .$fracBasesPossible."\t"                              # fraction of all possible bases
                      .$mean_hit_length."\t"                                # totalBPmapped/HitCount
                      .$mean_linear_hit_length."\t";                        # LinearLength/HitCount

        my @contigLengths = sort {$b <=> $a} keys %{ $giHisto{$gi} };

        # ----------------------------------
        # Calculate contig length statistics
        # ----------------------------------
        my ($totalCount, $currFreq, $mean, $stdev, $min, $max, $sum) = (0,0,0,0,0,0,0);
        $min = $contigLengths[-1];       # req's the array to be decr. sorted
        $max = $contigLengths[0];        # req's the array to be decr. sorted

        foreach (@contigLengths) {
            $currFreq    = $giHisto{$gi}->{$_}; #$hitTree->{$_[2]}->{$rankName}->{HISTO}->{$_};
            $totalCount += $currFreq;
            $sum        += ($_)*($currFreq);
        }
        
        # ----------------------------------
        # Calculate N50 stats
        # ----------------------------------
        my ($n50, $n60, $n70, $n80, $n90, $n95) = (0,0,0,0,0,0);     
        my $n50cutoff = (0.50)*$sum;
        my $n60cutoff = (0.60)*$sum;
        my $n70cutoff = (0.70)*$sum;
        my $n80cutoff = (0.80)*$sum;
        my $n90cutoff = (0.90)*$sum;
        my $n95cutoff = (0.95)*$sum;
        my $tmpSum    = 0;
        foreach (@contigLengths) {   
            $currFreq = $giHisto{$gi}->{$_}; #$hitTree->{$_[2]}->{$rankName}->{HISTO}->{$_};
            $tmpSum += ($_)*($currFreq);
            $n50 = $_ if(!$n50 && $tmpSum >= $n50cutoff);
            $n60 = $_ if(!$n60 && $tmpSum >= $n60cutoff);
            $n70 = $_ if(!$n70 && $tmpSum >= $n70cutoff);
            $n80 = $_ if(!$n80 && $tmpSum >= $n80cutoff);
            $n90 = $_ if(!$n90 && $tmpSum >= $n90cutoff);
            $n95 = $_ if(!$n95 && $tmpSum >= $n95cutoff);
        }

        $mean = $sum;
        $mean /= $totalCount;
        
        if($totalCount == 1) {
            $stdev = 0;
        }
        else {
            #$stdev += ($_ - $mean)**2 foreach (@contigLengths);
            foreach my $cLen (@contigLengths) {
                $stdev += ($cLen - $mean)**2 for (1..$giHisto{$gi}->{$cLen});
            }
            $stdev = sqrt($stdev/($totalCount-1));
        }
        $mean  = sprintf("%.1f", $mean);
        $stdev = sprintf("%.1f", $stdev);
        
        print $OUTFILE $totalCount."\t"
                      .$mean."\t"
                      .$stdev."\t"
                      .$min."\t"
                      .$max."\t";

        # Print out the histogram string:    LENGTH1:FREQ1,LENGTH2:FREQ2,...
        #my @contigLengths = sort {$a <=> $b} keys %{ $giHisto{$gi} };
        print $OUTFILE $contigLengths[0].":";
        for (1..$#contigLengths) {
            print $OUTFILE join(",",$giHisto{$gi}->{$contigLengths[$_-1]},$contigLengths[$_]).":";
        }
        print $OUTFILE $giHisto{$gi}->{$contigLengths[$#contigLengths]}."\t";

        # ----------------------------------
        # Print out the N50:N60:N70:N80:N90:N95 string
        # ----------------------------------
        print $OUTFILE $n50.":".$n60.":".$n70.":".$n80.":".$n90.":".$n95."\n";

    }
    
    close $OUTFILE;
    my $iter1 = new Benchmark;
    print "done. ".timestr(timediff($iter1,$iter0))."\n";

    return;
}
################################################################################
# Use "GI" from $nonOverlappedCoverage to get "ORG" from %genomeVitals
# Use "ORG" to get taxonomy from %taxTree by mapping to {S} or {SS}
#
#
# ARGS: \%nonOverlappedCoverage, $treeFile, \%genomeVitals
################################################################################
sub attachTaxTree {

    # Get list of GIs from uniques entries
    print "->Pulling replicon GIs from DB entries...";
    my $iter0 = new Benchmark;
    my @gis = keys %{ $_[0] };
    my $iter1 = new Benchmark;
    print "done. ".timestr(timediff($iter1,$iter0)),"\n";
    # --------------------------------------------------------------------------
    # Get scientific names of uniques entries associated with GIs 
    print "->Mapping replicon GIs to source organisms...";
    $iter0 = new Benchmark;
    my %org2gi = ();                                                                                                # <--------------- HERE
    # PROBLEM: !!!!!!!!!!!!
    #    There are multiple GIs corresponding to a single {ORG}, because the {ORG}
    #    may represent the chromsome sequence *and* the plasmid sequence.
    #    Need to record all the **RELEVANT** (i.e. replicon) GIs associated with
    #    the organism name rather than just map a 1:1.
    # SOLUTION: instead of simple scalar K->V pairs, make the VALUE an AREF of 
    #    relevant GIs.
    foreach my $gi (@gis) {
        my $org = $_[2]->{GI}->{$gi}->{ORG};                                                                        # <--------------- HERE
        push(@{ $org2gi{$org} }, $gi);
    }
    $iter1 = new Benchmark;
    print "done. ".timestr(timediff($iter1,$iter0)),"\n";
    # --------------------------------------------------------------------------
    # Retrieve $taxTree from disk
    print "->Retrieving Tax Tree from disk [".$_[1]."]...";                                                         # <--------------- HERE
    $iter0 = new Benchmark;
    my $taxTree = retrieve $_[1];                                                                                   # <--------------- HERE
    $iter1 = new Benchmark;
    print "done. ".timestr(timediff($iter1,$iter0)),"\n";
    # --------------------------------------------------------------------------
    # Lookup taxonomy tree in %taxTree using $org
    print "->Mapping source organism to its tax tree...";
    $iter0 = new Benchmark;

    my $numGIs     = scalar(@gis);
    my $origNumGIs = $numGIs;
    
    # WE'RE SEARCHING THROUGH THE $TAXTREE BECAUSE WE NEED TO IDENTIFY A TAXONOMY
    # TREE WITH EVERY GI PRESENT. IF AN ORGANISM HAS >1 GI ASSOCIATED WITH IT, ONCE
    # THE NAME IS FOUND ONCE, ALL ITS SUPPORTING GI'S WILL BE ASSIGNED THE SAME
    # TAXONOMY TREE.
    TAXID: foreach my $taxid (keys %{ $taxTree }) {                                                                 # <--------------- HERE
        my @variants = ();
        push(@variants, keys %{ $taxTree->{$taxid}->{SS} }) if(exists $taxTree->{$taxid}->{SS});                    # <--------------- HERE
        push(@variants, keys %{ $taxTree->{$taxid}->{S}  }) if(exists $taxTree->{$taxid}->{S} );                    # <--------------- HERE
        
        foreach (@variants) {
            if(exists $org2gi{$_}) {
                # Multiple GIs associated with organism?
                if(ref($org2gi{$_}) eq "ARRAY") {
                    my $currNumGIs = scalar(@{ $org2gi{$_} });              # no. of GIs assoc. w/org name
                    foreach my $gi (@{ $org2gi{$_} }) {
                        $_[0]->{$gi}->{TAXTREE} = $taxTree->{$taxid};       # add tree                              # <--------------- HERE
                    }
                    delete $org2gi{$_};                                     # delete key
                    $numGIs = $numGIs - $currNumGIs;                        # decrement by no. of GI's present
                }
                else {
                    $_[0]->{ $org2gi{$_} }->{TAXTREE} = $taxTree->{$taxid}; # add tree                              # <--------------- HERE
                    delete $org2gi{$_};                                     # delete key
                    $numGIs--;                                              # decrement
                }
                
            }
            #else {
            #    print "**WARNING**: Unrecognized organism \"$_\"!\n";
            #}
        } #@variants
        last TAXID if($numGIs == 0);
    }
    $iter1 = new Benchmark;
    print "done. ".timestr(timediff($iter1,$iter0)),"\n";
    # --------------------------------------------------------------------------
    # Any leftover keys in %org2gi are the unmatched org names
    my @unmappedOrgs = keys %org2gi;
    if(@unmappedOrgs) {
        print "->At least one organism name in \"$genomeVitalsFile\" is unrecognized in \"$treeFile\":\n";
        print "      [$_]\n" foreach (@unmappedOrgs);
        #die "Abort.\n";
        print "Continuing on...\n";
    }
    else {
        print "  Congratulations! All organisms have been identified!\n";
    }
    # --------------------------------------------------------------------------

    return;
}
################################################################################
#
# NOTE: No longer using $cRank; Forcing computation of ALL ranks from 
#       GI->S->SS->G->F->O->C->P
#
# ARGS: \%nonOverlappedCoverage, \%gi2mergedFrags, \%gi2mapFrags, \%genomeVitals
################################################################################
sub genTree {
    
    my %hitTree = ();

    foreach my $gi (keys %{ $_[0] }) {

        if(exists $_[0]->{$gi}->{TAXTREE}) {
        
            my $refsize       = $_[3]->{GI}->{$gi}->{SIZE};
            my $usize         = $_[0]->{$gi}->{USIZE};
            my $linlen        = $_[0]->{$gi}->{LINLEN};
            my $cov           = $_[0]->{$gi}->{COV};
            my $hitcount      = $_[1]->{$gi}->{HITCOUNT};
            my $fullReadHits  = 0; #$_[2]->{$gi}->{FULLREADHITS};
            my $totalBPmapped = $_[2]->{$gi}->{TOTALBPMAPPED};
        
            $hitTree{GI}->{$gi}->{REFSIZE}       = $refsize;
            $hitTree{GI}->{$gi}->{USIZE}         = $usize;
            $hitTree{GI}->{$gi}->{LINLEN}        = $linlen;
            $hitTree{GI}->{$gi}->{COV}           = $cov;
            $hitTree{GI}->{$gi}->{HITCOUNT}      = $hitcount;
            $hitTree{GI}->{$gi}->{FULLREADHITS}  = $fullReadHits;
            $hitTree{GI}->{$gi}->{MAPPED}        = $totalBPmapped;
            $hitTree{GI}->{$gi}->{HISTO}         = $giHisto{$gi};
            
            # Get {SS} and {S} SciNames
            my $repliconName = $_[0]->{$gi}->{REPLICON};
            my $strainName   = $_[0]->{$gi}->{STRAIN};
            my $speciesName = q{};
            my @rankNames = keys %{ $_[0]->{$gi}->{TAXTREE}->{S} };
            foreach (@rankNames) {
                $speciesName = $_ if($_[0]->{$gi}->{TAXTREE}->{S}->{$_} eq "scientific name");
                last if($speciesName);
            }
                                            
            my $genusName             = $_[0]->{$gi}->{TAXTREE}->{G};
            my $familyName            = $_[0]->{$gi}->{TAXTREE}->{F};
            my $orderName             = $_[0]->{$gi}->{TAXTREE}->{O};
            my $className             = $_[0]->{$gi}->{TAXTREE}->{C};
            my $phylumName            = $_[0]->{$gi}->{TAXTREE}->{P};
            $hitTree{GI}->{$gi}->{SS} = $strainName;
            $hitTree{GI}->{$gi}->{S}  = $speciesName;
            $hitTree{GI}->{$gi}->{G}  = $genusName;
            $hitTree{GI}->{$gi}->{F}  = $familyName;
            $hitTree{GI}->{$gi}->{O}  = $orderName;
            $hitTree{GI}->{$gi}->{C}  = $className;
            $hitTree{GI}->{$gi}->{P}  = $phylumName;
            
            # Transfer contig length:freq values to each tax level
            foreach my $contigLength (keys %{ $giHisto{$gi} }) {
                my $freq = $giHisto{$gi}->{$contigLength};
                $hitTree{SS}->{$strainName}->{HISTO}->{$contigLength} += $freq;
                $hitTree{S}->{$speciesName}->{HISTO}->{$contigLength} += $freq;
                $hitTree{G}->{$genusName}->{HISTO}->{$contigLength}   += $freq;
                $hitTree{F}->{$familyName}->{HISTO}->{$contigLength}  += $freq;
                $hitTree{O}->{$orderName}->{HISTO}->{$contigLength}   += $freq;
                $hitTree{C}->{$className}->{HISTO}->{$contigLength}   += $freq;
                $hitTree{P}->{$phylumName}->{HISTO}->{$contigLength}  += $freq;
            }
            
            # Populate {SS} level
            $hitTree{SS}->{$strainName}->{REFSIZE}       += $refsize;
            $hitTree{SS}->{$strainName}->{USIZE}         += $usize;
            $hitTree{SS}->{$strainName}->{LINLEN}        += $linlen;
            $hitTree{SS}->{$strainName}->{HITCOUNT}      += $hitcount;
            $hitTree{SS}->{$strainName}->{MAPPED}        += $totalBPmapped;
            $hitTree{SS}->{$strainName}->{FULLREADHITS}  += $fullReadHits;
            $hitTree{SS}->{$strainName}->{CHILDREN}->{NAMES}->{$gi} = ();

            # Populate {S} level
            $hitTree{S}->{$speciesName}->{REFSIZE}       += $refsize;
            $hitTree{S}->{$speciesName}->{USIZE}         += $usize;
            $hitTree{S}->{$speciesName}->{LINLEN}        += $linlen;
            $hitTree{S}->{$speciesName}->{HITCOUNT}      += $hitcount;
            $hitTree{S}->{$speciesName}->{MAPPED}        += $totalBPmapped;
            $hitTree{S}->{$speciesName}->{FULLREADHITS}  += $fullReadHits;
            $hitTree{S}->{$speciesName}->{CHILDREN}->{NAMES}->{$strainName} = ();
        
            # Populate {G} level
            $hitTree{G}->{$genusName}->{REFSIZE}       += $refsize;
            $hitTree{G}->{$genusName}->{USIZE}         += $usize;
            $hitTree{G}->{$genusName}->{LINLEN}        += $linlen;
            $hitTree{G}->{$genusName}->{HITCOUNT}      += $hitcount;
            $hitTree{G}->{$genusName}->{MAPPED}        += $totalBPmapped;
            $hitTree{G}->{$genusName}->{FULLREADHITS}  += $fullReadHits;
            $hitTree{G}->{$genusName}->{CHILDREN}->{NAMES}->{$speciesName} = ();
        
            # Populate {F} level
            $hitTree{F}->{$familyName}->{REFSIZE}       += $refsize;
            $hitTree{F}->{$familyName}->{USIZE}         += $usize;
            $hitTree{F}->{$familyName}->{LINLEN}        += $linlen;
            $hitTree{F}->{$familyName}->{HITCOUNT}      += $hitcount;
            $hitTree{F}->{$familyName}->{MAPPED}        += $totalBPmapped;
            $hitTree{F}->{$familyName}->{FULLREADHITS}  += $fullReadHits;
            $hitTree{F}->{$familyName}->{CHILDREN}->{NAMES}->{$genusName} = ();
        
            # Populate {O} level
            $hitTree{O}->{$orderName}->{REFSIZE}       += $refsize;
            $hitTree{O}->{$orderName}->{USIZE}         += $usize;
            $hitTree{O}->{$orderName}->{LINLEN}        += $linlen;
            $hitTree{O}->{$orderName}->{HITCOUNT}      += $hitcount;
            $hitTree{O}->{$orderName}->{MAPPED}        += $totalBPmapped;
            $hitTree{O}->{$orderName}->{FULLREADHITS}  += $fullReadHits;
            $hitTree{O}->{$orderName}->{CHILDREN}->{NAMES}->{$familyName} = ();
        
            # Populate {C} level
            $hitTree{C}->{$className}->{REFSIZE}       += $refsize;
            $hitTree{C}->{$className}->{USIZE}         += $usize;
            $hitTree{C}->{$className}->{LINLEN}        += $linlen;
            $hitTree{C}->{$className}->{HITCOUNT}      += $hitcount;
            $hitTree{C}->{$className}->{MAPPED}        += $totalBPmapped;
            $hitTree{C}->{$className}->{FULLREADHITS}  += $fullReadHits;
            $hitTree{C}->{$className}->{CHILDREN}->{NAMES}->{$orderName} = ();
        
            # Populate {P} level
            $hitTree{P}->{$phylumName}->{REFSIZE}       += $refsize;
            $hitTree{P}->{$phylumName}->{USIZE}         += $usize;
            $hitTree{P}->{$phylumName}->{LINLEN}        += $linlen;
            $hitTree{P}->{$phylumName}->{HITCOUNT}      += $hitcount;
            $hitTree{P}->{$phylumName}->{MAPPED}        += $totalBPmapped;
            $hitTree{P}->{$phylumName}->{FULLREADHITS}  += $fullReadHits;
            $hitTree{P}->{$phylumName}->{CHILDREN}->{NAMES}->{$className} = ();

        } #exists {TAXTREE}

        else {

            print "TAXTREE does not exist for GI \"$gi\"!\n";

        }
        
    } #gi
    
    # *!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
    #       UNDEF'ING UN-NEEDED LARGE VARIABLES TO FREE UP MEM
    # *!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
    # ARGS: \%nonOverlappedCoverage, \%gi2mergedFrags, \%gi2mapFrags, \%genomeVitals
    #
    undef $nonOverlappedCoverage;
    undef $gi2mergedFrags;
    undef $gi2mapFrags;
    undef $genomeVitals;
    # *!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!



    # Recursively (1) calculate the COV for each {SS},{S},{G},{F},{O},{C},{P}
    # and (2) store the BEST subrank for each hitTreeRank by their COV
    foreach my $rankIdx (1..$#hitTreeRanks) {
        my $hitTreeRank  = $hitTreeRanks[$rankIdx];
        my $prevTreeRank = $hitTreeRanks[$rankIdx-1];
        
        #print "rank = ".$hitTreeRank."\t".$prevTreeRank."\n";

        # Calculate COV for each entry
        foreach my $rankName (keys %{ $hitTree{$hitTreeRank} }) {

            $hitTree{$hitTreeRank}->{$rankName}->{COV} 
                = $hitTree{$hitTreeRank}->{$rankName}->{LINLEN} / $hitTree{$hitTreeRank}->{$rankName}->{USIZE};
            
            # Select BEST COV for each entry
            my %best = ();
            while(my($childName,undef) = each %{ $hitTree{$hitTreeRank}->{$rankName}->{CHILDREN}->{NAMES} }) {       # Get first entry and set as BEST (so we don't have to call exists in the loop below)
                $best{NAME} = $childName;
                $best{COV}  = $hitTree{$prevTreeRank}->{$childName}->{COV};
                last;
            }
            
            foreach my $prevRankName (keys %{ $hitTree{$hitTreeRank}->{$rankName}->{CHILDREN}->{NAMES} }) {
            
                #print "Looking for ".$hitTreeRank."\'s prevRankNameCOV in $prevTreeRank...";
                #print "FOUND!\n" if(exists $hitTree{$prevTreeRank}->{$prevRankName}->{COV});
                
                my $prevRankNameCOV = $hitTree{$prevTreeRank}->{$prevRankName}->{COV};
                if($prevRankNameCOV > $best{COV}) {
                    $best{COV}  = $prevRankNameCOV;
                    $best{NAME} = $prevRankName;
                }
            }
            
            # Set the BEST COV for each entry
            $hitTree{$hitTreeRank}->{$rankName}->{CHILDREN}->{BEST} = $best{NAME};
            
        }
        
    } #rankIdx

    # Recursively find the number of genome projects under each taxonomic rank, starting
    # from {S}->{G}->...-> {P}
    foreach my $rankIdx (2..$#hitTreeRanks) {
        my $currTreeRank = $hitTreeRanks[$rankIdx];
        my $prevTreeRank = $hitTreeRanks[$rankIdx-1];
        
        # Handle {S} separately since each child's name is essentially a genome project
        if($currTreeRank eq "S") {
            foreach my $rankName (keys %{ $hitTree{$currTreeRank} }) {
                $hitTree{$currTreeRank}->{$rankName}->{GPROJ} = 0 unless exists ($hitTree{$currTreeRank}->{$rankName}->{GPROJ});
                my $numGenomeProjects = scalar(keys %{ $hitTree{$currTreeRank}->{$rankName}->{CHILDREN}->{NAMES} });
                $hitTree{$currTreeRank}->{$rankName}->{GPROJ} = $numGenomeProjects;
            }
        }
        # Handle {G} on up...
        else {
            foreach my $rankName (keys %{ $hitTree{$currTreeRank} }) {
                $hitTree{$currTreeRank}->{$rankName}->{GPROJ} = 0 unless exists ($hitTree{$currTreeRank}->{$rankName}->{GPROJ});
                foreach my $childRankName (keys %{ $hitTree{$currTreeRank}->{$rankName}->{CHILDREN}->{NAMES} }) {
                    $hitTree{$currTreeRank}->{$rankName}->{GPROJ} += $hitTree{$prevTreeRank}->{$childRankName}->{GPROJ};
                }
            }
        }
    }
  
    store \%hitTree, $outdir."/hitTree.dmp";
    open my $HITTREEFILE, '>', $outdir."/hitTree.readable.dmp";
    print $HITTREEFILE Dump(\%hitTree)."\n";
    close $HITTREEFILE;
  
    return \%hitTree;
}
################################################################################
# ARGS: $parsedDB
################################################################################
sub retrieveParsedDB {
    
    STDOUT->autoflush(1);

    print "->Retrieving parsed DB from disk [".$_[0]."]...";
    my $iter0 = new Benchmark;
    my $gi2seqlen = retrieve $_[0];
    my $iter1 = new Benchmark;
    print "done. ".timestr(timediff($iter1,$iter0))."\n";
    
    return $gi2seqlen;
}
################################################################################
# FUNCTION: opens FASTA version of uniques DB, parses each FASTA sequence, and
#           determines the length of each fragment
# ARGS: DTAR $dbFile
################################################################################
sub calculateUniqueLengthsByGI {

    my $dtarDB = $_[0];
    
    # Parse FASTA file by GI, summing and storing lengths
    my %gi2seqlen = ();                 # KEY:gi    VAL:total_length_of_strain_in_uniqueDB
    print "->Parsing DB file [$dtarDB]...";
    my %fasta = %{ parseMultiFASTA($dtarDB) };
    #print $_," => ",length($fasta{$_}),"\n" foreach (keys %fasta);
    print "done.\n";
    
    my $iter0 = new Benchmark;
    print "->Associating GI with FASTA fragments...";
    foreach my $header (keys %fasta) {

        # *** NOTE: $header must include "$START|$STOP|$LENGTH|"
        #                      _gi_
        #                 ____________________giFRAG___________________________
        #if($header =~ m/^(gi\|(\d+)\|\S+\|(\w+(?:\.\d+))\|(\d+)\|(\d+)\|(\d+)\|)/) {
        if($header =~ m/^(gi\|(\d+)\|\S+\|(\S+)\|(\d+)\|(\d+)\|(\d+)\|)/) {
            my $giFrag = $1;            # get GI fragment (subset of GI entries)
            my $gi     = $2;            # get GI
            
            # Index
            #my $giFragLength = length($fasta{$header});
            #$gi2seqlen{$gi}->{FRAG}->{$giFrag} = $giFragLength; # length of unique frag
            #$gi2seqlen{$gi}->{LENGTH}         += $giFragLength; # total length of gi
            $gi2seqlen{$gi}->{FRAG}->{$giFrag} = $6; # length of unique frag
            $gi2seqlen{$gi}->{LENGTH}         += $6; # total length of gi
            $gi2seqlen{$gi}->{COUNT}++;              # no. of frags

        } #FIELDS[2]
        else {
            print "**Unrecognized FASTA entry!**    [$header]\n";
        }
        
    } #HEADER

    my $iter1 = new Benchmark;
    print "done. ".timestr(timediff($iter1,$iter0))."\n";

    # Store parsed DB
    my $outfile = $prefix ? "$outdir/$prefix.parsedGOTTCHA.dmp" : "$outdir/parsedGOTTCHA.dmp";
    $iter0 = new Benchmark;
    print "->Storing parsed version of DB to disk [$outfile]...";
    store \%gi2seqlen, $outfile;
    $iter1 = new Benchmark;
    print "done. ".timestr(timediff($iter1,$iter0))."\n";

    return \%gi2seqlen;    
}
################################################################################
#          $_[0]            $_[1]          $_[2]
# ARGS: \%giLengths, \%uniqueGIlengths, \%genomeVitals
################################################################################
sub calculateNonOverlappedCoverage {
    
    my %nonOverlappedCoverage = ();
    my $count = 0;

    STDOUT->autoflush(1);

    my $iter0 = new Benchmark;
    print "->Calculating non-overlapping coverage from mapping results...";
    foreach my $gi (keys %{ $_[0] }) {
        $nonOverlappedCoverage{$gi}->{LINLEN}   = $_[0]->{$gi};                         # LINEAR_LENGTH for GI
        $nonOverlappedCoverage{$gi}->{USIZE}    = $_[1]->{$gi}->{LENGTH};               # UNIQUE_REFDB_LINEAR_LENGTH for GI
        $nonOverlappedCoverage{$gi}->{COV}      = $_[0]->{$gi}/$_[1]->{$gi}->{LENGTH};  # COVERAGE for GI
        $nonOverlappedCoverage{$gi}->{REPLICON} = $_[2]->{GI}->{$gi}->{REPL};           # REPLICON NAME
        $nonOverlappedCoverage{$gi}->{STRAIN}   = $_[2]->{GI}->{$gi}->{ORG};            # STRAIN NAME
#            if(!$noVitalsFile && exists($_[2]->{GI}->{$gi}->{REPL}));
        if ($count % $lineBuffer == 0) {
            print "." if($verbose);
        }
        $count++;
    }
    my $iter1 = new Benchmark;
    print "done. ".timestr(timediff($iter1,$iter0))."\n";
    
    return \%nonOverlappedCoverage;
}
################################################################################
# FUNCTION: (1) Calls removeStacking(\%startStops, $giOriginalLength) on each
#               $giFrag entry to get its unique length.
#           (2) Writes contig start:stop coords to file
#           (3) Calculates histogram string by $giFrag
# ARGS: \%gi2mergedFrags, \%genomeVitals, $gi2mapFrags, \%giFragHisto, \%giHisto
################################################################################
sub sumMappedLengths {

    my %giLengths = ();
    my %giCoords  = ();     # $gi => { $giFrag => { $start => $stop } }
                            # Stored to disk and reloaded later for tax-dependent output

    my $coordsIter0 = new Benchmark;
    my $outFilename1 = $outdir."/".$prefix.".replicon.contig.coords.csv";
    print "->Storing coordinates to disk [".$outFilename1."]...";
    open my $COORDSFILE1, '>', $outFilename1;
    
    foreach my $gi (keys %{ $_[0] }) {
        
        # Make sure the GI is in the reference file so we can get its size
        if(!exists $_[1]->{GI}->{$gi}) {
            print "WARNING: GI $gi not found in GENOME VITALS file! Cannot profile this entry.\n";
            next;
        }
        my $giMappedLength = 0;
        my $giOrigLen = $_[1]->{GI}->{$gi}->{SIZE};

        while (my ($giFrag, $startStops) = each %{ $_[0]->{$gi}->{FRAG} }) {
            
            # Datastructure: \%contiguous = ( $start => $stop )
            my $contiguous = removeStacking($startStops, $giOrigLen);

#            print $COORDSFILE $_[2]->{$gi}->{FRAG}->{$giFrag};
            print $COORDSFILE1 $giFrag;

            # Length = Stop - Start + 1
            #$giMappedLength += $contiguous->{$_} - $_ + 1 foreach (keys %{ $contiguous});
            foreach my $contigStart (sort { $a <=> $b } keys %{ $contiguous }) {
                my $contigStop   = $contiguous->{$contigStart};
                my $contigLength = $contigStop - $contigStart + 1;
                $giMappedLength += $contigLength;
                $_[3]->{$giFrag}->{$contigLength}++;                    # %giFragHisto
                $_[4]->{$gi}->{$contigLength}++;                        # %giHisto
                $giCoords{$gi}->{$giFrag}->{$contigStart} = $contigStop;
                print $COORDSFILE1 ",".$contigStart.":".$contigStop;
            }
            print $COORDSFILE1 "\n";            
            
        }
        $giLengths{$gi} = $giMappedLength;
    }
    close $COORDSFILE1;
    my $coordsIter1 = new Benchmark;
    print "done. ".(timestr(timediff($coordsIter1, $coordsIter0)))."\n";

    # ARGS: \%hash, $filename, $description, $mode      $mode = "bin" or "ascii"
    dump2disk(\%giCoords, $outdir."/".$prefix.".giCoords.dmp", "GI coordinates", "bin");    

    #---------------------------------------------------------------------------
    # HISTOFORMAT:      GI,length1:freq1,length2:freq2,...
    dump2disk($_[3], $outdir."/".$prefix.".replicon.contig.HistoByEntry.dmp", "contig length histogram (by entry)", "bin");
    dump2disk($_[4], $outdir."/".$prefix.".replicon.contig.HistoByGI.dmp", "contig length histogram (by GI)", "bin");
    #---------------------------------------------------------------------------
    my $histoIter0 = new Benchmark;
    print "->Storing parseable contig length histogram(s) to disk...";

    # Store a string histogram of hits by giFrag (i.e. by each entry in the GOTTCHA DB)
    my $histoFilename1 = $outdir."/".$prefix.".replicon.contig.HistoByEntry.csv";
    open my $HISTOFILE1, '>', $histoFilename1;
    foreach my $giFrag (keys %{ $_[3] }) {
        print $HISTOFILE1 $giFrag;
        foreach my $contigLength (sort {$a <=> $b} keys %{ $_[3]->{$giFrag} }) {
            print $HISTOFILE1 ",".join(":", $contigLength, $_[3]->{$giFrag}->{$contigLength});
        }
        print $HISTOFILE1 "\n";
    }
    close $HISTOFILE1;
    #---------------------------------------------------------------------------
    # Store a string histogram of hits by GI only
    my $histoFilename2 = $outdir."/".$prefix.".replicon.contig.HistoByGI.csv";
    open my $HISTOFILE2, '>', $histoFilename2;
    foreach my $gi (keys %{ $_[4] }) {
        print $HISTOFILE2 $gi;
        foreach my $contigLength (sort {$a <=> $b} keys %{ $_[4]->{$gi} }) {
            print $HISTOFILE2 ",".join(":", $contigLength, $_[4]->{$gi}->{$contigLength});
        }
        print $HISTOFILE2 "\n";
    }
    close $HISTOFILE2;
 
    my $histoIter1 = new Benchmark;
    print "done. ".(timestr(timediff($histoIter1,$histoIter0)))."\n";

    return \%giLengths;     # ( $gi => $linearLength )
}
################################################################################
# Function: Each $startPos in %startStops represents a fragment that extends forward 
#           until $stopPos in a sequence that's $seqLen long. This sub() will
#           condense overlapping fragments into as many non-overlapping fragments
#           (each defined by their own START|STOP positions) as possible. If a 
#           fragment extends beyond the provided $seqLen, it is assumed that this
#           represents a circularize sequence and the fragment will wrap around
#           to the start of the sequence, and again merge anything that's overlapping.
#           The (START|STOP) coordinates for each contiguous fragment is returned
#           in a hash, where the KEYS are the START positions, and VALUES, the STOPs.
#
#           In the current application, these (START|STOP) pairs correspond to the
#           regions of a sequence that are common with at least 1 other sequence.
#           It is assumed that these (START|STOP) positions will be used to identify 
#           unique regions in the sequence (of length $seqLen). This is accomplished 
#           by determining the *complement* of the current START|STOP positions. This
#           function will be performed by some other subroutine.
#
# ARGS: $_[0] = \%{ $lmerCoords{$org} } = ( $start1 => $stop1, $start2 => $stop2, ... );
#       $_[1] = $seqLen
# Returns:  Nothing. The global var %contiguous is updated within the sub.
################################################################################
sub removeStacking {

    my $contiguous_href = ();
    my $seqLen          = $_[1];

    # ---------------------------------------------    
    # generate @unsortedStarts and @unsortedStops
    #     GENERAL case: stop positions explicitly specified in hash as VALUES
    # SPECIALIZED case: stop positions are all a fixed distance from starts
    # ---------------------------------------------
    my @unsortedStarts = keys %{ $_[0] };
    my @unsortedStops  = values %{ $_[0] };                                 # Generalized case
    #my @unsortedStops  = ();
    #push(@unsortedStops, $_ + $currWordSize) foreach (@unsortedStarts);    # Specialized case

    # ---------------------------------------------    
    # send to sub mergeStarts()
    # ---------------------------------------------    
    # Given:   \@unsortedStarts, \@unsortedStops
    # Returns: \%contiguous
    $contiguous_href = mergeStarts(\@unsortedStarts, \@unsortedStops);

    # ---------------------------------------------    
    # process circular ends
    # ---------------------------------------------    
    # Records all the start positions whose stop positions are out of bounds
    #print "Looking for all STOPS that are outside position $seqLen...\n";
    my @outOfBoundsStarts = grep { $contiguous_href->{$_} >= $seqLen } keys %{ $contiguous_href };
    my @outOfBoundsStops  = ();
    push(@outOfBoundsStops, $contiguous_href->{$_}) foreach (@outOfBoundsStarts);

    #print "THERE ARE ".(scalar(@outOfBoundsStarts))." out-of-bounds STARTS\n";

    # Remove the out-of-bounds entries
    delete $contiguous_href->{$_} foreach (@outOfBoundsStarts);

    # Re-initialize unsorted starts and stops
    @unsortedStarts = keys %{ $contiguous_href };
    @unsortedStops  = values %{ $contiguous_href };
        
    foreach my $idx (0..$#outOfBoundsStarts) {
        my $start = $outOfBoundsStarts[$idx];
        my $stop  = $outOfBoundsStops[$idx];
        
        # Circularize the out-of-bounds (START|STOP) pair by splitting it up into 
        # two separate pairs: 
        #     (START | seqLen-1)
        #     (0     | STOP-seqLen+1)
        unshift(@unsortedStarts,                0);     # Add "0" to front of StartArray
        unshift(@unsortedStops, ($stop-$seqLen-1));     # Add stopPos to same pos in StopArray
        push(@unsortedStarts,     $start);              # Add startPos to end of StartArray
        push(@unsortedStops, ($seqLen-1));              # Add stopPos to end of StopArray
     
        #print "Adding \"0\" to BEGINNING of unsortedSTARTS\n";
        #print "Adding \"".($stop-$seqLen-1)."\" to BEGINNING of unsortedSTOPS\n";
        #print "Adding \"$start\" to end of unsortedSTARTS\n";
        #print "Adding \"".($seqLen-1)."\" to end of unsortedSTOPS\n";
    }
    
    # ---------------------------------------------    
    # generate @unsortedStarts and @unsortedStops
    # ---------------------------------------------    
    # Given:   \@unsortedStarts, \@unsortedStops
    # Returns: \%contiguous
    $contiguous_href = mergeStarts(\@unsortedStarts, \@unsortedStops);

    # Crap In = Crap Out, so we check for blatant errors
    foreach (sort {$a <=> $b} keys %{ $contiguous_href }) {
        stat_log("**WARNING**: STOP position (".$contiguous_href->{$_}.") is >>OUTSIDE<< sequence string (Len=$seqLen)!") 
            if($contiguous_href->{$_} > ($seqLen-1));
    }
    
    return $contiguous_href;
}
################################################################################
# ARGS:   \@unsortedStarts, \@unsortedStops
# Return: contiguous start|stop positions
################################################################################
sub mergeStarts {

    # %contiguous: contiguous fragments of DNA demarcated by the contained START | STOP positions
    my %contiguous   = ();
    my @sortedStarts = ();
    my @sortedStops  = ();

    # Sort start positions in increasing order and store these indices in sorted order
    #my @sortedStartIndices = sort { $unsortedStarts[$a] <=> $unsortedStarts[$b] } 0..$#unsortedStarts; 
    #push(@sortedStarts, $unsortedStarts[$_]) foreach (@sortedStartIndices);
    #push(@sortedStops,  $unsortedStops[$_])  foreach (@sortedStartIndices);
    my @sortedStartIndices = sort { $_[0]->[$a] <=> $_[0]->[$b] } 0..(scalar(@{ $_[0] })-1); 
    push(@sortedStarts, $_[0]->[$_]) foreach (@sortedStartIndices);
    push(@sortedStops,  $_[1]->[$_]) foreach (@sortedStartIndices);

    my $bestStart = $sortedStarts[0];
    my $bestStop  = $sortedStops[0];
    
    START: foreach my $idx (0 .. (scalar(@sortedStarts)-1)) {
        my $start = $sortedStarts[$idx];
        my $stop  = $sortedStops[$idx];
        
        my $nextStart = defined($sortedStarts[$idx+1]) ? ($sortedStarts[$idx+1]) : (-1);
        my $nextStop  = defined($sortedStops[$idx+1])  ? ($sortedStops[$idx+1])  : (-1);
        
        if($nextStart == -1) {
            $contiguous{$bestStart} = $bestStop;
            last START;
        }
        else {
            if($nextStop <= $bestStop) {        # self-contained; skip this one
                next START;
            }
            else {                              # extension or disjointed; check start pos
                if($nextStart <= ($bestStop+1)) {   # extension: include immediately adjacet starts as well
                    $bestStop = $nextStop;
                    next START;
                }
                else {                          # disjointed
                    $contiguous{$bestStart} = $bestStop;    # store
                    $bestStart = $nextStart;                # update
                    $bestStop  = $nextStop;                 # update
                }
            }
        }
    }

    return \%contiguous;
    
}
################################################################################
# FUNCTION: Consolidates all the { start => { stop1 => count, stop2 => count, ... } }
#           entries into a single { start => stop } pair
# ARGS: \%gi2mapFrags
################################################################################
# Find all uniquely covered regions of each entry in the uniques DB
sub mergeOverlappingHits2 {                                                      # <-------------------- calculate total DEPTH in this sub()
    my %gi2mergedFrags = ();

    my $iter0 = new Benchmark;
    print "->Consolidating hits...";
    open my $DETAILFILE, '>', $detailFile;
    
    while (my ($gi,undef) = each %{ $_[0] }) {
#        print $DETAILFILE $gi;
        while (my ($frag,undef) = each %{ $_[0]->{$gi}->{FRAG} }) {
#            print $DETAILFILE "\t".$frag;
            print $DETAILFILE $frag;
            while (my ($start, $stophash) = each %{ $_[0]->{$gi}->{FRAG}->{$frag} }) {
                my $maxStop = 0;
                foreach (keys %{ $stophash }) {
                    $maxStop= $_ if($_ > $maxStop);
                }
                $_[0]->{$gi}->{FRAG}->{$frag}->{$start} = $maxStop;
                
                # Output the specific START|STOP coords to file
                # NOTE: these coords are relative to the START|STOP coords of the
                #       unique fragment as listed in $frag
#                print $DETAILFILE "\t".$start."|".$maxStop;
                print $DETAILFILE ",".$start.":".$maxStop;
            }
            print $DETAILFILE "\n";
        }
    }
    
    close $DETAILFILE;

    my $iter1 = new Benchmark;
    print "done. ".(timestr(timediff($iter1,$iter0)))."\n";
    
    return \%{ $_[0] };         # Return reference to input href
}
################################################################################
# FUNCTION: Consolidates all the { start => { stop1 => count, stop2 => count, ... } }
#           entries into a single { start => stop } pair
# ARGS: \%gi2mapFrags
################################################################################
# Find all uniquely covered regions of each entry in the uniques DB
sub mergeOverlappingHits {                                                      # <-------------------- calculate total DEPTH in this sub()

    my %gi2mergedFrags = ();
    
    # Foreach $giFrag, merge all the $rStartX and $rStopY's into as many contiguous START | STOP
    # pairs as possible
    my $counter = 0;
    
    foreach my $gi (keys %{ $_[0] }) {
        foreach my $giFrag (keys %{ $_[0]->{$gi}->{FRAG} }) {
        
            # Grab all starts for this current unique frag, sorted in ascending order
            my @starts = sort {$a <=> $b} keys %{ $_[0]->{$gi}->{FRAG}->{$giFrag} };
            if($debug) {
                $counter++;
                print "\n($counter)";
                print "$giFrag\tSTARTS|STOPS:\t"; # if(scalar(@starts) > 1);
                print $_."->".join(",",(keys %{ $_[0]->{$gi}->{FRAG}->{$giFrag}->{$_} })) 
                    foreach (@starts); 
                <STDIN>;
            }

            my $offset = 1;
            my $cStart = 0;           
            my $cStop  = 0;
            my $nStart = 0;
            my $nStop  = 0;

            for (my $idx = 0; $idx < scalar(@starts); $idx += $offset) {

                $offset = 1;                                                    # Reset
                $cStart = $starts[$idx];
                my @stops = keys %{ $_[0]->{$gi}->{FRAG}->{$giFrag}->{$cStart} };
                $cStop = getMax(\@stops);
                
                if($debug) {
                    print "cStart=$cStart, cStop=$cStop     [idx=$idx, offset=$offset]"; 
                    <STDIN>;
                }

                my $extend = 0;

                do {
                
                    # Check for addtional starts
                    if(defined $starts[$idx+$offset]) {
                        # get nStart|nStop
                        my $nStart = $starts[$idx+$offset];
                        my @nStops = keys %{ $_[0]->{$gi}->{FRAG}->{$giFrag}->{$nStart} };
                        my $nStop  = getMax(\@nStops);
                        if($debug) {
                            print "nStart=$nStart, nStop=$nStop     "
                                 ."(cStart=$cStart, cStop=$cStop)   "
                                 ."  [idx=$idx, offset=$offset, extend=$extend]"; 
                            <STDIN>;
                        }

                        # determine extension
                        if($nStart <= $cStop) {
                            # extend() if $extend
                            if($nStop > $cStop) {
                                # EXTENSION
                                $cStop  = $nStop;    # Update
                                $extend = 1;
                            }
                            #else {
                            #    # NON-EXTENSION
                            #}
                            $offset++;              # Look ahead
                        }
                        else {
                            # new fragment bc. nStart exists
                            print "-->STORING as NEW FRAGMENT: $cStart | $cStop\n" if($ddebug && ($gi == "126373828"));
                            $gi2mergedFrags{$gi}->{FRAG}->{$giFrag}->{$cStart} = $cStop;    # Store
                            $cStart = $nStart;                                              # Update 
                            $cStop  = $nStop;                                               # Update
                        }
                    } #IF DEFINED nStart
                    else {
                            # NO FURTHER EXTENSION OR NEW FRAGMENTS POSSIBLE; DUMP DATA TO HASH
                            print "-->STORING as LAST of EXT: $cStart | $cStop\n" if($ddebug && ($gi == "126373828"));
                            $gi2mergedFrags{$gi}->{FRAG}->{$giFrag}->{$cStart} = $cStop;
                            $extend = 0;
                    }
                    
                } while ($extend);
                
                # Store if last IDX
                print "-->STORING as LAST IDX: $cStart | $cStop\n" if($ddebug && ($gi == "126373828"));
                $gi2mergedFrags{$gi}->{FRAG}->{$giFrag}->{$cStart} = $cStop;

            } #FOR IDX

        } #GIFRAG
        
        # Record no. of actual hits to the reference GI
        $gi2mergedFrags{$gi}->{HITCOUNT} = $_[0]->{$gi}->{HITCOUNT};  
        
    } #GI
    
    return \%gi2mergedFrags;
}
################################################################################
# ARGS: \@integers
################################################################################
sub getMax {
    my $max = $_[0]->[0];           # assign first element
    foreach (@{ $_[0] }) {
        $max = $_ if($_ > $max);    # update with larger value
    }
    return $max;
}
################################################################################
#                (\w+(?:\.\d+))
#gi|253972022|gb|CP000819.1|790319|796620|Escherichia|Escherichia coli B str. REL606
#gi|218425442|emb|CU928162.2|2467207|2467240|Escherichia|Escherichia coli ED1a chromosome
#gi|257757386|dbj|AP010958.1|2879025|2879059|Escherichia|Escherichia coli O103:H2 str. 12009 DNA 
#
# ARGS: $mapFile (string)
################################################################################
sub parseMapFile {
    return my $href = $nucmer ? parseNucmer($_[0])
                    : $sam    ? parseSAM($_[0])
                    :           parseBLAST($_[0]);
}
################################################################################
#  0   1   2   3     4     5     6      7      8      9     10    11       12
# S1  E1  S2  E2  LEN1  LEN2  %IDY  LEN_R  LEN_Q  COV_R  COV_Q  TAGS
#                                                               [REF]   [QRY]
#
# S1    is the start position in the GI_frag found in the mapping (wrt GI_frag coords)
# E1    is the end   position in the GI_frag found in the mapping (wrt GI_frag coords)
# LEN1  is the mapped length of the GI_frag
# LEN_R is the entire length of the GI_Frag (in the uniques DB)
#
# *Record S1|E1 (out of LEN_R) for each hit found; condense
# 
#   (LEN_R)
#
# ARGS: $mapFile (string)
#-------------------------------------------------------------------------------
sub parseNucmer {

    my %gi2mapFrags = ();

    my $lineCount = 0;
    open my $NUCMER, '<', $_[0] || die "Cannot open DB file \"".$_[0]."\" for read!\n";
    
    my $iter0 = new Benchmark;
    print "->Parsing NUCMER mapping file [".$_[0]."] (each \'.\' = $lineBuffer lines)";
    
    LINE: while(my $line=<$NUCMER>) {
        $lineCount++;
        print "." if($lineCount % $lineBuffer == 0);
        
        #------------ Verify NUCMER coords file is actually given --------------
        if(($lineCount == 2) && ($line ne "NUCMER\n")) {
            print "*****File \"".$_[0]."\" does not appear to be a NUCMER coords file***** ...";
            for (my $idx=$countdown; $idx > 0; $idx--) {
                print "$idx..";
                sleep 1;
            }
            print "\n";
        }
        #-----------------------------------------------------------------------
        next LINE if($lineCount < 5);       # Data starts on line 5
        chomp $line;

        my @fields = split(/\t/,$line);
    
        if($fields[6] >= $minPctID) {
            # Do work here...

            # Get GI out of the REF TAG
            my $gi = $1 if($fields[11] =~ m/gi\|(\d+)\|/);

            #print "FIELD[11]=".$fields[11]."\n";

            # Get GI fragment out of the REF TAB
            if($fields[11] =~ m/^(gi\|(\d+)\|\S+\|(\w+(?:\.\d+))\|(\d+)\|(\d+)\|(\d+)\|)/) {
                my $giFrag     = $1;
                my $rStart     = $fields[0];
                my $rStop      = $fields[1];
                my $fragStart  = $4;
                my $fragStop   = $5;
                
                #print "GI: $2\t  GI_FRAG: ".$giFrag." starts\@$fragStart, stops\@$fragStop\n";
                #print "GI: $2\t  GI_FRAG: ".$giFrag." starts\@$rStart, stops\@$rStop\n";
                
                # Index the
                $gi2mapFrags{$gi}->{FRAG}->{$giFrag}->{$rStart}->{$rStop}++;
                $gi2mapFrags{$gi}->{HITCOUNT}++;

                # This is for the parseDB portion
                # $gi2mapFrags{$gi}->{$giFrag} = ($stop-$start+1);
                
                
            }
            #$rID2len{$gi}         += $fields[4];
            #$qID2len{$fields[12]} += $fields[4];
        }
    
    }
    close $NUCMER;
    
    my $iter1 = new Benchmark;
    print "done. ".timestr(timediff($iter1,$iter0))."\n";         # parsing SAM file

    return \%gi2mapFrags;
}
################################################################################
# ARGS: $mapFile (string)
################################################################################
sub parseSAM {

    my %gi2mapFrags = ();

    # -------------------------------------------    
    # Open list file
    # -------------------------------------------    
    my @samFiles = ();
    my @missing  = ();

    if( $_[0] ne "-" ){
        open my $LISTFILE, '<', $_[0] || die "Cannot open ".$_[0]." for reading!\n";
        while(my $line=<$LISTFILE>) {
            next if($line eq "\n");
            chomp $line;
            if(-e $line) {
                push(@samFiles, $line);
            }
            else {
                push(@missing, $line);
            }
        }
        close $LISTFILE;
        if(@missing) {
            print "*FATAL* The following SAM files do not exist:\n";
            print "=============================================\n";
            print "    $_\n" foreach(@missing);
            exit;
        }
    }
    else{
        @samFiles = ("-");
    }

    STDOUT->autoflush(1);

    # -------------------------------------------    
    # -------------------------------------------    
    my $mappedTrimmedReadsFilename   = $outdir."/".$prefix.".MAPPED.trimmed.fastq";
    my $unmappedTrimmedReadsFilename = $outdir."/".$prefix.".UNMAPPED.trimmed.fastq";

    my $MAPPED_TRIMMEDREADS_CNT   = 0;
    my $UNMAPPED_TRIMMEDREADS_CNT = 0; 
    my $MAPPED_TRIMMEDREADS   = q{};
    my $UNMAPPED_TRIMMEDREADS = q{};
    open $MAPPED_TRIMMEDREADS, '>', $mappedTrimmedReadsFilename unless ($noMappedFastq);
    open $UNMAPPED_TRIMMEDREADS, '>', $unmappedTrimmedReadsFilename unless ($noUnmappedFastq);
    
    # LOOP: Parse each SAM file (single-threaded)
    foreach(@samFiles) {
        my $SAMFILE;
        my $SAMFILE_is_STDIN=0;

        #open STDIN if samfile name is "-"
        if( $_ ne "-"){
            open $SAMFILE, '<', $_ || die "Cannot open SAM file \"".$_."\" for read!\n";
        } 
        else{
            $SAMFILE = *STDIN;
            $SAMFILE_is_STDIN = 1;
        }

        print "->Parsing SAM file [".$_."] "; #(each \'.\' = $lineBuffer lines)";
        my $iter0 = new Benchmark;

        ## ------------------------------------
        ## Skip all header lines, if they exist
        ## ------------------------------------
        #my $line = q{};
        #while($line=<$SAMFILE>) {
        #    last if($line !~ /^\@/);
        #}
        #seek($SAMFILE, -length($line), 1);      # Unread the line (req's chars = 8-bit)

        
        # ------------------------------------
        # Process all of the following entries
        # ------------------------------------
        LINE: while(my $line=<$SAMFILE>) {
            next LINE if($line =~ /^@/);        # no previous chomps
            next LINE if($line eq "\n");        # no previous chomps

            my @fields = split(/\t/,$line);
            if( scalar @fields < 11 ) {
               print "*FATAL* The following files are invalid SAM files:\n";
               print "==================================================\n";
               print "    $line\n";
               exit;
           }

            #print "LINE=".$line; <STDIN>;

            # ===================
            # Mapping entry found
            # ===================
            # The 2nd field in the alignment ($fields[1]) specifies the alignment code, i.e. whether or 
            # not an alignment has taken place. A set 3rd bit signifies an unmapped read.
            #next LINE if($fields[1] & hex(0b100));  # Checking if 3rd bit is set
            if($fields[1] & hex(0b100)) {
                #if($writeUnmapped) { print $UNMAPPED_TRIMMEDREADS "@".$fields[0]."\n".$fields[9]."\n"."+\n".$fields[10]."\n"; }
                print $UNMAPPED_TRIMMEDREADS "@".$fields[0]."\n".$fields[9]."\n"."+\n".$fields[10]."\n" unless ($noUnmappedFastq);
                $UNMAPPED_TRIMMEDREADS_CNT++;
                next LINE;
            }

###############################
# EDIT: OMITTING FULLREADHITS
###############################
#            my $lastMatch = q{};
#            my $lastRead  = q{};

            # Check field containing the number of matched bases
            #if($fields[5] =~ m/(\d+)M/) {
            ####################
            # EDIT: 2013-07-17
            ####################
            if($fields[5] =~ m/^(\d+)M$/) {
                my $lengthMatch = $1;
                $MAPPED_TRIMMEDREADS_CNT++;

                # Mapping entry conforms
                if($fields[2] =~ m/^(gi\|(\d+)\|\S+\|(\w+(?:\.\d+))\|(\d+)\|(\d+)\|)/) {
                #if($fields[2] =~ m/^(gi\|(\d+)\|\S+\|\S+\|(\d+)\|(\d+)\|(\d+)\|)/) {         # <------ Req's all parsedDTAR.dmp files to be re-parsed & stored using --make_dmp and --db=<DBFILE>*
                    my $giFrag = $1;            # get GI fragment (subset of GI entries)
                    my $gi     = $2;            # get GI
                    my $rStart = $fields[3];    # pos in REF where mapping starts

                    ##################################################################################################
                    # In order to track the no. of unique full-length reads that mapped,
                    # we have to track the unique header, minus the trailing "_".$counter."/".$end
                    # where $end = {1, 2}. 
                    #   @HISEQ:148:C0F5FACXX:1:1307:2735:88518_0/1 1:N:0:CGATGT 1:30      // Header line
                    #   @HISEQ:148:C0F5FACXX:1:1307:2735:88518_0/1                        // Trimmed read header
                    #   @HISEQ:148:C0F5FACXX:1:1307:2735:88518                            // Full-length read header
                    # If there are 300M reads, then tracking the entire full-length 
                    # header as a primary key becomes slow -- Perl hashes seem to slow down 
                    # around 350,000 primary keys. So we should split them into prefixes
                    # and suffixes. However, if a single FASTQ file contains only ONE
                    # Illumina lane, then the first 22 or so chars will be identical, meaning
                    # that 300M reads will be under a single hash key. For the worse case
                    # scenario, we assume only 1 lane 
                    ##################################################################################################

###############################
# EDIT: OMITTING FULLREADHITS
###############################
#                    if($fields[0] =~ m/(.+)_\d+\/\d+/) {
#                        # Assuming all frags from the same full-length read are clustered
#                        # together, then a new $1 means a new full-length read.
#                        if($1 ne $lastMatch) {
#                            $gi2mapFrags{$gi}->{FULLREADHITS}++;
#                            $lastMatch = $1;
#                        }
#                    }
#                    else {
#                        #print "*WARNING*: Unrecognized FASTQ header! \"".$fields[0]."\"\n";
#                        if($fields[0] ne $lastMatch) {
#                            $gi2mapFrags{$gi}->{FULLREADHITS}++;
#                            $lastMatch = $fields[0];
#                        }
#                    }
#
#                    # Store mapped trimmed reads to FASTQ file only if not yet accounted for
#                    if($fields[0] ne $lastRead) {
                     #if($writeMapped) { print $MAPPED_TRIMMEDREADS "@".$fields[0]."\n".$fields[9]."\n"."+\n".$fields[10]."\n"; }
                     print $MAPPED_TRIMMEDREADS "@".$fields[0]."\n".$fields[9]."\n"."+\n".$fields[10]."\n" unless ($noMappedFastq);
#                        $lastRead = $fields[0];
#                    }

                    # A match start at pos=1728 that is 119-bp long will bring the stop
                    # to pos=1846. Recall that length=stop-start+1, so stop=length+start-1
                    my $rStop  = $rStart+$lengthMatch-1;

                    #print "\"".$fields[2]."\"\n" if(!exists $uniqueGIlengths->{$gi}->{FRAG}->{$giFrag});
                    if(!exists $uniqueGIlengths->{$gi} || !exists $uniqueGIlengths->{$gi}->{FRAG}->{$giFrag}) {
                        print "\n  **WARNING**: GOTTCHA database inconsistency! Mapping entry has no match in reference [".$giFrag."]\n";
                    }                    

                    # ///////////////////////////////////////////
                    # // Doesn't this obviate local alignments?!
                    # ///////////////////////////////////////////
                    ## Omit 3'-overhanging reads
                    #next LINE if(($rStop+1) > $uniqueGIlengths->{$gi}->{FRAG}->{$giFrag});
                
                    # Index
                    $gi2mapFrags{$gi}->{FRAG}->{$giFrag}->{$rStart}->{$rStop}++;    
                    $gi2mapFrags{$gi}->{HITCOUNT}++;
                    $gi2mapFrags{$gi}->{TOTALBPMAPPED} += $lengthMatch;

                } #FIELDS[2]
                else {
                    print "  **Unrecognized mapping entry!**  [".$fields[2]."]\n";
                }

            } #FIELDS[5]

            ####################
            # EDIT: 2013-07-17
            ####################
            #else {
            #    print "  **Unrecognized SAM field!** [".$fields[5]."]\n";
            #}
 
        } #LINE
        close $SAMFILE unless $SAMFILE_is_STDIN;
        my $iter1 = new Benchmark;
        print "done. Mapped split-trimmed reads: $MAPPED_TRIMMEDREADS_CNT, Unmapped split-trimmed reads: $UNMAPPED_TRIMMEDREADS_CNT. ".timestr(timediff($iter1,$iter0))."\n";         # parsing SAM file
    } #@samFiles

    close $MAPPED_TRIMMEDREADS unless ($noMappedFastq);
    close $UNMAPPED_TRIMMEDREADS unless ($noUnmappedFastq);
    
    return \%gi2mapFrags;
}
################################################################################
sub parseBLAST {
    die "sorry... not yet implemented!\n";
}
################################################################################
#-------------------------------------------------------------------------------
sub parseMultiFASTA {
    my $filename   = shift;        # The passed filehandle
    my %filehash   = ();              # Hash of the non-emptyline headers from file
    my %dup_header = ();              # Hash holding seqs with duplicate headers
    my %unique     = ();              # Hash holding seqs with unique headers
    my $header     = q{};             # Holds the header extracted from the regex
    my $header_already_found = 0;     # Flag; TRUE when the first header is found
    my $dup_header_found     = 0;     # Flad; Cycles 'TRUE' when dup header found
    #share(%unique);
    #share($header);

    my $count = 0;
    *STDOUT->autoflush(1);

    open my $INFILE, '<', $filename || die "Cannot open file \'$filename\'!\n";

    while (my $line = <$INFILE>) {
        chomp $line;
        if ($line =~ m/^[^(\s||\n)]*/) {          # Line not empty

            # Does the line contain a FASTA header?
#            if ($line =~ m/^>(\S+).*$/) {          # Line contains a FASTA header
#            if ($line =~ m/^>(.+)\s.+$/) {         # Line contains a FASTA header
            if ($line =~ m/^>(\S+)\s+.+$/) {        # Line contains a FASTA header
                $header = $1;                     # ...store it
#               	$header =~ s/\|/__/g;               # Need to remove the "|" character
#                                                  # since it causes file access problems
#                print "HEADER = $header"; <STDIN>;                                                  

                # Has the header been encountered already?
                if (exists $filehash{$header}) {  # Header is a duplicate
                    $dup_header{$header} = q{};   # Save dup header;
                                                  # Q:should q{} be used?
                    $dup_header_found = 1;        # ... remember that dup was found
                } # EXISTS
                else {
                    $filehash{$header}    = q{};  # Create hash entry
                    $header_already_found = 1;    # ... remember that header found
                    $unique{$header}      = q{};  # ... save unique
                    $dup_header_found     = 0;    # ... reset dup found to FALSE
                } # ELSE

            } # IF
            elsif ($header_already_found) {         # Not a header, but one already found
                if ($dup_header_found) {
                    $line =~ m/^(\s*)(.*)$/;        # Capture the sequence info in $2
                    $dup_header{$header} .= $2; # ...need whitespace
                }
                else {
                    $line =~ m/^(\s*)(.*)$/;        # Capture the sequence info in $2
                    $filehash{$header} .= $2;   # Append seq to hash; whitespace
                    $unique{$header}   .= $2;   #   needed if module re-used for
                                                    #   importing QUALITY values.
                }
            } # ELSIF
            else {
#                print "RUBBISH LINE\n";
            }
        } # IF
        else {
#            print "EMPTY LINE\n";
        }
    if ($count % $lineBuffer == 0) {
        print "." if($verbose);
    }
    $count++;
    } # WHILE

    #print "\n" if($verbose);
    close ($INFILE);

    return \%unique;
}
################################################################################
#         $_[0]        $_[1]      $_[2]
# ARGS: \%hitTree, $taxRankAbbr, $taxRank
#-------------------------------------------------------------------------------
sub writeCoords {

    #my $writeIter0 = new Benchmark;
    my $outFilename = $outdir."/".$prefix.".".$_[2].".contig.coords.tsv";
    #print "->Storing ".$_[2]."-level coordinates to disk []...";

    my %giCoords = %{ retrieve $outdir."/".$prefix.".giCoords.dmp" };
    open my $COORDSFILE, '>', $outFilename;
    
    # $hitTree->{SS}->{$strainName}->{CHILDREN}->{NAMES}
    foreach my $rankName (keys %{ $_[0]->{$_[1]} }) {
        foreach my $gi (keys %{ $_[0]->{$_[1]}->{$rankName}->{CHILDREN}->{NAMES} }) {
            foreach my $giFrag (sort {$a cmp $b} keys %{ $giCoords{$gi} }) {
                print $COORDSFILE $rankName."\t".$giFrag;
                print $COORDSFILE ",".$_.":".$giCoords{$gi}->{$giFrag}->{$_} foreach (sort {$a <=> $b} keys %{ $giCoords{$gi}->{$giFrag} });
                print $COORDSFILE "\n";
            }
        }
    }

    close $COORDSFILE;
    #my $writeIter1 = new Benchmark;

    #print "done. ".(timestr(timediff($writeIter1, $writeIter0)))."\n";

}
################################################################################
#        $_[0]    $_[1]        $_[2]     $_[3]
# ARGS: \%hash, $filename, $description, $mode      $mode = "bin" or "ascii"
################################################################################
sub dump2disk {

    STDOUT->autoflush(1);

    my $time0 = new Benchmark;

    if($_[3] =~ m/^bin/i) {
        print "->Storing datastructure ".$_[2]." to disk in BINARY format as \"".$_[1]."\"...";
        #stat_log_("Storing datastructure ".$_[2]." to disk (BINARY) as \"".$_[1]."\"...");
        #store clone($_[0]), $_[1];
        store $_[0], $_[1];
    }
    elsif($_[3] =~ m/^ascii/i) {
        print "->Storing datastructure ".$_[2]." to disk in ASCII (text) format as \"".$_[1]."\"...";
        #stat_log_("Storing datastructure ".$_[2]." to disk (ASCII) as \"".$_[1]."\"...");
        open my $OUTFILE, '>', $_[1];
        #print $OUTFILE Dump clone($_[0]);
        print $OUTFILE Dump($_[0]);
        close $OUTFILE;
    }
    else {
        print "->Attempting to store datastructure \"".$_[2]."\" to \"".$_[1]."\" in an unknown format: skipping...\n";
        #stat_log("**Warning**: Attempting to store datastructure \"".$_[2]."\" to \"".$_[1]."\" in an unknown format: skipping...");
        return;
    }

    my $time1 = new Benchmark;
    my $timeTxt = timestr(timediff($time1,$time0));
       $timeTxt = $1 if($timeTxt =~ m/^\s*(\d+)\s+/);
    
    #my $timeTxt = q{};
    #my $timediff = timestr(timediff($time1,$time0));
    #if($timediff =~ m/^\s*(\S+)\s+/) {
    #    $benchRez = 5;
    #    $timeTxt = sciExpand($1, $benchRez);
    #}
    #elsif ($timediff =~ m/(\d+\.?\d+)\s+/) {
    #    $timeTxt = $1;
    #}
    #my $loadTimeText = sciExpand($1, $benchRez) if(timestr(timediff($time1,$time0)) =~ m/^\s*(\S+)\s+/);
    #_stat_log("done in [$timeTxt] wallsecs");
    print "done. ".$timeTxt." wallsecs\n";
    return;
}
################################################################################
# Objective: Parse out the "# of Reads" and "# of Bases" resulting from a FASTQ
#            trimming step. Expects the files to have the following (tab-
#            delimited) format:
# ******************************************************************************
#CMD: splitrim --inFile=/source/CDD1.fastq --ascii=33 --minL=40 --minQ=20 --maxN=1 --prefix=CDD1 --outPath=/home/traceyf/trim/ --outEncoding=33 --verbose 
#
#                 	RAW	SPLIT-TRIMMED
#                 	===	=============
#      # of Reads:	137440887	103703878	(75.45 %)
#      # of Bases:	13881529486	8127273766	(58.55 %)
#Mean Read Length:	101.00	78.00	(77.23 %)
#           stdev:	0.00	21.78
#---------------------------------------------
# Total Time (ms):	633822
#
# ******************************************************************************
# ARGS: \@trimStatsFiles
################################################################################
sub parseTrimStatsFiles {

    foreach my $trimStatsFile (@{ $_[0] }) {
        print "Parsing $trimStatsFile...";
        open my $INFILE, '<', $trimStatsFile;
        my $foundReads = 0;
        my $foundBases = 0;
        while(my $line=<$INFILE>) {
            if($line =~ m/\s+\# of Reads:\t\d+\t(\d+)/) {
                $foundReads = $1;
            }
            elsif($line =~ m/\s+\# of Bases:\t\d+\t(\d+)/) {
                $foundBases = $1;
            }
        }
        close $INFILE;
        $inputReads += $foundReads;
        $inputBases += $foundBases;
        print " found $foundReads reads and $foundBases bases.\n";
    }

}
################################################################################
# parsedDTAR
################################################################################
sub checkInput {

    # Need required files
    usage() if($help);
    usage() if(!$make_dmp && !$mapFile);    # if not making a *.dmp file, better provide a map file!

    createDir($outdir);

    # Both given
    if($dbFile && $parsedDTAR) {
        if((-e $dbFile && -e $parsedDTAR) || (-e $parsedDTAR)) {
            $dbFile = q{};
        }
        elsif(-e $dbFile) {
            $parsedDTAR = q{};
        }
        else {                      # none exist
            print "\n**Neither \"$dbFile\" nor \"$parsedDTAR\" exist. Abort**\n";
            usage();
        }
    }
    # Only dbFile given
    elsif($dbFile) {

        # Only want to parse GOTTCHA/DTAR FNA file and store the *.dmp file
        if($make_dmp) {
            if(!-e $dbFile) {
                print "\n**--db=$dbFile does not exist. Abort**\n";
                usage();
            }
            calculateUniqueLengthsByGI($dbFile);
            exit();
        }
        else {
            print "\n**--db=$dbFile does not exist. Abort**\n" if(!-e $dbFile);
            usage() if(!-e $dbFile);
            my $emptyFile = (-s $dbFile == 0) ? (1) : (0);
            usage() if($emptyFile);
        }
    }
    # Only parsedDTAR given
    elsif($parsedDTAR) {
        print "\n**--parsedDB=\"$parsedDTAR\" does not exist. Abort**\n" if(!-e $parsedDTAR);
        usage() if(!-e $parsedDTAR);
    }
    # None given
    else {
        print "\n**Neither --db=<DB_file> nor --parsedDB=<parsed_DB> specified. Abort**\n";
        usage();
    }
   
    die "\n**MAP file \"$mapFile\" does not exist**" if(!-e $mapFile && $mapFile ne "-");


    # both input file variations specified and exist; parsed version takes priority
    #$dbFile = q{} if($dbFile && $parsedDTAR && -e $dbFile && -e $parsedDTAR);
    #die "DTAR file \"$dbFile\" does not exist: $!" if($dbFile && !-e $dbFile && !$parsedDTAR);
    #die "Parsed DTAR file \"$parsedDTAR\" does not exist: $!" if($parsedDTAR && !-e $parsedDTAR && !$dbFile);

    if(!$nucmer && !$blast && !$sam) {                                    # Need format type
        print "**Please specify map file type. Abort**\n";
        usage();
    }

    #print "    dbFile=$dbFile\n";
    #print "parsedDTAR=$parsedDTAR\n";
    #print "   mapFile=$mapFile"; <STDIN>;
        
    $minPctID = (($minPctID <=  0) || ($minPctID > 100)) ? $MINPCTID
              :                                            $minPctID;
    print "\n";

    $prefix = $prefix ? $prefix : "profiled";
    $detailFile = $outdir."/".$prefix.".replicon.reads.coords.csv";
    
    # -----------------------------------------------------------------
    # Check for optional "genomeVitals.dmp" file mapping GI to strain name
    # -----------------------------------------------------------------
    $noVitalsFile = (!$genomeVitalsFile || !-e $genomeVitalsFile) ? (1) : (0);
    print "->Warning: genomeVitals file either does not exist or was not provided.\n" if($noVitalsFile);

    # -----------------------------------------------------------------
    # Check for trimmed stats file(s): provides total "# of Reads" and "# of Bases"
    # -----------------------------------------------------------------
    if($trimStats) {
        my @missingTrimStats = ();
        my @tmpFiles = split(/,/, $trimStats);
        foreach(@tmpFiles) {
            if(!-e $_) {
                push(@missingTrimStats, $_);
            }
            else {
                push(@trimStatsFiles, $_);
            }
        }
        if(@missingTrimStats) {
            foreach(@missingTrimStats) {
                print "**WARNING**: The following trim stats files do not exist!\n";
                print "             ============================================\n";
                print "             $_\n" foreach(@missingTrimStats);
            }
            print "Proceeding without trimmed input data...\n";
        }
        else {
            parseTrimStatsFiles(\@trimStatsFiles);
        }
    }
    else {
        print "**WARNING**: No trim stas file provided. Proceeding without trimmed input data...\n"
    }

    ######################
    # Abundance data
    ######################
    $method ||= 1;
    #if(!$method) {
    #    print "Please specify a method!\n" if(!$method);
    #    usage();
    #}   
    my @methods = split(/,/, $method);
    foreach my $currMethod (@methods) {
        if(!exists $methods{$currMethod}) {
            print "Invalid method [$currMethod]!\n";
            usage();
        }
        else {
            push(@wantedMethods, $currMethod);
        }
    }

    if($minLength < 1) {
        print "--minLen must be >= 1\n";
        usage();
    }
    if($minHits < 1) {
        print "--minHits must be >= 1\n";
        usage();
    }
    if($minCoverage < 0) {
        print "--minCov must be >= 0\n";
        usage();
    }
    if($criticalCoverage < 0) {
        print "--cCov must be >= 0\n";
        usage();
    }
    if($minMLHL < 0) {
        print "--minMLHL must be >= 0\n";
        usage();
    }
    $abundanceOptions{TOPHIT}     = $findTopHit;        # bool
    $abundanceOptions{ABU_METHOD} = $method;            # Comma-delimited string
    $abundanceOptions{FIELD_SEP}  = $fieldSep;          # "\t"
    $abundanceOptions{MIN_COV}    = $minCoverage;       # float
    $abundanceOptions{MIN_HITS}   = $minHits;           # int
    $abundanceOptions{MIN_LEN}    = $minLength;         # ulong
    $abundanceOptions{MIN_MLHL}   = $minMLHL;           # float
    $abundanceOptions{CCOV}       = $criticalCoverage;  # float

    if($noFastqOut) {
        $noUnmappedFastq = 1;
        $noMappedFastq   = 1;
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
# Explain the following data:
# --------------------------
#    $entries
#    $n_o_length
#    $uniqueDBlength
#    $coverage
#    $linearDOC
#    $uniqueDOC
#    $hit_count
#    $mean_hit_length
#    $mean_linear_hit_length
#    $totalBpMapped
#
sub usage {
    my $HELP = <<HELP;

This program profiles (organism ID + abundance) an alignment file (SAM, Nucmer, BLAST, etc.) containing reads 
that were mapped to a database of taxonomic rank-specific nucleotide signatures.
    
Usage:  $0 [REQUIRED1] [REQUIRED2] [OPTIONS]

    [REQUIRED1]
    --db=<filename>                     Filename of the signature (GOTTCHA) DB in multiFASTA format.
    or  
    --db=<filename>                     Filename of the signature (GOTTCHA) DB in multiFASTA format.
    --make_dmp                          Parse a multiFASTA formatted signature (GOTTCHA) DB, store in Perl's Storable format, and exit.
    or  
    --parsedDB=<filename>               Filename for the signature (GOTTCHA) DB in Perl's Storable format.

    [REQUIRED2: unless --make_dmp]
    --map=<filename>                    Mapping results file in one of the specified formats below.
    --sam, --nucmer, or --blast         The three supported map file format types.
    --genomeVitals=<genomeVitals file>  Perl Storable hash of genome vitals.
    --treeFile=<TaxTree file>           Perl Storable hash of taxonomic tree (includes GI).
    --trimStats=<file1,file2,...>       Comma-separated list of text files containing the no. of reads and bases of the
                                        input FASTQ files after trimming.

    [ID FILTER: Coverage-based]         PASS if values are greater than these thresholds
    --minCov=<FLOAT>                    minimum best_LINEAR_COVERAGE [0.005]
    --minHits=<INT>                     minimum best_HIT_COUNT [5]
    --minLen=<INT>                      minimum best_LINEAR_LENGTH [100]

    [ABUNDANCE FILTER: Stacking-based]  PASS if values are greater than these thresholds
    --method=<INT>                      Comma-separated integers corresponding to the method of abundance calculation, 
                                        being either 1 or 2 or 3 or 1,2 or 1,3 or 2,3 or 1,2,3 [1]
    --minMLHL=<FLOAT>                   mean Linear Hit Length [10]
    --cCov=<FLOAT>                      critical coverage [0.006]
    --field=<STRING>                    Field separator within table [TAB, "\\t"]

    [OPTIONS]
    --outdir                            Output directory [.]
    --prefix                            Prefix to add to output files
    --noMappedFastq                     Prevents writing of the *.MAPPED.trimmed.fastq file
    --noUnmappedFastq                   Prevents writing of the *.UNMAPPED.trimmed.fastq file
    --noFastqOut                        Same as --noMapped and --noUnmapped
    --topHit                            Identifies the ONE most probable organism in sample [<DISABLED>]
    --minID=<min \%ID>                   Filter OUT hits in --map with \%identity lower than this [98%] <DISABLED>
    --verbose                           Verbosity
    --debug                             Prints extra data for debugging purposes
    --ddebug                            Specialized debugging
    --help                              Prints this message

    *Note that entries in --genomeVitals and --treeFile must match those in --db or --parsedDB exactly. 
     If inconsistencies among entries exist:
     (1) recreate --genomeVitals and --treeFile with mkSpeciesTree_v##.pl with the version of NCBI's names/nodes/gi_taxid_nucl.dmp files
         that correspond to those entries specified in --db and/or --parsedDB;
     (2) recreate --genomeVitals and --treeFile with mkSpeciesTree_v##.pl with NCBI's most recent names/nodes/gi_taxid_nucl.dmp files and
         manually (!!) update the entries in --db to be compliant (DANGEROUS!);
     (3) recreate --db from scratch using NCBI's most recent names/nodes/gi_taxid_nucl.dmp files.


Example:

    1. Making a parsed version of GOTTCHA only
       perl $0 --make_dmp --db=DTAR.Ecoli_v9b.r93.s15.l35.N1.u30.n64.p1.m9.ALL.fna --prefix=ECOLI --outdir=/home/me/outdir
       
    2. Parsing GOTTCHA's FNA file from scratch and profiling a SAM alignment file
       perl $0 --db=DTAR.Ecoli_v9b.r93.s15.l35.N1.u30.n64.p1.m9.ALL.fna --map=samfiles.list --sam 
               --genomeVitals=/home/me/indir/genomeVitals2.dmp          --treeFile=/home/me/indir/speciesTreeGI2.dmp
               --trimStats=/home/me/indir/datasetXYZ_fixL30.stats.txt   --rank=all
               --prefix=ECOLI    --outdir=/home/me/outdir

    3. Loading a pre-parsed version of GOTTCHA and profiling a SAM alignment file
       perl $0 --parsedDB=DTAR.Ecoli_v9b.r93.s15.l35.N1.u30.n64.p1.m9.ALL.fna.dmp --map=samfiles.list --sam 
               --genomeVitals=/home/me/indir/genomeVitals2.dmp          --treeFile=/home/me/indir/speciesTreeGI2.dmp
               --trimStats=/home/me/indir/datasetXYZ_fixL30.stats.txt   --rank=all
               --prefix=ECOLI    --outdir=/home/me/outdir

HELP
print "$HELP\n";
exit;
}

################################################################################
# ------------------------------ DATASTRUCTURES --------------------------------
# $gi2mapFrags->{ $gi => { FRAG => {  $giFrag1 => {  $rStart1 => {  $rStop1 => $count,
#                                                                   $rStop2 => $count,
#                                                                      ...
#                                                                },
#                                                 }
#                                  }
#                          COUNT => $numHITS,
#                        }
#               };
#
#
#-------------------------------------------------------------------------------
# Holds START|STOP positions all uniquely mapped "unique" fragments
# $gi2mergedFrags->{ $gi => {  FRAG => { $giFrag1 => {  $rStart1 => $rStop1,      
#                                                       $rStart2 => $rStop2,
#                                                          ...
#                                                    },
#                                      },
#                              COUNT => $numHITS,
#                           }
#                     ...
#                  };
#-------------------------------------------------------------------------------
#    $uniqueGIlengths->{ $gi }->{ FRAGS  => $numFrags,
#                                 LENGTH => $length   }
################################################################################

=for documentation

GI => S => SS => G => F => O => C => P


%hitTree = (
    GI  => {
             $gi1 => {
                       REFSIZE      => $giRefLength,
                       USIZE        => $giUniqueRefLength,
                       LINLEN       => $giLinearLength,
                       COV          => $coverage,
                       HITCOUNT     => $hitCount,
                       FULLREADHITS => $fullHitCount,
                       MAPPED       => $totalBPmapped,
                     },
             $gi2 => { },
              ...
           },
    SS  => {
             $rankName1 => {
                             REFSIZE      => $strainRefLength,
                             USIZE        => $strainUniqueRefLength,
                             LINLEN       => $strainLinearLength,
                             COV          => $coverage,
                             HITCOUNT     => $hitCount,
                             FULLREADHITS => $fullHitCount,
                             MAPPED       => $totalBPmapped,
                             CHILDREN     => {
                                               NAMES => { $gi1 => (), $gi2 => (), ... },      <--- names of all N-1 rank names
                                               BEST  => $best_gi                              <--- rank name of best N-1 rank
                                         }
                             },
                           },
             $rankName2 => { ... },
                 ...
           },
    S   => {
             $rankName1 => {
                             REFSIZE      => $strainRefLength,
                             USIZE        => $strainUniqueRefLength,
                             LINLEN       => $strainLinearLength,
                             COV          => $coverage,
                             HITCOUNT     => $hitCount,
                             FULLREADHITS => $fullHitCount,
                             MAPPED       => $totalBPmapped,
                             GPROJ        => $numGenomeProjectsUnderSpecies,
                             CHILDREN     => {
                                               NAMES => { $strainRankName1 => (), $strainRankName2 => (), ... },    <--- names of all N-1 rank names
                                               BEST  => $bestStrainRankName                                         <--- rank name of best N-1 rank
                                             }
                             },
                           },
             $rankName2 => { ... },
                 ...
           },
    G   => {
             $rankName1 => {
                             REFSIZE      => $genusRefLength,
                             USIZE        => $genusUniqueRefLength,
                             LINLEN       => $genusLinearLength,
                             COV          => $coverage,
                             HITCOUNT     => $hitCount,
                             FULLREADHITS => $fullHitCount,
                             MAPPED       => $totalBPmapped,
                             GPROJ        => $numGenomeProjectsUnderGenus,
                             CHILDREN     => {
                                               NAMES => { $speciesRankName1 => (), $speciesRankName2 => (), ... },  <--- names of all N-1 rank names
                                               BEST  => $bestSpeciesRankName                                        <--- rank name of best N-1 rank
                                             }
                             },
                           },
             $rankName2 => { ... },
                 ...
           },
    F   => { },
    O   => { },
    C   => { },
    P   => { }
);


=cut
