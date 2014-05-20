#!/usr/bin/env perl
#
# gottchaFilterRollupAbundance4e.pl
#
use strict;
use warnings;
use Getopt::Long;
use Storable;
use YAML;
use Benchmark;
use File::Path;
#use feature 'switch';

my $strainLookupFile  = q{};
my $speciesLookupFile = q{};
my $genusLookupFile   = q{};
my $familyLookupFile  = q{};
my $orderLookupFile   = q{};
my $classLookupFile   = q{};
my $phylumLookupFile  = q{};
my $taxLookupFile     = q{};
my $indir             = ".";
my $prefix            = q{};    # prefix for the table files to be parsed
my $taxLevel          = q{};    # tax level that the GOTTCHA run was performed at
# --------------------------------
my $method            = q{};
#my @wantedMethods     = ();     # populated in sub checkInput()
my $fieldSep          = "\t";
my $minCoverage       = 0.005;
my $minHits           = 10;
my $minLength         = 100;
my $minMLHL           = 1;
my $criticalCoverage  = 0.006;  # The min. coverage req'd for MLHL filter to become active
#---------------------------------
my $abuField          = 'LINEAR_DOC'; # The field to calculate abundance
# --------------------------------
my $help              = 0;

my @taxAbbrs = ("SS", "S", "G", "F", "O", "C", "P");

my %taxAbbrLookup = ("SS" => "STRAIN",
                     "S"  => "SPECIES",
                     "G"  => "GENUS",
                     "F"  => "FAMILY",
                     "O"  => "ORDER",
                     "C"  => "CLASS",
                     "P"  => "PHYLUM",
                    );

my %taxAbbrRevLookup = ( "STRAIN"  => "SS",
                         "SPECIES" => "S",
                         "GENUS"   => "G",
                         "FAMILY"  => "F",
                         "ORDER"   => "O",
                         "CLASS"   => "C",
                         "PHYLUM"  => "P",
                       );

my %taxLevelConv = ("STRAIN"  => "species1",
                    "SPECIES" => "species0",
                    "GENUS"   => "genus0",
                    "FAMILY"  => "family0",
                    "ORDER"   => "order0",
                    "CLASS"   => "class0",
                    "PHYLUM"  => "phylum0"
                   );

# Load up variant Lookups: variantStrainLookup.dmp, variantSpeciesLookup.dmp, and taxLookupBySpecies.dmp
GetOptions(
    "strainLookup=s"    => \$strainLookupFile,
    "speciesLookup=s"   => \$speciesLookupFile,
    "genusLookup=s"     => \$genusLookupFile,
    "familyLookup=s"    => \$familyLookupFile,
    "orderLookup=s"     => \$orderLookupFile,
    "classLookup=s"     => \$classLookupFile,
    "phylumLookup=s"    => \$phylumLookupFile,
    "taxLookup=s"       => \$taxLookupFile,
    
    # If FILENAME="even-100genus1-100species-150Mpe-haserror_fixL30Q20_species0.strain.tsv"
    # then PREFIX="even-100genus1-100species-150Mpe-haserror_fixL30Q20"
    "dir=s"             => \$indir,                 # directory containing files
    "prefix=s"          => \$prefix,                # prefix for the table files to be parsed
    "taxLevel=s"        => \$taxLevel,              # tax level the GOTTCHA run was performed at
    # -------------------
    # FILTERING CRITERIA
    # -------------------
    "method=s"          => \$method,                # 1 or 2 or 3 or 1,2 or 1,3 or 2,3 or 1,2,3
    "field=s"           => \$fieldSep,              # field separator [default: tab]
    "minCov=f"          => \$minCoverage,           # FILTER1: minimum best_LINEAR_COVERAGE
    "minHits=i"         => \$minHits,               # FILTER2: minimum best_HIT_COUNT
    "minLen=i"          => \$minLength,             # FILTER3: minimum best_LINEAR_LENGTH
    "cCov=f"            => \$criticalCoverage,      # FILTER4a: critical coverage
    "minMLHL=f"         => \$minMLHL,               # FILTER4b: mean Linear Hit Length
    
    # -------------------
    # Lineage table
    # -------------------
    "abuField=s"        => \$abuField,              # field to represent abundance

    # -----------------------------------
    "help"              => \$help,                  # prints usage
);

usage() if($help);

my %inputOptions = ();
$inputOptions{STRAINFILE}    = $strainLookupFile;
$inputOptions{SPECIESFILE}   = $speciesLookupFile;
$inputOptions{GENUSFILE}     = $genusLookupFile;
$inputOptions{FAMILYFILE}    = $familyLookupFile;
$inputOptions{ORDERFILE}     = $orderLookupFile;
$inputOptions{CLASSFILE}     = $classLookupFile;
$inputOptions{PHYLUMFILE}    = $phylumLookupFile;
$inputOptions{TAXLOOKUPFILE} = $taxLookupFile;
$inputOptions{DIR}           = $indir;
$inputOptions{PREFIX}        = $prefix;
$inputOptions{TAXLEVEL}      = uc($taxLevel);
# -------------------------------------------------
$inputOptions{METHODS}       = $method;
$inputOptions{FIELDSEP}      = $fieldSep;
$inputOptions{MINCOV}        = $minCoverage;
$inputOptions{MINHITS}       = $minHits;
$inputOptions{MINLEN}        = $minLength;
$inputOptions{CCOV}          = $criticalCoverage;
$inputOptions{MINMLHL}       = $minMLHL;
# -------------------------------------------------
$inputOptions{ABUFIELD}      = $abuField;

# ------------------------------------------------------
# sub checkInput() adds the following to %inputOptions
# ------------------------------------------------------
# $inputOptions{TABLEFILE}    = $tableFilename;
# $inputOptions{STRAINLOOKUP} = retrieve $inputOptions{STRAINFILE};
# $inputOptions{SPECIESLOOKUP}= retrieve $inputOptions{SPECIESFILE};
# $inputOptions{TAXLOOKUP}    = retrieve $inputOptions{TAXLOOKUPFILE};
checkInput(\%inputOptions);

processTable(\%inputOptions);


################################################################################
# Given the STRAIN-level GOTTCHA, parse the *species1.strain.tsv file and roll up
# Given the SPECIES-level GOTTCHA, parse the *species0.strain.tsv file and roll up
# Given the GENUS-level GOTTCHA, parse the *genus0.strain.tsv file and roll up
# etc.
# Expected Format: ($PREFIX)_($TAXLEVEL).strain.tsv
#
# >>>> Parse the "strain.tsv" file and roll all the way up <<<
# ARGS: \%inputOptions
sub processTable {

    my %strainData = ();
    my %ancestry   = ();

    # 1. Load and parse the "strain.tsv" file
    parseTable($_[0], \%strainData);

    # 2. Filter with user-supplied criteria and roll-up each strain
    filterTable($_[0], \%strainData, \%ancestry);
    
    # 4. Write out each file as ($PREFIX)_($TAXLEVEL).strain.tsv.ABUX
    displayRollups(\%inputOptions, \%ancestry);

    # 5. Save to disk
    exportAbundances(\%inputOptions, \%ancestry, \%strainData);

    # 6. output lineage tsv 
    outputLineage(\%inputOptions, \%ancestry);
}
################################################################################
# ARGS: \%inputOptions, \%ancestry, \%strainData
sub exportAbundances {

#my %taxLevelConv = ("STRAIN"  => "species1",
#                    "SPECIES" => "species0",
#                    "GENUS"   => "genus0",
#                    "FAMILY"  => "family0",
#                    "ORDER"   => "order0",
#                    "CLASS"   => "class0",
#                    "PHYLUM"  => "phylum0"
#                   );

    # Name of output file is taken from: (i)  $_[0]->{PREFIX}
    #                                    (ii) $_[0]->{TAXLEVEL}
    my $currDir     = $_[0]->{DIR};
    my $prefix      = $_[0]->{PREFIX};
    my $taxLevel    = $taxLevelConv{ $_[0]->{TAXLEVEL} };
    my $taxLevelRes = lc( $_[0]->{TAXLEVEL} );
    my $tableFilename = $currDir."/".$prefix.".".$taxLevelRes.".tsv.ABUX";

    print "Exporting results to disk [".$tableFilename."]...";
    open my $OUTFILE, '>', $tableFilename;

    my $wantedRank = $taxAbbrRevLookup{ $_[0]->{TAXLEVEL} };
    #my @headers = @{ $_[0]->{HEADERS} };
    my @headers = ("LINEAR_LENGTH", "TOTAL_BP_MAPPED", "HIT_COUNT", "LINEAR_DOC", "NORM_COV", "NORM_COV_SCALED");
    
    # Headers
    print $OUTFILE $_[0]->{TAXLEVEL}."\t".join("\t", @headers)."\n";
    
    # Contents
    my %results = ();
    foreach my $org (sort keys %{ $_[1]->{$wantedRank} }) {
        my @headerValues = @{ $_[1]->{$wantedRank}->{$org} }{@headers};

        $results{$org}->{NORM_COV}        = $_[1]->{$wantedRank}->{$org}->{NORM_COV};
        $results{$org}->{NORM_COV_SCALED} = $_[1]->{$wantedRank}->{$org}->{NORM_COV_SCALED};
        
        print $OUTFILE $org."\t".join("\t", @headerValues)."\n";
    }
    close $OUTFILE;
    
    store \%results, $currDir."/".$prefix.".".$taxLevelRes.".tsv.ABUX.dmp";
    open my $READABLE, '>', $currDir."/".$prefix.".".$taxLevelRes.".tsv.ABUX.readable.dmp";
    print $READABLE Dump(\%results)."\n";
    close $READABLE;
    
    print "Done!\n";

}

################################################################################
# ARGS: \%inputOptions, \%ancestry
sub displayRollups {

    #print Dump($_[1])."\n"; exit;

#    my @taxAbbrs = ("SS", "S", "G", "F", "O", "C", "P");
#
#    my %taxAbbrLookup = ("SS" => "STRAIN",
#                         "S"  => "SPECIES",
#                         "G"  => "GENUS",
#                         "F"  => "FAMILY",
#                         "O"  => "ORDER",
#                         "C"  => "CLASS",
#                         "P"  => "PHYLUM",
#                        );

    my @headers = ("LINEAR_LENGTH", "TOTAL_BP_MAPPED", "HIT_COUNT", "LINEAR_DOC", "NORM_COV", "NORM_COV_SCALED");

    foreach my $taxAbbr (@taxAbbrs) {
        print "\n\n".$taxAbbrLookup{$taxAbbr}."\t".join("\t", @headers)."\n";
        foreach my $orgName (sort keys %{ $_[1]->{$taxAbbr} }) {
            my @headerValues = @{ $_[1]->{$taxAbbr}->{$orgName}}{@headers};
            print $orgName."\t".join("\t", @headerValues)."\n";
        }
        print "=================================================================\n";
    }
}
################################################################################
# ARGS: \%inputOptions, \%ancestry
sub outputLineage {

    my @headers = ("LINEAR_LENGTH", "TOTAL_BP_MAPPED", "HIT_COUNT", "LINEAR_DOC", "NORM_COV", "NORM_COV_SCALED");

    #my %taxLevel = ( "STRAIN"  => "90",
    #                 "SPECIES" => "80",
    #                 "GENUS"   => "70",
    #                 "FAMILY"  => "60",
    #                 "ORDER"   => "50",
    #                 "CLASS"   => "40",
    #                 "PHYLUM"  => "30"
    #               );

    # =================================
    # OUTPUT STRAIN DATA FOR KRONA PLOT
    # =================================
    my $abuField    = $_[0]->{ABUFIELD};
    my $currDir     = $_[0]->{DIR};
    my $prefix      = $_[0]->{PREFIX};
    my $taxLvl      = $_[0]->{TAXLEVEL};
    my $outfilename = $currDir."/".$prefix.".lineage.tsv";
    open my $OUTFILE, '>', $outfilename;
    my $taxAbbr = "SS";
    foreach my $orgName (sort keys %{ $_[1]->{$taxAbbr} }) {
        my @headerValues = @{ $_[1]->{$taxAbbr}->{$orgName}}{@headers};
		my ( $index ) = grep { $headers[$_] eq $abuField } 0..$#headers;
		my $kronaValue = $headerValues[$index];

        my $s = $_[0]->{STRAINLOOKUP}->{$orgName};
		$s = $_[0]->{SPECIESLOOKUP}->{$orgName} unless defined $s ;

        my @out;
		#my $k = $_[0]->{SUPERKINGDOM};
		push @out, $_[0]->{TAXLOOKUP}->{$s}->{P};
		push @out, $_[0]->{TAXLOOKUP}->{$s}->{C};
		push @out, $_[0]->{TAXLOOKUP}->{$s}->{O};
		push @out, $_[0]->{TAXLOOKUP}->{$s}->{F};
		push @out, $_[0]->{TAXLOOKUP}->{$s}->{G};
        push @out, $s;
        push @out, $orgName;

        my $text = join "\t", @out;
        print $OUTFILE "$kronaValue\t$text\n";
    }
    close $OUTFILE;
}
################################################################################
# ARGS: \%inputOptions, \%strainData, \%ancestry
sub filterTable {

    my $minLen  = $_[0]->{MINLEN};
    my $minHits = $_[0]->{MINHITS};
    my $minCov  = $_[0]->{MINCOV};
    my $cCov    = $_[0]->{CCOV};
    my $minMLHL = $_[0]->{MINMLHL};
    my %ancestry = ();

    ## Debug:
    ## ------
    #print "Starting orgs:\n";
    #print "    ".$_."\t".$_[1]->{$_}->{LINEAR_COV}."\t".$_[1]->{$_}->{HIT_COUNT}."\t".$_[1]->{$_}->{LINEAR_LENGTH}."\n" foreach (sort keys %{ $_[1] });
    #print "------------------------------------------------\n";
    #print "\n";
    #print "minCov  = $minCov\n";
    #print "minHits = $minHits\n";
    #print "minLen  = $minLen\n";

    # Filter out STRAINS
    foreach my $strain (keys %{ $_[1] }) {

        #my $test_org = "Hydrogenobacter thermophilus";
        #my $test_state = 0;
        #if($strain =~ m/$test_org/) {
        #    print "working on ... ".$test_org."\n";
        #    print "   LINEAR_COV = ".$_[1]->{$strain}->{LINEAR_COV}."\n";
        #    print "    HIT_COUNT = ".$_[1]->{$strain}->{HIT_COUNT}."\n";
        #    print "LINEAR_LENGTH = ".$_[1]->{$strain}->{LINEAR_LENGTH}."\n";
        #    $test_state = 1;
        #    <STDIN>;
        #}
    
        # Filter Stage1: Coverage Threshold
        delete $_[1]->{$strain} 
            if (
                ($_[1]->{$strain}->{LINEAR_COV}    < $minCov)
             || ($_[1]->{$strain}->{HIT_COUNT}     < $minHits)
             || ($_[1]->{$strain}->{LINEAR_LENGTH} < $minLen)
            );

        #if($test_state) {
        #    if(exists $_[1]->{$strain}) {
        #        print "strain is still here after STAGE #1 filter"; <STDIN>;
        #    }
        #    else {
        #        print $test_org." is GONE!\n"; <STDIN>;
        #    }
        #}

        # Filter Stage2: MLHL Threshold
        # A low MEAN_LINEAR_HIT_LENGTH with a low coverage suggests background or unknown bleed-through
        if(exists $_[1]->{$strain}) {
            delete $_[1]->{$strain}
                if (
                    ($_[1]->{$strain}->{MEAN_LINEAR_HIT_LENGTH} < $minMLHL)
                 && ($_[1]->{$strain}->{LINEAR_COV}             < $cCov)
                );
        }
    }
    # Remaining keys in %{ $_[1] } are valid strains

    ##Debug:
    ##------
    #print "Remaining orgs:\n";
    #print "    ".$_."\t".$_[1]->{$_}->{LINEAR_COV}."\t".$_[1]->{$_}->{HIT_COUNT}."\t".$_[1]->{$_}->{LINEAR_LENGTH}."\t".$_[1]->{$_}->{MEAN_LINEAR_HIT_LENGTH}."\n" foreach (sort keys %{ $_[1] });
    #exit;

    ##########################
    # STRAIN ROLL-UP
    ##########################
    # Calculate the LINEAR_DOC
    {
        my @strains = keys %{ $_[1] };
        foreach my $strain (@strains) {
            $ancestry{SS}->{$strain}->{TOTAL_BP_MAPPED} += $_[1]->{$strain}->{TOTAL_BP_MAPPED};
            $ancestry{SS}->{$strain}->{LINEAR_LENGTH}   += $_[1]->{$strain}->{LINEAR_LENGTH};
            $ancestry{SS}->{$strain}->{HIT_COUNT}       += $_[1]->{$strain}->{HIT_COUNT};
        }
        my $sum = 0;
        foreach my $strain (@strains) {
            my $linear_doc = $ancestry{SS}->{$strain}->{TOTAL_BP_MAPPED} / $ancestry{SS}->{$strain}->{LINEAR_LENGTH};
            $ancestry{SS}->{$strain}->{LINEAR_DOC} = $linear_doc;
            $sum += $linear_doc; 
        }
        my $max = 0;
        foreach my $strain (@strains) {
            my $norm_cov = $ancestry{SS}->{$strain}->{LINEAR_DOC} / $sum;
            $ancestry{SS}->{$strain}->{NORM_COV} = $norm_cov;
            
            $max = $norm_cov unless ($max > $norm_cov);
        }
        
        $ancestry{SS}->{$_}->{NORM_COV_SCALED} = $ancestry{SS}->{$_}->{NORM_COV} / $max foreach (@strains);
    }

    ############################################################################
    # Roll up each STRAIN to its SPECIES to get accurate LINEAR_DOC
    # NOTE: Not using the "best" values, since we need total values.
    #   $4  LINEAR_LENGTH
    #   $10 TOTAL_BP_MAPPED
    #   $11 LINEAR_DOC
    ############################################################################
    foreach my $org (keys %{ $_[1] }) {

        # =========================
        # Find SPECIES parent name
        # =========================
        my $species = q{};
        if(exists $_[0]->{STRAINLOOKUP}->{$org}) {
            $species = $_[0]->{STRAINLOOKUP}->{$org};
            #print "Placed STRAIN \"".$org."\" into SPECIES \"".$species."\" through STRAIN lookup.\n";
        }
        elsif(exists $_[0]->{SPECIESLOOKUP}->{$org}) {
            $species = $_[0]->{SPECIESLOOKUP}->{$org};
            #print "Placed STRAIN \"".$org."\" into SPECIES \"".$species."\" through SPECIES lookup.\n";
            
        }
        die "*FATAL*: Could not place STRAIN \"".$org."\" into parent SPECIES.\n" if(!$species);

        $ancestry{S}->{$species}->{TOTAL_BP_MAPPED} += $_[1]->{$org}->{TOTAL_BP_MAPPED};
        $ancestry{S}->{$species}->{LINEAR_LENGTH}   += $_[1]->{$org}->{LINEAR_LENGTH};
        $ancestry{S}->{$species}->{HIT_COUNT}       += $_[1]->{$org}->{HIT_COUNT};
        
        
    } #ORG

    ##########################
    # SPECIES ROLL-UP
    ##########################
    # Calculate the LINEAR_DOC
    {
        my $linear_doc_sum = 0;
        my @species = keys %{ $ancestry{S} };
        foreach my $species (@species) {
            my $linear_doc = $ancestry{S}->{$species}->{TOTAL_BP_MAPPED} / $ancestry{S}->{$species}->{LINEAR_LENGTH};
            $ancestry{S}->{$species}->{LINEAR_DOC} = $linear_doc;
            $linear_doc_sum += $linear_doc;
        }

        my $max_linear_doc_sum = 0;
        foreach my $species (@species) {
            my $norm_cov = $ancestry{S}->{$species}->{LINEAR_DOC} / $linear_doc_sum;
            $ancestry{S}->{$species}->{NORM_COV} = $norm_cov;

            # For scaling to largest {NORM_COV}
            $max_linear_doc_sum = $norm_cov unless ($max_linear_doc_sum > $norm_cov);
        }

        # Scale to largest {NORM_COV}
        $ancestry{S}->{$_}->{NORM_COV_SCALED} = $ancestry{S}->{$_}->{NORM_COV} / $max_linear_doc_sum foreach (@species);
    }

    #print Dump(\%ancestry)."\n"; <STDIN>;

    ##########################
    # HIGHER ROLL-UP
    ##########################
    
    
    {
        my @taxAbbr = ("S", "G", "F", "O", "C", "P");
        foreach my $taxAbbrIdx (1..5) {
            my $currTaxAbbr = $taxAbbr[$taxAbbrIdx];
            my $prevTaxAbbr = $taxAbbr[$taxAbbrIdx-1];
            
            my $lookup = ($currTaxAbbr eq "G") ? ($_[0]->{GENUSLOOKUP}) 
                       : ($currTaxAbbr eq "F") ? ($_[0]->{FAMILYLOOKUP})
                       : ($currTaxAbbr eq "O") ? ($_[0]->{ORDERLOOKUP})
                       : ($currTaxAbbr eq "C") ? ($_[0]->{CLASSLOOKUP})
                       : ($currTaxAbbr eq "P") ? ($_[0]->{PHYLUMLOOKUP})
                       : (q{});
            die "Could not access lookup table for \"".$currTaxAbbr."\"\n" if(!$lookup);

            # We roll up from the previous group of names            
            foreach my $prevRankName (keys %{ $ancestry{$prevTaxAbbr} }) {

                # ------------------
                # Get child values
                # ------------------
                if(exists $lookup->{$prevRankName}) {
                    my $currRankName = $lookup->{$prevRankName};
                    # E.g. $ancestry{G}->{ $_[0]->{TAXLOOKUP}->{$species}->{G} }->{TOTAL_BP_MAPPED} += $ancestry{S}->{$species}->{TOTAL_BP_MAPPED};
                    $ancestry{$currTaxAbbr}->{$currRankName}->{TOTAL_BP_MAPPED} += $ancestry{$prevTaxAbbr}->{$prevRankName}->{TOTAL_BP_MAPPED};
                    # E.g. $ancestry{G}->{ $_[0]->{TAXLOOKUP}->{$species}->{G} }->{LINEAR_LENGTH}   += $ancestry{S}->{$species}->{LINEAR_LENGTH};
                    $ancestry{$currTaxAbbr}->{$currRankName}->{LINEAR_LENGTH}   += $ancestry{$prevTaxAbbr}->{$prevRankName}->{LINEAR_LENGTH};
                    $ancestry{$currTaxAbbr}->{$currRankName}->{HIT_COUNT}   += $ancestry{$prevTaxAbbr}->{$prevRankName}->{HIT_COUNT};
                }
                else {
                    die "*FATAL*: Could not place $prevRankName [$prevTaxAbbr] into parent [$currTaxAbbr]\n";
                }
            } #prevRankName
            
            # ------------------
            # Calculate the roll-up values
            # ------------------
            my $sum = 0;
            my @names = keys %{ $ancestry{$currTaxAbbr} };
            foreach my $rankName (@names) {

                my $totalbpmapped = $ancestry{$currTaxAbbr}->{$rankName}->{TOTAL_BP_MAPPED};
                my $linLength     = $ancestry{$currTaxAbbr}->{$rankName}->{LINEAR_LENGTH};

                my $linear_doc = $ancestry{$currTaxAbbr}->{$rankName}->{TOTAL_BP_MAPPED} / $ancestry{$currTaxAbbr}->{$rankName}->{LINEAR_LENGTH};

                $ancestry{$currTaxAbbr}->{$rankName}->{LINEAR_DOC} = $linear_doc;
                $sum += $linear_doc;
            }
                
            my $max_sum = 0;
            foreach my $rankName (@names) {
                my $norm_cov = $ancestry{$currTaxAbbr}->{$rankName}->{LINEAR_DOC} / $sum;
                $ancestry{$currTaxAbbr}->{$rankName}->{NORM_COV} = $norm_cov;
                    
                # For scaling to largest {NORM_COV}
                $max_sum = $norm_cov unless ($max_sum > $norm_cov);
            }
                
            # Scale to largest {NORM_COV}
            $ancestry{$currTaxAbbr}->{$_}->{NORM_COV_SCALED} = $ancestry{$currTaxAbbr}->{$_}->{NORM_COV} / $max_sum foreach (@names);

        } #taxAbbrIdx
    }
    
    #print Dump(\%ancestry)."\n"; <STDIN>;
    
    %{ $_[2] } = %ancestry;
        
}
################################################################################
# ARGS: \%inputOptions, \%strainData
sub parseTable {    
    my $infilename = $_[0]->{TABLEFILE};
    my @headers = ();
    my $headerFound = 0;                    # Flag to ensure header was captured
    
    print "Parsing table \"".$infilename."\"...";
    open my $INFILE, '<', $infilename;
    HEADER: while(my $line=<$INFILE>) {
        chomp $line;
        next if($line eq "");
        my @fields = split(/\t/, $line);
        if($fields[0] =~ m/^STRAIN/) {
            @headers = @fields;             # Grab headers from STRAIN file
            $headerFound = 1;
            last HEADER;
        }
    }
    die "*FATAL*: Header could not be located!\n" if(!$headerFound);
    print "Done!\n";

    # -----------------------------------------------
    # $strainData{$org}->{$header1} = $header1value
    # -----------------------------------------------
    my $numHeaderFields = scalar(@headers);
    shift(@headers);

    $_[0]->{HEADERS} = \@headers;

    BODY: while(my $line=<$INFILE>) {
        chomp $line;
        next if($line eq "");
        my @fields = split(/\t/, $line);
        my $numDataFields = scalar(@fields);
        print "*WARNING*: No. of data fields [".$numDataFields."] does not match the no. of header fields [".$numHeaderFields."]!\n"
            if($numDataFields != $numHeaderFields);
        
        # Save data fields into hash
        # ... $strainData{} ...
        my $org = shift(@fields);           # remove the first element
        @{ $_[1]->{$org} }{@headers} = @fields;
    }
    close $INFILE;
}
################################################################################
# ARGS: \%inputOptions
sub checkInput {

    my %methods           = (1 => q{}, 2 => q{}, 3 => q{});

    # ---------------------------------------------------------------------
    # Scan for Lookup files    
    # ---------------------------------------------------------------------
    if(!$_[0]->{STRAINFILE}) {
        print "\nPlease specify strain lookup file.\n\n";
        usage();
    }
    if(!-e $_[0]->{STRAINFILE}) {
        print "\nStrain lookup file does not exist!\n\n";
        usage();
    }
    if(!-e $_[0]->{STRAINFILE}) {
        print "\nPlease specify species lookup file.\n\n";
        usage();
    }
    if(!-e $_[0]->{SPECIESFILE}) {
        print "\nSpecies lookup file does not exist!\n\n";
        usage();
    }
    if(!-e $_[0]->{GENUSFILE}) {
        print "\nGenus lookup file does not exist!\n\n";
        usage();
    }
    if(!-e $_[0]->{FAMILYFILE}) {
        print "\nFamily lookup file does not exist!\n\n";
        usage();
    }
    if(!-e $_[0]->{ORDERFILE}) {
        print "\nOrder lookup file does not exist!\n\n";
        usage();
    }
    if(!-e $_[0]->{CLASSFILE}) {
        print "\nClass lookup file does not exist!\n\n";
        usage();
    }
    if(!-e $_[0]->{PHYLUMFILE}) {
        print "\nPhylum lookup file does not exist!\n\n";
        usage();
    }
    if(!$_[0]->{TAXLOOKUPFILE}) {
        print "\nPlease specify tax lookup lookup file.\n\n";
        usage();
    }
    if(!-e $_[0]->{TAXLOOKUPFILE}) {
        print "\nTax lookup lookup file does not exist!\n\n";
        usage();
    }        

    # ---------------------------------------------------------------------
    # Check valid abundance calculation & filter options
    # ---------------------------------------------------------------------
    if(!$_[0]->{METHODS}) {
        print "Please specify a method!\n" if(!$_[0]->{METHODS});
        usage();
    }
    my @methods = split(/,/, $_[0]->{METHODS});
    foreach my $currMethod (@methods) {
        if(!exists $methods{$currMethod}) {
            print "Invalid method [$currMethod]!\n";
            usage();
        }
        else {
            push(@{ $_[0]->{WANTEDMETHODS} }, $currMethod);
        }
    }
    if($_[0]->{MINLEN} < 1) {
        print "--minLen must be >= 1\n";
        usage();
    }
    if($_[0]->{MINHITS} < 1) {
        print "--minHits must be >= 1\n";
        usage();
    }
    if($_[0]->{MINCOV} <= 0) {
        print "--minCov must be > 0\n";
        usage();
    }

    # ---------------------------------------------------------------------
    # Check tax level valid    
    # ---------------------------------------------------------------------
    my $taxLevel = $_[0]->{TAXLEVEL};
    if($taxLevel !~ m/(STRAIN|SPECIES|GENUS|FAMILY|ORDER|CLASS|PHYLUM)/) {
        print "\n--taxLevel must be one of: STRAIN, SPECIES, GENUS, FAMILY, ORDER, CLASS, or PHYLUM\n\n";
        usage();
    }
    $taxLevel = $taxLevelConv{$taxLevel};

    # ---------------------------------------------------------------------
    # Scan for Table files: looking for ($PREFIX)_($TAXLEVEL).strain.tsv
    # ---------------------------------------------------------------------
    my $currDir  = $_[0]->{DIR};
    my $prefix   = $_[0]->{PREFIX};
    my $tableFilename = $currDir."/".$prefix.".strain.tsv";
    if(!-e $tableFilename) {
        print "\nTable file \"".$tableFilename."\" does not exist!\n\n";
        usage();
    }
    $_[0]->{TABLEFILE} = $tableFilename;
    
    # ---------------------------------------------------------------------
    # Load Lookup files
    # ---------------------------------------------------------------------
    print "Loading STRAIN Lookup file...";
    $_[0]->{STRAINLOOKUP}  = retrieve $_[0]->{STRAINFILE};
    print "done!\n";
    print "Loading SPECIES Lookup file...";
    $_[0]->{SPECIESLOOKUP} = retrieve $_[0]->{SPECIESFILE};
    print "done!\n";
    print "Loading GENUS Lookup file...";
    $_[0]->{GENUSLOOKUP}   = retrieve $_[0]->{GENUSFILE};
    print "done!\n";
    print "Loading FAMILY Lookup file...";
    $_[0]->{FAMILYLOOKUP}  = retrieve $_[0]->{FAMILYFILE};
    print "done!\n";
    print "Loading ORDER Lookup file...";
    $_[0]->{ORDERLOOKUP}   = retrieve $_[0]->{ORDERFILE};
    print "done!\n";
    print "Loading CLASS Lookup file...";
    $_[0]->{CLASSLOOKUP}   = retrieve $_[0]->{CLASSFILE};
    print "done!\n";
    print "Loading PHYLUM Lookup file...";
    $_[0]->{PHYLUMLOOKUP}  = retrieve $_[0]->{PHYLUMFILE};
    print "done!\n";
    print "Loading TAX Lookup file...";
    $_[0]->{TAXLOOKUP}     = retrieve $_[0]->{TAXLOOKUPFILE};
    print "done!\n";
    
}
################################################################################
sub usage {
    my $HELP = <<HELP;

This program will parse a GOTTCHA *.tsv table file and generate accurate organism relatibve abundances
at all taxonomic levels.

Usage:  $0 [REQUIRED1] [REQUIRED2] [OPTIONS]

    --strainLookup=<STRING>     STRAIN lookup filename  (E.g. /lscratch/db/custom/variantStrainLookup.dmp)
    --speciesLookup=<STRING>    SPECIES lookup filename (E.g. /lscratch/db/custom/variantSpeciesLookup.dmp)
    --genusLookup=<STRING>      GENUS lookup filename   (E.g. /lscratch/db/custom/genusLookupBySpecies.dmp)
    --familyLookup=<STRING>     FAMILY lookup filename  (E.g. /lscratch/db/custom/familyLookupByGenus.dmp)
    --orderLookup=<STRING>      ORDER lookup filename   (E.g. /lscratch/db/custom/orderLookupByFamily.dmp)
    --classLookup=<STRING>      CLASS lookup filename   (E.g. /lscratch/db/custom/classLookupByOrder.dmp)
    --phylumLookup=<STRING>     PHYLUM lookup filename  (E.g. /lscratch/db/custom/phylumLookupByClass.dmp)
    --taxLookup=<STRING>        TAX lookup filename     (E.g. /lscratch/db/custom/taxLookupBySpecies.dmp)

    --taxLevel=<STRING>         Tax level that the GOTTCHA run was performed at. One of either: STRAIN, SPECIES, GENUS, FAMILY, ORDER, CLASS, or PHYLUM
    --dir=<STRING>              Directory name containing the GOTTCHA *.tsv tables
    --prefix=<STRING>           Prefix of the GOTTCHA *.tsv tables.
                                  If FILENAME="even-100genus1-100species-150Mpe-haserror_fixL30Q20.strain.tsv"
                                  Then PREFIX="even-100genus1-100species-150Mpe-haserror_fixL30Q20"

    --method=<1,2,3>            Method by which to calculate relative abundance. Multiple methods are comma-separated (no spaces).
                                1: Linear Depth-of-Coverage
                                2: 
                                3: 

    --field                     Field separator [default: <tab>]
    --minLen                    Minimum unique length to be considered valid in abundance calculation [100]
    --minCov                    Minimum linear coverage to be considered valid in abundance calculation [0.005]
    --minHits                   Minimum no. of hits to be considered valid in abundance calculation [10]
    --minMLHL                   Minimum Mean-Linear-Hit-Length to be considered valid in abundance calculation [1]
    --cCov                      Critical coverage below which --minMLHL will cause an organism to fail []

    --help                      Prints this message

HELP
print "$HELP\n";
exit 1;
}

