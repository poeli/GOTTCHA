#!/usr/bin/env perl
#
#
#   This is the UNSHARED version (of large variables) and is much FASTER but uses more RAM.
#
#
#

# 2014/10/28 Modified by Po-E (Paul) Li
#  - Use one-level-lower taxonomy name as unassigned level name

use strict;
use warnings;
use threads;
use threads::shared;
use IO::Handle;
use File::Basename;
use Getopt::Long;
use Tie::IxHash;
use Storable;
use Time::HiRes;
use Benchmark;
use XML::Simple;
use YAML::XS;

# GLOBALS
my $PREFIX       = "genomeVitals";
my $lineBuffer   = 250;
my $lineBuffer2  = 1000000;

my %dir2org     = ();       # KEY: original dirname;    VAL: removal of "_" and "uid#"
my %org2dir     = ();       # 
my %matchedOrgs = ();       # = ($dirname => {$rank => $sciName});

my %taxAbbr:shared;
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
   @revTaxAbbrExt{@ranksExtAbbr} = @ranksExt;
my $abbrTaxLevel:shared;

# INPUT OPTIONS
my $orgfile      = q{};
my $namesTaxFile = q{};      # names.ncbi (creating taxTree from scratch)
my $nodesTaxFile = q{};      # nodes.ncbi (creating taxTree from scratch)
my $parentTraceFile = q{};   # parentTrace.dmp file
my $taxTreeFile  = q{};      # Perl Storable format of species-level Tax Tree
my $loadVitals   = q{};
my $genbankList  = q{};
my $genomesList  = q{};
my $gi2taxidFile = q{};
my $prefix       = q{};
my $numCPUs      = 1;

# OTHER GLOBAL VARIABLES
my %taxTree     = ();               # HREF to %taxTree
my @unmatched   = ();
my @orgNames    = ();
my $taxOptions  = q{};
my $verbose     = 0;
my $vverbose    = 0;
my $help        = 0;

# INPUT:
#   [required]
#   1. Either --orgfile=<FILENAME>      Filename containing the full paths to the subdirs of interest where
#                                       the subdir names themselves are the organism names
#      or     --dir=<DIRECTORY>         Directory whose subdirectories are all the subdirs of interest where
#                                       the subdir names themselves are the organism names
#   2. --names
#   3. --nodes
#
#   [optional]
#   1. --genomes        Filename w/*.gbk filenames listed one per line from ftp://ftp.ncbi.nih.gov/genomes/Bacteria
#   2. --genbank        Filename w/*.gbk filenames listed one per line from ftp://ftp.ncbi.nih.gov/genbank/genomes/Bacteria


GetOptions(
    # ===============================    
    # CREATE FROM SCRATCH:
    # ===============================    
    "names=s"    =>  \$namesTaxFile,        # names.dmp
    "nodes=s"    =>  \$nodesTaxFile,        # nodes.dmp
    "pTrace=s"   =>  \$parentTraceFile,     # parentTrace.dmp
    # -------------------------------    
    # or LOAD:
    # -------------------------------    
    "loadTree=s" =>  \$taxTreeFile,         # Load %taxTree from disk
    # ===============================    
    # CREATE FROM SCRATCH:
    # ===============================    
    "genbank=s"  =>  \$genbankList,         # Filename of the GENBANK *.gbk filenames, listed one per line
    "genomes=s"  =>  \$genomesList,         # Filename of the GENOME/GENBANK *.gbk filenames, listed one per line
    "gi2taxid=s" =>  \$gi2taxidFile,        # gi_taxid_nucl.dmp
    "prefix=s"   =>  \$prefix,              # prefix used for the output file
    # -------------------------------    
    # or LOAD:
    # -------------------------------    
    "vitals=s"   =>  \$loadVitals,          # load %vitals instead of recomputing all the time (for debugging)
    # ===============================    
    # GENERAL
    # ===============================
    "threads=i"  =>  \$numCPUs,             # no. of threads to use in the computation
    "verbose"    =>  \$verbose,
    "vverbose"   =>  \$vverbose,
    "help"       =>  \$help,
);

STDOUT->autoflush(1);
my $iterX = new Benchmark;

##### GENERATE SPECIES TREE #####
my $iter0 = new Benchmark;
verifyInputOpts();
doTree();                                   # Generate TaxTree from scratch or Load
my $iter1 = new Benchmark;
print "-> Generation of SPECIES tree complete! [".timestr(timediff($iter1,$iter0))."]\n";

##### ADD GIs TO SPECIES TREE #####
$iter0 =  new Benchmark;
addVitals();                                # Add vitals to TaxTree
$iter1 =  new Benchmark;
print "-> Addition of GENOMIC VITALS complete! [".timestr(timediff($iter1,$iter0))."]\n";
print "\nAll done! [".timestr(timediff($iter1,$iterX))."]\n";


################################################################################
# Enable LOADTREE so GI and other vitals can be added
sub doTree {

    # Load %taxTree?
    if($taxOptions->{LOADTREE}) {
        my $taxTree = loadFromDisk($taxTreeFile, "Taxonomic Tree");

        print "-> Performing deep cloning of shared variable \"TAXONOMIC TREE\"...";
        my $iter0 = new Benchmark;
        $taxOptions->{TAXTREE}  = shared_clone($taxTree);
        my $iter1 = new Benchmark;
        print "done. [".timestr(timediff($iter1,$iter0))."]\n";
    }
    else {
        generateTaxTree($taxOptions);
    }
    
    return;
}    
################################################################################
sub addVitals {

    # Load %vitals?
    if($taxOptions->{LOADVITALS}) {
        my $vitals = loadFromDisk($loadVitals, "Genome Vitals");

        print "-> Performing deep cloning of shared variable \"GENOME VITALS\"...";
        my $iter0 = new Benchmark;
        $taxOptions->{VITALS} = shared_clone($vitals);
        my $iter1 = new Benchmark;
        print "done. [".timestr(timediff($iter1,$iter0))."]\n";
    }
    elsif($taxOptions->{GBKLISTFILE} || $taxOptions->{GENLISTFILE}) {

        my $gbkGIvitals = q{};                           #<-----------------------------HERE! BUILD UP "SHARED" %VITALS
        my $genGIvitals = q{};
        my %vitals :shared = ();

        # Generate %vitals
        if($taxOptions->{GBKLISTFILE}) {
            my $genbankFiles = getGBKfiles($taxOptions->{GBKLISTFILE});       # Build up array of files
            $gbkGIvitals     = parseGBK($genbankFiles, "GBK"); # Parse each *.gbk file and Populate $gbkGIvitals
        }
        if($taxOptions->{GENLISTFILE}) {
            my $genomesFiles = getGBKfiles($taxOptions->{GENLISTFILE});        # Build up array of files
            $genGIvitals     = parseGBK($genomesFiles, "GEN"); # Parse each *.gbk file and Populate $genGIvitals
        }

        if($gbkGIvitals && $genGIvitals) {
            %vitals = %{ mergeVitals($gbkGIvitals, $genGIvitals) };
        }
        else {
            %vitals = ($gbkGIvitals) ? (%{ $gbkGIvitals })
                    : ($genGIvitals) ? (%{ $genGIvitals })
                    : ();
        }

        $taxOptions->{VITALS} = \%vitals;
        addTaxID($taxOptions, $gi2taxidFile) if($taxOptions->{GI2TAXIDFILE});  # Parse gi_taxid_nucl.dmp for TAXIDs <--- THREAD
		changeOrgToSciName($taxOptions) if($taxOptions->{GI2TAXIDFILE});  # replace GBK organism to scientific name
		customStrainNameForSpeciesOrg($taxOptions);
		buildVitalOrgRepl($taxOptions);
		dump2disk($taxOptions->{VITALS}, $prefix.".dmp",          "Vitals", "bin");
        dump2disk($taxOptions->{VITALS}, $prefix.".readable.dmp", "Vitals", "ascii");
    }
    
    # Add %vitals to $taxOptions->{SPECIESTREE} if gi2taxid file was specified
    if($taxOptions->{GI2TAXIDFILE}) {
        appendGI2taxTree($taxOptions->{VITALS}, $taxOptions->{SPECIESTREE});
        dump2disk($taxOptions->{SPECIESTREE}, "speciesTreeGI.dmp", "Tax Tree (addition of GIs)", "bin");
        dump2disk($taxOptions->{SPECIESTREE}, "speciesTreeGI.readable.dmp", "Tax Tree (addition of GIs)", "ascii");
    }
    else {
        print "No GI or TAXID data to add to TAX TREE. Skipping... ";
        return;
    }

    return;    
}
################################################################################
sub verifyInputOpts {

    usage() if($help);

    $taxOptions = {
        REFRANKS => \@ranks,            # tax ranks, from genus --> superkingdom
        LOADTREE => q{},                 # create tree from scratch
    };
    
    # Create %taxTree from scratch or load?
    if( checkFile($taxTreeFile) == 1) {   # returns -1 if filename not given; 1 if given and exists; die() otherwise
        $taxOptions->{LOADTREE} = 1;      # Specify to LOADTREE instead of create from scratch
    }
    elsif( (checkFile($namesTaxFile) == 1) && (checkFile($nodesTaxFile) == 1) ) {
        $taxOptions->{NAMEFILE}   = $namesTaxFile;    # name of the NCBI "names.dmp" taxonomy file
        $taxOptions->{NODEFILE}   = $nodesTaxFile;    # name of the NCBI "nodes.dmp" taxonomy file
        $taxOptions->{PTRACEFILE} = $parentTraceFile if( checkFile($parentTraceFile) == 1);
    }
    else {
        usage();                                    # Must either load taxTree or create from scratch
    }

    ######### UNTESTED below...

    # Add genome vitals to %taxTree?    
    if( checkFile($loadVitals) == 1 ) {
        $taxOptions->{LOADVITALS} = 1;
    }
    elsif( $genbankList || $genomesList ) {
        $taxOptions->{LOADVITALS} = q{};
        $taxOptions->{GBKLISTFILE} = $genbankList if(checkFile($genbankList) == 1);
        $taxOptions->{GENLISTFILE} = $genomesList if(checkFile($genomesList) == 1);
    }

    $taxOptions->{GI2TAXIDFILE} = q{};
    $taxOptions->{GI2TAXIDFILE} = $gi2taxidFile if(checkFile($gi2taxidFile) == 1);

    # prefix
    $prefix = $PREFIX unless $prefix;
    
    # threads
    $numCPUs = ($numCPUs < 1) ? (1) : ($numCPUs);

    $vverbose = 1 if($verbose);
        
    return;   
}
################################################################################
# ARGS: $taxOptions
#       populates globally shared %taxTree
#-------------------------------------------------------------------------------
sub generateTaxTree {
    
    STDOUT->autoflush(1);

    # ////////////////////////////////////
    #   Parse out nodes
    # ////////////////////////////////////
    my $nodes = parseTaxNodesNCBI($_[0]);           # $nodes = { $node => {PNODE => $parentNode, RANK  => $rank } }

    #print "-> Performing deep cloning of shared variable \"NODES\"...";
    #my $iter0 = new Benchmark;
    $_[0]->{NODEHASH} = $nodes;                     # Add HREF to %options
    #my $iter1 = new Benchmark;
    #print "done. [".timestr(timediff($iter1,$iter0))."]\n";

    # ////////////////////////////////////
    # Parse out names
    # ////////////////////////////////////
    my $names = parseTaxNamesNCBI($_[0]);           # $names = { $node => { $name1 => $type, $name2 => $type } }

    #print "-> Performing deep cloning of shared variable \"NAMES\"...";
    #$iter0 = new Benchmark;
    $_[0]->{NAMEHASH} = $names;                     # Add HREF to %options
    #$iter1 = new Benchmark;
    #print "done. [".timestr(timediff($iter1,$iter0))."]\n";

    # ////////////////////////////////////
    # Generate parent trace for each node
    # ////////////////////////////////////
    if($_[0]->{PTRACEFILE}) {
        my $pTrace = loadFromDisk($_[0]->{PTRACEFILE}, "Parent Trace");

        print "-> Performing deep cloning of shared variable \"PARENT TRACE\"...";
        my $iter0 = new Benchmark;
        $_[0]->{PTRACE} = shared_clone($pTrace);    # Add it to the \%options
        my $iter1 = new Benchmark;
        print "done. [".timestr(timediff($iter1,$iter0))."]\n";
    }
    else {
        my $parentTrace = generateParentTrace($_[0]);

        print "-> Performing deep cloning of shared variable \"PARENT TRACE\"...";
        my $iter0 = new Benchmark;
        $_[0]->{PTRACE} = shared_clone($parentTrace);
        my $iter1 = new Benchmark;
        print "done. [".timestr(timediff($iter1,$iter0))."]\n";

        dump2disk($parentTrace, "parentTrace.dmp", "Parent Trace", "bin");
        #dump2disk($parentTrace, "parentTrace.readable.dmp", "Parent Trace", "ascii");
    }    

    # ////////////////////////////////////
    #   DEPRECATED! DEPRECATED! DEPRECATED! 
    # Construct tree from %nodes and %names and add to $options
    #   DEPRECATED! DEPRECATED! DEPRECATED!     
    # ////////////////////////////////////
    #my $taxTree = buildTree($_[0]);
    #my $taxTree = retrieve "taxTree.dmp";
    #$_[0]->{TAXTREE} = $taxTree;                   # Add tree to %options

    # ////////////////////////////////////
    # Split tree into species-grouped nodes, ready for export
    # ////////////////////////////////////
    reorganizeTree($_[0]);

    #print "Exiting..."; exit;

    return;
}

################################################################################
# ARGS: \%options
################################################################################
sub reorganizeTree {

    STDOUT->autoflush(1);

    # -------------------------------  SPECIES  --------------------------------

    print "-> Looking for all SPECIES nodes...";
    my $iter0 = new Benchmark;

    my %speciesNodes = ();      # Tracks all nodes in %parentTrace that are SPECIES
    my @speciesNodes = ();
    foreach my $currNode (keys %{ $_[0]->{PTRACE} }) {
        #print join("|", (keys %{ $_[0]->{NAMEHASH}->{$currNode} }))."\n" 
        #    if($_[0]->{NODEHASH}->{$currNode}->{RANK} eq "species");
        push(@speciesNodes, $currNode) if($_[0]->{NODEHASH}->{$currNode}->{RANK} eq "species");
        
    }

    my $iter1 = new Benchmark;
    print "done. [".timestr(timediff($iter1,$iter0))."]\n";
    #my $numSpecies = scalar(@speciesNodes);
    #print "Total of $numSpecies\n";
    @speciesNodes{ @speciesNodes } = ();        # Create lookup for valid species nodes
    #dump2disk(\%speciesNodes,"speciesNodes.readable.dmp","SPECIES NODES","ascii");

    # -------------------------------  SUBSPECIES  -----------------------------

    # Look for each subspecies (and variants thereof) for each species by looking for
    # non-standard ranks in %parentTrace (i.e. *not* SK,P,C,O,F,G and *not* S) 
    # whose nodes contain a node in @speciesNodes.
    print "-> Looking for all SUBSPECIES nodes...";
    $iter0 = new Benchmark;
    my %screenedRanks = %taxAbbr;   
    my @toScreen = ("subkingdom","superphylum","subphylum","superclass","subclass",
                    "infraclass","forma","subtribe","tribe","varietas","superorder",
                    "suborder","infraorder","parvorder","superfamily","subfamily",
                    "supergenus","subgenus","species");
    @screenedRanks{ @toScreen } = ();

    # Loop through %parentTrace for potential subspecies nodes: 
    # %subspeciesNodes = { $speciesNode => { $subspeciesNode1 => (),
    #                                        $subspeciesNode2 => (), ... }
    #                    };
    # Tracks all nodes in %parentTrace that are SUBSPECIES, indexed by the species node
    # that it belongs to.
    my %subspeciesNodes = ();
    NODE: foreach my $currNode (keys %{ $_[0]->{PTRACE} }) {
        next NODE if(exists $screenedRanks{ $_[0]->{NODEHASH}->{$currNode}->{RANK} });
        
        # Any nodes making it here have the potential for being a subspecies node.
        my @ancestors = @{ $_[0]->{PTRACE}->{$currNode} };
        
        # Any $ancestor nodes to the $currNode that are found in %speciesNodes are
        # actual subspecies nodes.
        foreach my $ancestor (@ancestors) {
            if(exists $speciesNodes{$ancestor}) {
                $subspeciesNodes{$ancestor}->{$currNode} = ();
                next NODE;
            }
        } #ancestor
    } #NODE
    $iter1 = new Benchmark;
    print "done. [".timestr(timediff($iter1,$iter0))."]\n";
    #my $numSpeciesWithSubspecies = scalar(keys %subspeciesNodes);
    #print "Total of $numSpeciesWithSubspecies species that have subspecies.\n";
    #dump2disk(\%subspeciesNodes,"subspeciesNodes.readable.dmp","SUBSPECIES NODES","ascii");


    # ---------------------- BUILD SPECIES TREE --------------------------------
    
    # Now build up the speciesTree and attach all subspecies to it and voila!
    # ...
    # Foreach speciesNode, use the parentTrace to get its history
    #   Verify that each ancestor in the parentTrace is in %taxAbbr, else do not add to speciesTree
    #   Add speciesNode to speciesTree

	print "-> Reconstructing the SPECIES Taxonomic Tree...";
    #my %speciesTree = ();
    my %speciesTree :shared = ();
    $iter1 = new Benchmark;
    my $count = 0;
    foreach my $speciesNode (@speciesNodes) {

        # Keep list of ranks found in each species' parentTrace        
        my %validRanks = ();
        @validRanks{ (keys %taxAbbr) } = ();

        # Get the species' parent trace
        my @ancestors = @{ $_[0]->{PTRACE}->{$speciesNode} };
        # Loop through the ancestral nodes, saving any that have a valid rank
        foreach my $ancestor (@ancestors) {
            my $ancestorRank = $_[0]->{NODEHASH}->{$ancestor}->{RANK};
            if(exists $taxAbbr{$ancestorRank}) {            # Valid rank?

                # Returns sciName if found; q{} otherwise
                my $sciName = getSciName($_[0]->{NODEHASH}->{$ancestor}->{RANK}, 
                                         $_[0]->{NAMEHASH}->{$ancestor});
                die "Fatal Exception: no \"scientific name\" found!\n"
                   ."  NODE:  $ancestor\n"
                   ."  RANK:  ".$_[0]->{NODEHASH}->{$ancestor}->{RANK}
                   ."  NAMES: ".join(" | ", (keys %{ $_[0]->{NAMEHASH}->{$ancestor} }))."\n"
                   if(!$sciName);
                   
                # Add sciName for SK,P,C,O,F,G ranks
                $speciesTree{$speciesNode} = &share({}) unless (exists $speciesTree{$speciesNode});
                $speciesTree{$speciesNode}->{ $taxAbbr{$ancestorRank} } = $sciName;
                
                # Remove the discovered rank
                delete $validRanks{ $ancestorRank };  
            }
        } #ancestor
       
        # Any remaining ranks in %validRanks went undiscovered, and are labeled as "Unassigned"
        my @missingRanks = keys %validRanks;
        if(@missingRanks) {
			my %rank_oneLvlDown = (
				"S" => "SS",
				"G" => "S",
				"F" => "G",
				"O" => "F",
				"C" => "O",
				"P" => "C"
			);
			my %tax_fulfilling_order = (
				"species" => 1,
				"genus"   => 2,
				"family"  => 3,
				"order"   => 4,
				"class"   => 5,
				"phylum"  => 6,
				"superkingdom"  => 7
			);
            $speciesTree{$speciesNode} = &share({}) unless (exists $speciesTree{$speciesNode});

			# UPDATE: 2014-10-28
			# Appending one-level-lower taxonomy name to unassigned level name

			my $speSciName = getSciName("genus",$_[0]->{NAMEHASH}->{$speciesNode});

			foreach my $rank ( sort {$tax_fulfilling_order{$a}<=>$tax_fulfilling_order{$b}} @missingRanks ){
				last if $rank eq "superkingdom";
				my $usename = $rank eq "genus" ? $speSciName : $speciesTree{$speciesNode}->{ $rank_oneLvlDown{$taxAbbr{$rank}} };
            	$speciesTree{$speciesNode}->{ $taxAbbr{$rank} } = "Unassigned $rank - $usename";
			}
		} #missingRanks

        # Incrementally add the SPECIES names to SPECIES TREE
        $speciesTree{$speciesNode}->{ $taxAbbrExt{"species"} } = &share({}) 
            unless exists ($speciesTree{$speciesNode}->{ $taxAbbrExt{"species"} });
        $speciesTree{$speciesNode}->{ $taxAbbrExt{"species"} }->{$_}
            = $_[0]->{NAMEHASH}->{$speciesNode}->{$_} foreach (keys %{ $_[0]->{NAMEHASH}->{$speciesNode} });
                
        # Add the SUBSPECIES names, if any
        if(exists $subspeciesNodes{$speciesNode}) {

            # Grab subspecies node
            foreach my $subspPNode (keys %{ $subspeciesNodes{$speciesNode} }) {

                # Incrementally add the names associated with the subspecies node
                $speciesTree{$speciesNode}->{ $taxAbbrExt{"subspecies"} } = &share({}) 
                    unless (exists $speciesTree{$speciesNode}->{ $taxAbbrExt{"subspecies"} });
                $speciesTree{$speciesNode}->{ $taxAbbrExt{"subspecies"} }->{$_} 
                    = $_[0]->{NAMEHASH}->{$subspPNode}->{$_} foreach (keys %{ $_[0]->{NAMEHASH}->{$subspPNode} });
                    
                #$speciesTree{$speciesNode}->{ $taxAbbrExt{"subspecies"} } = $subspNames;

            } #subspNode

        } #exists
        
        $count++;
        print "." if($count % 30000 == 0);
    }
    $iter1 = new Benchmark;
    print "done. [".timestr(timediff($iter1,$iter0))."]\n";

    # %speciesTree is complete! Add it to %options
    $_[0]->{SPECIESTREE} = \%speciesTree;

    # Save these as-is ONLY IF no GIs are going to be added
    if(!$taxOptions->{GBKLISTFILE} && !$taxOptions->{GENLISTFILE}) {
        dump2disk($_[0]->{SPECIESTREE}, "speciesTree.dmp", "Species Taxonomic Tree", "bin");
        dump2disk($_[0]->{SPECIESTREE}, "speciesTree.readable.dmp", "Species Taxonomic Tree", "ascii");
    }
    
    return;
}
################################################################################
# ARGS: reads $taxOptions, \%names
#       updates original %names in memory
#-------------------------------------------------------------------------------
sub parseTaxNodesNCBI {

    STDOUT->autoflush(1);

    my %nodes :shared = ();

    # Parse "nodes.dmp" file
    print "-> Parsing NODES file \"".$_[0]->{NODEFILE}."\"...";
    my $iter0 = new Benchmark;
    open my $NODES, '<', $_[0]->{NODEFILE} 
        || die "Cannot open file \"".$_[0]->{NODEFILE}."\"!\n";

    NODE: while(<$NODES>) {
        chomp; next if /\A\#/;
        my ($node, $parent, $rank, @others) = split /\t\|\t/;
        
        $nodes{$node} = &share({}) unless (exists $nodes{$node});
        $nodes{$node}->{PNODE} = $parent;
        $nodes{$node}->{RANK}  = $rank;
    }
    close($NODES);
    
    my $iter1 = new Benchmark;
    print "done. [".timestr(timediff($iter1,$iter0))."]\n";
    
    return \%nodes;
}
################################################################################
# ARGS: $taxOptions, \%names
#       writes/updates original %names in memory
#-------------------------------------------------------------------------------
sub parseTaxNamesNCBI {

    STDOUT->autoflush(1);

    my %names :shared = ();

    # Parse "names.dmp" file
    print "-> Parsing NAMES file \"".$_[0]->{NAMEFILE}."\"...";
    my $iter0 = new Benchmark;
    open my $NAMES, '<', $_[0]->{NAMEFILE} 
        || die "Cannot open file \"".$_[0]->{NAMEFILE}."\"!\n";
    while(<$NAMES>) {
        chomp; next if /\A\#/;
        my ($node, $name, $nada, $class) = split /\t\|\t/;
        my ($class1, @others) = split("\\s+\\|", $class);

        $names{$node} = &share({}) unless (exists $names{$node});
        $names{$node}->{$name} = $class1;
    }
    close($NAMES);
    
    my $iter1 = new Benchmark;
    print "done. [".timestr(timediff($iter1,$iter0))."]\n";
    
    return \%names;
}    
################################################################################
# ARGS:    \%options, where
#          $options{NODEHASH} = { $node => { PNODE => $parentNode, RANK => $rank } }
# Returns: \%parentTrace = { $node => [$pNodeA, $pNodeB, ... $pNodeZ] }
#          where $pnodeA is the most recent ancestor, and $pNodeZ is the root
# Function: generate an array representation for each node's parents, with most 
#           recent ancestor listed first:
#               %parentTrace = { $node => [$parentA, $parentB, ..., $parentZ ] }
################################################################################
sub generateParentTrace {
    
    my %parentTrace = ();
    
    print "-> Generating Parent Trace...";
    my $iter0 = new Benchmark;
    NODE: foreach my $currNode (keys %{ $_[0]->{NODEHASH} }) {
        my $pnode = $_[0]->{NODEHASH}->{$currNode}->{PNODE};
        
        PARENT: while(exists $_[0]->{NODEHASH}->{$pnode}) {
            push(@{ $parentTrace{$currNode} }, $pnode);
            
            # Get pnode's PNODE
            if(exists $_[0]->{NODEHASH}->{$pnode}->{PNODE}) {
                last PARENT if($pnode == 1);
                $pnode = $_[0]->{NODEHASH}->{$pnode}->{PNODE};
            }
            else {
                last PARENT;
            }
        } #PARENT
        
    } #currNode
    my $iter1 = new Benchmark;
    print "done. [".timestr(timediff($iter1,$iter0))."]\n";

    return \%parentTrace;
}
################################################################################
# ARGS: $nodes->{$ancestor}->{RANK}, $names->{$ancestor}
# Returns the scientific name if one is found; otherwise, returns q{}.
################################################################################
sub getSciName {

    # filter out name type "scientific name" from standard ranks: SK,P,C,O,F,G
    my $sciName = q{};
    if(exists $taxAbbr{ $_[0] }) {
        my @names     = keys   %{ $_[1] };          # NAME
        my @nameTypes = values %{ $_[1] };          # TYPE ("scientific name", "misspelling", ...)
        IDX: for my $idx (0 .. (scalar(@nameTypes)-1)) {
            return ($names[$idx]) if($nameTypes[$idx] eq "scientific name");
        } #IDX
    } #exists
    return $sciName;
}            
################################################################################
# Get all *.gbk filenames from input file
# ARGS: $gbklist
################################################################################
sub getGBKfiles {

    my $iter0 = new Benchmark;

    my @gbkfiles   = ();
    my @missingGBK = ();

    print "-> Acquiring list of Genbank files from \"".$_[0]."\"...";
    # Track inconsistent *.gbk files (present in $gbklist, but not where it says it is);
    open my $GBKLIST, '<', $_[0] || die "Fatal: Cannot open ".$_[0]." for read. Abort.\n";
    LINE: while(my $line=<$GBKLIST>) {
        chomp $line;
        next LINE if($line eq "");
        if(-e $line) {
            push(@gbkfiles, $line);
        }
        else {
            push(@missingGBK, $line);
        }
        
    }
    close $GBKLIST;
    my $iter1 = new Benchmark;
    print "done. ".timestr(timediff($iter1,$iter0))."\n";

    # Notify of missing *.gbk files
    print "The following files in ".$_[0]." were not found on disk:\n" if(@missingGBK);
    print "    ".$_."\n" foreach (@missingGBK);
    
    return \@gbkfiles;
}
################################################################################
# Parse each *.gbk file for GI, NAME(Org), REPL, ACC, DATE, and SIZE
# ARGS: \%giVitals, $gbk_type ("GBK" for GENBANK, or "GEN" for GENOMES/GENBANK)
################################################################################
sub parseGBK {

    my $iter0 = new Benchmark;
    
    my %giVitals :shared = ();
    $giVitals{GI}   = &share({});
    $giVitals{ORG}  = &share({});
    $giVitals{REPL} = &share({});
    $giVitals{ACC}  = &share({});        
    print "-> Parsing Genbank files for vital data (each \'.\' = $lineBuffer records)";
  
    # Parse *.gbk files for required info
    #foreach my $gbk (@{ $_[0] }) {
    for my $idx (0..(scalar(@{ $_[0]})-1)) {
        my $gbk = $_[0]->[$idx];
        # 1. open file
		#
		die "Genbank file \"$gbk\" is empty.\n" unless -s $gbk;
        open my $GBK, '<', $gbk || die "Cannot open Genbank file \"$gbk\" for read. Abort.\n";
        my $size     = q{}; 
        my $nType    = q{};
        my $topology = q{};     # "circular" or "linear"
        my $date     = q{};
        my $acc      = q{};
        my $gi       = q{};
        my $org      = q{};
        my $repl     = q{};
        my $defFound = 0;
        my $stopDef  = 0;
        my $orgFound = 0;
        my $stopOrg  = 0;

        LINE: while(my $line=<$GBK>) {
            chomp $line;
            # 2. skip until ^LOCUS; parse
            # LOCUS       NC_006946               7494 bp    RNA     linear   VRL 21-APR-2009
            # LOCUS       NC_011646              34952 bp    DNA     circular VRL 10-DEC-2008
            #if($line =~ m/^LOCUS\s+\w+\s+(\d+)\sbp\s+DNA\s+(\w+)\s+\w+\s+(\S+)/) {

            # Removing "DNA" requirement because there are RNA viruses
            if($line =~ m/^LOCUS\s+\w+\s+(\d+)\sbp\s+(\S+)\s+(\w+)\s+\w+\s+(\S+)/) {
                $size     = $1;
                $nType    = $2;
                $topology = $3;
                $date     = $4;
                #print "SIZE = $size\n";
                #print "DATE = $date\n";
            }
            # 3. skip until ^DEFINITION; parse multiple lines until ^ACCESSION
            elsif($line =~ m/^DEFINITION\s+(.+)/) {
                $defFound = 1;
                $repl = $1;
                #print "DEF: $repl\n";
            }
            # 4. skip until ^VERSION; parse
            elsif($line =~ m/^VERSION\s+(\w+\.\d+)\s+GI\:(\d+)/) {
                $acc = $1;
                $gi  = $2;
                $stopDef = 1;
                #print "ACC = $acc\n";
                #print "GI  = $gi\n";
            }
            # 5. skip until ^SOURCE; parse
            elsif($line =~ m/^SOURCE\s+(.+)/) {
            #elsif($line =~ m/^\s\sORGANISM\s+(.+)/) {
                $orgFound = 1;
                $org = $1;
                #print "SOURCE = $org\n";
                #last LINE;  # Last piece of data found; terminate loop
            }
            # Continuation of DEFINITION
            elsif(!$stopDef && $defFound && ($line =~ m/^\s+(.+)/)) {
                #print "ORIG=\"$repl\"\n";
                $repl .= " ".$1;
                #print " NEW=\"$repl\"\n";
            }
            elsif(!$stopOrg && $orgFound && ($line =~ m/^\s+(.+)/)) {
                #print "\nORIG_ORG=$org\n";
                my $text = $1;
                if($text =~ m/ORGANISM/) {
                    $stopOrg = 1;
                    if($org =~ m/(.+)\s\(.+[^\(]\)$/) {
                        #my $oldOrg = $org;
                        $org = $1;
                        #print "OLD_ORG=$oldOrg\n";
                        #print "NEW_ORG=$org"; <STDIN>;
                    }
                }
                else {
                    $org .= " ".$text;
                }
                #print " NEW_ORG=$org\n";
                last LINE if($stopOrg);  # Last piece of data found; terminate loop
            }
        } #LINE
        close $GBK;

        #if($org =~ m/^Baumannia/) {
        #    print "\nORG=\"$org\"";
        #    <STDIN>;
        #}
            
        $repl =~ s/\.$//;                               # Remove tailing periods
        $giVitals{GI}->{$gi} = &share({}) unless (exists $giVitals{GI}->{$gi});
        $giVitals{GI}->{$gi}->{ORG}     = $org;         # NAME
        $giVitals{GI}->{$gi}->{REPL}    = $repl;        # Replicon Name
        $giVitals{GI}->{$gi}->{DATE}    = $date;        # GenBank published date
        $giVitals{GI}->{$gi}->{NTYPE}   = $nType;       # Nucleic Acid Type (DNA/RNA/etc.)
        $giVitals{GI}->{$gi}->{TOPO}    = $topology;    # "linear" or "circular"
        $giVitals{GI}->{$gi}->{ACC}     = $acc;         # Accession.Version
        $giVitals{GI}->{$gi}->{SIZE}    = $size;        # Length of genome
        $giVitals{GI}->{$gi}->{SOURCE}  = $_[1];        # GBK or GEN
        
#        # $_[1] = GEN or GBK
#        $giVitals{ORG}->{$org}                                   = &share({}) 
#            unless (exists $giVitals{ORG}->{$org});
#        $giVitals{ORG}->{$org}->{REPL}                           = &share({}) 
#            unless (exists $giVitals{ORG}->{$org}->{REPL});
#        $giVitals{ORG}->{$org}->{REPL}->{$repl}                  = &share({}) 
#            unless (exists $giVitals{ORG}->{$org}->{REPL}->{$repl});
#        $giVitals{ORG}->{$org}->{REPL}->{$repl}->{$_[1]}         = &share({}) 
#            unless (exists $giVitals{ORG}->{$org}->{REPL}->{$repl}->{$_[1]});
#        $giVitals{ORG}->{$org}->{REPL}->{$repl}->{$_[1]}->{$gi}  = &share({}) 
#            unless (exists $giVitals{ORG}->{$org}->{REPL}->{$repl}->{$_[1]}->{$gi});
#        $giVitals{ORG}->{$org}->{REPL}->{$repl}->{$_[1]}->{$gi}->{SIZE} = $size;
#        $giVitals{ORG}->{$org}->{REPL}->{$repl}->{$_[1]}->{$gi}->{DATE} = $date;
#        $giVitals{ORG}->{$org}->{REPL}->{$repl}->{$_[1]}->{$gi}->{NTYPE}= $nType;
#        $giVitals{ORG}->{$org}->{REPL}->{$repl}->{$_[1]}->{$gi}->{TOPO} = $topology;
#        $giVitals{ORG}->{$org}->{REPL}->{$repl}->{$_[1]}->{$gi}->{ACC}  = $acc;
#        #$giVitals{ORG}->{$org}->{REPL}->{$repl}->{$_[1]}->{GI}->{$gi} = ();
#
#
#        $giVitals{REPL}->{$repl}                        = &share({}) 
#            unless (exists $giVitals{REPL}->{$repl});
#        $giVitals{REPL}->{$repl}->{$_[1]}               = &share({}) 
#            unless (exists $giVitals{REPL}->{$repl}->{$_[1]});
#        $giVitals{REPL}->{$repl}->{$_[1]}->{$gi}        = &share({}) 
#            unless (exists $giVitals{REPL}->{$repl}->{$_[1]}->{$gi});
#        $giVitals{REPL}->{$repl}->{$_[1]}->{$gi}->{ORG}  = $org;
#        $giVitals{REPL}->{$repl}->{$_[1]}->{$gi}->{SIZE} = $size;
#        $giVitals{REPL}->{$repl}->{$_[1]}->{$gi}->{DATE} = $date;
#        $giVitals{REPL}->{$repl}->{$_[1]}->{$gi}->{NTYPE}= $nType;
#        $giVitals{REPL}->{$repl}->{$_[1]}->{$gi}->{TOPO} = $topology;
#        $giVitals{REPL}->{$repl}->{$_[1]}->{$gi}->{ACC}  = $acc;
#        #$giVitals{REPL}->{$repl}->{$_[1]}->{GI}->{$gi} = ();

        $giVitals{ACC}->{$acc} = $gi;
        print "." if($idx % $lineBuffer == 0);
    } #FOREACH

    my $iter1 = new Benchmark;
    print "done. ".timestr(timediff($iter1,$iter0))."\n";

    return \%giVitals;
}    
################################################################################
# ARGS: \%gbkGIvitals, \%genGIvitals
# KEYS TO MERGE:    "ACC", "GI", "ORG", "REPL"
################################################################################
sub mergeVitals {
    my %vitals :shared = ();
    
    # ACC: unique among both hashes, so just append
    $vitals{ACC} = mergeACC($_[0]->{ACC}, $_[1]->{ACC});

    # GI: unique among both hashes, so just append
    $vitals{GI} = mergeGI($_[0]->{GI}, $_[1]->{GI});

    # ORG: 
    $vitals{ORG} = mergeORG($_[0]->{ORG}, $_[1]->{ORG});
    
    # REPL:
    $vitals{REPL} = mergeREPL($_[0]->{REPL}, $_[1]->{REPL});
    
    return \%vitals;
}
################################################################################
# [0] = gbk
# [1] = gen
################################################################################
sub mergeACC {
    my %newACC :shared = (%{ $_[0] }, %{ $_[1] });
    return \%newACC;
}
################################################################################
# [0] = gbk
# [1] = gen
################################################################################
sub mergeGI {
    my %newGI :shared = (%{ $_[0]}, %{ $_[1] });
    return \%newGI;
}
################################################################################
# [0] = gbk
# [1] = gen
################################################################################
sub mergeORG {
    my %newOrg :shared = %{ $_[0] };            # Taken care of the GBK GI's; now for the GEN GI's

    # Append or add new to %newOrg from [1]
    my @orgs = keys %{ $_[1] };
    foreach my $org (@orgs) {
        if(exists $newOrg{$org}) {
            # append
            my @replicons = keys %{ $_[1]->{$org}->{REPL} };
            foreach my $replicon (@replicons) {
                $newOrg{$org}->{REPL}->{$replicon}        = &share({}) 
                    unless (exists $newOrg{$org}->{REPL}->{$replicon});
                $newOrg{$org}->{REPL}->{$replicon}->{GEN} = &share({}) 
                    unless (exists $newOrg{$org}->{REPL}->{$replicon}->{GEN});
                foreach my $gi (keys %{ $_[1]->{$org}->{REPL}->{$replicon}->{GEN} }) {
                    $newOrg{$org}->{REPL}->{$replicon}->{GEN}->{$gi} = 
                        $_[1]->{$org}->{REPL}->{$replicon}->{GEN}->{$gi};
                } #gi
            } #replicon
        }
        else {
            # add new
            $newOrg{$org} = &share({}) unless (exists $newOrg{$org});
            $newOrg{$org} = $_[1]->{$org};
        }
    } 
    
    return \%newOrg;
    
}
################################################################################
# [0] = gbk
# [1] = gen
################################################################################
sub mergeREPL {
    my %newRepl :shared = %{ $_[0] };
    
    my @replicons = keys %{ $_[1] };
    foreach my $replicon (@replicons) {
        $newRepl{$replicon}        = &share({}) unless (exists $newRepl{$replicon});
        $newRepl{$replicon}->{GEN} = &share({}) unless (exists $newRepl{$replicon}->{GEN});
        foreach my $gi (keys %{ $_[1]->{$replicon}->{GEN} }) {
            $newRepl{$replicon}->{GEN}->{$gi} = $_[1]->{$replicon}->{GEN}->{$gi};
        }
    }
    
    return \%newRepl;
   
}
################################################################################
# ARGS: \%taxOptions, $gi2taxidFile
################################################################################
sub addTaxID {

    my $splitPrefix = "split";
    $_[0]->{VITALS}         = &share({}) unless (exists $_[0]->{VITALS});
    $_[0]->{VITALS}->{GI}   = &share({}) unless (exists $_[0]->{VITALS}->{GI});
    $_[0]->{VITALS}->{TAXID}= &share({}) unless (exists $_[0]->{VITALS}->{TAXID});

    STDOUT->autoflush(1);
    
    print "-> Processing $gi2taxidFile...\n";
    my $iterX = new Benchmark;

#    print "     ...Skipping split & using existing files...\n";    
    # 1. COUNT LINES: 
    print "     ...Determining size of $gi2taxidFile...";
    my $totalLines = (`wc -l $gi2taxidFile 2>&1 | awk '{print \$1}'`);
    print "done.\n";

    # 2. LOADBALANCE: 
    my $linesPerCPU = int($totalLines/$numCPUs)+1;

    # 3. SPLIT: determine unique temp directory name to store split files
    print "     ...Creating temporary directory for split files...";
    my $tempdir    = q{};
    {
        my $counter    = 0;                 # Counter to ensure uniqueness of subdir
        my $tempPrefix = "tmp";  
        while ($counter < 999) {
            $counter = pad_zeroes(++$counter, 3);
            $tempdir = "$tempPrefix$counter";
            last if(!-d $tempdir);   # exit loop when $tempdir !found; or max iterations
            $tempdir = q{};          # flag that max iterations met and no $tempdir found
        }
        die "Fatal: addTaxID(). Could not create temporary directory. "
           ."Please remove all \'$tempPrefix\*\' subdirectories and ensure "
           ."current directory is writeable.\n" 
           if(!$tempdir);

        # create the subdirectory
        my $mkdirCmd= `mkdir $tempdir 2>&1`;
        chomp $mkdirCmd;
        die "\nFatal: Unable to create directory \"$tempdir\"!\n" unless $mkdirCmd eq "";
        print "done.\n";
    }

    print "     ...Splitting $gi2taxidFile into $numCPUs partitions [$tempdir/]...";
    # <------- need to put in unique tmpdir so no foreign files can be mixed in
    my $splitCmd    = (`split -d --lines=$linesPerCPU $gi2taxidFile $tempdir/$splitPrefix`); 
    my @splitFiles  = (`ls -x1 $tempdir/$splitPrefix* 2>&1`);
       chomp $_ foreach(@splitFiles);
#    #print "splitFiles: ".join(",",@splitFiles); <STDIN>;
    print "done.\n";

    # 4. ASSIGN EACH $FILE TO A THREAD: 
    my @threads = ();
    my $iter0 = new Benchmark;

    print "     ...Processing the partitions...";

    my @GIs = sort {$a <=> $b} keys %{ $_[0]->{VITALS}->{GI} };

    ############################################################################
    # gi2taxidFile contains the GIs in numerical order, SO... each SPLITFILE is in order, too
    # if the the GI reference sent to each thread contains ONLY those GIs that could possible hit,
    # then we could shorten the lookup time on each thread.
    # Just need to sort the @GIs and then create individual lookups that correspond to the
    # range of GIs in the splitfiles
    # 
    # 1. open each $splitfile and get the first $gi in each.
    # 2. sort the GI's in increasing order -- this is the range of the GIs
    # 3. partition the GIs in @GIs according to these ranges
    # 4. send to own thread

    my %gi2splitfile = ();      # KEY=first GI of splitfile     VAL=splitfileName
    my @ranges = ();
    foreach my $splitFile (@splitFiles) {
        open my $INFILE, '<', $splitFile;
        while(my $line = <$INFILE>) {
            my ($gi, $taxid) = split(/\t/, $line);
            next if($gi !~ m/\d+/);
            push(@ranges, $gi);
            $gi2splitfile{$gi} = $splitFile;
            #print "\nUsing GI=$gi as splitFile=$splitFile";
            last;
        }    
        close $INFILE;
    }
    #print "\n";
    # contains the sorted GIs (use %gi2splitfile to lookup the corresponding splitfileName)
    my @sorted_ranges = sort { $a <=> $b } @ranges;     

    my %partitionedGIs = ();    # KEY=$idx  VAL=\@GI_subset

    # Determine GI borders within the SPLITFILES
    foreach my $idx (0..(scalar(@sorted_ranges)-1)) {
        # lookup the corrresponding $splitfileName with $gi2splitfile{$splitStart}
        my $splitStart = $sorted_ranges[$idx];         
        #print "Idx=$idx\tSTART($splitStart)\t";
        my $splitStop  = 0;
        if(exists $sorted_ranges[$idx+1]) {
            $splitStop = $sorted_ranges[$idx+1];
        }
        else {
            $splitStop = -1;                     # flag; last boundary, so accept all remaining GIs
        }
        
        #print "STOP($splitStop)\n";
        
        if($splitStop == -1) {
            #print "At splitStop = -1\n";
            #print "$_\n" foreach (@GIs);
            #print "--------\n";
            push(@{ $partitionedGIs{$splitStart} }, @GIs);   # add reamining GIs
        }
        else {
            # Partition the GI's according to the SPLITFILE GI's
            GI: while (1) {
                my $currGI = shift (@GIs);  # Grab first GI in array
                
                # pull out GIs from @GIs until the $currGI is >= $splitStop
                if($currGI < $splitStop) {
                    push(@{ $partitionedGIs{$splitStart} }, $currGI);
                }
                else {
                    # reached end of boundary, put GI back and quit
                    unshift(@GIs, $currGI);         # put it back
                    #print "Found 313011998!\n";
                    last GI;                        # exit loop
                }
            } #GI
            
        } #else
    } #idx

    # Pair each splitfile with its corresponding GI Loookup hash and send to thread
    foreach my $splitStart (keys %partitionedGIs) {         # Splitstart is the first GI in the splitFile
        my %giLookup = ();
        @giLookup{ @{$partitionedGIs{$splitStart}} } = ();
        my $splitFile = $gi2splitfile{$splitStart};
        my $totalGIs = scalar(keys %giLookup);
        
        my $opts = { SPLITFILE => $splitFile,
                     GILOOKUP  => \%giLookup,
                     TOTALGIS  => $totalGIs,
                   };
        my $t = threads->new(\&addTaxID_THREADED, $opts);
        
        #print "Adding GI=$splitStart    and    SplitFile=$splitFile to thread!\n";
        
        push(@threads, $t);
    }

    
    foreach (@threads) {
    	#$_->detach;
    	my $vitals_href = $_->join();       # {GI}->{$gi}->{TAXID} = $taxid
    	                                    # {TAXID}->{$taxid}    = $gi
    	
        # The correct TAXID for the given strain is: $_[0]->{VITALS}->{GI}->{$gi}->{TAXID}
    	foreach my $gi (keys %{ $vitals_href->{GI} }) {
    	    #$vitals{GI}->{$gi}->{TAXID} = $vitals_href->{GI}->{$gi}->{TAXID};
    	    $_[0]->{VITALS}->{GI}->{$gi}          = &share({}) unless (exists $_[0]->{VITALS}->{GI}->{$gi});
    	    $_[0]->{VITALS}->{GI}->{$gi}->{TAXID} = $vitals_href->{GI}->{$gi}->{TAXID};
    	}

    	my $iter1 = new Benchmark;
    	my $timeText = timestr(timediff($iter1,$iter0));
        $timeText  = $1 if($timeText  =~ m/^\s*(\d+)\s+/);
        print "\n   thread ".$_->tid().": done in $timeText wallsecs";
    }
    my $iterY = new Benchmark;
    print "\ndone. ".timestr(timediff($iterY,$iterX))."\n";

    # TAXID will likely have many GIs associated with it, so attach the GI href
    print "Indexing GIs by TAXID...";
    $iterX = new Benchmark;
    foreach my $gi (keys %{ $_[0]->{VITALS}->{GI} }) {
        my $taxid = q{};
        $taxid = $_[0]->{VITALS}->{GI}->{$gi}->{TAXID} if(exists $_[0]->{VITALS}->{GI}->{$gi}->{TAXID});

        if(!$taxid) {
            print "\n    **Warning: GI \"$gi\" has no TAXID! Skipping this GI!**";
            delete $_[0]->{VITALS}->{GI}->{$gi};
        }
		else{
        	$_[0]->{VITALS}->{TAXID}->{$taxid}        = &share({}) unless (exists $_[0]->{VITALS}->{TAXID}->{$taxid});
        	$_[0]->{VITALS}->{TAXID}->{$taxid}->{$gi} = ();
		}
    } #gi
    $iterY = new Benchmark;
    print "\n    done. ".timestr(timediff($iterY,$iterX))."\n";

    return;
}
################################################################################
# ARGS: \%opts = (SPLITFILE => $splitfile,
#                 GIVITALS  => \%giVitals,
#                );
################################################################################
sub addTaxID_THREADED {

#        my $opts = { SPLITFILE => $splitFile,
#                     GILOOKUP  => \%giLookup,
#                     TOTALGIS  => $totalGIs,
#                   };
    
    my $originalGIs = $_[0]->{TOTALGIS};
    my $lineCount = 0;
    
    my %localVitals :shared = ();
    $localVitals{GI}    = &share({});
    $localVitals{TAXID} = &share({});

    # 5. PROCESS & COLLECT (in each thread):
    open my $INFILE, '<', $_[0]->{SPLITFILE} || die "Cannot open splitfile \"".$_[0]->{SPLITFILE}."\"\n";
    LINE: while(my $line=<$INFILE>) {
        chomp $line; 
        $lineCount++; 
        print "." if($lineCount % $lineBuffer2 == 0);
        my ($gi, $taxid) = split(/\t/, $line);
        #print "Split GI($gi) and TAXID($taxid)\n";
        next LINE if($gi !~ m/\d+/);
        next LINE unless (exists $_[0]->{GILOOKUP}->{$gi});

        #print "GI: $gi      TAXID: $taxid       ORG: ".$vitals{GI}->{$gi}->{ORG}."\n";

        $localVitals{GI}->{$gi}       = &share({}) unless (exists $localVitals{GI}->{$gi});
        $localVitals{TAXID}->{$taxid} = &share([]);

        $localVitals{GI}->{$gi}->{TAXID} = $taxid;
        push(@{ $localVitals{TAXID}->{$taxid} }, $gi);

        delete $_[0]->{GILOOKUP}->{$gi};
        $_[0]->{TOTALGIS}--;
        last LINE if($_[0]->{TOTALGIS} == 0);
    }
    close $INFILE;

    return \%localVitals;
}
################################################################################
# ARGS: \%vitals, \%taxTree
################################################################################
sub appendGI2taxTree {

    # foreach $vitals->{ORG}
    #   find {ORG} in $taxTree->{S} or $taxTree->{SS}
    #   get all {REPL} in %vitals: $vitals->{ORG}->{$org}->{REPL}
    #     determine if {REPL} is (p|P)lasmid or not
    #   add to %taxTree under GI=>{ PLASMID => { $gi1 => $replicon1,
    #                                            $gi2 => $replicon2, ... },
    #                               CHR     => { $gi1 => $replicon1,
    #                                            $gi2 => $replicon2, ... }, }
    
    my @unmatched = ();     # organisms in genomeVitals that couldn't not be
                            # placed in the taxTree
    
    # --------------------------------------------------
    # Build up a hash of {S} and {SS} from %taxTree, each indexing their $node in %taxTree
    # --------------------------------------------------
    #   %speciesHash = { $(sub)species1a => $node1,
    #                    $(sub)species1b => $node1,
    #                          ...       => $node1,
    #                    $(sub)species2a => $node2,
    #                    $(sub)species2b => $node2,
    #                          ...
    #              };
    print "-> Generating (sub)species lookup hash from $taxTreeFile...";
    my $iter0 = new Benchmark;
    my %speciesHash :shared = ();
    foreach my $node (keys %{ $_[1] }) {
        foreach my $species (keys %{ $_[1]->{$node}->{S} }) {
            $speciesHash{$species} = $node;
        } #species
        
        # Subspecies may not exist for every node
        if(exists $_[1]->{$node}->{SS}) {
            foreach my $subspecies (keys %{ $_[1]->{$node}->{SS} }) {
                $speciesHash{$subspecies} = $node;
            } #subspecies
        } #exists
    } #node
    my $iter1 = new Benchmark;
    print "done. ".timestr(timediff($iter1,$iter0))."\n";
    
    # --------------------------------------------------
    # Extract the GIs for each REPLICON of each ORGANISM
    # --------------------------------------------------
    print "-> Appending GIs to Taxonomic Tree...";
    $iter0 = new Benchmark;
    foreach my $org (keys %{ $_[0]->{ORG} }) {
        my %branch = ();
           $branch{PLASMID}        = &share({}) unless (exists $branch{PLASMID});
           $branch{PLASMID}->{GBK} = &share({}) unless (exists $branch{PLASMID}->{GBK});
           $branch{PLASMID}->{GEN} = &share({}) unless (exists $branch{PLASMID}->{GEN});
           $branch{CHR}            = &share({}) unless (exists $branch{CHR});
           $branch{CHR}->{GBK}     = &share({}) unless (exists $branch{CHR}->{GBK});
           $branch{CHR}->{GEN}     = &share({}) unless (exists $branch{CHR}->{GEN});
        my @replicons = keys %{ $_[0]->{ORG}->{$org}->{REPL} };
        foreach my $replicon (@replicons) {

            # Check for *GENBANK/GENOMES* GI: each replicon will have 1 or 2 entries,
            # indexed by {GEN} or {GBK}.
            if(exists $_[0]->{ORG}->{$org}->{REPL}->{$replicon}->{GEN}) {

                # For Reference:
                #        $giVitals{ORG}->{$org}->{REPL}->{$repl}->{$_[1]}->{$gi}->{SIZE} = $size;
                #        $giVitals{ORG}->{$org}->{REPL}->{$repl}->{$_[1]}->{$gi}->{DATE} = $date;
                #        $giVitals{ORG}->{$org}->{REPL}->{$repl}->{$_[1]}->{$gi}->{ACC}  = $acc;

                my @gis = keys %{ $_[0]->{ORG}->{$org}->{REPL}->{$replicon}->{GEN} };
                foreach (@gis) {
                    if($replicon =~ m/plasmid/i) {
                        $branch{PLASMID}->{GEN}->{$_} = $org;
                    }
                    else {
                        $branch{CHR}->{GEN}->{$_} = $org;
                    }
                } #gis
            }

            # Check for *GENOMES* GI: each replicon will have 1 or 2 entries,
            # indexed by {GEN} or {GBK}.
            if(exists $_[0]->{ORG}->{$org}->{REPL}->{$replicon}->{GBK}) {
                my @gis = keys %{ $_[0]->{ORG}->{$org}->{REPL}->{$replicon}->{GBK} };
                foreach (@gis) {
                    if($replicon =~ m/plasmid/i) {
                        $branch{PLASMID}->{GBK}->{$_} = $org;
                    }
                    else {
                        $branch{CHR}->{GBK}->{$_} = $org;
                    }
                } #gis
            }            
        } #replicon
        
        # Locate $org in %taxTree
        if(exists $speciesHash{$org}) {         # Locate ORG in %taxTree

            my $node = $speciesHash{$org};      # Get its node to place in %taxTree
            delete $speciesHash{$org};          # Reduce search space

            foreach my $sourceType (keys %branch) {                     # PLASMID or CHR
                foreach my $gbkType (keys %{ $branch{$sourceType} }) {  # GEN or GBK
                    my @gis = keys %{ $branch{$sourceType}->{$gbkType} };
                    
                    # Place the GIs in the %taxTree

                    # [Existing GIs are strings, not HREFs, so remove them or use a treeFile that does not have GIs yet...]
                    #$_[1]->{$node}->{GI} = &share({}) if(exists $_[1]->{$node}->{GI} and ref($_[1]->{$node}->{GI}) ne "HASH");
                    
                    $_[1]->{$node}->{GI}                                  = &share({}) 
                        unless (exists $_[1]->{$node}->{GI});
                    $_[1]->{$node}->{GI}->{$sourceType}                   = &share({}) 
                        unless (exists $_[1]->{$node}->{GI}->{$sourceType});
                    $_[1]->{$node}->{GI}->{$sourceType}->{$gbkType}       = &share({}) 
                        unless (exists $_[1]->{$node}->{GI}->{$sourceType}->{$gbkType});
                    $_[1]->{$node}->{GI}->{$sourceType}->{$gbkType}->{$_} = $branch{$sourceType}->{$gbkType}->{$_} 
                        foreach (@gis);
                    
                } #gbkType
            } #sourceType
        } #exists
        else {
            my $matchFound = 0;
            # Check if any of the REPL names match
            foreach my $replicon (@replicons) {
            
                # Remove anything past the parentheses
                if($replicon =~ m/^(.+\(.+\))/) {
                    $replicon = $1;
                }
                elsif($replicon =~ m/^(.[^,]+),(\S[^,]+)$/) {
                    $replicon = $1;
                }
            
                #if($replicon =~ m/Methanococcus maripaludis/) {
                #    print "Looking for \"$replicon\"..."; 
                #    <STDIN>;
                #}
                
                if(exists $speciesHash{$replicon}) {
                    $matchFound = 1;
                    my $node = $speciesHash{$replicon};
                    delete $speciesHash{$replicon};
                    foreach my $sourceType (keys %branch) {
                        foreach my $gbkType (keys %{ $branch{$sourceType} }) {
                            my @gis = keys %{ $branch{$sourceType}->{$gbkType} };
                            
                            # Place the GIs in the %taxTree
                            $_[1]->{$node}->{GI}                                  = &share({}) 
                                unless (exists $_[1]->{$node}->{GI});
                            $_[1]->{$node}->{GI}->{$sourceType}                   = &share({}) 
                                unless (exists $_[1]->{$node}->{GI}->{$sourceType});
                            $_[1]->{$node}->{GI}->{$sourceType}->{$gbkType}       = &share({}) 
                                unless (exists $_[1]->{$node}->{GI}->{$sourceType}->{$gbkType});
                            $_[1]->{$node}->{GI}->{$sourceType}->{$gbkType}->{$_} = $branch{$sourceType}->{$gbkType}->{$_} 
                                foreach (@gis);
                            
                        } #gbkType
                    } #sourceType
                } #exists
            } #replicon

            # Save all organism names that couldn't be matched to the speciesTree            
            push(@unmatched, $org) if(!$matchFound);
            
        } #else
        
    } #org
    $iter1 = new Benchmark;
    print "done. ".timestr(timediff($iter1,$iter0))."\n";

    if(@unmatched) {
        print "GIs for the following organism names (parsed from *.gbk files) could not be placed in $taxTreeFile:\n";
        print "---------------------------------------------------------------------------------------------------------\n"; 
        print "    $_\n" foreach (@unmatched);
    }
    
    return;
}
################################################################################
#
# !!!!!!!!!!!!! NOT IMPLEMENTED !!!!!!!!!!!!!!!!
#
# ARGS: reads $taxOptions
#       update $taxOptions
#
# Note: If this takes too long, can use the unix "split" command to split the file
#       into a set no. of pieces (< $numCPUs) and then thread the parser
#-------------------------------------------------------------------------------
sub addGI2taxTree {

    #my %gi2taxID = ();

    STDOUT->autoflush(1);

    # Collect current taxids and index them in lookup hash
    my %taxids = ();
    my @taxids = keys %{ $_[0]->{TAXTREE} };
    foreach (@taxids) {
        my $index = substr($_,0,1);
        $taxids{$index}->{$_} = ();
    }

    my $iter0 = new Benchmark;
    print "-> Parsing \"".$_[0]->{GIFILE}."\" for GIs and TAXIDs...";
    # Load $giTaxFile into temporary hash %gi2taxID: KEY=node (or taxid), VAL=GI
    open my $GIFILE, '<', $_[0]->{GIFILE} || die "Cannot open GI file \"".$_[0]->{GIFILE}."\"!\n";
    LINE: while(my $line=<$GIFILE>) {
        chomp $line;
        next LINE if($line eq "");
        # GI\tNODE
        my ($gi, $taxid) = split("\t",$line);
        my $index = substr($taxid,0,1);
        
        if(exists $taxids{$index}->{$taxid}) {      # Only retain entries whose taxids are in TAXTREE
            $_[0]->{TAXTREE}->{$taxid}->{GI}->{$gi} = ();
        }
    }
    close $GIFILE;
    my $iter1 = new Benchmark;
    print "done. [".timestr(timediff($iter1,$iter0))."]\n";

    return;    

}
################################################################################
# ARGS: $filename, $description
################################################################################
sub loadFromDisk {

    die "Attempted to load non-existent file \"".$_[0]."\" from disk. Abort\n" if(!-e $_[0]);

    print "-> Loading ".$_[1]." from disk [".$_[0]."]...";
    my $iter0 = new Benchmark;
    my $href = retrieve $_[0];
    my $iter1 = new Benchmark;
    print "done. [".timestr(timediff($iter1,$iter0))."]\n";
    
    return $href;
}
################################################################################
#        $_[0]    $_[1]        $_[2]     $_[3]
# ARGS: \%hash, $filename, $description, $mode      $mode = "bin" or "ascii"
################################################################################
sub dump2disk {

    STDOUT->autoflush(1);

    my $iter0 = new Benchmark;

    if($_[3] =~ m/^bin/i) {
        print "-> Storing ".$_[2]." to disk in BINARY format as \"".$_[1]."\"...";
        store clone($_[0]), $_[1];
    }
    elsif($_[3] =~ m/^ascii/i) {
        print "-> Storing ".$_[2]." to disk in ASCII (text) format as \"".$_[1]."\"...";
        open my $OUTFILE, '>', $_[1];
        print $OUTFILE Dump clone($_[0]);
        close $OUTFILE;
    }
    else {
        print "-> Attempting to store hash \"".$_[2]."\" to \"".$_[1]."\" in an unknown format: skipping...\n";
        return;
    }

    my $iter1 = new Benchmark;
    print "done. [".timestr(timediff($iter1,$iter0))."]\n";
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
sub loadBalanceStatic {
    my ($totalPoints, $numCPUs) = @_;
    
    my %dataPoints = ();       # KEY: dataPoints; VAL: NumNodes
                               # Spread dataPoints across $numCPUs NumNodes times

    my $n = int($totalPoints/$numCPUs); # no. of times $numCPUs are fully loaded 
    my $r = $totalPoints % $numCPUs;    # leftover datapoints

    my $numFullCycles = $n;                     # const, regardless of $r
    my $numPointsPerFullCycle = $n;             # const, regardless of $r
    my $numPointsPerUnbalancedNode = $n;
    my $numPointsPerBalancedNode = $n+1;
    my $numBalancedNodes = $r;
    my $numUnbalancedNodes = $numCPUs - $r;

    if($r == 0) {
        $numPointsPerBalancedNode   = $n;
        $numBalancedNodes           = $numCPUs; 
        $numUnbalancedNodes         = 0;
        $numPointsPerUnbalancedNode = 0;
    }
    $numUnbalancedNodes = 0 if($numPointsPerUnbalancedNode == 0);

    return ($numBalancedNodes,                  #   $r           or $numCPUs
            $numPointsPerBalancedNode,          #   $n+1         or $n
            $numUnbalancedNodes,                #   $numCPUs-$r  or  0
            $numPointsPerUnbalancedNode);       #   $n           or  0
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
sub changeOrgToSciName {
	foreach my $gi ( keys %{$_[0]->{VITALS}->{GI}} ) {
		my $taxid = $_[0]->{VITALS}->{GI}->{$gi}->{TAXID};
		if($taxid){
			my $sciname = getSciName( "genus", $_[0]->{NAMEHASH}->{$taxid}); # use "genus" to pass check
			$_[0]->{VITALS}->{GI}->{$gi}->{ORG} = $sciname if $sciname;
		}
		else{
			delete $_[0]->{VITALS}->{GI}->{$gi};
		}
	}
}

sub customStrainNameForSpeciesOrg {
	foreach my $gi ( keys %{$_[0]->{VITALS}->{GI}} ) {
		my $taxid = $_[0]->{VITALS}->{GI}->{$gi}->{TAXID};
		my $rank  = $_[0]->{NODEHASH}->{$taxid}->{RANK};
		if( $rank ne "no rank" && $rank ne "subspecies" ){
			my $sciname = $_[0]->{VITALS}->{GI}->{$gi}->{REPL};
			$sciname =~ s/,[^,]+$//;
			unless ( exists $_[0]->{SPECIESTREE}->{$taxid} ){
				print "WARNING: Skipped taxid: $taxid. Taxonomy speciestree not found.\n";
				next;
			}
			$_[0]->{SPECIESTREE}->{$taxid}->{SS} = &share({}) unless (exists $_[0]->{SPECIESTREE}->{$taxid}->{SS} );
			$_[0]->{SPECIESTREE}->{$taxid}->{SS}->{$sciname} = "scientific name - custom";
			$_[0]->{VITALS}->{GI}->{$gi}->{ORG} = $sciname;
		}
	}
}

sub buildVitalOrgRepl {
	foreach my $gi ( keys %{$_[0]->{VITALS}->{GI}} ) {
		my $org   = $_[0]->{VITALS}->{GI}->{$gi}->{ORG};
		my $src   = $_[0]->{VITALS}->{GI}->{$gi}->{SOURCE};
		my $repl  = $_[0]->{VITALS}->{GI}->{$gi}->{REPL};
		my $size  = $_[0]->{VITALS}->{GI}->{$gi}->{SIZE};
		my $date  = $_[0]->{VITALS}->{GI}->{$gi}->{DATE};
		my $nType = $_[0]->{VITALS}->{GI}->{$gi}->{NTYPE};
		my $topo  = $_[0]->{VITALS}->{GI}->{$gi}->{TOPO};
		my $acc   = $_[0]->{VITALS}->{GI}->{$gi}->{ACC};

		my $giVitals = $_[0]->{VITALS};

		$giVitals->{ORG}->{$org}                                         = &share({})
			unless (exists $giVitals->{ORG}->{$org});
		$giVitals->{ORG}->{$org}->{REPL}                                 = &share({})
			unless (exists $giVitals->{ORG}->{$org}->{REPL});
		$giVitals->{ORG}->{$org}->{REPL}->{$repl}                        = &share({})
			unless (exists $giVitals->{ORG}->{$org}->{REPL}->{$repl}); 
        $giVitals->{ORG}->{$org}->{REPL}->{$repl}->{$src}                = &share({})
			unless (exists $giVitals->{ORG}->{$org}->{REPL}->{$repl}->{$src});
        $giVitals->{ORG}->{$org}->{REPL}->{$repl}->{$src}->{$gi}         = &share({})
			unless (exists $giVitals->{ORG}->{$org}->{REPL}->{$repl}->{$src}->{$gi});
        $giVitals->{ORG}->{$org}->{REPL}->{$repl}->{$src}->{$gi}->{SIZE} = $size;
        $giVitals->{ORG}->{$org}->{REPL}->{$repl}->{$src}->{$gi}->{DATE} = $date;
        $giVitals->{ORG}->{$org}->{REPL}->{$repl}->{$src}->{$gi}->{NTYPE}= $nType;
        $giVitals->{ORG}->{$org}->{REPL}->{$repl}->{$src}->{$gi}->{TOPO} = $topo;
        $giVitals->{ORG}->{$org}->{REPL}->{$repl}->{$src}->{$gi}->{ACC}  = $acc;
        #$giVitals{ORG}->{$org}->{REPL}->{$repl}->{$_[1]}->{GI}->{$gi} = ();

		$giVitals->{REPL}->{$repl}                        = &share({})
			unless (exists $giVitals->{REPL}->{$repl});
        $giVitals->{REPL}->{$repl}->{$src}                = &share({})
			unless (exists $giVitals->{REPL}->{$repl}->{$src});
        $giVitals->{REPL}->{$repl}->{$src}->{$gi}         = &share({})
			unless (exists $giVitals->{REPL}->{$repl}->{$src}->{$gi});
        $giVitals->{REPL}->{$repl}->{$src}->{$gi}->{ORG}  = $org;
        $giVitals->{REPL}->{$repl}->{$src}->{$gi}->{SIZE} = $size;
        $giVitals->{REPL}->{$repl}->{$src}->{$gi}->{DATE} = $date;
        $giVitals->{REPL}->{$repl}->{$src}->{$gi}->{NTYPE}= $nType;
        $giVitals->{REPL}->{$repl}->{$src}->{$gi}->{TOPO} = $topo;
        $giVitals->{REPL}->{$repl}->{$src}->{$gi}->{ACC}  = $acc;
        #$giVitals{REPL}->{$repl}->{$_[1]}->{GI}->{$gi} = ();
	}
}


#-------------------------------------------------------------------------------
sub usage {
    print <<END;

This program (1) creates a Taxonomic Tree file based on NCBI's taxonomic dump files "names.dmp" and "nodes.dmp" and
(2) optionally parses a set of completed genome GenBank files and updates the Taxonomic Tree with the respective GIs.
Since two GIs exist for any completed genome (one for ~/genbank/genomes/Bacteria and one for ~/genomes/Bacteria),
both GIs are included and distinguished in the updated tree.

Usage: $0 [REQUIRED] [OPTIONS]

  [REQUIRED]
  
    CREATE DENOVO TAXTREE:
        --names     Name of NCBI's taxonomic dump file "names.dmp"
        --nodes     Name of NCBI's taxonomic dump file "nodes.dmp"
       [--pTrace]   Filename of pre-computed parent trace in Perl Storable format
    or
    LOAD PRE-EXISTING TAXTREE:
        --loadTree  Filename of pre-existing species-level tax tree, in Perl Storable format
                    Used for adding GIs to "speciesTree.dmp" after-the-fact. [**UNDER CONSTRUCTION**]

  [OPTIONAL]
    --genbank       Filename containing GenBank (*.gbk) files to process from NCBI's genbank/genomes, listed one per line
    --genomes       Filename containing GenBank (*.gbk) files to process from NCBI's genomes directory, listed one per line
    --gi2taxid      NCBI's dmp file mapping GI to TAXID [gi_taxid_nucl.dmp]
    --threads       No. of threads to process data with   
    --help          Print this message

Example:
1) Create TREE from scratch and add GIs parsed out of both ~/genomes/Bacteria and ~/genbank/genomes/Bacteria *.gbk files:

    $0 --names=/db/taxonomy/work/names.dmp --nodes=/db/taxonomy/working/nodes.dmp --genbank=genbankList.txt --genomes=genomesList.txt --gi2taxid=/db/taxonomy/working/gi_taxid_nucl.dmp --threads=8

2) Load pre-computed Tax Tree and add GIs parsed out of ~/genbank/genomes/Bacteria *.gbk files:

    $0 --loadTree=speciesTree.dmp --genbank=genbankList.txt --gi2taxid=/db/taxonomy/working/gi_taxid_nucl.dmp --threads=8

END
exit;
}

