#!/usr/bin/env perl
#
#
##################################################################################################################################
#
#   To Add
#   ======
#   -Subdirs of "/db/ftp.ncbi.nih.gov/genbank/genomes/Viruses" do not all have their *.fna files downloaded from NCBI's FTP site.
#
#   Actually, some viruses do NOT have a single, contiguous genomic sequence, but rather discrete FASTA sequences. 
#   Script needs to be updated to ACCEPT MULTIFASTA FILES ASSOCIATED WITH A SINGLE REPLICON OF A SINGLE ORGANISM.
#       -these FASTA sequences are not linear, so wrap-around l-mers should NOT be generated (MAKE THIS AN OPTION 
#        IN sub ripRepliconsWithStarts)
#       -both strands (like single FASTA sequences) will be l-mer extracted
#       -need to differentiate the following:
#           ...multiple contigs/genes from a single replicon in ONE file       (common)     <--- virus *.ffn file; SUPPORT
#           ...multiple contigs/genes from a single replicon in MULTIPLE files (not likely) <--- do NOT support
#       -Need hierarchy of file extension search:
#           ... specify "fna" > "ffn";
#           ... if an organism's subdir does NOT have *.fna files, then *.ffn files will be searched
#           ... **POTENTIAL PROBLEM** if chromosome is complete (*.fna) but plasmid is not (*.ffn), and both are
#               in the directory, as this will only read in the chromosome data and ignore the plasmid data;
#
#
##################################################################################################################################
#
#   README
#   ======
#
#   Two files must be created in a separate program (mkSpeciesTree.pl) prior to running this one.
#      -speciesTreeGI.dmp
#      -genomeVitals.dmp
#   
#   Though this program creates an XML of all the metadata surrounding the FASTA files, it does not store the FASTA sequence.
#   This could lead to a data integrity issue if the FASTA files from which the XML file is created changes within the time 
#   the XML file is created and the FASTA files are actually used in the generation of the unique sequences. To avoid this,
#   ALWAYS CREATE THE XML FILE JUST BEFORE GENERATING THE UNIQUE SEQUENCES.
#
#   Some viruses (like ss-RNA viruses) do not have a typical genomic FASTA file (*.fna), but one is needed for proper processing.
#   So unless directed otherwise, if the *.fna file is missing in an organism's directory, this program will extract the genomic
#   sequence from the *.gbk file and create a new *.fna file from it.
#
#   If the addition of a GI from any FASTA file attempts to clobber an already pre-existing GI, a descriptive WARNING will be 
#   thrown and the new GI will not be added.
#
#   If any genome projects have duplicate SPECIES or STRAIN names, then they will be de-replicated (i.e. made unique) by the
#   addition of the NCBI's uid number to their scientific name. Consequently, the corresponding XML file will have to be updated
#   with this new information, as well as the original data in the TAXTREE file, "speciesTreeGI.dmp," making it "speciesTreeGI2.dmp".
#   Within the new "speciesTreeGI2.dmp" file, the new SPECIES or STRAIN names were appended to the corresponding {S} or {SS} entry.
#   In addition, all GI's in the {GI} entry are likewise updated, since correspondence with downstream results using profileGOTTCHA.pl
#   will likely result in an unmatched organism had these GIs not been updated. THERE MAY BE AN ISSUE WITH genomeVitals.dmp, AND
#   THE SPECIES/STRAIN NAMES MAY NEED TO BE UPDATED, TOO, HOWEVER THIS HAS NOT YET BEEN EXAMINED.
#
#
#
#   <root>  = path to the subdirectories containing all vital sequence data 
#             (Eg. "/home/traceyf/tmp/ftp.ncbi.nih.gov/genbank/genomes/Bacteria/").
#   <dir>   = subdirectory under "root" which isolates all genomic data for a single organism
#   <name>  = the official NCBI taxonomic name assigned at the NCBI tax rank identified in <rank>.
#   <rank>  = taxonomic rank assigned to organism, which will be either species
#             (S) or subspecies (SS). Subspecies encompasses everything below 
#             the NCBI species level, whether or not an official rank was assigned.
#   <taxid> = the taxonomic ID (node) associated with the Species or SubSpecies
#             assigned.
#   <gi>    = the GI of the associated GenBank file. Since there are two versions
#             (the ~/genomes and the ~/genbank/genomes) for Bacteria, only one is
#             listed here (~/genbank/genomes).
#   <date>  = the sequence submission date parsed out of the ^LOCUS line from the
#             corresponding GenBank file.
#   <source>= either "PLASMID" or "CHR", depending on the FASTA file under <dir> 
#             that was first parsed for the GI number. Eventually, all FASTA file 
#             data both will be listed here.
#   <stype> = either "GEN" (for ~/genomes/Bacteria) or "GBK" (for ~/genbank/genomes/Bacteria),
#             depending on which GenBank file the <gi> parsed out of the FASTA file was found.
#   <acc>   = the ACCESSION number associated with the <gi>.
#
#
#
use strict;
use warnings;
use threads;
use threads::shared;
use IO::Handle;
use File::Basename;
use Getopt::Long;
#use Tie::IxHash;
use Storable;
use Time::HiRes;
use Benchmark;
use XML::Simple;
use YAML::XS;
use Data::Dumper;

#------------------------------------------------------------------------------------
# %filehash = ( DRAFTS    => {
#                            },
#               COMPLETED => {
#                              $kingdom1 => {  $rootdir1 => { $dirname1 => $basename,
#                                                             $dirname2 => $basename,
#                                                                 ...
#                                                           },
#                                              $rootdir2 => { ... },
#                                                 ...
#                                           },
#                              $kingdom2 => { ... },
#                                 ...
#                            }
#             );
#------------------------------------------------------------------------------------
my %filehash    = ();       # main struct for ID'ing draft/completed input files and their kingdoms
my %xmlTemplate = ();       # main struct for XML output
my %allGIs      = ();       # {gi}=>{KINGDOM}; {gi}=>{STATUS}. All GIs found in local organisms
my %org2dir     = ();       # 
my %orgName2TaxData = ();   # $sciName->{NODE} = $taxNode (in speciesTreeGI.dmp)
                            # $sciName->{RANK} = $taxRank (in speciesTreeGI.dmp)

my %taxAbbr:shared;
#tie %taxAbbr, "Tie::IxHash";
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

my %validKingdoms = ("VIRUSES"     => 1, 
                     "PROKARYOTES" => 1,
                     "EUKARYOTES"  => 1,
                     "OTHERS"      => 1
                    );
my %statusAbbr    = ("drafts"    => "I",    # (I)ncomplete
                     "completed" => "C",    # (C)ompleted
                    );
my %topoAbbr      = ("linear" => "L", "circular" => "C");

my $PREFIX      = "matched";
my $XML_OUTFILE = "matched.xml";
my %taxTree     = ();               # HREF to %taxTree
my @fastaexts   = ("fna","gbk");    # Input file search order: fna > gbk
my @unmatched   = ();
my @dirNames    = ();
my @validDirNames=();               # Org dirs that have FNA files in them
my $xmlExt      = "xml";
my $xml_outfile = q{};
my $lineBuffer  = 100;

# INPUT OPTIONS
my $draftsFile    = q{};
my $completedFile = q{};
my $taxTreeFile = q{};
my $vitalsFile  = q{};
my $prefix      = q{};
my $skipOrgNoGI   = 0;
my $continueNoHostFound = 0;
my $noDenovoFNA = 0;
my $help        = 0;
my $verbose     = 0;
my $vverbose    = 0;

GetOptions(
    "drafts=s"        =>  \$draftsFile,          # file listing the draft organisms, w/kingdom headers
    "completed=s"     =>  \$completedFile,       # file listing the completed organisms, w/kingdom headers
    "taxTree=s"       =>  \$taxTreeFile,         # denovoTaxTree.dmp
    "vitals=s"        =>  \$vitalsFile,          # genomeVitals.dmp
    "prefix=s"        =>  \$prefix,              # Prefix for output Tax Registry File in XML format
    "skipOrgNoGI"     =>  \$skipOrgNoGI,         # Flag; If 1, skips org if no GI could be parsed, else die; [0]
    "continueNoMatch" =>  \$continueNoHostFound, # Flag; If 1, continues even if no TAXID associated w/GI [0]
    "noDenovoFNA"     =>  \$noDenovoFNA,         # If *.fna files missing, do NOT create *.fna files from *.gbk files
    "help"            =>  \$help,
    "verbose"         =>  \$verbose,
    "vverbose"        =>  \$vverbose,
);

STDOUT->autoflush(1);

verifyInputOpts();
#--------------------------------
my $taxOptions = {
        TAXFILE    => $taxTreeFile,       # name of the input (Storable) Tax file
        VITALSFILE => $vitalsFile,        # genomeVitals.dmp
        TAXABBR    => \%taxAbbr,          # maps full rank names to their abbreviations
        REFRANKS   => \@ranks,            # tax ranks, from genus --> superkingdom
        TAXTREE    => \%taxTree,
};
#--------------------------------
loadOrganismNames($draftsFile,"drafts") if($draftsFile);            # Populates  @dirNames () and %root2orgDir
loadOrganismNames($completedFile,"completed") if($completedFile);   # Appends to @dirNames () and %root2orgDir
parseOrgGIs(\@dirNames);                        # Populates %allGIs
loadTaxTree($taxOptions);
parseReplTypes($taxOptions);                    # Populates %replTypeLookupByGI = ($gi => $replType);
loadVitals($taxOptions);
matchOrgs2TaxTreeByGI($taxOptions);
dereplicateOrgNames(\%xmlTemplate, $taxOptions);
exportXML(\%xmlTemplate, $taxOptions);

# ARGS: \%hash, $filename, $description, $mode      $mode = "bin" or "ascii"
dump2disk(XMLin($xml_outfile), "$prefix.xmlIn.readable.dmp", "xmlIN", "ascii");

print "All done!\n";

################################################################################
sub verifyInputOpts {

    usage() if($help);

    #if(!$infile && !$taxTreeFile && !$vitalsFile) {
    if(!$draftsFile && !$completedFile && !$taxTreeFile && !$vitalsFile) {
        usage();
    }

    $prefix      = $PREFIX unless $prefix;
    $xml_outfile = $prefix.".".$xmlExt if($prefix);
    $xml_outfile = $XML_OUTFILE unless $xml_outfile;

    $vverbose = 1 if($verbose);

    die "File \"$taxTreeFile\" does not exist!\n"   if(!-e $taxTreeFile);
    die "File \"$vitalsFile\" does not exist!\n"    if(!-e $vitalsFile);
    die "File \"$draftsFile\" does not exist!\n"    if($draftsFile && (!-e $draftsFile));
    die "File \"$completedFile\" does not exist!\n" if($completedFile && (!-e $completedFile));

    return;
        
}
################################################################################
# $infile contains the list of all the organism subdirs, listed one per line
# Example:
# |__________________________________________________________|_____________________________________________
# /home/traceyf/tmp/ftp.ncbi.nih.gov/genbank/genomes/Bacteria/Acetobacter_pasteurianus_IFO_3283_12_uid32203
# /home/traceyf/tmp/ftp.ncbi.nih.gov/genbank/genomes/Bacteria/Acholeplasma_laidlawii_PG_8A_uid19259
# /home/traceyf/tmp/ftp.ncbi.nih.gov/genbank/genomes/Bacteria/Acidithiobacillus_ferrooxidans_ATCC_23270_uid53
#
# ROOTDIR=
#
# ARGS: $infile, $COMPLETION_STATUS ("DRAFTS" or "COMPLETED")
################################################################################
sub loadOrganismNames {

    STDOUT->autoflush(1);

    # 1. Look for [$header], where $header = VIRUSES, PROKARYOTES, EUKARYOTES, or OTHERS
    #    a. save $header as current $kingdom
    #    b. read until [$header] or EOF; these non-empty lines are the org names
    #    c. repeat
    # 2. Look for 


    print "Identifying source directories in \"".$_[0]."\"...";
    my $iter0 = new Benchmark;
    open my $INFILE, '<', $_[0];
    
    my $headerFound = 0;
    my $validHeader = q{};
    
    LINE: while(my $line = <$INFILE>) {
        chomp $line;

        # Check for HEADER line
        if($line =~ m/^\[(\w+)\]/) {
            my $currHeader = uc($1);
            if(!exists $validKingdoms{$currHeader}) {
                die "Fatal: Unrecognized kingdom header in \"".$_[0]."\"!\n";
            }
            else {
                $validHeader = $currHeader;
            }
            next LINE;
        }
        
        next LINE unless $validHeader;
        next LINE if($line eq "");
        
        # ---- Reach here if current HEADER is VALID, and subsequent line(s) are non-empty ---- #
        
        # Remove basename -- this becomes the ROOTDIR. Everything else is the filename
        # ***** Directory (i.e. $line) cannot end with a "/" *****
        my ($basename, $rootdir) = fileparse($line);    # basename = organism dir;  rootdir = path to basename
        #print "DIR      = $dir\n";
        #print "BASENAME = $basename\n";
        #<STDIN>;
        my $dirname = $basename;
        $basename =~ tr/_/ /;                    # Remove underscores from name
        if($basename =~ m/^(.+)\s+(uid\d+)/) {   # Remove the "uid###" from name
            $basename = $1;
            
            $filehash{$_[1]}->{$validHeader}->{$rootdir}->{$dirname} = $basename;
            $org2dir{$rootdir}->{$basename} = $dirname;
            
            push(@dirNames, $basename);
        }
        else {
            print "Unrecognized organism directory format! Try removing trailing \"/\" from input dir names\n";
        }
    }
    close $INFILE;
    my $iter1 = new Benchmark;
    print "done. [".timestr(timediff($iter1,$iter0))."]\n";

    if($validHeader) {
        if($verbose) {
            print "organisms:\n";
            print "    \"$_\"\n" foreach (@dirNames);
        }
        print "Total of ".scalar(@dirNames)." organisms in database.\n";

        #print Dump(\%filehash); 
        #print "\n";
        #print "Press ENTER..."; <STDIN>;
    }
    else {
        print "No valid header found in file \"".$_[0]."\"!\n";
        print "Abort.\n";
        exit;
    }


    #print "searching for organisms in tax tree...\n";
    return;
}
################################################################################
# Function: Extracts the genomic sequence under the ORIGIN tag within *.gbk files
#           and saves them under new *.fna files written to the same directory.
# Input:    $gbkFilename    (fully-qualified filename)
# Returns:  $fastaFilename
################################################################################
sub gbk2fna {
    
    # Generate a *.fna-compatible header from the *.gbk file. Example:
    # GBK tags:     ___GI____     __VERSION__  _______________________DEFINITION__________________________
    #           >gi|116326069|ref|NC_008520.1| Anticarsia gemmatalis nucleopolyhedrovirus, complete genome    
    #

    # Verify the *.gbk file is still present: warn & skip if not
    if(!-e $_[0]) {
        print "WARNING: GenBank file \"".$_[0]."\" doesn't exist!\n";
        return q{};
    }

    # 1. Open the *.gbk file
    open my $GBKFILE, '<', $_[0] || die "GenBank file \"".$_[0]."\" cannot be opened!";
    
    # 2. Parse the DEFINITION line (DEF), VERSION (GI and ACC.VER), and ORIGIN
    my $repl     = q{};
    my $gi       = q{};
    my $seq      = q{};
    my $acc      = q{};
    my $nuclType = q{};
    my $defFound = 0;
    my $stopDef  = 0;
    my $originFound = 0;
    my $stopOrigin  = 0;
    LINE: while(my $line=<$GBKFILE>) {
        chomp $line;
        
        if($originFound) {
            # Sequence Extraction
            if($line =~ /^\/\//) {      # Termination of GBK record
                $stopOrigin = 1;
            }
            else {
                $seq .= $line;
            }
            last LINE if($stopOrigin);
        } #originFound
        else {
            # Process everything else
            #if($line =~ m/^LOCUS\s+\w+\s+(\d+)\sbp\s+\S+\s+(\w+)\s+\w+\s+(\S+)/) {
            #                     _ACC_  _SIZE_   DNA/RNA linear     _DATE_
            if($line =~ m/^LOCUS\s+\w+\s+\d+\sbp\s+(\S+)\s+\w+\s+\w+\s+\S+/x) {
                $nuclType = "$1";
            }
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
            elsif($line =~ m/^ORIGIN/) {
                $originFound = 1;
            }
            # Continuation of DEFINITION
            elsif(!$stopDef && $defFound && ($line =~ m/^\s+(.+)/)) {
                #print "ORIG=\"$repl\"\n";
                $repl .= " ".$1;
                #print " NEW=\"$repl\"\n";
            }
        } # !$originFound
        
    } #LINE
    close $GBKFILE;    

    $seq  =~ tr/a-z0-9 /A-Z/d;          # Remove all digits and spaces
    $repl =~ s/\.$//;                   # Remove tailing periods

    my $seqLen  = length($seq);         # Check valid sequence length (> 0)
    return q{} if($seqLen == 0);        # Empty sequence = FAILURE

    # ----------------------
    # Write FNA file to disk
    # ----------------------
    
    my @exts = (".gbk");
    my ($filename, $directory, $suffix) = fileparse($_[0], @exts);  # Deconstruct GBK file into components
    my $fnaFilename = $directory."$filename.fna";                   # Create FNA filename
    my $header      = ">gi|$gi|custom|$acc| $repl ($nuclType)";     # Generate FASTA header
    my $lineLen = 70;                                               # No. of bases to write per line

    # Write FNA file to disk
    open my $FNAFILE, '>', $fnaFilename || die "Fatal: Cannot write \"".$fnaFilename."\" to disk!\n";
    print $FNAFILE $header."\n";
    for(my $offset=0; $offset < $seqLen; $offset=$offset+$lineLen) {
        print $FNAFILE substr($seq, $offset, $lineLen)."\n";
    }
    close $FNAFILE;
            
    return $fnaFilename;
    
}
################################################################################
# ARGS: $status, $kingdom, $orgDir, $noDenovoFNA
################################################################################
sub displayFNAwarning {
    print "\n  **WARNING**: [".$_[1]."] ".$_[0]." project ".$_[2]." has no genome (*.fna)"
         ."files.\n";
    print "               Skipping de novo FNA file creation..." if($_[3]);
    print "               Attempting to extract from GenBank file..." if(!$_[3]);     
    return;
}
################################################################################
# Function: populates %allGIs
# ARGS: \@dirNames
################################################################################
sub parseOrgGIs {
    
    my @noFastaExts = ();   # Organism directories (fully qualified) that do NOT contain any valid FASTAEXTs
    my @noGIs       = ();   # Filenames (fully qualified) where no valid GI could be parsed
    my $count       = 0; 
    
    print "Locating FASTA files (".join(" > ",@fastaexts).") and parsing for GIs [each \".\" equals $lineBuffer files]...";
    #print "Parsing GIs from all *.$fastaext files (each \".\" equals 100 files) ...";
    my $iter0 = new Benchmark;
    
    # Loop through input files and parse the GIs in from each
    STATUS: foreach my $status (keys %filehash) {                       # DRAFTS or COMPLETED
        KINGDOM: foreach my $kingdom (keys %{ $filehash{$status} }) {   # VIRUSES, PROKARYOTES, ...
            ROOTDIR: foreach my $rootDir (keys %{ $filehash{$status}->{$kingdom} }) {
                SUBDIR: foreach my $subDir (keys %{ $filehash{$status}->{$kingdom}->{$rootDir} }) {

                    my $orgDir = $rootDir.$subDir;
            
                    # Get list of all *.$fastaext files under $orgDir
                    my $extFound   = 0;
                    my @inputFiles = ();

                    # Some organisms (viruses) on the NCBI FTP site do not have *.fna files in their 
                    # subdirectory, but so far, all have *.gbk files. So if there is no *.fna file, 
                    # we look for *.gbk file(s) and generate *.fna files out of them... 
                    # ... unless $noDenovoFNA is specified
                    FASTAEXT: foreach my $fastaext (@fastaexts) {
                        @inputFiles  = (`ls -x1 $orgDir/*.$fastaext 2>&1`);
                        if((scalar(@inputFiles) == 1) && ($inputFiles[0] =~ m/No such file or directory/)) {
                            @inputFiles = ();
                            next FASTAEXT;
                        }
                        else {
                            $extFound = 1;
                            chomp $_ foreach(@inputFiles);

                            # If GBK file(s) found, extract sequence to FNA file; return new FNA filename
                            if($fastaext eq "gbk") {
                                displayFNAwarning($status, $kingdom, $orgDir, $noDenovoFNA);
                                if($noDenovoFNA) {
                                    @inputFiles = ();
                                    last FASTAEXT;
                                }
                                my @fnaFiles = ();
                                foreach (@inputFiles) {
                                    my $fnaFile = gbk2fna($_);                     # Converts *.gbk file to *.fna file
                                    if($fnaFile ne "") {                           # Check for conversion failure
                                        push(@fnaFiles, $fnaFile);
                                        print "success!\n";
                                    } #if
                                } #foreach
                                @inputFiles = @fnaFiles;
                            } #if GBK
                            
                            # If FFN, extract unique ID from ">gi|########|ref|XXXXXXXX|:START-STOP"
                            # ... deferred ...
                            
                            last FASTAEXT;
                        }
                    } #FASTAEXT
                    
                    if(!$extFound) {
                        push(@noFastaExts, $orgDir);
                        next SUBDIR;
                    }
            
                    # Record those org subdirs that have *.fna files inside
                    push(@validDirNames, $subDir) if(@inputFiles);
            
                    # Load each fasta file (handle multiple ">" entries per file)
                    foreach (@inputFiles) {

                        my $filename = fileparse($_);       # Remove path info

                        # UPDATED: add support for contigs (2012-01-25)
                        # parse multiFASTA file: return header, GI, and description
                        my ($gi_aref, $contigName_aref, $contigDesc_aref) = parseHeader($_);
                        my @gis         = @{ $gi_aref };
                        my @contigNames = @{ $contigName_aref };
                        my @contigDescs = @{ $contigDesc_aref };

                        #print "GIs = ".join(",",@gis)."\n";
                        #print "Names = ".join(",",@contigNames)."\n";
                        #print "Descs = ".join(",",@contigDescs)."\n";

                        if(@gis) {
                            for(my $idx = 0; $idx < (scalar(@gis)); $idx++) {
                                my $gi         = $gis[$idx];
                                my $contigName = $contigNames[$idx];
                                my $contigDesc = $contigDescs[$idx];
                                if(exists $allGIs{$gi}) {
                                    print "\n  **WARNING**: Attempted to clobber pre-existing GI ($gi)!\n";
                                    print "           OLD = ".$allGIs{$gi}->{FILE}." \"".$allGIs{$gi}->{NAME}." ".$allGIs{$gi}->{DESC}."\"\n";
                                    print "                 [".$allGIs{$gi}->{ROOTDIR}.$allGIs{$gi}->{ORGDIR}."]\n";
                                    print "           NEW = $filename \"$contigName $contigDesc\"\n";
                                    print "                 [$rootDir$subDir]\n";
                                    print "           SKIPPING new GI...\n";
                                }
                                else {
                                    $allGIs{$gi}->{STATUS}  = $status;
                                    $allGIs{$gi}->{KINGDOM} = $kingdom;
                                    $allGIs{$gi}->{ROOTDIR} = $rootDir;
                                    $allGIs{$gi}->{ORGDIR}  = $subDir;
                                    $allGIs{$gi}->{FILE}    = $filename;
                                    $allGIs{$gi}->{NAME}    = $contigName;
                                    $allGIs{$gi}->{DESC}    = $contigDesc;
                                    #print "GI=$gi\n";
                                    
                                    # UPDATE: 2012-03-04
                                    #$orgLookup{$subDir}
                                }
                            } #for
                        }
                        else {
                            push(@noGIs, $_);
                        }
                        print "." if($count % $lineBuffer == 0);
                        $count++;
                    } #inputFiles
                } #subDir
            } #rootDir
        } #KINGDOM
    } #STATUS
    

    my $iter1 = new Benchmark;
    print "done. [".timestr(timediff($iter1,$iter0))."]\n";

    if(@noFastaExts) {
        print "The following organism directories did not contain any valid FASTA extensions (".join(" > ",@fastaexts)."):\n";
        print "    $_\n" foreach (@noFastaExts);
    }

    if(@noGIs) {
        print "GIs could not be parsed out of the following FASTA files:\n";
        print "    $_\n" foreach (@noGIs);
        unless ($skipOrgNoGI) {
            print "Abort.\n";
            exit;
        }
    }
    
    return;
}
################################################################################
# ARGS: Fully qualified filename
# Returns: \@gis, \@contigNames, \@contigDescriptions
################################################################################
sub parseHeader {

    my @gis          = ();
    my @contigNames  = ();
    my @contigDescriptions = ();

    open my $INFILE, '<', $_[0] || die "Cannot open FASTA file $_[0]!\n";
    LINE: while(my $line=<$INFILE>) {
        chomp $line;
        next LINE if($line !~ /^\>/);

        # ____________NAME____________ _______________DESCRIPTION_____________________
        # >gi|193001753|gb|CP001047.1| Mycoplasma arthritidis 158L3-1, complete genome
        if($line =~ m/^\>(gi\|(\d+)\|\S+\|\S+)\s+(.+)$/) {
            push(@contigNames, $1);
            push(@gis, $2);
            push(@contigDescriptions, $3);
            next LINE;
        }

    } #while
    close $INFILE;

    return (\@gis, \@contigNames, \@contigDescriptions);

}
################################################################################
# Input:    (multi)FASTA filename
# Function: Parses the GI out of a FASTA header. Supports single multiFASTA file
# Returns:  array of valid GIs or q{}
# ARGS: $fastaFilename
################################################################################
sub parseGIfromHeader {

    my @gis = ();
    #print "opening $_..."; <STDIN>;
    
    open my $INFILE, '<', $_[0] || die "Cannot open FASTA file $_[0]!\n";
    LINE: while(my $line=<$INFILE>) {
        chomp $line;
        next LINE if($line !~ /^\>/);

        # >gi|193001753|gb|CP001047.1| Mycoplasma arthritidis 158L3-1, complete genome
        if($line =~ m/\>gi\|(\d+)\|/) {
            push(@gis, $1);
            next LINE;
        }

    } #while
    close $INFILE;
    
    return @gis;
}
################################################################################
# ARGS: \%taxTreeByGI, \%speciesTreeGI
################################################################################
sub indexTaxTreeByGI {

    my @sources = ("CHR", "PLASMID");
    my @sTypes  = ("GEN", "GBK");
    my %giHash = ();

    # grab all GIs from speciesTree
    NODE: foreach my $node (keys %{ $_[1] }) {
        next NODE unless exists ($_[1]->{$node}->{GI});
        
        # search through CHR/PLASMID first
        SOURCE: foreach my $source (@sources) {

            # now search through GEN/GBK
            STYPE: foreach my $sType (@sTypes) {
            
                # grab underlying GIs (and ORGs) for SOURCE TYPE
                my @gis    = keys   %{ $_[1]->{$node}->{GI}->{$source}->{$sType} };
                next STYPE if(!@gis);
                my @values = values %{ $_[1]->{$node}->{GI}->{$source}->{$sType} };
            
                foreach my $idx (0..(scalar(@gis)-1)) {
                    my $gi  = $gis[$idx];
                    my $org = $values[$idx];                    
                    $_[0]->{$gi}->{ORG}    = $org;      # Organism name
                    $_[0]->{$gi}->{SOURCE} = $source;   # CHR or PLASMID
                    $_[0]->{$gi}->{STYPE}  = $sType;    # GBK or GEN
                    $_[0]->{$gi}->{NODE}   = $node;     # matched Node
                } #gi

            } #sType
        } #source
    } #node

    #dump2disk($_[0], "taxTreeByGI.readable.dmp", "TAX TREE by GI", "ascii");
    #exit;
    
    return;
}
################################################################################
# ARGS: $taxOptions, $node, $orgName
################################################################################
sub determineRank {

    my ($options, $node, $orgName) = @_;
    my $rank = q{};

    # Rank is obtained by locating the valid tax name $orgName in either the
    # {SS} or {S} entries, searched in that order.
    my %subspecies = ();
    my @subspecies = keys %{ $options->{TAXTREE}->{$node}->{ $taxAbbrExt{"subspecies"} } };
       @subspecies{@subspecies} = ();
            
    if(exists $subspecies{$orgName}) {
        $rank = $taxAbbrExt{"subspecies"};      # i.e. "SS"
    }
    else {
        # Verify in species {S}
        my %species = ();
        my @species = keys %{ $options->{TAXTREE}->{$node}->{ $taxAbbrExt{"species"} } };
           @species{@species} = ();
        if(exists $species{$orgName}) {
            $rank = $taxAbbrExt{"species"};     # i.e. "S"
        }
        else {
            $rank = "UNKNOWN";
        }
    }
            
    return $rank;
}
################################################################################
# Input: Either   (i) a list of directories of organism names
#           or   (ii) the root directory of the list of directories of organism names
#           or  (iii) the fully qualified path of the *.fna files
#
# Function: If   (i), finds all *.$fastaext files in each directory. Loads only ONE from
#                     each directory, parses out the GI, then assigns that directory the
#                     taxonomy.
#           If  (ii), finds all subdirs that have the *uid$ in them, and repeats (i).
#           If (iii), separates the basename from the FASTA filename proper, and foreach
#                     unique subdir, loads a FASTA file, parses the GI, and assigns the
#                     taxonomy to that subdir.
#
# ARGS: $taxOptions
################################################################################
sub matchOrgs2TaxTreeByGI {

    my %taxTreeByGI = ();       # %taxTreeByGI = ( $gi => { ORG     => $organismName,  [fully valid tax name]
                                #                           SOURCE  => $source,        [CHR or PLASMID]
                                #                           STYPE   => $sourceType,    [GBK or GEN]
                                #                         },
                                #                );

    # GIs in query that can't be found in the tree
    my @unmatched = ();
                                
    # Populates %taxTreeByGI
    indexTaxTreeByGI(\%taxTreeByGI, $_[0]->{TAXTREE});
    
    my $doRC = 1;   # current default for all
    
    # Loop through all GIs and match to tax tree
    foreach my $gi (keys %allGIs) {
      my $status     = $allGIs{$gi}->{STATUS};      # "drafts" or "completed"
      my $kingdom    = $allGIs{$gi}->{KINGDOM};     # "VIRUSES", "PROKARYOTES", ...
      my $rootdir    = $allGIs{$gi}->{ROOTDIR};     # full path to org subdir
      my $orgdir     = $allGIs{$gi}->{ORGDIR};      # org subdir
      my $filename   = $allGIs{$gi}->{FILE};        # file in org subdir containing GI
      my $contigName = $allGIs{$gi}->{NAME};        # Contig name: gi|193001753|gb|CP001047.1|
      my $contigDesc = $allGIs{$gi}->{DESC};        # Contig description: Mycoplasma arthritidis 158L3-1, complete genome

      # UPDATE: 2012-03-04
      # ...
      
      # Segregate those files whose taxonomic links could be established (LINKED)
      # with those that could not (UNLINKED).
      if(exists $taxTreeByGI{$gi}) {
        # reorganized speciesTree.dmp
        my $orgName = $taxTreeByGI{$gi}->{ORG};             # scientific name
        my $sType   = $taxTreeByGI{$gi}->{STYPE};           # GEN or GBK        or $_[0]->{VITALS}->{GI}->{$gi}->{SOURCE}   <--------- VERIFY THIS
        my $node    = $taxTreeByGI{$gi}->{NODE};            # node
        my $rank    = determineRank($_[0], $node, $orgName);# S or SS
        
        # Track for each reverse lookup in sub dereplicateOrgNames()
        $orgName2TaxData{$orgName}->{NODE} = $node;
        $orgName2TaxData{$orgName}->{RANK} = $rank;
        
        # genomeVitals.dmp
        my $taxid   = $_[0]->{VITALS}->{GI}->{$gi}->{TAXID};# NCBI taxonomic ID
        my $date    = $_[0]->{VITALS}->{GI}->{$gi}->{DATE}; # published date in GEN/GBK's *.gbk file
        my $nType   = $_[0]->{VITALS}->{GI}->{$gi}->{NTYPE};# Nucleic Acid Type: DNA, RNA, ss-RNA, ...
        my $topology= $topoAbbr{ $_[0]->{VITALS}->{GI}->{$gi}->{TOPO} }; # "L" or "C"
        my $acc     = $_[0]->{VITALS}->{GI}->{$gi}->{ACC};  # accession # in GEN/GBK's *.gbk file
        # local functions
        my $replType= $_[0]->{REPLTYPELOOKUP}->{$gi};       # "CHR" or "PLASMID"
        #my $topology= topologyLookup($status, $kingdom, $replType); # "L" or "C"

        # Support multiFASTA input files (contigs, discontiguous viral segments, etc.) 
        # by associating multiple GIs to a single file.
        $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$orgdir}->{NAME}=$orgName;
        $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$orgdir}->{TAXID}=$taxid;
        $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$orgdir}->{RANK}=$rank;
        $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$orgdir}->{FILE}
                            ->{$filename}->{STYPE}    = $sType;
        $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$orgdir}->{FILE}
                            ->{$filename}->{STATUS}   = $status;
        $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$orgdir}->{FILE}
                            ->{$filename}->{REPLTYPE} = $replType;
        $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$orgdir}->{FILE}
                            ->{$filename}->{TOPO}     = $topology;
        $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$orgdir}->{FILE}
                            ->{$filename}->{NTYPE}= $nType;
        $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$orgdir}->{FILE}
                            ->{$filename}->{RC} = $doRC;
        $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$orgdir}->{FILE}
                            ->{$filename}->{CONTIGS}->{$gi}->{ACC} = $acc;
        $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$orgdir}->{FILE}
                            ->{$filename}->{CONTIGS}->{$gi}->{DATE} = $date;
        $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$orgdir}->{FILE}
                            ->{$filename}->{CONTIGS}->{$gi}->{NAME} = $contigName;
        $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$orgdir}->{FILE}
                            ->{$filename}->{CONTIGS}->{$gi}->{DESC} = $contigDesc;
      }
      else {
        $xmlTemplate{UNLINKED}->{$status}->{$kingdom}->{$rootdir}->{$orgdir}->{FILE}->{$filename}->{STYPE}   = q{};
        $xmlTemplate{UNLINKED}->{$status}->{$kingdom}->{$rootdir}->{$orgdir}->{FILE}->{$filename}->{STATUS}  = $status;
        $xmlTemplate{UNLINKED}->{$status}->{$kingdom}->{$rootdir}->{$orgdir}->{FILE}->{$filename}->{REPLTYPE}= q{};
        $xmlTemplate{UNLINKED}->{$status}->{$kingdom}->{$rootdir}->{$orgdir}->{FILE}->{$filename}->{TOPO}    = q{};
        $xmlTemplate{UNLINKED}->{$status}->{$kingdom}->{$rootdir}->{$orgdir}->{FILE}->{$filename}->{NTYPE}   =q{};
        $xmlTemplate{UNLINKED}->{$status}->{$kingdom}->{$rootdir}->{$orgdir}->{FILE}->{$filename}->{RC}      =$doRC;
        $xmlTemplate{UNLINKED}->{$status}->{$kingdom}->{$rootdir}->{$orgdir}->{FILE}->{$filename}->{CONTIGS}->{$gi}->{ACC} =q{};
        $xmlTemplate{UNLINKED}->{$status}->{$kingdom}->{$rootdir}->{$orgdir}->{FILE}->{$filename}->{CONTIGS}->{$gi}->{DATE}=q{};
        $xmlTemplate{UNLINKED}->{$status}->{$kingdom}->{$rootdir}->{$orgdir}->{FILE}->{$filename}->{CONTIGS}->{$gi}->{NAME}=q{};
        $xmlTemplate{UNLINKED}->{$status}->{$kingdom}->{$rootdir}->{$orgdir}->{FILE}->{$filename}->{CONTIGS}->{$gi}->{DESC}=q{};
        push(@unmatched, $gi);
      }
    }
    return;
    
}
################################################################################
# ARGS: \%xmlTemplate, $taxOptions
################################################################################
sub exportXML {

    # ------------------------------------------
    # PRINT OUT XML FILE OF TAX-(UN)MATCHED ORGS
    # ------------------------------------------
    my $xml_fh = init_file_write($xml_outfile);
    print {$xml_fh} "<taxidx>\n";
    appendUnlinked($xml_fh);    # write to XML file those orgs w/GIs not matching tax tree
    appendLinked($xml_fh);      # write to XML file those orgs w/GIs matching tax tree
    print {$xml_fh} "</taxidx>\n";
    close $xml_fh;
    print "Created file \"$xml_outfile\".\n";

    # ------------------------------------------------
    # Report Unmatched (orgs w/o TAXID)
    # ------------------------------------------------
    if(@unmatched) {
        print "The following ".(scalar(@unmatched))." organisms went unmatched:\n";
        print "========================================\n";
        foreach (@unmatched) {
            print "    $_\n";
        }
        print "STOPPING.\n" unless ($continueNoHostFound);
        die "\n"              unless ($continueNoHostFound);
        print "CONTINUING...\n";
    }
    else {
        print "CONGRATULATIONS! All ".(scalar(@validDirNames))." organisms were assigned "
                ."taxonomic identities at the species level.\n";
    }
    
    return;
}
################################################################################
# ARGS: $xml_fh             filehandle to the XML outfile
################################################################################
sub appendUnlinked {

    print "Exporting UNLINKED organisms...";
    #------------ UNLINKED --------------
    print {$_[0]} "  <unlinked>\n";
    foreach my $status (keys %{ $xmlTemplate{UNLINKED} }) {
      print {$_[0]} "    <$status>\n";
      foreach my $kingdom (sort {$a cmp $b} keys %{ $xmlTemplate{UNLINKED}->{$status} }) {
        print {$_[0]} "      <$kingdom>\n";
        foreach my $rootdir (sort {$a cmp $b} keys %{ $xmlTemplate{UNLINKED}->{$status}->{$kingdom} }) {
          print {$_[0]} "        <root>\n";
          print {$_[0]} "          <name>$rootdir</name>\n";
          foreach my $dirname (sort {$a cmp $b} keys %{ $xmlTemplate{UNLINKED}->{$status}->{$kingdom}->{$rootdir} }) {
            print {$_[0]} "          <org>\n";
            print {$_[0]} "            <dir>$dirname</dir>\n";
            print {$_[0]} "                  <name></name>\n";
            print {$_[0]} "                  <rank></rank>\n";
            print {$_[0]} "                  <taxid></taxid>\n";
            foreach my $filename (keys %{ $xmlTemplate{UNLINKED}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{FILE} }) {
              my $rc = $xmlTemplate{UNLINKED}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{FILE}->{$filename}->{RC};
              print {$_[0]} "                  <file>\n";
              print {$_[0]} "                    <name>$filename</name>\n";                    # OK
              print {$_[0]} "                    <stype></stype>\n";                           # OK
              print {$_[0]} "                    <status>".$statusAbbr{$status}."</status>\n"; # OK
              print {$_[0]} "                    <repltype></repltype>\n";                     # OK
              print {$_[0]} "                    <topo></topo>\n";                             # OK
              print {$_[0]} "                    <ntype></ntype>\n";                           # OK
              print {$_[0]} "                    <rc>$rc</rc>\n";                              # OK
              foreach my $gi ($xmlTemplate{UNLINKED}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{FILE}->{$filename}->{CONTIGS}) {
                print {$_[0]} "                    <contig>\n";                                # OK
                print {$_[0]} "                      <gi></gi>\n";                             # OK
                print {$_[0]} "                      <name></name>\n";                             # OK
                print {$_[0]} "                      <desc></desc>\n";                             # OK
                print {$_[0]} "                      <acc></acc>\n";                           # OK
                print {$_[0]} "                      <date></date>\n";                         # OK
                print {$_[0]} "                    </contig>\n";                               # OK
                print {$_[0]} "                  </file>\n";                                   # OK
              }
            } #FILENAME
            print {$_[0]} "          </org>\n";
          } #dirname
          print {$_[0]} "        </root>\n";
        } #rootdir      
        print {$_[0]} "      </$kingdom>\n";
      } #kingdom
      print {$_[0]} "    </$status>\n";
    } #status
    print {$_[0]} "  </unlinked>\n";
    # ...
    print "done.\n";

    return;
}
################################################################################
# ARGS: $xml_fh             filehandle to the XML outfile
################################################################################
sub appendLinked {

#    $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$orgdir}->{FILE}->{$filename}->{CONTIGS}->{$gi}->{RC} = $doRC;

    #------------ LINKED --------------
    print "Exporting LINKED organisms...";
    print {$_[0]} "  <linked>\n";

    foreach my $status (keys %{ $xmlTemplate{LINKED} }) {                   # DRAFTS or COMPLETED
      print {$_[0]} "    <$status>\n";
        foreach my $kingdom (sort {$a cmp $b} keys %{ $xmlTemplate{LINKED}->{$status} }) {
          print {$_[0]} "      <$kingdom>\n";
          foreach my $rootdir (sort {$a cmp $b} keys %{ $xmlTemplate{LINKED}->{$status}->{$kingdom} }) {
            print {$_[0]} "        <root>\n";
            print {$_[0]} "          <name>$rootdir</name>\n";
            foreach my $dirname (sort {$a cmp $b} keys %{ $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir} }) {
              my $sciname = $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{NAME};
              my $rank    = $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{RANK};
              my $taxid   = defined($xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{TAXID})
                          ? $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{TAXID}
                          : "n/a";
              print {$_[0]} "          <org>\n";
              print {$_[0]} "            <dir>$dirname</dir>\n";
              print {$_[0]} "            <name>$sciname</name>\n";
              print {$_[0]} "            <rank>$rank</rank>\n";
              print {$_[0]} "            <taxid>$taxid</taxid>\n";
              foreach my $filename (keys %{ $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{FILE} }) {
                my $sType    = $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{FILE}->{$filename}->{STYPE};
                my $replType = $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{FILE}->{$filename}->{REPLTYPE};
                my $topology = $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{FILE}->{$filename}->{TOPO};
                my $nType    = $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{FILE}->{$filename}->{NTYPE};
                my $doRC     = $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{FILE}->{$filename}->{RC};
                print {$_[0]} "            <file>\n";
                print {$_[0]} "              <name>$filename</name>\n";                        # OK
                print {$_[0]} "              <stype>$sType</stype>\n";                         # OK
                print {$_[0]} "              <status>".$statusAbbr{$status}."</status>\n";     # OK
                print {$_[0]} "              <repltype>$replType</repltype>\n";                # <---- extract from "speciesTreeGI.dmp"
                print {$_[0]} "              <topo>$topology</topo>\n";                        # <---- verify
                print {$_[0]} "              <ntype>$nType</ntype>\n";                         # OK
                print {$_[0]} "              <rc>$doRC</rc>\n";                                # OK
                foreach my $gi (keys %{ $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{FILE}->{$filename}->{CONTIGS} }) {
                  my $acc  = $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{FILE}->{$filename}->{CONTIGS}->{$gi}->{ACC};
                  my $date = $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{FILE}->{$filename}->{CONTIGS}->{$gi}->{DATE};
                  my $name = $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{FILE}->{$filename}->{CONTIGS}->{$gi}->{NAME};
                  my $desc = $xmlTemplate{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{FILE}->{$filename}->{CONTIGS}->{$gi}->{DESC};
                  print {$_[0]} "              <contig>\n";                                    #
                  print {$_[0]} "                  <gi>$gi</gi>\n";                            # OK
                  print {$_[0]} "                  <name>$name</name>\n";
                  print {$_[0]} "                  <desc>$desc</desc>\n";
                  print {$_[0]} "                  <acc>$acc</acc>\n";                         # OK
                  print {$_[0]} "                  <date>$date</date>\n";                      # OK
                  print {$_[0]} "              </contig>\n";                                   #
                } #gi
                print {$_[0]} "            </file>\n";
              } #filename
              print {$_[0]} "          </org>\n";
            } #dirname
            print {$_[0]} "        </root>\n";
          } #rootdir
          print {$_[0]} "      </$kingdom>\n";
        } #kingdom
      print {$_[0]} "    </$status>\n";
    } #status
    print {$_[0]} "  </linked>\n";
    print "done.\n";

    return;
}
################################################################################
# FUTURE IMPLEMENTATION ARGS:  $gi
#                                $_[0]    $_[1]     $_[2]
# CURRENT IMPLEMENTATION ARGS: $status, $kingdom, $replType
################################################################################
sub topologyLookup {
    
    # FUTURE IMPLEMENTATION: Search a pre-computed hash that indexes GI with 
    #                        replicon type: (L)inear or (C)ircular.
    # ...
    
    # CURRENT IMPLEMENTATION: Everything returns linear unless its a finished
    #                         Prokaryotic chromosome.
    my $topology  = (($_[0] eq "COMPLETED") && ($_[1] eq "PROKARYOTES") && ($_[2] eq "CHR"))
                  ? "C"
                  : "L";
    
    return $topology;
}
################################################################################
# ARGS: \%xmlTemplate, $taxOptions
################################################################################
sub dereplicateOrgNames {

    my %sciNames        = ();
    my $duplicatesFound = 0;

    # Cycle through all $status, $kingdom, $rootdir, $dirname, $sciname and build a hash of encountered $scinames
    foreach my $status (keys %{ $_[0]->{LINKED} }) {                   # DRAFTS or COMPLETED
        foreach my $kingdom (sort {$a cmp $b} keys %{ $_[0]->{LINKED}->{$status} }) {
            foreach my $rootdir (sort {$a cmp $b} keys %{ $_[0]->{LINKED}->{$status}->{$kingdom} }) {
                foreach my $dirname (sort {$a cmp $b} keys %{ $_[0]->{LINKED}->{$status}->{$kingdom}->{$rootdir} }) {
                    my $sciName = $_[0]->{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{NAME};
                    
                    # Store no. of encounters, and their sources
                    $sciNames{$sciName}->{COUNT}++;
                    $sciNames{$sciName}->{SOURCE}->{$status}->{$kingdom}->{$rootdir}->{$dirname} = ();
                    $sciNames{$sciName}->{NODE} = $orgName2TaxData{$sciName}->{NODE};
                    $sciNames{$sciName}->{RANK} = $orgName2TaxData{$sciName}->{RANK};

                    foreach my $filename (keys %{ $_[0]->{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{FILE} }) {
                        my $sType    = $_[0]->{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{FILE}->{$filename}->{STYPE};
                        my $replType = $_[0]->{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{FILE}->{$filename}->{REPLTYPE};

                        foreach my $gi (keys %{ $_[0]->{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{FILE}->{$filename}->{CONTIGS} }) {
                            $sciNames{$sciName}->{SOURCE}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{GI}->{$replType}->{$sType}->{$gi} = ();
                        } #GI

                    } #FILENAME
                } #DIRNAME
            } #ROOTDIR
        } #KINGDOM
    } #STATUS

    NAME: foreach my $name (keys %sciNames) {
        next NAME unless $sciNames{$name}->{COUNT} > 1;

        $duplicatesFound = 1;

        print "Duplicate organism name found! [\"$name\"]\n";
        print "    Reformatting as unique using UID; Updating XML file and Tax Tree File...";
        foreach my $status (keys %{ $sciNames{$name}->{SOURCE} }) {
            foreach my $kingdom (keys %{ $sciNames{$name}->{SOURCE}->{$status} }) {
                foreach my $rootdir (keys %{ $sciNames{$name}->{SOURCE}->{$status}->{$kingdom} }) {
                    DIRNAME: foreach my $dirname (keys %{ $sciNames{$name}->{SOURCE}->{$status}->{$kingdom}->{$rootdir} }) {
                        my $uid = q{};
                        if($dirname =~ /_(uid\d+)$/) {
                            $uid = $1;
                        }
                        else {
                            die "**FATAL**: Unrecognized directory name \"$dirname\"!\n";
                        }

                        # 1: Update $xmlTemplate
                        my $uniqueName = $name." ".$uid;
                        $_[0]->{LINKED}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{NAME} = $uniqueName;
                        
                        # 2a: Update "speciesTreeGI.dmp" NAME
                        my $taxNode = $sciNames{$name}->{NODE};
                        my $taxRank = $sciNames{$name}->{RANK};
                        
                        ## Update "genomeVitals.dmp" ORG
                        #my $href = $_[1]->{VITALS}->{ORG}->{$name};
                        #if(exists $_[1]->{VITALS}->{ORG}->{$uniqueName}) {
                        #    print "\n**WARNING**: Attempting to update organism name [$name] with new name [$uniqueName],\n"
                        #         ."             but new name already exists! Skipping...\n";
                        #}
                        #else {
                        #    $_[1]->{VITALS}->{ORG}->{$uniqueName} = $href;  # Add new
                        #    delete $_[1]->{VITALS}->{ORG}->{$name};         # Delete old
                        #}
                        
                        # Check if name already exists at rank; add if not
                        $_[1]->{TAXTREE}->{$taxNode}->{$taxRank}->{$uniqueName} = "scientific name - custom" unless (exists $_[1]->{TAXTREE}->{$taxNode}->{$taxRank}->{$uniqueName});
                        
                        # Still have to update all references to $name with $uniqueName in {GI}
                        # $sciNames{$sciName}->{GI}->{$replType}->{$sType}->{$gi} = ();
                        
                        foreach my $replType (keys %{ $sciNames{$name}->{SOURCE}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{GI} }) {                                   # CHR/PLASMID
                            foreach my $sType (keys %{ $sciNames{$name}->{SOURCE}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{GI}->{$replType} }) {                     # GEN/GBK
                                foreach my $gi (keys %{ $sciNames{$name}->{SOURCE}->{$status}->{$kingdom}->{$rootdir}->{$dirname}->{GI}->{$replType}->{$sType} }) {

                                    # 2b: Update "speciesTreeGI.dmp" GIs
                                    $_[1]->{TAXTREE}->{$taxNode}->{GI}->{$replType}->{$sType}->{$gi} = $uniqueName;

                                    # 3a: Update {GI}->{$gi}->{ORG}
                                    $_[1]->{VITALS}->{GI}->{$gi}->{ORG} = $uniqueName;

                                    # 3b: Update {REPL}->{$stype}->{$gi}->{ORG}
                                    my $repl = $_[1]->{VITALS}->{GI}->{$gi}->{REPL};
                                    $_[1]->{VITALS}->{REPL}->{$repl}->{$sType}->{$gi}->{ORG} = $uniqueName;                                    
                                    
                                    ## 3c: Update {ORG}->{$org}->{REPL}->...
                                    $_[1]->{VITALS}->{ORG}->{$uniqueName}->{REPL}->{$repl}->{$sType}->{$gi}->{ACC}   = $_[1]->{VITALS}->{GI}->{$gi}->{ACC};
                                    $_[1]->{VITALS}->{ORG}->{$uniqueName}->{REPL}->{$repl}->{$sType}->{$gi}->{DATE}  = $_[1]->{VITALS}->{GI}->{$gi}->{DATE};
                                    $_[1]->{VITALS}->{ORG}->{$uniqueName}->{REPL}->{$repl}->{$sType}->{$gi}->{NTYPE} = $_[1]->{VITALS}->{GI}->{$gi}->{NTYPE};
                                    $_[1]->{VITALS}->{ORG}->{$uniqueName}->{REPL}->{$repl}->{$sType}->{$gi}->{SIZE}  = $_[1]->{VITALS}->{GI}->{$gi}->{SIZE};
                                    $_[1]->{VITALS}->{ORG}->{$uniqueName}->{REPL}->{$repl}->{$sType}->{$gi}->{TOPO}  = $_[1]->{VITALS}->{GI}->{$gi}->{TOPO};

                                } #GI
                            } #STYPE
                        } #REPLTYPE

                        delete $_[1]->{VITALS}->{ORG}->{$name};

                    } #DIRNAME
                } #ROOTDIR
            } #KINGDOM
        } #STATUS

        print "done.\n";

    } #NAME

    # Re-write $_[1]->{TAXTREE} if duplicates found
    if($duplicatesFound) {
   
        my ($filename, $dir, $ext) = fileparse($_[1]->{TAXFILE}, qr/\.[^.]*/);
        my $binName   = $dir.$filename."2".$ext;
        my $asciiName = $dir.$filename."2.readable".$ext;

        #        $_[0]    $_[1]        $_[2]     $_[3]
        # ARGS: \%hash, $filename, $description, $mode      $mode = "bin" or "ascii"
        dump2disk($_[1]->{TAXTREE}, $binName, "Updated Species Tree (w/GI)", "bin");
        dump2disk($_[1]->{TAXTREE}, $asciiName, "Updated Species Tree (w/GI)", "ascii");
        
        my ($filename2, $dir2, $ext2) = fileparse($_[1]->{VITALSFILE}, qr/\.[^.]*/);
        my $binName2   = $dir2.$filename2."2".$ext2;
        my $asciiName2 = $dir2.$filename2."2.readable".$ext2;

        #        $_[0]    $_[1]        $_[2]     $_[3]
        # ARGS: \%hash, $filename, $description, $mode      $mode = "bin" or "ascii"
        dump2disk($_[1]->{VITALS}, $binName2, "Updated Genome Vitals", "bin");
        dump2disk($_[1]->{VITALS}, $asciiName2, "Updated Genome Vitals", "ascii");
        
        print "\n**WARNING**: Duplicate strains were found in your dataset.\n"
             ."             Inclusion of these in your generation of a unique DB will likely result\n"
             ."             in ~100% removal of their host DNA, effectively eliminating both organisms\n"
             ."             from the dataset.\n"
             ."             Please consider retaining only ONE representative from each duplicate listed above.\n\n";
    }

    return;
    
}
################################################################################
# ARGS: $taxOptions
# Populates %taxTree = ( $species1 => { 
#                                      G  => $genus,
#                                      F  => $family,
#                                      O  => $order,
#                                      C  => $class,
#                                      P  => $phylum,
#                                      SK => $superkingdom
#                                      S  => { $species1 => $details, ...
#                                      SS => { $strain1 => $details, ...};
#                                     }
#                           ...
#                       );
# where $species1 is the TAXID for the species
#-------------------------------------------------------------------------------
sub loadTaxTree {

    $_[0]->{TAXTREE} = loadFromDisk($_[0]->{TAXFILE}, "TAX TREE");

    return;
}
################################################################################
# Reorganize $_[0]->{TAXTREE} into hash: 
#   KEY => VALUE
#   $gi => $repliconType [CHR or PLASMID]
#
# ARGS: $taxOptions
################################################################################
sub parseReplTypes {

    my %replTypeLookupByGI = ();
    $_[0]->{REPLTYPELOOKUP} = \%replTypeLookupByGI;


    NODE: foreach my $node (keys %{ $_[0]->{TAXTREE} }) {
        next NODE unless (exists $_[0]->{TAXTREE}->{$node}->{GI});
        foreach my $replType (keys %{ $_[0]->{TAXTREE}->{$node}->{GI} }) {           # CHR or PLASMID
            foreach my $source (keys %{ $_[0]->{TAXTREE}->{$node}->{GI}->{$replType} }) {
                foreach my $gi (keys %{ $_[0]->{TAXTREE}->{$node}->{GI}->{$replType}->{$source} }) {
                    $replTypeLookupByGI{$gi} = $replType;
                } #GI
            } #SOURCE
        } #REPLTYPE
    } #NODE
    
    return;
}
################################################################################
# ARGS: \%taxOptions
################################################################################
sub loadVitals {

    $_[0]->{VITALS} = loadFromDisk($_[0]->{VITALSFILE},"GENOME VITALS");

    return;
}
################################################################################
# Function: Used to recursively traverse a HoH by returning a reference to the
#           deeper hash.
# ARGS: \%hash, $key1
################################################################################
sub getHref {
    
    $_[0]->{$_[1]} = {} unless (exists $_[0]->{$_[1]});
    return $_[0]->{$_[1]};
    
}
################################################################################
# ARGS: $filename, $description
################################################################################
sub loadFromDisk {

    die "Attempted to load non-existent file \"".$_[0]."\" from disk. Abort\n" if(!-e $_[0]);

    print "Loading ".$_[1]." from disk [".$_[0]."]...";
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
        print "Storing ".$_[2]." to disk in BINARY format as \"".$_[1]."\"...";
        store $_[0], $_[1];
    }
    elsif($_[3] =~ m/^ascii/i) {
        print "Storing ".$_[2]." to disk in ASCII (text) format as \"".$_[1]."\"...";
        open my $OUTFILE, '>', $_[1];
        print $OUTFILE YAML::XS::Dump($_[0]);
        close $OUTFILE;
    }
    else {
        print "Attempting to store hash \"".$_[2]."\" to \"".$_[1]."\" in an unknown format: skipping...\n";
        return;
    }

    my $iter1 = new Benchmark;
    print "done. [".timestr(timediff($iter1,$iter0))."]\n";
    return;
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
# Given:   filename
# Returns: filehandle
#-------------------------------------------------------------------------------
sub init_file_write {
  my $fname = shift;
  open my $fh, '>', $fname || die "Fatal: cannot open \'$fname\' for writing.\n";
  return \*$fh;
}
################################################################################
#-------------------------------------------------------------------------------
sub usage {
    print <<END;

  This program creates a Taxonomic Index XML file based on a list of organisms whose 
  names are the subdirectories itself, and within these subdirs reside their corresponding 
  FASTA files.

  Usage: $0 --drafts=<FILE> --completed=<FILE> --taxTree=<FILE> --vitals=<FILE> [OPTIONS]

  [REQUIRED]
    *AT LEAST ONE of* --drafts or --completed must be specified.
    Each specifies a file containing fully-qualified directory names -- separated by valid
    kingdom name headers and listed one per line -- of the organisms intended for match against
    the tax file. 
               
    Kingdom headers are enclosed within square brackets, reside on their own line, 
    and are in no particular order. Valid kingdoms: [VIRUSES], [PROKARYOTES], 
    [EUKARYOTES], and [OTHERS]. Kingdom [OTHERS] is for nucleic acid sequences from 
    an as-yet undetermined source or for situations where the kingdom distinction is 
    undesireable.
        
    --drafts            File of directory names containing draft genome sequences
    --completed         File of directory names containing completed genome sequences
                        Example contents of the --completed file:
                            [VIRUSES]
                            /db/genomes/Viruses/White_clover_mosaic_virus_uid15069
                            /db/genomes/Viruses/Whitewater_Arroyo_virus_uid29833
                            [PROKARYOTES]
                            /db/genbank/genomes/Bacillus_anthracis_Ames_uid309
                            /db/genbank/genomes/Bacillus_anthracis_CDC_684_uid31329
                            /db/genbank/genomes/Bacillus_anthracis_CI_uid36309
                            [EUKARYOTES]
                            [OTHERS]
                        **Note: There is a bug in parsing the directory format. If you receive
                        "Unrecognized directory format", then remove any trailing "/" in the dir name. 
    
    --taxTree           Species-level taxonomy tree file in Perl's Storable format. 
                        Usually "speciesTreeGI.dmp".

    --vitals            Genome vitals file in Perl's Storable format. 
                        Usually "genomeVitals.dmp".

  [OPTIONS]
    --prefix            Name of the output XML file. [$XML_OUTFILE]
        
    --skipOrgNoGI       If specified, the program will proceed, but will simply exclude 
                        (and report) all organism's whose FASTA files do not contain (valid) 
                        GIs. [abort]
    --continueNoMatch   If specified, the program will proceed, but will simply exclude
                        (and report) all organism's whose GI has no associated TAXID. [abort]
    --noDenovoFNA       If *.fna files missing in an organism's subdir, do NOT create *.fna 
                        files from corresponding *.gbk files. [0]
        
END
exit;
}

