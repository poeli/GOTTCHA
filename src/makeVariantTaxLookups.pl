#!/usr/bin/perl 
#===============================================================================
#
#         FILE: makeSpeciesLookup.pl
#
#        USAGE: ./makeSpeciesLookup.pl  
#
#  DESCRIPTION: Reads in speciesTreeGI.dmp file, and creates three hash lookups:
#               (1) %variantStrainLookup  => { $strain_variant_name => $species_SciName }
#               (2) %variantSpeciesLookup => { $species_variant_name => $species_SciName }
#               (3) %taxLookupBySpecies   => { $species_SciName      => { G => $genus_SciName,
#                                                                         F => $family_SciName,
#                                                                         O => $order_SciName,
#                                                                         C => $class_SciName,
#                                                                         P => $phylum_SciName
#                                                                       }
#                                            }
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: YOUR NAME (), 
# ORGANIZATION: 
#      VERSION: 1.0
#      CREATED: 06/26/2013 04:18:33 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use Storable;
use YAML::XS;
use IO::Handle;
use Benchmark;

die "Please specify the location of the \"speciesTreeGI.dmp\" file\n" if(!$ARGV[0]);
die "The file does not exist!\n" if(!-e $ARGV[0]);

my %speciesTree = %{ retrieve $ARGV[0] };
my %variantStrainLookup  = ();
my %variantSpeciesLookup = ();
my %taxLookupBySpecies   = ();
my %genusLookupBySpecies = ();
my %familyLookupByGenus  = ();
my %orderLookupByFamily  = ();
my %classLookupByOrder   = ();
my %phylumLookupByClass  = ();

reshuffleTree(\%speciesTree, \%variantStrainLookup, \%variantSpeciesLookup, \%taxLookupBySpecies, 
              \%genusLookupBySpecies, \%familyLookupByGenus, \%orderLookupByFamily, \%classLookupByOrder, \%phylumLookupByClass);
dump2disk(\%variantStrainLookup, "variantStrainLookup.dmp", "VARIANT STRAIN LOOKUP", "bin");
dump2disk(\%variantStrainLookup, "variantStrainLookup.readable.dmp", "VARIANT STRAIN LOOKUP", "ascii");
dump2disk(\%variantSpeciesLookup, "variantSpeciesLookup.dmp", "VARIANT SPECIES LOOKUP", "bin");
dump2disk(\%variantSpeciesLookup, "variantSpeciesLookup.readable.dmp", "VARIANT SPECIES LOOKUP", "ascii");
dump2disk(\%taxLookupBySpecies, "taxLookupBySpecies.dmp", "TAX LOOKUP BY SPECIES", "bin");
dump2disk(\%taxLookupBySpecies, "taxLookupBySpecies.readable.dmp", "TAX LOOKUP BY SPECIES", "ascii");

dump2disk(\%genusLookupBySpecies, "genusLookupBySpecies.dmp", "GENUS LOOKUP BY SPECIES", "bin");
dump2disk(\%genusLookupBySpecies, "genusLookupBySpecies.readable.dmp", "GENUS LOOKUP BY SPECIES", "ascii");
dump2disk(\%familyLookupByGenus, "familyLookupByGenus.dmp", "FAMILY LOOKUP BY GENUS", "bin");
dump2disk(\%familyLookupByGenus, "familyLookupByGenus.readable.dmp", "FAMILY LOOKUP BY GENUS", "ascii");
dump2disk(\%orderLookupByFamily, "orderLookupByFamily.dmp", "ORDER LOOKUP BY FAMILY", "bin");
dump2disk(\%orderLookupByFamily, "orderLookupByFamily.readable.dmp", "ORDER LOOKUP BY FAMILY", "ascii");
dump2disk(\%classLookupByOrder, "classLookupByOrder.dmp", "CLASS LOOKUP BY ORDER", "bin");
dump2disk(\%classLookupByOrder, "classLookupByOrder.readable.dmp", "CLASS LOOKUP BY ORDER", "ascii");
dump2disk(\%phylumLookupByClass, "phylumLookupByClass.dmp", "PHYLUM LOOKUP BY CLASS", "bin");
dump2disk(\%phylumLookupByClass, "phylumLookupByClass.readable.dmp", "PHYLUM LOOKUP BY CLASS", "ascii");

################################################################################
#           $_[0]             $_[1]                    $_[2]                $_[3]               $_[4]                       $_[5]                   $_[6]               $_[7]                   $_[8]
# ARGS: \%speciesTree, \%variantStrainLookup, \%variantSpeciesLookup, \%taxLookupBySpecies, \%genusLookupBySpecies, \%familyLookupByGenus, \%orderLookupByFamily, \%classLookupByOrder, \%phylumLookupByClass
sub reshuffleTree {

    foreach my $node (keys %{ $_[0] }) {
        
        # Get STRAINs
        my @strainNames = keys %{ $_[0]->{$node}->{SS} };

        # Get SPECIES SciName
        my $speciesSciNameIdx = -1;
        my @speciesNames = keys %{ $_[0]->{$node}->{S} };
        my @speciesNameTypes = values %{ $_[0]->{$node}->{S} };
        IDX: for (0..$#speciesNameTypes) {
            if($speciesNameTypes[$_] eq "scientific name") {
                $speciesSciNameIdx = $_;
                last IDX;
            }
        }
        my $speciesSciName = ($speciesSciNameIdx == -1) 
                           ? ($speciesNames[0])                     # Default: return first one if no SciName found
                           : ($speciesNames[$speciesSciNameIdx]);
        print "No Scientific Name found for SPECIES \"".$speciesSciName."\"!\n" if($speciesSciNameIdx == -1);

        # Populate Variant STRAIN Tree
        $_[1]->{$_} = $speciesSciName foreach (@strainNames);
        
        # Populate Variant SPECIES Tree and TaxLookupBySpecies Tree
        foreach my $species (@speciesNames) {
            $_[2]->{$species} = $speciesSciName;
            $_[3]->{$speciesSciName}->{G} = $_[0]->{$node}->{G};
            $_[3]->{$speciesSciName}->{F} = $_[0]->{$node}->{F};
            $_[3]->{$speciesSciName}->{O} = $_[0]->{$node}->{O};
            $_[3]->{$speciesSciName}->{C} = $_[0]->{$node}->{C};
            $_[3]->{$speciesSciName}->{P} = $_[0]->{$node}->{P};
            
            $_[4]->{$species}               = $_[0]->{$node}->{G};        # genusLookup
            $_[5]->{ $_[0]->{$node}->{G} }  = $_[0]->{$node}->{F};        # familyLookup
            $_[6]->{ $_[0]->{$node}->{F} }  = $_[0]->{$node}->{O};        # orderLookup
            $_[7]->{ $_[0]->{$node}->{O} }  = $_[0]->{$node}->{C};        # classLookup
            $_[8]->{ $_[0]->{$node}->{C} }  = $_[0]->{$node}->{P};        # phylumLookup
        }

    }

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
        store $_[0], $_[1];
    }
    elsif($_[3] =~ m/^ascii/i) {
        print "->Storing datastructure ".$_[2]." to disk in ASCII (text) format as \"".$_[1]."\"...";
        open my $OUTFILE, '>', $_[1];
        print $OUTFILE Dump($_[0]);
        close $OUTFILE;
    }
    else {
        print "->Attempting to store datastructure \"".$_[2]."\" to \"".$_[1]."\" in an unknown format: skipping...\n";
        return;
    }

    my $time1 = new Benchmark;
    my $timeTxt = timestr(timediff($time1,$time0));
       $timeTxt = $1 if($timeTxt =~ m/^\s*(\d+)\s+/);
    
    print "done. ".$timeTxt." wallsecs\n";
    return;
}





