#! /usr/bin/env perl
use strict;
use Getopt::Long;

my ($abu_col,$sort_col) = ("LINEAR_DOC", "LINEAR_DOC");
my $sig_lvl = "replicon";

GetOptions(
    'sig_lvl=s',  \$sig_lvl,
    'abu_col=s',  \$abu_col,
    'sort_col=s', \$sort_col,
);

$sig_lvl = lc $sig_lvl;
$abu_col = uc $abu_col;
$sort_col = uc $sort_col;

$|=1;

my %major_level = (
	'superkingdom' => 10,
	'phylum'       => 20,
	'class'        => 30,
	'order'        => 40,
	'family'       => 50,
	'genus'        => 60,
	'species'      => 70,
	'strain'       => 80,
	'replicon'     => 90
);

my $taxa;

my @files = @ARGV;
my @headers;
my $rank;
my $abu_sum;
foreach my $fileName ( @files ) {
        my $index;
	open TAXA, $fileName || die "Can't open $fileName\n";
	while(<TAXA>)
	{
		chomp;
		next if /^#/;
		next if /^$/;
	
		my @fields = split /\t/, $_;
		next if scalar @fields < 3;
		if ( $_ =~ /\tTOTAL_BP_MAPPED\t/ ){
			@headers = @fields[1..$#fields];
		        ( $index ) = grep { $headers[$_] eq $abu_col } 0..$#headers;
                        die "Can't find $abu_col.\n" unless defined $index;
			$rank = lc $fields[0];
			next;
		}

		my $name = $fields[0];
		$abu_sum->{$rank}                   += $fields[$index+1];
		$taxa->{$rank}->{$name}->{ABU}       = $fields[$index+1];
		@{$taxa->{$rank}->{$name}->{RESULT}} = @fields[1..$#fields];
	}
	close TAXA;
}

my $header = join "\t", @headers;
print "LEVEL\tTAXA\tREL_ABUNDANCE\t$header\n";
#print "LEVEL\tTAXA\t$header\n";
foreach my $rank ( sort {$major_level{$a}<=>$major_level{$b}} keys %$taxa )
{
	last if $major_level{$rank} > $major_level{$sig_lvl};
	
	my ( $index ) = grep { $headers[$_] eq $sort_col } 0..$#headers;
        die "Can't find $sort_col.\n" unless defined $index;
	
	foreach my $name ( sort { $taxa->{$rank}->{$b}->{RESULT}[$index] <=> $taxa->{$rank}->{$a}->{RESULT}[$index] } keys %{$taxa->{$rank}}  ) {
		next unless $abu_sum->{$rank};
		printf "%s\t%s\t%.4f\t%s\n",
			$rank,
			$name,
			$taxa->{$rank}->{$name}->{ABU}/$abu_sum->{$rank},
			join "\t", @{$taxa->{$rank}->{$name}->{RESULT}};
	}
}

print STDERR "Done converting ".(scalar @files)." files to list.\n"
