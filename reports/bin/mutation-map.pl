#!/usr/bin/perl -w

################################################################################
# The effect of base modification on RNA polymerase and reverse transcriptase fidelity
# Copyright (C) 2017 New England Biolabs, Inc.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

use strict;
use Getopt::Long;
use Text::ParseWords;
use List::Util qw(min);

### command-line options
my $o_mode = "FSE";

GetOptions(
    "mode=s" => \$o_mode,
    );

### command-line usage
if( @ARGV == 0 )
{
    print "usage: [options] $0 samples.csv sampledir\n\n";
    print "options:\n";
    print " --mode\t FSE | SSE\n";
    print "\n";
    exit;
}

### command-line arguments
my $file = shift @ARGV;
my $sdir = shift @ARGV;

### load list of samples
my ($samples,$head) = read_csv($file);

my %data = ();
my %base = ();

for my $entry ( @$samples )
{
    my $logdir = sprintf("%s/%05i/summary", $sdir, $$entry{"SampleID"});
    
    ### find log files
    opendir(DIR,$logdir) || die "Can't open '$logdir'";
    my @files = grep(/summary.log./,readdir(DIR));
    closedir(DIR);

    ### accumulate stats
    for my $file ( @files )
    {
	open(IN,"$logdir/$file") || die "Can't open '$file'";

	while( my $line = <IN> )
	{
	    if( $line =~ m/^(FSE|SSE) (fine|rare|complex) (.*)$/ )
	    {
		my $strand = $1;
		my $class = $2;
		my $stats = $3;

		### parse stats
		my ($movie,$zmw,$rname,$pos,$rb,$type,$t1,$i1,$t2,$i2,$c1,$c2) = split(/,/,$stats);
		
		### accumulate by type / position / mutation
		if( $strand eq $o_mode && $class eq "fine" )
		{
		    my $pos = sprintf("%s:%i", $$entry{"Amplicon"}, $pos);
		    
		    if( $type eq "SNP" )
		    {
			$data{"Substitution"}{$pos}{$rb}++;
		    }
		    else
		    {
			$data{$t2}{$pos}{$i2}++;
		    }
		}
		
		$base{$pos} = $rb;
	    }
	}

	close(IN);
    }
}

### combine evrything into a single string
my @values = ();

for my $type ( "Substitution" , "Deletion" , "Insertion" )
{
    my $total = 0;
    my %counts = ( "1" => 0, "2" => 0, "3" => 0 );
    my %position = ();

    for my $pos ( keys %{$data{$type}} )
    {
	for my $mutation ( keys %{$data{$type}{$pos}} )
	{
	    $total += $data{$type}{$pos}{$mutation};
	    $counts{length($mutation)} += $data{$type}{$pos}{$mutation};
	    $position{$pos.":".$mutation} += $data{$type}{$pos}{$mutation};
	}
    }
    
    my @sorted = sort { $position{$b} <=> $position{$a} } keys %position;

    my @subset = ();

    for my $m ( @sorted[0..min(2,scalar(@sorted)-1)] )
    {
	if( $position{$m} > 1 )
	{
	    push @subset, $m;
	}
    }

    my @spectrum = map { sprintf("%s=%i",$_,$position{$_}) } @subset;
    
    push
	@values,
	$total,
	$counts{"1"},
	$counts{"2"},
	$counts{"3"},
	$total - $counts{"1"} - $counts{"2"} - $counts{"3"},
	join(" ",@spectrum);
}

print join(",", @values), "\n";

################################################################################
#                             SUBROUTINES                                      #
################################################################################

sub read_csv {
    my ($file) = @_;

    my @head = ();
    my @array = ();

    open(CSV,$file) || die "Can't open '$file'";

    while( my $line = <CSV> )
    {
	chomp($line);

	my @tokens = quotewords(',',0,$line);
	
	if( @head == 0 )
	{
	    @head = @tokens;
	}
	else
	{
	    my %entry = ();

	    for( my $i = 0; $i < @tokens; $i++ )
	    {
		$entry{$head[$i]} = $tokens[$i];
	    }

	    push @array, \%entry;
	}
    }

    close(CSV);

    return (\@array,\@head);
}
