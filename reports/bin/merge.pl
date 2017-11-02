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
use Text::ParseWords;
use Getopt::Long;

### command-line options
my $o_mode = "";

GetOptions( "mode=s" => \$o_mode );

### command-line usage
if( @ARGV == 0 )
{
    print "usage: $0 [options] path/to/summary.csv.*\n\n";
    exit;
}

### command-line arguments
my $nfiles = scalar @ARGV;

### accumulate data here
my %data = ();
my @head = ();

for my $file ( @ARGV )
{
    ### load data
    my ($data,$head) = read_csv($file);

    ### save for future use
    if( @head == 0 ) { @head = @$head; }

    ### split data
    if( $o_mode eq "fse" )
    {
	$data = [ @{$data}[0] ];
    }
    elsif( $o_mode eq "sse" )
    {
	$data = [ @{$data}[1] ];
    }

    ### accumulate data
    for my $entry ( @$data )
    {
	for my $col ( keys %$entry )
	{
	    if( index($col,"Spectrum") != -1 )
	    {
		### combine spectrum
		for my $pair ( split(/;/,$$entry{$col}) )
		{
		    my ($type,$count) = split(/=/,$pair);
		    
		    $data{$col}{$type} += $count;
		}
	    }
	    else
	    {
		### or, just sum up values
		$data{$col} += $$entry{$col};
	    }
	}
    }
}

### print CSV headers
print join(",",@head), "\n";

### CSV values go here
my @values = ();

for my $col ( @head )
{
    if( index($col,"Spectrum") != -1 )
    {
	push @values, spectrum($data{$col});
    }
    else
    {
	if( exists $data{$col} )
	{
	    if( $col eq "NP" || $col eq "QV" )
	    {
		push @values, $data{$col} / $nfiles;
	    }
	    else
	    {
		push @values, $data{$col};
	    }
	}
	else
	{
	    push @values, 0;
	}
    }
}

### print CSV values
print join(",",@values), "\n";

sub spectrum {
    my ($data) = @_;

    if( defined $data )
    {
	my @sorted = sort { $$data{$b} <=> $$data{$a} } keys %$data;
	my @spectrum = map { sprintf("%s=%i",$_,$$data{$_}) } @sorted;

	return join(";",@spectrum);
    }
    else
    {
	return "";
    }
}

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
