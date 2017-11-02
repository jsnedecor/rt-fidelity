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

### command-line options
my $o_qv  = 0;
my $o_np  = 1;
my $o_mapq = 254;
my $o_rla = 0;
my $o_rlp = 0.0;
my $o_rlen_min_a = -1;
my $o_rlen_max_a = -1;
my $o_rlen_min_f = -1;
my $o_rlen_max_f = -1;
my $o_acov_min_a = -1;
my $o_acov_max_a = -1;
my $o_acov_min_f = -1;
my $o_acov_max_f = -1;
my $o_zm  = -1;
my $o_sna = 1e6;
my $o_snc = 1e6;
my $o_sng = 1e6;
my $o_snt = 20.0;
my $o_lb  = -1;
my $o_ub  = -1;
my $o_log = "";
my @o_strand = ();
my $o_range = "";
my $o_aln_check = "";
my $o_mutation_dir = "mutation";
my $o_ctxfile = "";

GetOptions(
    "qv=f"  => \$o_qv,
    "np=f"  => \$o_np,
    "mapq=f"  => \$o_mapq,
    "rla=f" => \$o_rla,
    "rlp=f" => \$o_rlp,
    "rlen-min-a=f" => \$o_rlen_min_a,
    "rlen-max-a=f" => \$o_rlen_max_a,
    "rlen-min-f=f" => \$o_rlen_min_f,
    "rlen-max-f=f" => \$o_rlen_max_f,
    "acov-min-a=f" => \$o_acov_min_a,
    "acov-max-a=f" => \$o_acov_max_a,
    "acov-min-f=f" => \$o_acov_min_f,
    "acov-max-f=f" => \$o_acov_max_f,
    "zm=f"  => \$o_zm,
    "sna=f" => \$o_sna,
    "snc=f" => \$o_snc,
    "sng=f" => \$o_sng,
    "snt=f" => \$o_snt,
    "lb=f"  => \$o_lb,
    "ub=f"  => \$o_ub,
    "log"   => \$o_log,
    "strand=s" => \@o_strand,
    "range=s"  => \$o_range,
    "aln-check=s" => \$o_aln_check,
    "mutation-dir=s" => \$o_mutation_dir,
    "ctxfile=s" => \$o_ctxfile,
    );

### post-process command line options
@o_strand = split(/,/,join(",",@o_strand));
if( ! @o_strand ) { @o_strand = qw/fwd rev/; }

### command-line arguments
if( @ARGV == 0 )
{
    print "usage: $0 [options] reference.fasta sampleDir\n";
    print "\n";
    print "options:\n";
    print " --qv\t\tminimum quality value ($o_qv)\n";
    print " --np\t\tminimum number of passes ($o_np)\n";
    print " --rlp\t\tminimum read length ($o_rlp% of expected length)\n";
    print " --rla\t\tminimum read length ($o_rla nucleotides)\n";
    print " --zm\t\tanalyse a specific ZMW ($o_zm)\n";
    print " --sna\t\tmaximum SNR-A ($o_sna)\n";
    print " --snc\t\tmaximum SNR-C ($o_snc)\n";
    print " --sng\t\tmaximum SNR-G ($o_sng)\n";
    print " --snt\t\tmaximum SNR-T ($o_snt)\n";
    print " --lb\t\tskip first N bases ($o_lb)\n";
    print " --ub\t\tskip last N bases ($o_ub)\n";
    print " --log\t\tprint log info ($o_log)\n";
    print " --strand\tsummarize data for 'fwd', 'rev' or both strands (@o_strand)\n";
    print "\n";
    exit;
}

my $reference = shift @ARGV;
my $sampleDir = shift @ARGV;

### load reference data
my $refs = load_references($reference);

### accumulate data for fwd and rev strands
my %data = ();

for my $strand ( @o_strand )
{
    my $variants_file = sprintf( "%s/%s/%s/variants.csv.7z",       $sampleDir, $o_mutation_dir, $strand );
    my $zmws_file     = sprintf( "%s/%s/%s/zmws.csv.7z",           $sampleDir, $o_mutation_dir, $strand );
    my $alns_file     = sprintf( "%s/%s/%s/aln.csv.7z",            $sampleDir, $o_mutation_dir, $strand );
    my $chimeric_file = sprintf( "%s/chimeric/%s/chimeric.csv.7z", $sampleDir, $strand );
    
    if( $o_log )
    {
	print STDERR $variants_file, "\n";
	print STDERR $chimeric_file, "\n";
	print STDERR $zmws_file, "\n";
	print STDERR $alns_file, "\n";
    }

    ### load pacbio read stats
    my $zmws = load_zmw_data($zmws_file);
    
    ### load alignment stats
    my $alns = load_zmw_data($alns_file);

    ### load chimeric reads
    my $blacklist = load_blacklist($chimeric_file);
    
    ### filter variants
    $data{$strand} = filter($variants_file, $zmws, $alns, $blacklist, $strand);
}

sse(\%data);

sub filter {
    my ($variants_file,$zmws,$alns,$blacklist,$strand) = @_;

    my %data = ();
    my %warn = ();
    
    my @head = ();
    
    my ($b,$e) = ( $o_range ne "" ? split(/-/,$o_range) : (-1e9,+1e9) );

    my %indel = ();
    my $prev_qname = "";

    open(VARIANTS,"7za e -so $variants_file |") || die "Can't open '$variants_file': $!";
    
    while( my $line = <VARIANTS> )
    {
	chomp($line);
	
	### parse line
	my @tokens = split(/,/,$line);

	my %entry = ();
	
	if( @head == 0 )
	{
	    @head = @tokens;
	    next;
	}
	else
	{
	    for( my $i = 0; $i < @tokens; $i++ )
	    {
		$entry{$head[$i]} = $tokens[$i];
	    }
	}
	
	if( ! exists $$refs{$entry{"Ref"}} )
	{
	    print STDERR "error :: missing reference '$entry{'Ref'}'\n\n";
	    exit;
	}
	
	### filter
	my $refLength = length($$refs{$entry{"Ref"}});
	
	### frequently used vars
	my $movie = $entry{"Movie"};
	my $zmw = $entry{"ZMW"};
	my $pos = $entry{"Pos"};
	my $type = $entry{"Type"};
	my $qname = "$movie/$zmw";

	### track indel info (before any filter applied)
	if( $o_aln_check ne "" )
	{
	    ### keep indel positions
	    if( $type eq "INDEL" )
	    {
		if( $entry{"IndelType"} eq "Deletion" )
		{
		    ### store deletion anchor for convenience
		    $indel{$pos}{$entry{"IndelType"}."Start"} = 1;

		    for( my $i = 1; $i <= $entry{"Length"}; $i++ )
		    {
			$indel{$pos+$i}{$entry{"IndelType"}} = 1;
		    }
		}
		elsif( $entry{"IndelType"} eq "Insertion" )
		{
		    ### store insertion anchor for convenience
		    $indel{$pos}{$entry{"IndelType"}."Start"} = 1;

		    $indel{$pos}{$entry{"IndelType"}} = $entry{"Length"};
		}
	    }

	    ### keep track of individual reads
	    if( $prev_qname eq "" )
	    {
		$prev_qname = $qname;
	    }
	    
	    ### start of a new read
	    if( $qname ne $prev_qname || eof(VARIANTS) )
	    {
		# print "end of '$prev_qname'\n";
		
		### inject data for a previous read
		inject_aln_data(\%data,\%indel,$prev_qname,$refLength);
		
		## reset
		%indel = ();
		$prev_qname = $qname;
	    }
	}
	
	### filter by ZMW range
	next if( $zmw < $b || $zmw > $e );
	
	### filter by number of passes
	next if( $$zmws{$movie}{$zmw}{"NP"} < $o_np );
	
	### filter by mapping quality
	next if( $o_mapq != -1 && $$alns{$movie}{$zmw}{"MAPQ"} < $o_mapq );

	### filter by read length
	next if( $o_rlen_min_a != -1 && $$zmws{$movie}{$zmw}{"ReadLength"} < $refLength - $o_rlen_min_a );
	next if( $o_rlen_max_a != -1 && $$zmws{$movie}{$zmw}{"ReadLength"} > $refLength + $o_rlen_max_a );
	next if( $o_rlen_min_f != -1 && $$zmws{$movie}{$zmw}{"ReadLength"} < $refLength * $o_rlen_min_f );
	next if( $o_rlen_max_f != -1 && $$zmws{$movie}{$zmw}{"ReadLength"} > $refLength * $o_rlen_max_f );

	### filter by alignment coverage
	next if( $o_acov_min_a != -1 && $$alns{$movie}{$zmw}{"AlnCoverage"} < $refLength - $o_acov_min_a );
	next if( $o_acov_max_a != -1 && $$alns{$movie}{$zmw}{"AlnCoverage"} > $refLength + $o_acov_max_a );
	next if( $o_acov_min_f != -1 && $$alns{$movie}{$zmw}{"AlnCoverage"} < $refLength * $o_acov_min_f );
	next if( $o_acov_max_f != -1 && $$alns{$movie}{$zmw}{"AlnCoverage"} > $refLength * $o_acov_max_f );
	
	### make sure read covers the entire region of interest (in other words, avoid truncated reads)
	next if( $o_lb != -1 && $$alns{$movie}{$zmw}{"AlnStart"} > $o_lb );
	next if( $o_ub != -1 && $$alns{$movie}{$zmw}{"AlnEnd"} < $refLength - $o_ub );

	### filter by SNR
	next if( $$zmws{$movie}{$zmw}{"SnrA"} > $o_sna );
	next if( $$zmws{$movie}{$zmw}{"SnrC"} > $o_snc );
	next if( $$zmws{$movie}{$zmw}{"SnrG"} > $o_sng );
	next if( $$zmws{$movie}{$zmw}{"SnrT"} > $o_snt );
	
	### skip bases at the begining and at the end of reference (primer sites)
	next if( $o_lb != -1 && $entry{"Pos"} < $o_lb );
	next if( $o_ub != -1 && $entry{"Pos"} > $refLength - $o_ub );
	
	### exclude blacklisted reads
	next if( exists $$blacklist{$movie}{$zmw} );
	
	### ensure correct strand directionality
	next if( $strand eq "fwd" && $$alns{$movie}{$zmw}{"Flag"} != 0 );
	next if( $strand eq "rev" && $$alns{$movie}{$zmw}{"Flag"} != 16 );

	### store substitutions, deletions, insertions
	$data{$qname}{$pos}{$type} = \%entry;

	$warn{$qname}{$pos}{$type}++;
	$warn{$qname}{$pos}{$entry{"IndelType"}}++;

	### print out data lines (optional)
	if( $o_log )
	{
	    print STDERR $line, "\n";
	}
    }
    
    close(VARIANTS);

    check_warnings(\%warn);

    return \%data;
}

sub load_references {
    my ($file) = @_;
    
    my %refs = ();
    
    my $name = "";
    my $seq = "";
    
    open(FA,$file) || die "Can't open '$file': $!";
    
    while( my $line = <FA> )
    {
	chomp($line);
	
	if( substr($line,0,1) eq ">" )
	{
	    if( $name ne "" )
	    {
		$refs{$name} = uc($seq);
	    }
	    
	    $name = substr($line,1);
	    $seq = "";
	}
	else
	{
	    $seq .= $line;
	}
    }
    
    if( $name ne "" )
    {
	$refs{$name} = uc($seq);
    }
    
    close(FA);
    
    return \%refs;
}

sub load_blacklist {
    my ($file) = @_;

    my %data = ();

    open(SR,"7za e -so $file |") || die "Can't open '$file': $!";

    while( my $line = <SR> )
    {
	chomp($line);

	my ($movie,$zmw) = split(/,/,$line);

	$data{$movie}{$zmw} = 1;
    }

    close(SR);

    return \%data;
}

sub load_zmw_data {
    my ($file) = @_;
    
    my %zmws = ();
    my @head = ();

    open(IN,"7za x -so $file |") || die "Can't open '$file': $!";
    
    while( my $line = <IN> )
    {
	chomp($line);
	
	### parse line
	my @tokens = split(/,/,$line);
	
	if( @head == 0 )
	{
	    ### extract column names
	    @head = @tokens;
	}
	else
	{
	    ### store data fields
	    my %entry = ();
	    
	    for( my $i = 0; $i < @head; $i++ )
	    {
		$entry{$head[$i]} = $tokens[$i];
	    }

	    $zmws{$entry{"Movie"}}{$entry{"ZMW"}} = \%entry;
	}
    }

    close(IN);

    return \%zmws;
}

sub complement {
    my ($seq) = @_;

    my %bases = (
	"A" => "T",
	"C" => "G",
	"G" => "C",
	"T" => "A",
	"a" => "t",
	"c" => "g",
	"g" => "c",
	"t" => "a",
	"-" => "-",
	);

    my $complement = join("",map { $bases{$_} } split(//,$seq));

    return $complement;
}

sub check_warnings {
    my ($warn) = @_;
    
    for my $qname ( sort_reads(keys %$warn) )
    {
	for my $pos ( sort keys %{$$warn{$qname}} )
	{
	    for my $type ( keys %{$$warn{$qname}{$pos}} )
	    {
		if( $$warn{$qname}{$pos}{$type} > 1 )
		{
		    printf( STDERR "WARN :: multiple mutations assigned to the same position: %s,%i,%s,%i\n",
			    $qname,
			    $pos,
			    $type,
			    $$warn{$qname}{$pos}{$type} );
		}
	    }
	}
    }
}

sub sort_reads {
    my (@array) = @_;
    
    my @mov = map { substr($_,0,index($_,"/")) } @array;
    my @zmw = map { substr($_,index($_,"/")+1) } @array;

    my @sorted = @array[ sort { $mov[$a] cmp $mov[$b] || $zmw[$a] <=> $zmw[$b] } 0..$#array ];

    return @sorted;
}

sub sse {
    my ($data) = @_;

    ### list of CCS reads
    my @sorted = sort_reads(keys %{$$data{"fwd"}});

    my %context = ();

    ### ~~~~~ first-strand & second-strand error ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    my %fse_fine = ();
    my %fse_rare = ();
    my %fse_cplx = ();

    my %sse_fine = ();
    my %sse_rare = ();
    my %sse_cplx = ();

    for my $qname ( @sorted )
    {
	### both strands must be present
	next if( ! exists $$data{"rev"}{$qname} );
	
	for my $pos ( sort { $a <=> $b } keys %{$$data{"fwd"}{$qname}} )
	{
	    ### mutational data for FWD and REV strands
	    my $f = $$data{"fwd"}{$qname}{$pos};
	    my $s = $$data{"rev"}{$qname}{$pos};

	    ### both 'anchors' must be present
	    next if( ! exists $$f{"SNP"} || ! exists $$s{"SNP"} );
	    
	    ### alignment check (optional)
	    next if( $o_aln_check ne "" && $$f{"SNP"}{"Aln"} ne "pass" );
	    next if( $o_aln_check ne "" && $$s{"SNP"}{"Aln"} ne "pass" );
	    
	    ### determine reference base at the current position
	    my $rb = $$f{"SNP"}{"RefBP"};
	    my $s1 = $$f{"SNP"}{"AltBP"};
	    my $s2 = $$s{"SNP"}{"AltBP"};
	    my $q1 = $$f{"SNP"}{"QV"};
	    my $q2 = $$s{"SNP"}{"QV"};
            my $c1 = $$f{"SNP"}{"Context"};
            my $c2 = $$s{"SNP"}{"Context"};

	    my $movie = $$f{"SNP"}{"Movie"};
	    my $zmw   = $$f{"SNP"}{"ZMW"};
	    my $rname = $$f{"SNP"}{"Ref"};

	    ### ~~~~~ dealing with SNPs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	    
	    next if( $q1 < $o_qv );
	    next if( $q2 < $o_qv );

	    $context{"fse"}{substr($c1,6,3)}++;
	    $context{"sse"}{substr($c2,6,3)}++;

	    if( $s1 eq $s2 )
	    {
		### R     1     2   Status
		### ======================
		### A  >  A  >  A = OK
		### A  >  T  >  T = OK

		$sse_fine{$s1.$s2}++;
		#(do not log matches)

		$fse_fine{$rb.$s1}++;
		printf( STDERR "FSE fine %s,%i,%s,%i,%s,%s,%s,%s,%s,%s,%s,%s\n",$movie,$zmw,$rname,$pos,$rb,"SNP",$s1,"NA",$s2,"NA",$c1,$c2) if ( $rb ne $s1 );
	    }
	    else
	    {
		if( $rb eq $s1 )
		{
		    ### R     1     2   Status
		    ### ======================
		    ### A  >  A  >  T = OK

		    $sse_fine{$s1.$s2}++;
		    printf( STDERR "SSE fine %s,%i,%s,%i,%s,%s,%s,%s,%s,%s,%s,%s\n",$movie,$zmw,$rname,$pos,$rb,"SNP",$s1,"NA",$s2,"NA",$c1,$c2);

		    $fse_fine{$rb.$s1}++;
		    #(do not log matches)
		}
		else
		{
		    ### R     1     2   Status
		    ### ======================
		    ### A  >  T  >  A = RARE

		    $sse_rare{$s1.$s2}++;
		    printf( STDERR "SSE rare %s,%i,%s,%i,%s,%s,%s,%s,%s,%s,%s,%s\n",$movie,$zmw,$rname,$pos,$rb,"SNP",$s1,"NA",$s2,"NA",$c1,$c2);

		    $fse_rare{$rb.$s1}++;
		    printf( STDERR "FSE rare %s,%i,%s,%i,%s,%s,%s,%s,%s,%s,%s,%s\n",$movie,$zmw,$rname,$pos,$rb,"SNP",$s1,"NA",$s2,"NA",$c1,$c2);
		}
	    }

	    ### ~~~~~ dealing with INDELs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	    if( exists $$f{"INDEL"} && ! exists $$s{"INDEL"} )
	    {
		### R   1   2   Status
		### ==========================
		### . > D > . = rare (1>2 ins)
		### . > I > . = rare (1>2 del)

		my $i1 = $$f{"INDEL"}{"Bases"};
		my $t1 = $$f{"INDEL"}{"IndelType"};
		my $q1 = $$f{"INDEL"}{"QV"};
                my $c1 = $$f{"INDEL"}{"Context"};
		
		next if( $q1 < $o_qv );

		if( $t1 eq "Deletion" )
		{
		    ### count total number of inserted bases
		    $sse_rare{"Insertion"} += length($i1);
		    $fse_rare{"Deletion"}  += length($i1);

		    ### insertion spectrum
		    $sse_rare{"InsSpectrum"}{$i1}++;
		    $fse_rare{"DelSpectrum"}{$i1}++;

		    if( length($i1) == 1 )
		    {
			### count single-base insertions
			$sse_rare{"I1"}++;
			$fse_rare{"D1"}++;
		    }
		    else
		    {
			### count multi-base insertions
			$sse_rare{"IX"}++;
			$fse_rare{"DX"}++;
		    }
		}
		elsif( $t1 eq "Insertion" )
		{
		    ### count total number of deleted bases
		    $sse_rare{"Deletion"}  += length($i1);
		    $fse_rare{"Insertion"} += length($i1);

		    ### deletion spectrum
		    $sse_rare{"DelSpectrum"}{$i1}++;
		    $fse_rare{"InsSpectrum"}{$i1}++;

		    if( length($i1) == 1 )
		    {
			### count single-base deletions
			$sse_rare{"D1"}++;
			$fse_rare{"I1"}++;
		    }
		    else
		    {
			### count multi-base deletions
			$sse_rare{"DX"}++;
			$fse_rare{"IX"}++;
		    }
		}

		printf( STDERR "SSE rare %s,%i,%s,%i,%s,%s,%s,%s,%s,%s,%s,%s\n",$movie,$zmw,$rname,$pos,$rb,"INDEL",$t1,$i1,"NA","NA",$c1,"NA");
		printf( STDERR "FSE rare %s,%i,%s,%i,%s,%s,%s,%s,%s,%s,%s,%s\n",$movie,$zmw,$rname,$pos,$rb,"INDEL",$t1,$i1,"NA","NA",$c1,"NA");
	    }
	    elsif( ! exists $$f{"INDEL"} && exists $$s{"INDEL"} )
	    {
		### R   1   2   Status
		### ========================
		### . > . > D = ok (1>2 del)
		### . > . > I = ok (1>2 ins)
		
		my $i2 = $$s{"INDEL"}{"Bases"};
		my $t2 = $$s{"INDEL"}{"IndelType"};
		my $q2 = $$s{"INDEL"}{"QV"};
                my $c2 = $$s{"INDEL"}{"Context"};

		next if( $q2 < $o_qv );

		if( $t2 eq "Deletion" )
		{
		    ### ~~~~~ second-strand error ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		    ### count total number of deleted bases
		    $sse_fine{"Deletion"} += length($i2);

		    ### deletion spectrum
		    $sse_fine{"DelSpectrum"}{$i2}++;

		    if( length($i2) == 1 )
		    {
			### count single-base deletions
			$sse_fine{"D1"}++;
		    }
		    else
		    {
			### count multi-base deletions
			$sse_fine{"DX"}++;
		    }

		    ### ~~~~~ first-strand error ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		    
		    ### nothing to do
		}
		elsif( $t2 eq "Insertion" )
		{
		    ### ~~~~~ second-strand error ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		    ### count total number of inserted bases
		    $sse_fine{"Insertion"} += length($i2);

		    ### insertion spectrum
		    $sse_fine{"InsSpectrum"}{$i2}++;

		    if( length($i2) == 1 )
		    {
			### count single-base insertions
			$sse_fine{"I1"}++;
		    }
		    else
		    {
			### count multi-base insertions
			$sse_fine{"IX"}++;
		    }

		    ### ~~~~~ first-strand error ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		    
		    ### nothing to do
		}

		printf( STDERR "SSE fine %s,%i,%s,%i,%s,%s,%s,%s,%s,%s,%s,%s\n",$movie,$zmw,$rname,$pos,$rb,"INDEL","NA","NA",$t2,$i2,"NA",$c2);
	    }
	    elsif( exists $$f{"INDEL"} && exists $$s{"INDEL"} )
	    {
		my $i1 = $$f{"INDEL"}{"Bases"};
		my $t1 = $$f{"INDEL"}{"IndelType"};
		my $q1 = $$f{"INDEL"}{"QV"};
                my $c1 = $$f{"INDEL"}{"Context"};
		
		my $i2 = $$s{"INDEL"}{"Bases"};
		my $t2 = $$s{"INDEL"}{"IndelType"};
		my $q2 = $$s{"INDEL"}{"QV"};
                my $c2 = $$s{"INDEL"}{"Context"};

		next if( $q1 < $o_qv );
		next if( $q2 < $o_qv );

		if( $i1 ne $i2 || $t1 ne $t2 )
		{
		    ### R   1   2   Status
		    ### ==============================
		    ### . > D > D = complex (D1 <> D2)
		    ### . > I > I = complex (I1 <> I2)
		    ### or
		    ### . > I > D = complex (I <> D)
		    ### . > D > I = complex (D <> I)

		    ### ~~~~~ second-strand error ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		    $sse_cplx{"Complex"}++;
		    printf( STDERR "SSE complex %s,%i,%s,%i,%s,%s,%s,%s,%s,%s,%s,%s\n",$movie,$zmw,$rname,$pos,$rb,"INDEL",$t1,$i1,$t2,$i2,$c1,$c2);
		    
		    ### ~~~~~ first-strand error ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		    $fse_cplx{"Complex"}++;
		    printf( STDERR "FSE complex %s,%i,%s,%i,%s,%s,%s,%s,%s,%s,%s,%s\n",$movie,$zmw,$rname,$pos,$rb,"INDEL",$t1,$i1,$t2,$i2,$c1,$c2);
		}
		elsif( $t1 eq "Deletion" )
		{
		    ### R   1   2   Status
		    ### ===================================
		    ### . > D > D = ok (D1 == D2; 1>2 none)

		    ### ~~~~~ second-strand error ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		    
		    ### nothing to do

		    ### ~~~~~ first-strand error ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		    ### count total number of deleted bases
		    $fse_fine{"Deletion"} += length($i1);

		    ### deletion spectrum
		    $fse_fine{"DelSpectrum"}{$i1}++;

		    if( length($i1) == 1 )
		    {
			### count single-base deletions
			$fse_fine{"D1"}++;
		    }
		    else
		    {
			### count multi-base deletions
			$fse_fine{"DX"}++;
		    }

		    printf( STDERR "FSE fine %s,%i,%s,%i,%s,%s,%s,%s,%s,%s,%s,%s\n",$movie,$zmw,$rname,$pos,$rb,"INDEL",$t1,$i1,$t2,$i2,$c1,$c2);
		}
		elsif( $t1 eq "Insertion" )
		{
		    ### R   1   2   Status
		    ### ======================================
		    ### . > I > I = ok (I1 == I2; 1>2 correct)

		    ### ~~~~~ second-strand error ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		    ### we might have an insertion longer than 1 nt
		    my $len = length($i1);
		    
		    for( my $k = 0; $k < $len; $k++ )
		    {
			my $b1 = substr($i1,$k,1);
			my $b2 = substr($i2,$k,1);
			
			$sse_fine{$b1.$b2}++;
		    }

		    ### ~~~~~ first-strand error ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		    ### count total number of deleted bases
		    $fse_fine{"Insertion"} += length($i1);

		    ### deletion spectrum
		    $fse_fine{"InsSpectrum"}{$i1}++;

		    if( length($i1) == 1 )
		    {
			### count single-base deletions
			$fse_fine{"I1"}++;
		    }
		    else
		    {
			### count multi-base deletions
			$fse_fine{"IX"}++;
		    }

		    printf( STDERR "FSE fine %s,%i,%s,%i,%s,%s,%s,%s,%s,%s,%s,%s\n",$movie,$zmw,$rname,$pos,$rb,"INDEL",$t1,$i1,$t2,$i2,$c1,$c2);
		}
	    } ### indel
	} ### pos
    } ### qname

    my @type = qw(AA AC AT AG CA CC CT CG TA TC TT TG GA GC GT GG Deletion Insertion D1 DX I1 IX);
    
    ### print headers
    print join(",",
	       "NP",
	       "QV",
	       @type,
	       "DelSpectrum",
	       "InsSpectrum",
	       @type,
	       "DelSpectrum",
	       "InsSpectrum",
	       "Complex",
	), "\n";
    
    ### print FSE counts
    print join(",",
	       $o_np,
	       $o_qv,
	       (map { (exists $fse_fine{$_} ? $fse_fine{$_} : 0) } @type),
	       spectrum($fse_fine{"DelSpectrum"}),
	       spectrum($fse_fine{"InsSpectrum"}),
	       (map { (exists $fse_rare{$_} ? $fse_rare{$_} : 0) } @type),
	       spectrum($fse_rare{"DelSpectrum"}),
	       spectrum($fse_rare{"InsSpectrum"}),
	       (exists  $fse_cplx{"Complex"} ? $fse_cplx{"Complex"} : 0),
	), "\n";

    ### print SSE counts
    print join(",",
	       $o_np,
	       $o_qv,
	       (map { (exists $sse_fine{$_} ? $sse_fine{$_} : 0) } @type),
	       spectrum($sse_fine{"DelSpectrum"}),
	       spectrum($sse_fine{"InsSpectrum"}),
	       (map { (exists $sse_rare{$_} ? $sse_rare{$_} : 0) } @type),
	       spectrum($sse_rare{"DelSpectrum"}),
	       spectrum($sse_rare{"InsSpectrum"}),
	       (exists  $sse_cplx{"Complex"} ? $sse_cplx{"Complex"} : 0),
	), "\n";

    if( $o_ctxfile )
    {
	open(CTX,">",$o_ctxfile) || die "Can't open '$o_ctxfile'";

	print CTX join(",", "Context", "FSE", "SSE"), "\n";

	for my $m ( qw/A C G T/ )
	{
	    for my $l ( qw/A C G T/ )
	    {
		for my $r ( qw/A C G T/ )
		{
		    my @values = ( "$l$m$r" );

		    for my $type ( "fse" , "sse" )
		    {
			if( exists $context{$type}{"$l$m$r"} )
			{
			    push @values, $context{$type}{"$l$m$r"};
			}
			else
			{
			    push @values, 0;
			}
		    }

		    print CTX join(",",@values), "\n";
		}
	    }
	}

	close(CTX);
    }
}

sub inject_aln_data {
    my ($data,$indel,$qname,$reflen) = @_;
    
    ### cumulative number of indels
    my %csum = ();
    my %nsum = ();
    my $number_of_indels = 0;
    my $number_of_distinct_indels = 0;

    for( my $pos = 1; $pos <= $reflen; $pos++ )
    {
	if( exists $$indel{$pos} )
	{
	    for my $type ( keys %{$$indel{$pos}} )
	    {
		if( index($type,"Start") == -1 )
		{
		    $number_of_indels += $$indel{$pos}{$type};
		}
		else
		{
		    $number_of_distinct_indels++;
		}
	    }
	}
	
	$csum{$pos} = $number_of_indels;
	$nsum{$pos} = $number_of_distinct_indels;
    }

    ### window size and maximum number of indels
    my ($wsize,$limit) = split(/,/,$o_aln_check );

    ### insert alignment check data
    for my $pos ( keys %{$$data{$qname}} )
    {
	next if( ! exists $$data{$qname}{$pos}{"SNP"} );
	
	### default value
	my $status = "fail";

	### ensure that N bases present around a given position
	### 
	### NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
	###      <---------*--------->
	###     i-N        i        i+N
	### 
	
	if( exists $$data{$qname}{$pos-$wsize} && exists $$data{$qname}{$pos+$wsize} )
	{
	    my $number_of_indels = $csum{$pos+$wsize} - $csum{$pos-$wsize};
	    my $number_of_distinct_indels = $nsum{$pos+$wsize} - $nsum{$pos-$wsize};
	    
	    if( $number_of_distinct_indels <= $limit )
	    {
	    	$status = "pass";
	    }
	}
	
	$$data{$qname}{$pos}{"SNP"}{"Aln"} = $status;
    }
}

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
