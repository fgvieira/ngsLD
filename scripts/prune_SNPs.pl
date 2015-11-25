#!/usr/bin/perl -w
#
# Cared for by Filipe G. Vieira <>
#
# Copyright Filipe G. Vieira
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

    prune_SNPs.pl v1.0.5

=head1 SYNOPSIS

    perl prune_SNPs.pl [--in_file /path/to/input/file] [--subset /path/to/subset/file] [--max_dist +Inf] [--min_LD 0.7] [--field 4] [--method 0]

    --in_file       = File with input network (default STDIN)
    --subset        = File with node IDs to include (one per line)
    --max_kb_dist   = Maximum distance between nodes (input file 3rd column) to assume they are connected (in Kb)
    --min_LD        = Minimum level of LD to consider to SNPs as linked (as in --field option)
    --field         = Which column from input file to use
    --method        = Prunning method to use

=head1 DESCRIPTION

    This script will prune SNPs to LD and distance between them, ending up in selecting only representative SNPs. Input file should be tab-sepparated into: node1, node2, dist, LD.

=head1 AUTHOR

    Filipe G. Vieira - fgarrettvieira _at_ gmail _dot_ com

=head1 CONTRIBUTORS

    Additional contributors names and emails here

=cut


# Let the code begin...



use strict;
use Getopt::Long;
use Math::BigInt;
use IO::Zlib;

my ($in_file, $subset_file, $max_kb_dist, $min_LD, $field, $method, $debug);
my ($cnt, %subset);
my (%keep, %prunned);   # Only used on first method
my ($prev_id, $search); # Only used on second method

$in_file = '-';
$max_kb_dist = '+inf';
$min_LD = 0.7;
$field = 4;
$method = 0;
$debug = 0;
$search = 0;

# Parse args
GetOptions('h|help'             => sub { exec('perldoc',$0); exit; },
           'i|in_file=s'        => \$in_file,
	   's|subset:s'         => \$subset_file,
	   'd|max_kb_dist:s'    => \$max_kb_dist,
           'r|min_LD:s'         => \$min_LD,
           'f|field:s'          => \$field,
	   'method:i'           => \$method,
	   'debug!'             => \$debug,
    );


print(STDERR "### Assuming input first two columns are sorted by genomic coordinates and that there was no random sampling!\n");

print(STDERR "WARNING: \"--method 1\" is not guaranteed to work, e.g. if an unlinked site has no SNPs less than \"--max_kb_dist\" away the remaining comparisons are ignored.\n") if($method == 1);


# Parse subset file
if($subset_file){
    open(SUBSET, $subset_file) || die("ERROR: cannot open subset file!");
    %subset = map {chomp($_); $_ => 1} <SUBSET>;
    close(SUBSET);
}


# Open file
my $FILE = new IO::Zlib;
$FILE->open($in_file, "r") || die("ERROR: cannot open input file!");


# Loop through all comparisons
while(<$FILE>){
    $cnt++;
    my @interact = split(m/\t/, $_);
    chomp();

    # Skip if NaN, +Inf, -Inf, ...
    my $x = Math::BigInt->new($interact[$field-1]);
    next if($x->is_nan() || $x->is_inf('+') || $x->is_inf('-'));

    # Skip if SNP is not in the subset file
    next if($subset_file && !defined($subset{$interact[0]}));

    ## If first comparison, print first SNP and fill in $prev_id variable
    if(!defined($prev_id)){
	print($interact[0]."\n");
	$keep{$interact[0]}++;
	$prev_id = $interact[0];
    }

    if($method == 0){
	# Algorithm that "marks" SNPs as linked and excludes them from future comparisons
	if( !defined($keep{$interact[0]}) && !defined($prunned{$interact[0]}) ){
	    $keep{$interact[0]}++;
	    print($interact[0]."\n");
	}

	if( defined($keep{$interact[0]}) ){
	    if($interact[2] <= $max_kb_dist*1000 && $interact[$field-1] > $min_LD){
		$prunned{$interact[1]}++;
	    }
	}
    }elsif($method == 1){
	# Algorithm that "jumps" to next unlinked SNP
	next if($search && $interact[0] ne $prev_id);
	$search = 0;

	if($interact[0] eq $prev_id){
	    if($interact[2] >= $max_kb_dist*1000 || $interact[$field-1] < $min_LD){
		print($interact[1]."\n");
		$prev_id = $interact[1];
		$search = 1;
	    }
	}else{
	    print($interact[0]."\n");
	    $prev_id = $interact[0];
	}
    }else{
	die("ERROR: wrong method specified!");
    }
}


$FILE->close;
exit(0);
