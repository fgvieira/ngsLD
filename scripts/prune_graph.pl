#!/usr/bin/perl -w
#
# Cared for by Filipe G. Vieira <>
#
# Copyright Filipe G. Vieira
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

    prune_graph.pl v1.0.6

=head1 SYNOPSIS

    perl prune_graph.pl --in_file /path/to/input/file [--subset /path/to/subset/file] [--max_dist +Inf] [--min_weight 0] [--weight_field 4] [--print_excl /path/to/subset/out_excl.gz]

    --in_file       = File with input network (default STDIN)
    --subset        = File with node IDs to include (one per line)
    --max_dist      = Maximum distance between nodes (input file 3rd column) to assume they are connected
    --min_weight    = Minimum weight to consider an edge (in --weight_field)
    --weight_field  = Which column from input file has the weight
    --print_excl    = File to dump the name of excluded nodes

=head1 DESCRIPTION

    Assuming a network of linked nodes, this script will return a set of unlinked nodes. For that, 
    it will iteratively prune the until there are only unlinked nodes. Only nodes at less than 
    '--max_dist' and with more than '--min_weight' weight are considered and, on each round, the 
    most connected node (as in the node with the higher sum of edge weights) is removed. Input file 
    should be tab-sepparated into: node1, node2, dist, weight. The method assumes an undirected graph 
    with, at most, one edge between any two nodes (only the highest weight edge is used).

    Examples of uses for this script might be LD prunning, where each SNP is a node and the edges 
    are linkage levels between them. Or to obtain unique sequences from a all-against-all BLAST 
    analyses, where each sequence is a node and the edges are BLAST hits.

=head1 AUTHOR

    Filipe G. Vieira - fgarrettvieira _at_ gmail _dot_ com

=head1 CONTRIBUTORS

    Additional contributors names and emails here

=cut

# Let the code begin...



use strict;
use Getopt::Long;
use Graph::Easy;
use Math::BigInt;
use IO::Zlib;

my ($in_file, $subset_file, $max_dist, $min_weight, $weight_field, $print_excl, $debug);
my ($cnt, $graph, @excl, %subset);

$in_file = '-';
$max_dist = '+inf';
$min_weight = 0;
$weight_field = 4;
$debug = 0;
$cnt = 0;

# Parse args
GetOptions('h|help'             => sub { exec('perldoc',$0); exit; },
           'i|in_file=s'        => \$in_file,
	   's|subset:s'         => \$subset_file,
	   'd|max_dist:s'       => \$max_dist,
           'r|min_weight:s'     => \$min_weight,
           'f|weight_field:s'   => \$weight_field,
	   'x|print_excl:s'     => \$print_excl,
	   'debug!'             => \$debug,
    );


# Parse subset file
if($subset_file){
    my $SUBSET = new IO::Zlib;
    $SUBSET->open($subset_file, "r") || die("ERROR: cannot open subset file ".$subset_file.": ".$!."!");
    %subset = map {chomp($_); $_ => 1} <$SUBSET>;
    $SUBSET->close;
}


# Create graph
$graph = Graph::Easy->new();
$graph->timeout(60);


# Read file and build network
my $FILE = new IO::Zlib;
$FILE->open($in_file, "r") || die("ERROR: cannot open input file ".$in_file.": ".$!."!");
while(<$FILE>){
    $cnt++;
    my @interact = split(m/\t/, $_);
    die("ERROR: invalid entry in input file ".$in_file." line ".$cnt.": ".$!.".\n") if($#interact < 1);

    $graph->add_node($interact[0]) if(!$subset_file || defined($subset{$interact[0]}));
    $graph->add_node($interact[1]) if(!$subset_file || defined($subset{$interact[1]}));

    # Skip SNP if distance more than $max_dist
    next if( defined($interact[2]) && $interact[2] >= $max_dist );
    # Skip SNP if weight less than $min_weight
    next if( defined($interact[$weight_field-1]) && $interact[$weight_field-1] <= $min_weight );
    # Skip if NaN, +Inf, -Inf, ...
    my $x = Math::BigInt->new($interact[$weight_field-1]);
    next if($x->is_nan() || $x->is_inf('+') || $x->is_inf('-'));
    # Skip if not on the subset file
    next unless( !$subset_file || (defined($subset{$interact[0]}) &&  defined($subset{$interact[1]})) );

    # Check if edge already exists
    my $edge;
    if( ($edge = $graph->edge($interact[0], $interact[1])) ||
        ($edge = $graph->edge($interact[1], $interact[0])) ){

	my $weight = $edge->label;
	print(STDERR "### Edge between nodes ".$interact[0]." and ".$interact[1]." already exists with weight = ".$weight.".\n");
	# Check 'dist' attribute
	my $dist = $edge->get_attribute('comment');
	die("ERROR: distance between nodes ".$interact[0]." and ".$interact[1]." is not consistent across records (".$interact[2]." != ".$dist.")!") if($interact[2] != $dist);

	# Update 'label' attribute to store weight
	if($interact[$weight_field-1] > $weight){
	    print(STDERR "# Edge will be replaced by another with weight = ".$interact[$weight_field-1].".\n");
	    $edge->set_attribute('label', $interact[$weight_field-1]);
	}
    }else{
	$edge = $graph->add_edge($interact[0], $interact[1]);
	$edge->undirected(1);

	# Use 'label' attribute to store weight
	$edge->set_attribute('label', $interact[$weight_field-1]);
	# Use 'comment' attribute to store distance
	$edge->set_attribute('comment', $interact[2]);
    }
}
$FILE->close if($in_file ne '-');


# Print stats for read nodes and edges
print(STDERR "### Parsed ".$cnt." interactions between ".$graph->nodes()." nodes.\n");

# Parse unconnected nodes
foreach my $node ( $graph->nodes() ){
    if($node->edges() == 0){
	print($node->label."\n");
	$graph->del_node($node);
    }
}

# Print stats for actually connected nodes and edges
print(STDERR "### Kept ".$graph->edges()." interactions between ".$graph->nodes()." nodes.\n");


# DEBUG
if($debug){
    print(STDERR "### Initial network\n");
    print(STDERR $graph->as_ascii());
}


# Open file to store EXCLUDED nodes
open(EXCL, ">".$print_excl) || die("ERROR: cannot open EXCL file ".$print_excl.": ".$!."!") if($print_excl);
#my $EXCL = new IO::Zlib;
#$EXCL->open($print_excl, (substr($print_excl, -3) eq '.gz' ? "wb9" : "wT")) || die("ERROR: cannot open EXCL file ".$print_excl.": ".$!."!") if($print_excl);


if(0){
    # Slow algorithm
    @excl = &prune_graph($graph, $debug);
}else{
    @excl = &prune_graph_idx($graph, $debug);
}


if($print_excl) {
    # Print EXCL nodes
    print(EXCL $_."\n") for(@excl);
    # Close file-handle for 
    close(EXCL);
    #$EXCL->close;
}


# DEBUG
if($debug){
    print(STDERR "### Final network!\n");
    print(STDERR $graph->as_ascii());
}


# Print pruned nodes
for my $node ($graph->nodes()){
    print($node->label."\n");
}


exit();



sub prune_graph($$){
    my ($graph, $debug) = @_;
    my ($node, $edge, $edges_weight, $max_edges_weight, $most_heavy_node, @excl);

    while(1){
	# Find top-connected node
	$max_edges_weight = 0;
	foreach $node ( $graph->nodes() ){
	    $edges_weight = 0;
	    foreach $edge ( $node->edges() ){
		$edges_weight += $edge->label;
	    }
	    if($edges_weight > $max_edges_weight){
		$most_heavy_node = $node->name();
		$max_edges_weight = $edges_weight;
	    }
	}

	# If top node has no edge, break!
	last if($max_edges_weight <= 0);
	# Delete node from graph
	$graph->del_node($most_heavy_node);
	# Add node to the list of EXCL
	push(@excl, $most_heavy_node);

	if($debug){
	    print(STDERR "### Pruned network after most \"heavy\" node (".$most_heavy_node.") removed!\n");
	    print(STDERR $graph->as_ascii());
	}
    }

    return @excl;
}

sub prune_graph_idx($$){
    my ($graph, $debug) = @_;
    my ($node, $edge, $most_heavy_node, @excl, %nodes_idx);

    # Fill in Index table
    foreach $node ( $graph->nodes() ){
	$nodes_idx{$node->name}{'weight'} = 0;
	$nodes_idx{$node->name}{'n_edges'} = 0;
	foreach $edge ( $node->edges() ){
	    $nodes_idx{$node->name}{'weight'} += $edge->label;
	    $nodes_idx{$node->name}{'n_edges'}++;
	}
    }

    while(1){
	# Get ID of most "heavy" node
	$most_heavy_node = ( sort { $nodes_idx{$b}{'weight'} <=> $nodes_idx{$a}{'weight'} } grep {$nodes_idx{$_}{'n_edges'} > 0} keys(%nodes_idx) )[0];

	# If top node has no edges, break!
	last if(!defined($most_heavy_node));

	# Update Index table
	$node = $graph->node($most_heavy_node);
	foreach $edge ( $node->edges() ){
	    foreach $node ( $edge->nodes() ){
		$nodes_idx{$node->name}{'weight'} -= $edge->label;
		$nodes_idx{$node->name}{'n_edges'}--;
	    }
	}

	# Delete node from graph
	$graph->del_node($most_heavy_node);
	# Delete node from index
	delete($nodes_idx{$most_heavy_node});
	# Add node to the list of EXCL
	push(@excl, $most_heavy_node);

	if($debug){
	    print(STDERR "### Pruned network after most \"heavy\" node (".$most_heavy_node.") removed!\n");
	    print(STDERR $graph->as_ascii());
	}
    }

    if($debug){
	foreach $node (keys(%nodes_idx)){
	    print(STDERR $nodes_idx{$node}{'weight'}."\t".$nodes_idx{$node}{'n_edges'}."\n");
	}
    }

    return @excl;
}
