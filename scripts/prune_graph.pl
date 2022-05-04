#!/usr/bin/env perl
#
# Cared for by Filipe G. Vieira <>
#
# Copyright Filipe G. Vieira
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

    prune_graph.pl v1.2.2

=head1 SYNOPSIS

    perl prune_graph.pl [OPTIONS]

    --in_file       = File with input network [STDIN]
    --subset        = File with node IDs to include (one per line)
    --field_dist    = Column from input file with distances [3]
    --max_kb_dist   = Maximum distance (in kb) between nodes (in `--field_dist` column) to assume they are connected
    --field_weight  = Column from input file with weights [7]
    --min_weight    = Minimum weight (in `--field_weight` column) of an edge to assume nodes are connected
    --weight_type   = How to calculate most connected node: sum of (a)bsolute edges' weight [default], sum of (e)dges' weight, or (n)umber of connections
    --keep_heavy    = Keep 'heaviest' nodes, instead of removing them (default)
    --print_excl    = File to dump excluded nodes
    --out           = Path to output file [STDOUT]

=head1 DESCRIPTION

    Assuming a network of linked nodes, this script will return a set of unlinked nodes. For that, 
    it will iteratively prune the graph until there are only unlinked nodes. Only nodes at less than 
    '--max_kb_dist' and with more than '--min_weight' weight are considered. On each round, the most 
    'heavy' node (either the one with more connections or with the higher sum of edge weights) is 
    either kept or removed. Input file should be tab-sepparated into: node1, node2, dist, weight. 
    The method assumes an undirected graph with, at most, one edge between any two nodes (only the 
    highest weight edge is used).

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
use Math::BigFloat;
use List::AllUtils qw(max);
use Scalar::Util::Numeric qw(isint isfloat);
use IO::Zlib;

my ($in_file, $subset_file, $max_kb_dist, $field_dist, $min_weight, $field_weight, $weight_type, $weight_precision, $print_excl, $keep_heavy, $out_file, $debug);
my ($graph, @excl, %subset);

$in_file = '-';
$field_dist = 3;
$max_kb_dist = '+inf';
$field_weight = 7;
$min_weight = 0;
$weight_type = 'a';
$weight_precision = 4;
$out_file = '-';
$debug = 0;

# Parse args
GetOptions('h|help'             => sub { exec('perldoc',$0); exit; },
           'i|in_file=s'        => \$in_file,
           's|subset:s'         => \$subset_file,
           'field_dist:s'       => \$field_dist,
           'd|max_kb_dist:s'    => \$max_kb_dist,
           'field_weight:s'     => \$field_weight,
           'r|min_weight:s'     => \$min_weight,
           'weight_type:s'      => \$weight_type,
           'weight_precision:i' => \$weight_precision,
           'x|print_excl:s'     => \$print_excl,
           'keep_heavy!'        => \$keep_heavy,
           'o|out:s'            => \$out_file,
           'debug!'             => \$debug,
    );


# Convert --max_kb_dist to base pairs
my $max_dist = $max_kb_dist * 1000;
# Convert precision
$weight_precision = 10 ** $weight_precision;


# Parse subset file
if($subset_file){
    my $SUBSET = new IO::Zlib;
    $SUBSET->open($subset_file, "r") || die("ERROR: cannot open subset file ".$subset_file.": ".$!."!");
    %subset = map {chomp($_); $_ => 1} <$SUBSET>;
    $SUBSET->close;
}


# Read file
my $FILE = new IO::Zlib;
$FILE->open($in_file, "r") || die("ERROR: cannot open input file ".$in_file.": ".$!."!");


###############
# Build graph #
###############
my $cnt = 0;
my $edge;
$graph = Graph::Easy->new();
$graph->timeout(60);
print(STDERR "### Reading data from ".$in_file."\n");
while(<$FILE>){
    print(STDERR "# Found ".$cnt." edges (".$graph->edges." passed filters) between ".$graph->nodes." nodes so far...\n") unless(!$cnt || $cnt % 1e6);
    $cnt++;
    my @interact = split(m/\t/, $_);
    my $weight = $interact[$field_weight-1];
    my $dist = $interact[$field_dist-1];
    die("ERROR: invalid entry in input file ".$in_file." line ".$cnt.": \n".$_."\n") if($#interact < max($field_weight, $field_dist) - 1 ||
                                                                                        !defined($weight) || !isfloat($weight) ||
                                                                                        !defined($dist) || !isint($dist)
        );

    $graph->add_node($interact[0]) if(!$subset_file || defined($subset{$interact[0]}));
    $graph->add_node($interact[1]) if(!$subset_file || defined($subset{$interact[1]}));

    # Check if edge already exists
    die("ERROR: Edge between nodes ".$interact[0]." and ".$interact[1]." already exists!\n") if( $graph->edge($interact[0], $interact[1]) || $graph->edge($interact[1], $interact[0]) );

    # Skip if NaN, +Inf, -Inf, ...
    my $x = Math::BigFloat->new($weight);
    next if($x->is_nan() || $x->is_inf('+') || $x->is_inf('-'));
    # Skip SNP if distance more than $max_dist
    next if($dist > $max_dist);
    # Use the absolute weight value
    $weight = abs($weight) if($weight_type eq 'a');
    # Skip SNP if weight less than $min_weight
    next if($weight < $min_weight);
    # Skip if not on the subset file
    next unless(!$subset_file || (defined($subset{$interact[0]}) &&  defined($subset{$interact[1]})));

    # Use weights or just the number of connections?
    $weight = 1 if($weight_type eq 'n');
    # Add edge to graph
    $edge = $graph->add_edge($interact[0], $interact[1]);
    $edge->undirected(1);
    # Use 'label' attribute to store weight as integer (to avoid precision issues)
    $edge->set_attribute('label', int($weight * $weight_precision));
    # Use 'comment' attribute to store distance
    $edge->set_attribute('comment', $dist);
}
$FILE->close if($in_file ne '-');

# Print stats for read nodes and edges
print(STDERR "### Parsed a total of ".$cnt." edges (".$graph->edges." kept) between ".$graph->nodes." nodes\n");

# Open output filehandle
open(OUT, ">".$out_file) || die("ERROR: cannot open OUTPUT file ".$out_file.": ".$!."!");

# Parse unconnected nodes
my $cnt = 0;
foreach my $node ($graph->nodes){
    if($node->edges == 0){
        $cnt++;
        print(OUT $node->label."\n");
        $graph->del_node($node);
    }
}
print(STDERR "### Found ".$cnt." unlinked nodes; these will be printed and removed from graph now!\n") if($cnt);



###############
# Prune graph #
###############
# Print stats for actually connected nodes and edges
print(STDERR "### Pruning graph (".$graph->edges." edges between ".$graph->nodes." nodes)\n");


# DEBUG
if($debug){
    print(STDERR "### Initial network\n");
    print(STDERR $graph->as_ascii);
}


@excl = &prune_graph_idx($graph, $keep_heavy, $debug);


# Print EXCLUDED nodes to file
if($print_excl){
    # Open file to store EXCLUDED nodes
    if(substr($print_excl, -3) eq '.gz'){
        my $EXCL = IO::Zlib->new($print_excl, "wb9") || die("ERROR: cannot open EXCL file ".$print_excl.": ".$!."!");
        print($EXCL $_."\n") for(@excl);
        $EXCL->close;
    }else{
        open(EXCL, ">".$print_excl) || die("ERROR: cannot open EXCL file ".$print_excl.": ".$!."!");
        print(EXCL $_."\n") for(@excl);
        close(EXCL);
    }
}


# DEBUG
if($debug){
    print(STDERR "### Final network!\n");
    print(STDERR $graph->as_ascii);
}


# Print remaining nodes
for my $node ($graph->nodes){
    print(OUT $node->label."\n");
}
close(OUT);

print(STDERR "### Run finished (final graph with ".$graph->nodes." nodes)\n");
exit(0);



# prune graph (keep/remove most connected node) by searching the graph - SLOW!
sub prune_graph($$$){
    my ($graph, $keep_heavy, $debug) = @_;
    my ($cnt, $most_heavy_node, $most_heavy_node_weight, @excl);

    while(1){
        # Find top-connected node
        $most_heavy_node_weight = 0;
        foreach my $node ( sort{lc($a->name) cmp lc($b->name)} $graph->nodes ){
            my $node_weight = 0;
            foreach my $edge ($node->edges){
                $node_weight += $edge->label;
            }
            if($node_weight > $most_heavy_node_weight){
                $most_heavy_node = $node->name;
                $most_heavy_node_weight = $node_weight;
            }
        }

        # If top node has no edges, break!
        last if($most_heavy_node_weight <= 0);

        if($keep_heavy){
            foreach my $edge ($graph->node($most_heavy_node)->edges){
                foreach my $node_child ($edge->nodes){
                    next if($most_heavy_node eq $node_child->name);
                    # Delete node from graph
                    $graph->del_node($node_child);
                    # Add node to the list of EXCL
                    push(@excl, $node_child->name);
                }
            }
        }else{
            # Delete node from graph
            $graph->del_node($most_heavy_node);
            # Add node to the list of EXCL
            push(@excl, $most_heavy_node);
        }

        if($debug){
            print(STDERR "### Pruned network after processing most \"heavy\" node ".$most_heavy_node." (weight ".$most_heavy_node_weight."):\n");
            print(STDERR $graph->as_ascii);
        }

        # Print report to STDERR
        print(STDERR "# ".$cnt." nodes pruned...\n") unless($cnt == 0 || $cnt % 1e3);
        $cnt++;
    }

    return @excl;
}

# prune graph (keep/remove most connected node) using an index to speed up - FAST!
sub prune_graph_idx($$$){
    my ($graph, $keep_heavy, $debug) = @_;
    my ($cnt, @linked_nodes, $most_heavy_node, $most_heavy_node_weight, @excl, %nodes_idx);

    ### Fill in Index table
    # For each node, stores how many connected nodes (edges) and its total weight
    foreach my $node ($graph->nodes){
        $nodes_idx{$node->name}{'weight'} = 0;
        $nodes_idx{$node->name}{'n_edges'} = 0;
        foreach my $edge ($node->edges){
            $nodes_idx{$node->name}{'weight'} += $edge->label;
            $nodes_idx{$node->name}{'n_edges'}++;
        }
    }

    while(1){
        # Get ID of most "heavy" node
        $most_heavy_node = ( sort{$nodes_idx{$b}{'weight'} <=> $nodes_idx{$a}{'weight'} || lc($a) cmp lc($b)} keys(%nodes_idx) )[0];
        $most_heavy_node_weight = $nodes_idx{$most_heavy_node}{'weight'};

        # If top node has no edges, break!
        last if($most_heavy_node_weight <= 0);

        if($keep_heavy){
            foreach my $edge ($graph->node($most_heavy_node)->edges){
                foreach my $node_child ($edge->nodes){
                    next if($most_heavy_node eq $node_child->name);
                    # Delete node
                    remove_node_idx($graph, \%nodes_idx, $node_child);
                    # Add node to the list of EXCL
                    push(@excl, $node_child->name);
                }
            }
        }else{
            # Delete node
            remove_node_idx($graph, \%nodes_idx, $most_heavy_node);
            # Add node to the list of EXCL
            push(@excl, $most_heavy_node);
        }

        if($debug){
            print(STDERR "### Pruned network after processing most \"heavy\" node ".$most_heavy_node." (weight ".$most_heavy_node_weight."):\n");
            print(STDERR $graph->as_ascii);
        }

        # Print report to STDERR
        @linked_nodes = grep {$nodes_idx{$_}{'n_edges'} > 0} keys(%nodes_idx);
        print(STDERR "# ".$cnt." nodes pruned; ".($#linked_nodes + 1)." linked nodes remaining ...\n") unless($cnt == 0 || $cnt % 1e3);
        $cnt++;
    }

    # Sanity check if all nodes have been checked/removed
    foreach my $node (keys(%nodes_idx)){
        if($nodes_idx{$node}{'n_edges'} != 0 ||
           abs($nodes_idx{$node}{'weight'}) > 0){
            warn("ERROR: prunning done but node ".$node." still connected (".($nodes_idx{$node}{'weight'}/$weight_precision)." weight over ".$nodes_idx{$node}{'n_edges'}." edges)");
        }
    }
    foreach my $node ($graph->nodes){
        if($node->edges != 0){
            warn("ERROR: prunning done but node ".$node." still connected (".$node->edges." edges)");
        }
    }

    return @excl;
}



sub remove_node_idx($$$){
    my ($graph, $nodes_idx, $node) = @_;

    # Delete node from index
    foreach my $edge ($graph->node($node)->edges){
        foreach my $child_node ($edge->nodes){
            next if($node eq $child_node);
            $nodes_idx->{$child_node->name}{'weight'} -= $edge->label;
            $nodes_idx->{$child_node->name}{'n_edges'}--;
        }
    }

    # Delete node from graph
    $graph->del_node($node);
    # Delete node from index
    delete($nodes_idx->{$node});

    return $nodes_idx
}
