#!/usr/bin/env python3
#
# Maintained by Zachary J. Nolen
#
# Copyright Zachary J. Nolen
#
# This script aims to prune SNPs from pairwise linkage disequilibrium 
# estimates output by ngsLD. Not multithreaded, so only assign 
# multiple cores if required to increase memory allowance on your 
# cluster. Memory usage seems to require RAM equivalent to ~60-90% of 
# uncompressed input file size. Runtime depends on number of pairwise 
# comparisons, so it may be best to break up input into linkage groups, 
# prune per group, and merge afterwards. Assumes input first two fields 
# of input are the two positions under comparison in each line.

####### Housekeeping #######

# import modules needed to show help info
import argparse
import sys

# parse input arguments
parser = argparse.ArgumentParser(description='Prunes SNPs from ngsLD output to produce a list of sites in linkage equilibrium.')

parser.add_argument("--input",
					help="The .ld output file from ngsLD to be pruned. Can also be gzipped. [STDIN]",
					default=sys.stdin)
parser.add_argument("--output",
					help="The file to output pruned SNPs to. [STDOUT]",
					default=sys.stdout)
parser.add_argument("--field_dist", help="Field from input with distances. [3]", default=3)
parser.add_argument("--field_weight", help="Field from input with weights. [7]", default=7)
parser.add_argument("--max_dist", help="Maximum distance in bp between nodes to assume they are connected.", required=True)
parser.add_argument("--min_weight", help="Minimum weight of an edge to assume nodes are connected.", required=True)
parser.add_argument("--weight_type", help="How to calculate most connected node: sum of (a)bsolute edges' weight [default], sum of (e)dges' weight, or (n)umber of connections.", default="a")
parser.add_argument("--keep_heavy", help="Keep 'heaviest' nodes, instead of removing them (default)", action='store_true')
args = parser.parse_args()

# import remaining necessary modules
import datetime
import csv
import pandas as pd
from graph_tool.all import *
import gzip

# log start time
begin_time = datetime.datetime.now()

####### Read in data #######

# check if input is compressed or plain text
print("Checking if file is gzipped...", file=sys.stderr)
def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

# find out how many columns are in input
if is_gz_file(args.input):
	with gzip.open(args.input, mode="rt") as f:
		reader = csv.reader(f,delimiter="\t")
		ncol = len(next(reader))
else:
	with open(args.input) as f:
		reader = csv.reader(f,delimiter="\t")
		ncol = len(next(reader))
		
# set up properties based on input size (forces unused columns
# to smaller object type to reduce memory usage)
eprop_type_list = ["bool"] * (ncol-2)
eprop_name_list = ["na"] * (ncol-2)

eprop_type_list[int(args.field_dist)-3] = "int32_t"
eprop_name_list[int(args.field_dist)-3] = "dist"
eprop_type_list[int(args.field_weight)-3] = "double"
eprop_name_list[int(args.field_weight)-3] = "weight"

# read input into graph object, whether compressed or not
if is_gz_file(args.input):
	with gzip.open(args.input, mode="rt") as f:
		print("Reading in data...", file=sys.stderr)
		G = load_graph_from_csv(f, directed = False, 
			eprop_types = eprop_type_list,
			eprop_names = eprop_name_list, hashed = True, 
			csv_options = {'delimiter': '\t'})
else:
	print("Reading in data...", file=sys.stderr)
	G = load_graph_from_csv(args.input, directed = False, 
		eprop_types = eprop_type_list,
		eprop_names = eprop_name_list, hashed = True, 
		csv_options = {'delimiter': '\t'})

####### Filter graph #######

# get rid of unused properties (other weight measures) from graph
del G.ep["na"]

# scale weights by preferred weight method from arguments
if args.weight_type == "a":
	map_property_values(G.ep["weight"], G.ep["weight"], lambda x: abs(x))
elif args.weight_type == "n":
	map_property_values(G.ep["weight"], G.ep["weight"], lambda x: 1)

# create properties needed to filter out edges where dist > threshold and
# weight < threshold
drop_dist = G.new_edge_property("bool")
drop_weight = G.new_edge_property("bool")
weight = G.new_vertex_property("double")

# determine which edges to filter by dist and weight based on input arguments
if args.max_dist:
	print("Filtering edges by distance...", file=sys.stderr)
	map_property_values(G.ep["dist"], drop_dist, lambda x: x > int(args.max_dist))
	G.set_edge_filter(drop_dist, inverted=True)
	G.purge_edges()
	G.clear_filters()
if args.min_weight:
	print("Filtering edges by weight...", file=sys.stderr)
	map_property_values(G.ep["weight"], drop_weight, lambda x: x < float(args.min_weight))
	G.set_edge_filter(drop_weight, inverted=True)
	G.purge_edges()
	G.clear_filters()

####### Prune graph #######

# print out messages at start of pruning describing starting point and method
if args.keep_heavy:
	print("Beginning pruning by dropping neighbors of heaviest position...", 
			file=sys.stderr)
else:
	print("Beginning pruning by dropping heaviest position...",
		file=sys.stderr)
print("Starting with "+str(G.num_vertices())+" positions with "+
	str(G.num_edges())+" edges between them...", file=sys.stderr)

# prune while edges > 0; uncomment print lines for some debugging
while True:
	nodes = G.num_vertices()
	edges = G.num_edges()
	if edges == 0:
		break
	if (nodes % 1000 == 0) and nodes > 0:
		print(str(nodes)+" positions remaining with "+str(edges)+
			" edges between them...", file=sys.stderr)
	#print("Nodes remaining = "+str(nodes), file=sys.stderr)
	#print("Edges remaining = "+str(edges), file=sys.stderr)
	incident_edges_op(G, "out", "sum", G.ep["weight"], weight)
	max_weight = max(weight)
	#print("Heaviest node weight = "+str(max_weight), file=sys.stderr)
	heavy = find_vertex(G, weight, max_weight)
	if args.keep_heavy:
		#print("Dropping heavy neighbors...", file=sys.stderr)
		heavy_neighbors = G.get_out_neighbors(heavy[0])
		G.remove_vertex(heavy_neighbors, fast = True)
	else:
		#print("Dropping heavy...", file=sys.stderr)
		G.remove_vertex(heavy[0], fast = True)

####### Output creation #######

print("Pruning complete! Pruned graph contains "+str(nodes)+" positions.",
	file=sys.stderr)
print("Exporting kept positions to file...", file=sys.stderr)
pruned_df = pd.DataFrame([G.vp["name"][v] for v in G.get_vertices()])
pruned_df = pruned_df[0].str.split(pat=":", expand = True)
pruned_df.columns = ['chr','pos']
pruned_df.chr = pruned_df.chr.astype('string')
pruned_df.pos = pruned_df.pos.astype('int')
pruned_df = pruned_df.sort_values('pos')
pruned_df.to_csv(args.output, sep="_", quoting=csv.QUOTE_NONE, 
	header = False, index = False)

print("Total runtime: "+str(datetime.datetime.now() - begin_time),
	file=sys.stderr)