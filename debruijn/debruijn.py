#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    fq=open(fastq_file, "r")
    for line in fq :
        if line.startswith("@") or line.startswith("+") or line.startswith("J"):
            continue
        else :
            yield line.strip()
    fq.close()


def cut_kmer(read, kmer_size):
    i = 0
    while i < len(read):
        if i+(kmer_size-1) >= len(read):
            yield (read[i:])
            i = len(read)
        else :
            yield (read[i:i+kmer_size])
            i = i+1


def build_kmer_dict(fastq_file, kmer_size):
    kmer_dict = {}
    for Seq in read_fastq(fastq_file):
        for read in cut_kmer(Seq, kmer_size):
            if read not in kmer_dict:
                kmer_dict[read] = 1
            else:
                kmer_dict[read] += 1
    return kmer_dict

def build_graph(kmer_dict):
    kmer_list = list(kmer_dict.keys())
    graph = nx.DiGraph()
    for i in range(0,len(kmer_list),1):
        prefixe = len(kmer_list[i])-1
        suffixe = len(kmer_list[i])
        val_1 = kmer_list[i][0:prefixe]
        val_2 = kmer_list[i][1:suffixe]
        weight = kmer_dict[kmer_list[i]]
        graph.add_edge(val_1, val_2, weight = weight)
    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        sink = len(path)-1
        if delete_entry_node= True and delete_sink_node= True:
	    graph.remove_nodes_from(path)
        if delete_entry_node= True and delete_sink_node= False:
    	    graph.remove_nodes_from(path[0:sink])
        if delete_entry_node= False and delete_sink_node= True:
	    graph.remove_nodes_from(path[1:])
        if delete_entry_node= False and delete_sink_node= False:
	    graph.remove_nodes_from(path[1:sink])
    return graph

def std(data):
    st_dev = statistics.pstdev(data)
    return st_dev


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    if std(weight_avg_list) > 0:
	best_path = path_list[weight_avg_list.index(max(weight_avg_list))]
    elif std(path_length) > 0:
	best_path = path_list[path_length.index(max(path_length))]
    else:
	best_path = path_list[randint(0,len(path_list)-1)]
    path_list.remove(best_path)
    graph = remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
    return graph


def path_average_weight(graph, path):
    avg_weight = statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])
    return avg_weight 


def solve_bubble(graph, ancestor_node, descendant_node):
    path_list = []
    path_length = []
    weight_avg_list = []
    for path in nx.all_simple_paths(graph, source = ancestor_node, target = descendant_node):
	path_list.append(path)
	path_length.append(len(path))
	weight_avg_list.append(path_average_weight(graph, path))
    graph = select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False)
    return graph


def simplify_bubbles(graph):
    bubble = False
    graph_nodes = list(graph)
    for node in graph_nodes:
	if node in graph_nodes:
	    list_pred = list(graph.predecessors(node))
	    if len(list_pred) > 1:
		for i in range(0 ,len(list_pred), 1):
		    ancest_node = nx.lowest_common_ancestor(graph, list_pred[i], list_pred[i+1])
		    if ancest_node != None:
			bubble = True
			break
	if bubble:
	    break
    if bubble:
	graph = simplify_bubbles(solve_bubble(graph, ancest_node, node))
    return graph


def solve_entry_tips(graph, starting_nodes):
    graph_nodes = list(graph)
    for node in graph_nodes:
	path_list = []
	path_length = []
	weight_avg_list = []
	list_pred = list(graph.predecessors(node))
	if len(list_pred) > 1:
	    for start in starting_nodes:
		try:
		    for path in nx.all_simple_paths(graph, source = start, target = node):
			path_list.append(path)
			path_length.append(len(path))
			weight_avg_list.append(path_average_weight(graph, path))
		except:
		    pass
	graph = select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=True, delete_sink_node=False)
    return graph


def solve_out_tips(graph, ending_nodes):
    graph_nodes = list(graph)
    for node in graph_nodes:
	path_list = []
	path_length = []
	weight_avg_list = []
	list_pred = list(graph.successors(node))
	if len(list_pred) > 1:
	    for end in ending_nodes:
		try:
		    for path in nx.all_simple_paths(graph, source = node, target = end):
			path_list.append(path)
			path_length.append(len(path))
			weight_avg_list.append(path_average_weight(graph, path))
		except:
		    pass
	graph = select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=True)
    return graph


def get_starting_nodes(graph):
    starting_nodes = []
    node_list = list(graph.nodes)
    for node in node_list:
        node_pred = list(graph.predecessors(node))
        if len(node_pred) == 0:
            starting_nodes.append(node)
    return starting_nodes


def get_sink_nodes(graph):
    ending_nodes = []
    node_list = list(graph.nodes)
    for node in node_list:
        node_succ = list(graph.successors(node))
        if len(node_pred) == 0:
            ending_nodes.append(node)
    return ending_nodes


def get_contigs(graph, starting_nodes, ending_nodes):
    all_path = []
    tuples = []
    for start in starting_nodes:
	for ending in ending_nodes:
	    paths = (nx.all_simple_paths(graph, source = start, target = ending))
	    all_path.append(paths)

    for singles in all_path:
        for path in singles:
            contig = ""
            for i in range(0, len(path), 1):
		if i == 0:
                    contig = contig + path[i]
		elif path[i][0] == path[i-1][-1]:
		    contig = contig + path[i][-1]
            tuples.append([contig, len(contig)])
    return tuples


def save_contigs(contigs_list, output_file):
    save = open(output_file, "a")
    for i in range(0, len(contigs_list), 1):
        save.write(">contig_{} len={}".format(i, contigs_list[i][1]))
        save.write(fill(contigs_list[i][0]))
    save.close()


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """                                    
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
            pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit
    # graphe

    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(kmer_dict)

    # Résolution des bulles
    graph = simplify_bubbles(graph)

    # Résolution des pointes d’entrée et de sortie
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    graph = solve_entry_tips(graph, starting_nodes)
    graph = solve_out_tips(graph, ending_nodes)
    
    # Plot the graph
    if args.graphimg_file:
        draw_graph(graph, args.graphimg_file)

    # Ecriture du/des contigs
    if args.output_file:
	contigs_list = get_contigs(graph, starting_nodes, ending_nodes)
	save_contigs(contigs_list, args.output_file)

    # Save the graph in file
    # if args.graph_file:
    #    save_graph(graph, args.graph_file)


if __name__ == '__main__':
    main()
