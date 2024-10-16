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

from typing import Iterator, Dict, List
import numpy as np
import matplotlib.pyplot as plt
import textwrap
import statistics
from random import randint
import argparse
import os
import sys
from pathlib import Path
from networkx import (
    DiGraph,
    all_simple_paths,
    lowest_common_ancestor,
    has_path,
    random_layout,
    draw,
    spring_layout,
)
import matplotlib
from operator import itemgetter
import random
import re

random.seed(9001)

matplotlib.use("Agg")

__author__ = "Siann Chalvin"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Siann Chalvin"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Siann Chalvin"
__email__ = "siann@chalvin.org"
__status__ = "Developpement"


def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments():  # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(
        description=__doc__, usage="{0} -h".format(sys.argv[0])
    )
    parser.add_argument(
        "-i", dest="fastq_file", type=isfile, required=True, help="Fastq file"
    )
    parser.add_argument(
        "-k", dest="kmer_size", type=int, default=22, help="k-mer size (default 22)"
    )
    parser.add_argument(
        "-o",
        dest="output_file",
        type=Path,
        default=Path(os.curdir + os.sep + "contigs.fasta"),
        help="Output contigs in fasta file (default contigs.fasta)",
    )
    parser.add_argument(
        "-f", dest="graphimg_file", type=Path, help="Save graph as an image (png)"
    )
    return parser.parse_args()


def read_fastq(fastq_file: Path) -> Iterator[str]:
    """Extract reads from fastq files.

    :param fastq_file: (Path) Path to the fastq file.
    :return: A generator object that iterate the read sequences.
    """
    with open(fastq_file, "r") as f:
        for line in f:
            sequence = next(f).strip()
            next(f)
            next(f)
            yield sequence


def cut_kmer(read: str, kmer_size: int) -> Iterator[str]:
    """Cut read into kmers of size kmer_size.

    :param read: (str) Sequence of a read.
    :return: A generator object that provides the kmers (str) of size kmer_size.
    """
    k_over = kmer_size + (len(read) %
                          kmer_size)  # length of sequence end to remove

    for i in range(0, k_over+1):
        yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file: Path, kmer_size: int) -> Dict[str, int]:
    """Build a dictionnary object of all kmer occurrences in the fastq file

    :param fastq_file: (str) Path to the fastq file.
    :return: A dictionnary object that identify all kmer occurrences.
    """
    kmer_dict = {}
    for seq in read_fastq(fastq_file):
        kmers = cut_kmer(seq, kmer_size)
        for kmer in kmers:
            if kmer in kmer_dict:
                kmer_dict[kmer] = kmer_dict[kmer]+1
            else:
                kmer_dict[kmer] = 1
    return kmer_dict


def build_graph(kmer_dict: Dict[str, int]) -> DiGraph:
    """Build the debruijn graph

    :param kmer_dict: A dictionnary object that identify all kmer occurrences.
    :return: A directed graph (nx) of all kmer substring and weight (occurrence).
    """
    graph = DiGraph()
    for seq, w in kmer_dict.items():
        prefix = seq[:-1]
        suffix = seq[1:]

        graph.add_edge(prefix, suffix, weight=w)
    return graph


def remove_paths(
    graph: DiGraph,
    path_list: List[List[str]],
    delete_entry_node: bool,
    delete_sink_node: bool,
) -> DiGraph:
    """Remove a list of path in a graph. A path is set of connected node in
    the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param delete_entry_node: (boolean) True->We remove the first node of a path
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    for path in path_list:
        if (delete_entry_node is True) & (delete_sink_node is True):
            graph.remove_nodes_from(path)
        elif delete_entry_node is True:
            graph.remove_nodes_from(path[:-1])
        elif delete_sink_node is True:
            graph.remove_nodes_from(path[1:])
        else:
            graph.remove_nodes_from(path[1:-1])
    return graph


def select_best_path(
    graph: DiGraph,
    path_list: List[List[str]],
    path_length: List[int],
    weight_avg_list: List[float],
    delete_entry_node: bool = False,
    delete_sink_node: bool = False,
) -> DiGraph:
    """Select the best path between different paths

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param path_length_list: (list) A list of length of each path
    :param weight_avg_list: (list) A list of average weight of each path
    :param delete_entry_node: (boolean) True->We remove the first node of a path
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    if len(weight_avg_list) > 1:
        stdev_weight = statistics.stdev(weight_avg_list)
    else:
        stdev_weight = 0

    if len(weight_avg_list) > 1:
        stdev_length = statistics.stdev(path_length)
    else:
        stdev_length = 0

    if stdev_weight != 0:
        index = np.array(weight_avg_list).argmax()
    elif stdev_length != 0:
        index = np.array(path_length).argmax()
    else:
        index = randint(0, len(path_list))

    path_list.pop(index)
    graph = remove_paths(graph, path_list, delete_entry_node, delete_sink_node)

    return graph


def path_average_weight(graph: DiGraph, path: List[str]) -> float:
    """Compute the weight of a path

    :param graph: (nx.DiGraph) A directed graph object
    :param path: (list) A path consist of a list of nodes
    :return: (float) The average weight of a path
    """
    return statistics.mean(
        [d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)]
    )


def solve_bubble(graph: DiGraph, ancestor_node: str, descendant_node: str) -> DiGraph:
    """Explore and solve bubble issue

    :param graph: (nx.DiGraph) A directed graph object
    :param ancestor_node: (str) An upstream node in the graph
    :param descendant_node: (str) A downstream node in the graph
    :return: (nx.DiGraph) A directed graph object
    """
    path_list = []
    len_list = []
    weight_avg_list = []
    for path in all_simple_paths(graph, ancestor_node, descendant_node):
        len_path = len(path)
        path_list.append(path)
        len_list.append(len_path)
        avg_weight = path_average_weight(graph, path)
        weight_avg_list.append(avg_weight)
    path = select_best_path(graph, path_list, len_list, weight_avg_list)
    return path


def simplify_bubbles(graph: DiGraph) -> DiGraph:
    """Detect and explode bubbles

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    bubble = False
    for node in graph.nodes():
        predecessors = list(graph.predecessors(node))
        if len(predecessors) > 1:
            for i in range(0, len(predecessors)):
                for j in range(i+1, len(predecessors)):
                    ancestor_node = lowest_common_ancestor(
                        graph, predecessors[i], predecessors[j])
                    if ancestor_node is not None:
                        bubble = True
                        break
    if bubble is True:
        print("simplified")
        graph = simplify_bubbles(solve_bubble(
            graph, ancestor_node, node))
    return graph


def solve_entry_tips(graph: DiGraph, starting_nodes: List[str]) -> DiGraph:
    """Remove entry tips

    :param graph: (nx.DiGraph) A directed graph object
    :param starting_nodes: (list) A list of starting nodes
    :return: (nx.DiGraph) A directed graph object
    """
    multi_tips = False
    for node in graph.nodes():
        predecessors = list(graph.predecessors(node))
        print(predecessors)
        tip_predecessors = list(set(predecessors).intersection(starting_nodes))
        if len(tip_predecessors) > 1:
            print(1)
            path_list = []
            len_list = []
            weight_avg_list = []
            paths = [list(all_simple_paths(graph, tip_pred, node))
                     for tip_pred in tip_predecessors]
            print(paths)
            for path in paths:
                len_path = len(path[0])
                path_list.append(path[0])
                len_list.append(len_path)
                avg_weight = path_average_weight(graph, path[0])
                weight_avg_list.append(avg_weight)
            best_path = select_best_path(
                graph, path_list, len_list, weight_avg_list)

            print(paths)
            paths.pop(best_path)
            graph = remove_paths(
                graph, paths, delete_entry_node=True, delete_sink_node=False)
            multi_tips = False

    return graph


def solve_out_tips(graph: DiGraph, ending_nodes: List[str]) -> DiGraph:
    """Remove out tips

    :param graph: (nx.DiGraph) A directed graph object
    :param ending_nodes: (list) A list of ending nodes
    :return: (nx.DiGraph) A directed graph object
    """
    pass


def get_starting_nodes(graph: DiGraph) -> List[str]:
    """Get nodes without predecessors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without predecessors
    """
    no_pred = []
    for node in graph.nodes():
        if any(True for _ in graph.predecessors(node)):
            continue
        else:
            no_pred.append(node)
    return no_pred


def get_sink_nodes(graph: DiGraph) -> List[str]:
    """Get nodes without successors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without successors
    """
    no_succ = []
    for node in graph.nodes():
        if any(True for _ in graph.successors(node)):
            continue
        else:
            no_succ.append(node)
    return no_succ


def get_contigs(
    graph: DiGraph, starting_nodes: List[str], ending_nodes: List[str]
) -> List:
    """Extract the contigs from the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param starting_nodes: (list) A list of nodes without predecessors
    :param ending_nodes: (list) A list of nodes without successors
    :return: (list) List of [contiguous sequence and their length]
    """
    paths = ()
    contigs = ()
    for start in starting_nodes:
        for end in ending_nodes:
            if has_path(graph, start, end):
                for path in all_simple_paths(graph, start, end):
                    paths += (path,)
    for path in paths:
        contig = path[0]
        for i in range(1, len(path)):
            contig += path[i][-1]
        contigs += ((contig, len(contig)),)
    return contigs


def save_contigs(contigs_list: List[str], output_file: Path) -> None:
    """Write all contigs in fasta format

    :param contig_list: (list) List of [contiguous sequence and their length]
    :param output_file: (Path) Path to the output file
    """
    with open(output_file, mode='w') as f_out:
        contig_n = 0
        for contig in contigs_list:
            f_out.write(f">contig_{contig_n} len={contig[1]}\n")
            f_out.write(textwrap.fill(contig[0], width=80))
            f_out.write(f"\n")
            contig_n += 1


def draw_graph(graph: DiGraph, graphimg_file: Path) -> None:  # pragma: no cover
    """Draw the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param graphimg_file: (Path) Path to the output file
    """
    fig, ax = plt.subplots()
    elarge = [(u, v)
              for (u, v, d) in graph.edges(data=True) if d["weight"] > 3]
    # print(elarge)
    esmall = [(u, v)
              for (u, v, d) in graph.edges(data=True) if d["weight"] <= 3]
    # print(elarge)
    # Draw the graph with networkx
    # pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(
        graph, pos, edgelist=esmall, width=6, alpha=0.5, edge_color="b", style="dashed"
    )
    # nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file.resolve())


# ==============================================================
# Main program
# ==============================================================
def main() -> None:  # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    file_in = args.fastq_file
    kmer_size = args.kmer_size
    file_out = args.output_file
    kmer_dict = build_kmer_dict(file_in, kmer_size)
    graph = build_graph(kmer_dict)
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    contigs = get_contigs(graph, starting_nodes, ending_nodes)
    save_contigs(contigs, file_out)

    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)


if __name__ == "__main__":  # pragma: no cover
    main()
