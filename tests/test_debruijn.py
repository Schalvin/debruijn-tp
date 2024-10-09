"""Tests for graph characteristic"""
import pytest
import os
import networkx as nx
import hashlib
from pathlib import Path
from .test_fixtures import global_data
from .context import debruijn
from debruijn import read_fastq
from debruijn import cut_kmer
from debruijn import build_kmer_dict
from debruijn import build_graph
from debruijn import get_starting_nodes
from debruijn import get_sink_nodes
from debruijn import get_contigs
from debruijn import save_contigs
from debruijn import path_average_weight
from debruijn import remove_paths
from debruijn import select_best_path
from debruijn import solve_bubble
from debruijn import simplify_bubbles
from debruijn import solve_entry_tips
from debruijn import solve_out_tips


def test_read_fastq(global_data):
    """Test fastq reading"""
    fastq_reader = read_fastq(Path(__file__).parent / "test_two_reads.fq")
    assert (
        next(fastq_reader)
        == "TCAGAGCTCTAGAGTTGGTTCTGAGAGAGATCGGTTACTCGGAGGAGGCTGTGTCACTCATAGAAGGGATCAATCACACCCACCACGTGTACCGAAACAA"
    )
    assert (
        next(fastq_reader)
        == "TTTGAATTACAACATCCATATGTTCTTGATGCTGGAATTCCAATATCTCAGTTGACAGTGTGCCCTCACCAGTGGATCAATTTACGAACCAACAATTGTG"
    )
    global_data.grade += 1


def test_cut_kmer(global_data):
    """Test cut_kmer"""
    kmer_reader = cut_kmer("TCAGA", 3)
    assert next(kmer_reader) == "TCA"
    assert next(kmer_reader) == "CAG"
    assert next(kmer_reader) == "AGA"
    global_data.grade += 1


def test_build_kmer_dict(global_data):
    """Test kmer dict"""
    kmer_dict = build_kmer_dict(Path(__file__).parent / "test_build.fq", 3)
    assert len(kmer_dict.keys()) == 4
    assert "TCA" in kmer_dict
    assert "CAG" in kmer_dict
    assert "AGA" in kmer_dict
    assert "GAG" in kmer_dict
    assert kmer_dict["AGA"] == 2
    global_data.grade += 2


def test_build_graph(global_data):
    """Test build graph"""
    kmer_dict = {"GAG": 1, "CAG": 1, "AGA": 2, "TCA": 1}
    graph = build_graph(kmer_dict)
    # TCAGAGA
    # TCA TC CA
    # CAG CA AG
    # AGA AG GA
    # GAG GA AG
    # AGA AG GA
    assert graph.number_of_nodes() == 4
    assert graph.number_of_edges() == 4
    assert "AG" in graph
    assert "GA" in graph
    assert graph.edges["AG", "GA"]["weight"] == 2
    global_data.grade += 4


def test_get_starting_nodes(global_data):
    graph = nx.DiGraph()
    graph.add_edges_from([(1, 2), (3, 2), (2, 4), (4, 5), (5, 6), (5, 7)])
    nodes = get_starting_nodes(graph)
    assert len(nodes) == 2
    assert 1 in nodes
    assert 3 in nodes
    global grade
    global_data.grade += 1


def test_get_sink_nodes(global_data):
    graph = nx.DiGraph()
    graph.add_edges_from([(1, 2), (3, 2), (2, 4), (4, 5), (5, 6), (5, 7)])
    nodes = get_sink_nodes(graph)
    assert len(nodes) == 2
    assert 6 in nodes
    assert 7 in nodes
    global_data.grade += 1


def test_get_contigs(global_data):
    graph = nx.DiGraph()
    graph.add_edges_from(
        [
            ("TC", "CA"),
            ("AC", "CA"),
            ("CA", "AG"),
            ("AG", "GC"),
            ("GC", "CG"),
            ("CG", "GA"),
            ("GA", "AT"),
            ("GA", "AA"),
        ]
    )
    contig_list = get_contigs(graph, ["TC", "AC"], ["AT", "AA"])
    results = ["TCAGCGAT", "TCAGCGAA", "ACAGCGAT", "ACAGCGAA"]
    assert len(contig_list) == 4
    for contig in contig_list:
        assert contig[0] in results
        assert contig[1] == 8
    global_data.grade += 3


def test_save_contigs(global_data):
    test_file = Path(__file__).parent / "test.fna"
    contig = [("TCAGCGAT", 8), ("TCAGCGAA", 8),
              ("ACAGCGAT", 8), ("ACAGCGAA", 8)]
    save_contigs(contig, test_file)
    with test_file.open("rb") as contig_test:
        assert (
            hashlib.md5(contig_test.read()).hexdigest()
            == "ca84dfeb5d58eca107e34de09b3cc997"
        )
    global_data.grade += 2


def test_path_weight(global_data):
    graph = nx.DiGraph()
    graph.add_weighted_edges_from(
        [(1, 2, 5), (3, 2, 10), (2, 4, 10), (4, 5, 3), (5, 6, 10), (5, 7, 10)]
    )
    assert path_average_weight(graph, [1, 2, 4, 5]) == 6.0
    global_data.grade += 1


def test_remove_paths(global_data):
    graph_1 = nx.DiGraph()
    graph_2 = nx.DiGraph()
    graph_3 = nx.DiGraph()
    graph_4 = nx.DiGraph()
    graph_1.add_edges_from([(1, 2), (3, 2), (2, 4), (4, 5), (5, 6), (5, 7)])
    graph_2.add_edges_from([(1, 2), (3, 2), (2, 4), (4, 5), (5, 6), (5, 7)])
    graph_3.add_edges_from([(1, 2), (3, 2), (2, 4), (4, 5), (5, 6), (5, 7)])
    graph_4.add_edges_from([(1, 2), (3, 2), (2, 4), (4, 5), (5, 6), (5, 7)])
    graph_1 = remove_paths(graph_1, [(1, 2)], True, False)
    graph_2 = remove_paths(graph_2, [(5, 7)], False, True)
    graph_3 = remove_paths(graph_3, [(2, 4, 5)], False, False)
    graph_4 = remove_paths(graph_4, [(2, 4, 5)], True, True)
    assert (1, 2) not in graph_1.edges()
    assert (3, 2) in graph_1.edges()
    assert (5, 7) not in graph_2.edges()
    assert (5, 6) in graph_2.edges()
    assert 4 not in graph_3.nodes()
    assert (2, 4) not in graph_4.edges()
    assert (4, 5) not in graph_4.edges()
    assert 2 not in graph_4.nodes()
    assert 4 not in graph_4.nodes()
    assert 5 not in graph_4.nodes()
    global_data.grade += 3


def test_select_best_path(global_data):
    graph_1 = nx.DiGraph()
    graph_1.add_edges_from([(1, 2), (3, 2), (2, 4), (4, 5), (5, 6), (5, 7)])
    graph_1 = select_best_path(
        graph_1, [[1, 2], [3, 2]], [1, 1], [5, 10], delete_entry_node=True
    )
    assert (1, 2) not in graph_1.edges()
    assert (3, 2) in graph_1.edges()
    assert 1 not in graph_1.nodes()
    graph_2 = nx.DiGraph()
    graph_2.add_edges_from(
        [(1, 2), (3, 2), (2, 4), (4, 5), (5, 6), (5, 7), (7, 8)])
    graph_2 = select_best_path(
        graph_2, [[5, 6], [5, 7, 8]], [1, 2], [13, 10], delete_sink_node=True
    )
    assert (5, 7) not in graph_2.edges()
    assert (7, 8) not in graph_2.edges()
    assert (5, 6) in graph_2.edges()
    assert 7 not in graph_2.nodes()
    assert 8 not in graph_2.nodes()
    # Select heavier
    graph_3 = nx.DiGraph()
    graph_3.add_edges_from(
        [(1, 2), (3, 2), (2, 4), (4, 5), (2, 8), (8, 9), (9, 5), (5, 6), (5, 7)]
    )
    graph_3 = select_best_path(
        graph_3, [[2, 4, 5], [2, 8, 9, 5]], [1, 4], [13, 10])
    assert (2, 8) not in graph_3.edges()
    assert (8, 9) not in graph_3.edges()
    assert (9, 5) not in graph_3.edges()
    assert (2, 4) in graph_3.edges()
    assert (4, 5) in graph_3.edges()
    assert 8 not in graph_3.nodes()
    assert 9 not in graph_3.nodes()
    assert 2 in graph_3.nodes()
    assert 5 in graph_3.nodes()
    # Select longest
    graph_4 = nx.DiGraph()
    graph_4.add_edges_from(
        [(1, 2), (3, 2), (2, 4), (4, 5), (2, 8), (8, 9), (9, 5), (5, 6), (5, 7)]
    )
    graph_4 = select_best_path(
        graph_4, [[2, 4, 5], [2, 8, 9, 5]], [1, 4], [10, 10])
    assert (2, 4) not in graph_4.edges()
    assert (4, 5) not in graph_4.edges()
    assert (2, 8) in graph_4.edges()
    assert (8, 9) in graph_4.edges()
    assert (9, 5) in graph_4.edges()
    # Select random
    graph_5 = nx.DiGraph()
    graph_5.add_edges_from(
        [(1, 2), (3, 2), (2, 4), (4, 5), (2, 8), (8, 9), (9, 5), (5, 6), (5, 7)]
    )
    graph_5 = select_best_path(
        graph_5, [[2, 4, 5], [2, 8, 9, 5]], [1, 4], [10, 10])
    global_data.grade += 4


def test_solve_bubble(global_data):
    graph_1 = nx.DiGraph()
    graph_1.add_weighted_edges_from(
        [
            (1, 2, 10),
            (3, 2, 10),
            (2, 4, 15),
            (4, 5, 15),
            (2, 10, 10),
            (10, 5, 10),
            (2, 8, 3),
            (8, 9, 3),
            (9, 5, 3),
            (5, 6, 10),
            (5, 7, 10),
        ]
    )
    graph_1 = solve_bubble(graph_1, 2, 5)
    assert (2, 8) not in graph_1.edges()
    assert (8, 9) not in graph_1.edges()
    assert (9, 5) not in graph_1.edges()
    assert (2, 10) not in graph_1.edges()
    assert (10, 5) not in graph_1.edges()
    assert (2, 4) in graph_1.edges()
    assert (4, 5) in graph_1.edges()
    assert 8 not in graph_1.nodes()
    assert 9 not in graph_1.nodes()
    assert 10 not in graph_1.nodes()
    assert 2 in graph_1.nodes()
    assert 5 in graph_1.nodes()
    graph_2 = nx.DiGraph()
    graph_2.add_weighted_edges_from(
        [
            (1, 2, 10),
            (3, 2, 10),
            (2, 4, 10),
            (4, 5, 10),
            (2, 10, 10),
            (10, 5, 10),
            (2, 8, 10),
            (8, 9, 10),
            (9, 5, 10),
            (5, 6, 10),
            (5, 7, 10),
        ]
    )
    graph_2 = solve_bubble(graph_2, 2, 5)
    assert (2, 4) not in graph_2.edges()
    assert (4, 5) not in graph_2.edges()
    assert (2, 10) not in graph_2.edges()
    assert (10, 5) not in graph_2.edges()
    assert (2, 8) in graph_2.edges()
    assert (8, 9) in graph_2.edges()
    assert (9, 5) in graph_2.edges()
    global_data.grade += 2


def test_simplify_bubbles(global_data):
    graph_1 = nx.DiGraph()
    graph_1.add_weighted_edges_from(
        [
            (3, 2, 10),
            (2, 4, 15),
            (4, 5, 15),
            (2, 10, 10),
            (10, 5, 10),
            (2, 8, 3),
            (8, 9, 3),
            (9, 5, 3),
            (5, 6, 10),
            (5, 7, 10),
        ]
    )
    graph_1 = simplify_bubbles(graph_1)
    assert (2, 8) not in graph_1.edges()
    assert (8, 9) not in graph_1.edges()
    assert (9, 5) not in graph_1.edges()
    assert (2, 10) not in graph_1.edges()
    assert (10, 5) not in graph_1.edges()
    global_data.grade += 1


def test_solve_entry_tips(global_data):
    graph_1 = nx.DiGraph()
    graph_1.add_weighted_edges_from(
        [(1, 2, 10), (3, 2, 2), (2, 4, 15), (4, 5, 15)])
    graph_1 = solve_entry_tips(graph_1, [1, 3])
    assert (3, 2) not in graph_1.edges()
    assert (1, 2) in graph_1.edges()
    graph_2 = nx.DiGraph()
    graph_2.add_weighted_edges_from(
        [(1, 2, 2), (6, 3, 2), (3, 2, 2), (2, 4, 15), (4, 5, 15)]
    )
    graph_2 = solve_entry_tips(graph_2, [1, 6])
    assert (1, 2) not in graph_2.edges()
    assert (6, 3) in graph_2.edges()
    assert (3, 2) in graph_2.edges()
    global_data.grade += 4


def test_solve_out_tips(global_data):
    graph_1 = nx.DiGraph()
    graph_1.add_weighted_edges_from(
        [(1, 2, 15), (2, 3, 15), (3, 4, 15), (4, 5, 15), (4, 6, 2)]
    )
    graph_1 = solve_out_tips(graph_1, [5, 6])
    assert (4, 6) not in graph_1.edges()
    assert (4, 5) in graph_1.edges()
    graph_2 = nx.DiGraph()
    graph_2.add_weighted_edges_from(
        [(1, 2, 15), (2, 3, 15), (3, 4, 15), (4, 5, 2), (4, 6, 2), (6, 7, 2)]
    )
    graph_2 = solve_out_tips(graph_2, [5, 7])
    assert (4, 5) not in graph_2.edges()
    assert (6, 7) in graph_2.edges()
    global_data.grade += 4


# Explicitly call finalize function after all tests
def pytest_sessionfinish(session, exitstatus):
    global_data_instance = global_data(None)  # Pass None as a dummy request
    global_data_instance.grade = 0
    global_data_instance.finalize()
