# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 12:00:30 2023

@author: user
"""

def read_fasta_file(filename):
    """
    Reads a FASTA file and returns a dictionary with the sequence ID as key and the sequence as value.
    """
    sequences = {}
    with open(filename, "r") as file:
        sequence_id = None
        sequence = []
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if sequence_id:
                    sequences[sequence_id] = "".join(sequence)
                sequence_id = line[1:]
                sequence = []
            else:
                sequence.append(line)
        if sequence_id:
            sequences[sequence_id] = "".join(sequence)
    return sequences

import subprocess

def run_clustalo(input_file, output_file):
    """
    Runs Clustal Omega to perform multiple sequence alignment on the sequences in input_file and writes the result to output_file.
    """
    command = f"clustalo --infile={input_file} --outfile={output_file} --outfmt=clu"
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise Exception(f"Clustal Omega failed with error:\n{stderr.decode()}")

import numpy as np

def jukes_cantor_distance(s1, s2, p=0.01):
    """
    Computes the Jukes-Cantor distance between two DNA sequences s1 and s2, assuming a nucleotide substitution rate of p per site.
    """
    n = len(s1)
    assert n == len(s2), "Sequences must have the same length"
    diffs = sum(1 for i in range(n) if s1[i] != s2[i])
    return -1 * (3/4) * np.log(1 - (4/3) * (diffs/n) * (1/p))

def pairwise_distance_matrix(sequences, p=0.01):
    """
    Computes the pairwise distance matrix between a list of DNA sequences using the Jukes-Cantor model with a nucleotide substitution rate of p per site.
    """
    n = len(sequences)
    distances = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            distance = jukes_cantor_distance(sequences[i], sequences[j], p=p)
            distances[i,j] = distance
            distances[j,i] = distance
    return distances
