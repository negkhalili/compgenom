# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 12:02:37 2023

@author: user
"""

import pandas as pd

def read_expression_file(filename):
    """
    Reads an RNA-seq expression file and returns a pandas DataFrame with gene expression values for each sample.
    """
    df = pd.read_csv(filename, sep="\t", index_col=0)
    return df

def filter_genes_by_expression(df, min_mean_expression=0, min_fraction_samples=0):
    """
    Filters genes based on mean expression and the fraction of samples in which the gene is expressed.
    """
    mean_expression = df.mean(axis=1)
    fraction_samples_expressed = (df > 0).mean(axis=1)
    mask = (mean_expression >= min_mean_expression) & (fraction_samples_expressed >= min_fraction_samples)
    return df.loc[mask]

import pybedtools
import dexseq_prepare_annotation as dpa
import dexseq_count as dc
import dexseq_compare as dcomp
import subprocess

def prepare_annotation(annotation_file, output_dir):
    """
    Runs dexseq_prepare_annotation.py to prepare the annotation file for use with DEXSeq.
    """
    command = f"dexseq_prepare_annotation.py {annotation_file} {output_dir}"
    subprocess.run(command, shell=True)

def count_reads(bam_file, annotation_file, output_dir):
    """
    Runs dexseq_count.py to count reads overlapping exons in the annotation file.
    """
    bam = pybedtools.BedTool(bam_file)
    exons = dpa.AnnotationExons(annotation_file)
    counts = dc.count_reads(bam, exons, output_dir)
    return counts

def compare_conditions(counts1, counts2, output_dir, design_file):
    """
    Runs dexseq_compare.py to identify differentially spliced exons between two conditions.
    """
    result = dcomp.DexSeqComparison(counts1, counts2, output_dir, design_file)
    return result



"""import subprocess"""

def run_rmats(bam_file1, bam_file2, output_dir):
    """
    Runs rMATS to identify alternative splicing events between two samples.
    """
    command = f"python rMATS/rMATS-turbo-Linux-UCS4/rmats.py " \
              f"--b1 {bam_file1} " \
              f"--b2 {bam_file2} " \
              f"--gtf rMATS/gencode.v38.annotation.gtf " \
              f"--od {output_dir}"
    subprocess.run(command, shell=True)

def plot_rmats_event(event_id, output_file):
    """
    Plots the inclusion levels of alternative splicing events identified by rMATS.
    """
    command = f"python rMATS/rMATS-turbo-Linux-UCS4/rmats2sashimiplot.py " \
              f"-b1 rMATS/{event_id}.bam " \
              f"-b2 rMATS/{event_id}.bam " \
              f"-t {event_id} " \
              f"-e rMATS/{event_id}.txt " \
              f"-o {output_file}"
    subprocess.run(command, shell=True)


import pysam
"""import pybedtools"""

def calculate_psi(bam_file, annotation_file):
    """
    Calculates the percent spliced in (PSI) value for each alternative splicing event in a BAM file.
    """
    # Load the annotation file into a BedTool object
    annotation = pybedtools.BedTool(annotation_file)
    
    # Open the BAM file using pysam
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    # Loop over each alternative splicing event in the annotation file
    for feature in annotation:
        # Extract the reads overlapping the event
        reads = bam.fetch(feature.chrom, feature.start, feature.end)
        
        # Count the number of reads supporting each splice variant
        left_count = right_count = 0
        for read in reads:
            if read.is_unmapped or read.is_secondary:
                continue
            if read.is_reverse:
                if read.pos < feature.start:
                    right_count += 1
                elif read.aend > feature.end:
                    left_count += 1
            else:
                if read.pos < feature.start:
                    left_count += 1
                elif read.aend > feature.end:
                    right_count += 1
        
        # Calculate the PSI value
        total_count = left_count + right_count
        if total_count == 0:
            psi = 0
        else:
            psi = 100.0 * right_count / total_count
        
        # Print the result for this event
        print(f"{feature.name}\t{psi}")
