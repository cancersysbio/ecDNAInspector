# ecDNAInspector_functions
# This file provides the functionality of ecDNAInspector.
# Last update: 07/16/2025

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import csv
import ast
import ecDNAInspector.consensusClustering as consensusClustering
from sklearn.cluster import KMeans
from functools import partial

# FILE CONVERSIONS

"""
NOTE: File conversions are currently automated only for AmpliconArchitect/CoRAL input. 
To use output of other tools, please manually create the standard documentation provided 
for a cycle file (see sample converted cycle file for example).
"""


def convert_cycle_file_set(sample_amp_dict, input_cycles_path, output_cycles_path, copy_ct_threshold, cycle_type):
    """
    cycle_file_set()

    Inputs:
    -- sample_amp_dict (dictionary) --> key = sample (string), value = list of valid amplicon numbers
    -- input_cycles_path (file path) --> path to directory where raw cycle files are stored
    -- output_cycles_path (file path) --> path to direc tory where converted cycle files should be stored
    -- copy_ct_threshold (float) --> minimum copy number for cycle to be included in downstream analysis (recommended = 4)
    -- cycle_type (string) --> types of cycles to include ("circular", "linear", "all")

    Output:
    -- folder of converted files
    """
    for sample in sample_amp_dict:
        for amp in sample_amp_dict[sample]:
            input_cycle_file = input_cycles_path + sample + "_amplicon" + amp + "_cycles.txt"
            convert_single_cycle_file(input_cycle_file, output_cycles_path, copy_ct_threshold,
                                         cycle_type, sample, amp)
    return


def convert_single_cycle_file(cycle_file, output_cycle_directory, copy_ct_threshold, cycle_type, sample_name, amp):
    """
    convert_single_cycle_file()

    Inputs:
    -- cycle_file (string) --> path to file containing AA cycle file
    -- output_cycle_directory (string) --> path to directory to save new file with cycle information into
    -- copy_ct_threshold (float) --> minimum copy count level to include cycle
    -- cycle_type (string) --> types of cycles to include ("circular", "linear", "all")
    -- sample_name (string) --> name of sample for cycle
    -- amp (string) --> amplicon number of cycle

    Outputs:
    -- converted, stored file

    #####################################################
    Important notes about the AA/CoRAL cycle file output format:

    0's at either end of the cycle segment lists means the cycle is a linear contig type. Otherwise, the cycle
    is a circular contig type.

    AA/CoRAL provides cycle information such that positively-oriented strands go in a counter-clockwise
    direction. A "-" after the segment number means the segment order (smaller base pair, larger base pair) is oriented
    clockwise, while a "+" after the segment number means the segment order is oriented counter-clockwise. This affects
    how paired ends are determined.

    1-, 2-, 3+

    then the segments should be arranged:

    (1large, 1small) <--> (2large, 2small) <--> (3small, 3large)

    This ordering is automatically done when converting.
    """
    # read and store information from cycle file
    file = open(cycle_file, 'r')
    all_segments_dict = {}
    all_cycles_dict = {}
    for line in file:
        line = line.strip()
        # identify segment info lines
        if line[0] == 'S':
            segment_info = line.split("\t")
            segment_num = segment_info[1]
            segment_chromosome = segment_info[2]
            segment_start = int(segment_info[-2])
            segment_end = int(segment_info[-1])
            # add this segment into the segment dictionary
            all_segments_dict[segment_num] = [segment_chromosome, segment_start, segment_end]
        # identify cycle info lines
        if line[0] == 'C':
            cycle_info = line.split(";")
            cycle_num = cycle_info[0][cycle_info[0].index("=") + 1:]
            # check the copy count of the cycle and if it meets our criteria
            copy_ct = float(cycle_info[1][cycle_info[1].index("=") + 1:])
            if copy_ct >= copy_ct_threshold:
                segments_list = cycle_info[2][cycle_info[2].index("=") + 1:].split(",")
                # zero_start = 1 if the first segment is 0, otherwise zero_start = 0
                zero_start = not bool(int(segments_list[0][0]))
                # check if this is a circular or linear contig cycle, and if that meets our criteria
                if (cycle_type == "all") or (not zero_start and cycle_type == "circular") or (zero_start and cycle_type == "linear"):
                    if cycle_num not in all_cycles_dict:
                        all_cycles_dict[cycle_num] = [zero_start]
                    for segment in segments_list:
                        # get the orientation
                        orient = segment[-1]
                        segment_num = segment[:-1]
                        if segment_num == "0": continue
                        segment_chrom = all_segments_dict[segment_num][0]
                        segment_start = all_segments_dict[segment_num][1]
                        segment_end = all_segments_dict[segment_num][2]
                        if orient == "+":
                            cycle_tuple = (segment_chrom, segment_start, segment_end)
                        else:
                            cycle_tuple = (segment_chrom, segment_end, segment_start)
                        all_cycles_dict[cycle_num].append(cycle_tuple)

    # write information to output file
    os.makedirs(output_cycle_directory, exist_ok=True)
    output_file_name = output_cycle_directory + "/" + sample_name + "_processed_amplicon" + amp + "_cycles.txt"
    output_file = open(output_file_name, 'w')
    output_file.write("CycleNum\tCycleType\tSegmentsTuples(chr,start,end)\n")
    for cycle_num, cycle_segments_list in all_cycles_dict.items():
        if cycle_segments_list[0]:
            cycle_type = "linear"
        else: cycle_type = "circular"
        line = [cycle_num, cycle_type]
        for cycle_segment in cycle_segments_list[1:]:
            line.append(f"{cycle_segment[0]},{cycle_segment[1]},{cycle_segment[2]}")
        output_file.write("\t".join(line) + "\n")
    return


# PARSING FUNCTIONS


def parse_amp_info_file_for_sample_amp_dict(amp_file):
    """
    parse_amp_info_file_for_sample_amp_dict()

    This function returns a dictionary of samples with amplicons matching the requested amplicon type.

    Inputs:
    -- amp_file (file path) --> path to amplicon information file

    Outputs:
    -- sample_amp_dict (dictionary) --> {sample1: [amp1, amp2, ...], sample2: [amp1, amp2, ...], ...}
    """
    sample_amp_dict = {}
    f = open(amp_file)
    f.readline()
    for line in f:
        line = line.strip()
        line = line.split(',')
        sample = line[0]
        amp = line[1]
        if sample not in sample_amp_dict:
            sample_amp_dict[sample] = []
        sample_amp_dict[sample].append(amp)
    return sample_amp_dict


def parse_cycle_file_for_cycle_nums(cycle_file):
    """
    parse_cycle_file_for_cycle_nums()

    This function returns a list of the cycle numbers from a processed cycle file.

    Inputs:
    -- cycle_file (file path) --> path processed cycle file

    Outputs:
    -- cycles (list)
    """
    cycle_nums = []
    with open(cycle_file) as f:
        next(f)
        for line in f:
            line = line.strip()
            line = line.split('\t')
            cycle_num = line[0]
            cycle_nums.append(cycle_num)
    return cycle_nums


def parse_cycle_file_for_cycle_type(cycle_file, cycle_num):
    """
    parse_cycle_file_for_cycle_type()

    This function returns the cycle type for a specific cycle from a processed cycle file.

    Inputs:
    -- cycle_file (file path) --> path to processed cycle file
    -- cycle_num (string)

    Outputs:
    -- cycle_type (string) --> type of cycle ("circular", "linear", "all")
    """
    cycle_type = ''
    f = open(cycle_file)
    for line in f:
        line = line.strip()
        line = line.split('\t')
        if line[0] == cycle_num:
            cycle_type = line[1]
    return cycle_type


def parse_cycle_file_for_cycle_region_dict(cycle_file, cycle_num):
    """
    parse_cycle_file_for_cycle_region_dict()

    This function gets the region dictionary for a specific cycle from a processed cycle file.

    Inputs:
    -- cycle_file (file path) --> path to processed cycle file
    -- cycle_num (string)

    Outputs:
    -- region_dict (dict) --> {chrom1: [(start1, end1), (start2, end2), ...], chrom2: [...], ...}
    """
    region_dict = {}
    with open(cycle_file) as f:
        next(f)
        for line in f:
            line = line.strip()
            line = line.split('\t')
            if line[0] == cycle_num:
                for segment in line[2:]:
                    segment = segment.split(',')
                    chrom = segment[0]
                    if chrom[0] == "c":
                        chrom = chrom[3:]
                    start = int(segment[1])
                    end = int(segment[2])
                    if chrom not in region_dict:
                        region_dict[chrom] = []
                    region_dict[chrom].append((start, end))

    # condense any continuous regions
    for chrom in region_dict:
        condensed_regions = []
        regions = region_dict[chrom]
        cur_region = regions[0]
        if len(regions) > 1:
            for next_region in regions[1:]:
                # forward direction first
                if cur_region[0] < cur_region[1] and cur_region[1] + 1 == next_region[0] and next_region[0] < \
                        next_region[1]:
                    cur_region = (cur_region[0], next_region[1])
                elif cur_region[1] < cur_region[0] and cur_region[1] - 1 == next_region[0] and next_region[1] < \
                        next_region[0]:
                    cur_region = (cur_region[0], next_region[1])
                else:
                    condensed_regions.append(cur_region)
                    cur_region = next_region
        condensed_regions.append(cur_region)
        region_dict[chrom] = condensed_regions

    return region_dict


def parse_cycle_file_for_cycle_region_list(cycle_file, cycle_num):
    """
    parse_cycle_file_for_cycle_region_list()

    This function gets the region list for a specific cycle from a processed cycle file.

    Inputs:
    -- cycle_file (file path) --> path to processed cycle file
    -- cycle_num (string)

    Outputs:
    -- region_list (list) --> [(chrom1, start1, end1), (chrom1, start2, end2), ..., (chrom2, start1, end1), ...]
    """
    region_list = []
    region_dict = parse_cycle_file_for_cycle_region_dict(cycle_file, cycle_num)
    for chrom in region_dict:
        for start, end in region_dict[chrom]:
            region_list.append((chrom, start, end))
    return region_list


def parse_cycle_file_for_cycle_paired_ends_list(cycle_file, cycle_num, condense=False):
    """
    This function gets the paired ends list for a specific cycle from a processed cycle file.

    Inputs:
    -- cycle_file (file path) --> path to processed cycle file
    -- cycle_num (string)
    -- condense (Boolean) --> if True, exactly consecutive cycle segments will be condensed (suggested for long-read based predictions)

    Outputs:
    -- paired_ends_list (list) --> [(chrom1, be1, chrom2, be2), ...]

    """

    paired_ends_list = []
    f = open(cycle_file)
    for line in f:
        line = line.strip()
        line = line.split('\t')
        if line[0] == cycle_num:
            first_segment_info = line[2].split(',')
            first_chrom = first_segment_info[0]
            first_bp1 = int(first_segment_info[1])

            cur_chrom = first_chrom
            cur_bp1 = first_bp1
            cur_bp2 = int(first_segment_info[2])

            # if there is only one segment, connect it to itself
            # if you want this only for the "circular" type, include this: if line[1] == "circular"
            if len(line[2:]) == 1:
                paired_ends_list.append((cur_chrom, cur_bp2, first_chrom, first_bp1))

            # if there are multiple segments, continue
            elif len(line[2:]) > 1:
                for segment in line[3:]:
                    next_segment_info = segment.split(',')
                    next_chrom = next_segment_info[0]
                    next_bp1 = int(next_segment_info[1])
                    next_bp2 = int(next_segment_info[2])
                    if cur_chrom == next_chrom:
                        if cur_bp1 < cur_bp2 and cur_bp2 + 1 == next_bp1 and next_bp1 < next_bp2:
                            cur_bp2 = next_bp2
                        elif cur_bp1 > cur_bp2 and cur_bp2 - 1 == next_bp1 and next_bp1 > next_bp2:
                            cur_bp2 = next_bp2
                        else:
                            paired_ends_list.append((cur_chrom, cur_bp2, next_chrom, next_bp1))
                            cur_chrom = next_chrom
                            cur_bp1 = next_bp1
                            cur_bp2 = next_bp2
                    else:
                        paired_ends_list.append((cur_chrom, cur_bp2, next_chrom, next_bp1))
                        cur_chrom = next_chrom
                        cur_bp1 = next_bp1
                        cur_bp2 = next_bp2

                # connect the last region to the first region
                # if you want this only for the "circular" type, include this: if line[1] == "circular"
                paired_ends_list.append((next_chrom, next_bp2, first_chrom, first_bp1))

    # suggested for long-read based predictions
    if condense:
        paired_ends_list = condense_cycle_paired_ends_list(paired_ends_list)

    return paired_ends_list


def condense_cycle_paired_ends_list(paired_ends_list):
    """
    condense_cycle_paired_ends_list()

    This function is often helpful for long read sequencing-based predictions, which may span multiple consecutive reads of slightly different copy number
    (but likely do not represent genuine breakpoints).

    Inputs:
    -- paired_ends_list (list) --> [(chrom1, be1, chrom2, be2), (chrom1, be1, chrom2, be2), ..., (chrom1, be1, chrom2, be2), ...]

    Outputs:
    -- condensed_paired_ends_list (list) --> [(chrom1, be1, chrom2, be2), (chrom1, be1, chrom2, be2), ..., (chrom1, be1, chrom2, be2), ...]
    """

    condensed_paired_ends_list = []

    for paired_end in paired_ends_list:
        chr1 = paired_end[0]
        be1 = paired_end[1]
        chr2 = paired_end[2]
        be2 = paired_end[3]
        if chr1 == chr2:
            if be2 > be1 and be2 == be1 + 1:
                continue
            elif be2 < be1 and be2 == be1 - 1:
                continue
            else:
                condensed_paired_ends_list.append(paired_end)
        else:
            condensed_paired_ends_list.append(paired_end)

    return condensed_paired_ends_list


def parse_genome_dict_for_features_in_cycle(basepairs_by_chrom_dict, genome_dict, intersect_prop_buffer=0.5):
    """
    parse_genome_dict_for_features_in_cycle()

    This function returns a dictionary and list of genes present in the cycle.
    Gene information comes from appropriate reference gencode annotation file.

    Inputs:
    -- basepairs_by_chrom_dict (dict) --> {chr1: [(start1, end1), (start2, end2), ...], chr2: [...], ...}
    -- genome_dict (dict) --> {"genes": {1: {gene1: (pos1, pos2), gene2: (pos1, pos2), ...}, ...}, "transcripts": {}, ...}
    -- intersect_prop_buffer (float) --> proportion of the gene that must intersect with the cycle to be considered included in cycle (if 1.0, entire gene length must be included)

    Outputs:
    -- genes_dict (dict) --> {chr1: {gene1: {"type": type1, "loc": (pos1, pos2)}, gene2: {"type": type2, "loc": (pos1, pos2)}, ...}, ...}
    -- genes_list (list) --> [gene1, gene2, gene3, ...]
    """
    genes_dict = {}
    genes_list = []

    for chrom in basepairs_by_chrom_dict:
        basepair_segments = basepairs_by_chrom_dict[chrom]
        gene_info = genome_dict["gene"]
        if chrom in gene_info:
            genes_chrom_dict = genome_dict["gene"][chrom]
            for gene in genes_chrom_dict:
                gene_segment = genes_chrom_dict[gene]["loc"]
                for segment in basepair_segments:
                    intersection = calc_basepair_intersect([gene_segment], [segment])
                    feature_length = int(gene_segment[1]) - int(gene_segment[0])
                    # checks that the intersection is at least a certain proportion of the feature itself (e.g. at least X% of it intersects with the cycle region)
                    if intersection >= intersect_prop_buffer * feature_length:
                        if chrom not in genes_dict:
                            genes_dict[chrom] = {}
                        gene_type = genes_chrom_dict[gene]["gene_type"]
                        genes_dict[chrom][gene] = {"type": gene_type, "loc": gene_segment}
                        genes_list.append(gene)

        return genes_dict, genes_list


def parse_cosmic_dict_for_genes_in_cycle(basepairs_by_chrom_dict, cosmic_dict, intersect_prop_buffer=0.5):
    """
    parse_cosmic_dict_for_genes_in_cycle()

    This function returns a dict of cancer genes present in the cycle.
    Cancer gene information comes from appropriate reference COSMIC gene file.

    Inputs:
    -- basepairs_by_chrom_dict (dict) --> {chr1: [(start1, end1), (start2, end2), ...], chr2: [...], ...}
    -- cosmic_dict (dict) --> {1: {gene1: {"type": type1, "loc": (pos1, pos2)}, gene2: {"type": type2, "loc": (pos1, pos2)}, ...}, ...}
    -- intersect_prop_buffer (float) --> proportion of the gene that must intersect with the cycle to be considered included in cycle (if 1.0, entire gene length must be included)

    Outputs:
    -- cancer_genes_dict (dict) --> {chr1: {gene1: {"type": type1, "loc": (pos1, pos2)}, gene2: {"type": type2, "loc": (pos1, pos2)}, ...}, ...}
    -- cancer_genes_list (list) --> [gene1, gene2, gene3, ...]
    """
    cancer_genes_dict = {}
    cancer_genes_list = []
    for chrom in basepairs_by_chrom_dict:
        basepair_segments = basepairs_by_chrom_dict[chrom]
        if chrom in cosmic_dict:
            cosmic_chrom_dict = cosmic_dict[chrom]
            for gene in cosmic_chrom_dict:
                gene_segment = cosmic_chrom_dict[gene]["loc"]
                if gene_segment[0] == '' or gene_segment[1] == '':
                    continue
                for segment in basepair_segments:
                    intersection = calc_basepair_intersect([gene_segment], [segment])
                    segment_length = int(gene_segment[1]) - int(gene_segment[0])
                    if intersection >= intersect_prop_buffer * segment_length:
                        cancer_genes_list.append(gene)
                        if chrom not in cancer_genes_dict:
                            cancer_genes_dict[chrom] = {}
                        cancer_genes_dict[chrom][gene] = cosmic_dict[chrom][gene]
    return cancer_genes_dict, cancer_genes_list


def parse_reference_for_genome_dict(gene_file):
    """
    parse_reference_for_genome_dict()

    This function returns a dictionary of all genes and their locations, from the appropriate reference.
    Built for gencode.v[REFERENCE#].annotation.txt files.

    Inputs:
    -- reference_gene_file (file) --> path to file with gene names and locations

    Outputs:
    -- genome_dict (dict) --> {"genes": {1: {gene1: (pos1, pos2), gene2: (pos1, pos2), ...}, ...}, "transcripts": {}, ...}
    """
    ref_gene_file = open(gene_file)
    # skip first 5 lines
    ref_gene_file.readline()
    ref_gene_file.readline()
    ref_gene_file.readline()
    ref_gene_file.readline()
    ref_gene_file.readline()
    # begin storing gene info
    genome_dict = {}
    for line in ref_gene_file.readlines():
        line = line.split('\t')
        chrom = line[0][3:]
        pos1 = line[3]
        pos2 = line[4]
        # e.g. gene, transcript, exon, UTR, CDS, start codon, ...
        region_type = line[2]
        if region_type not in genome_dict:
            genome_dict[region_type] = {}
        if chrom not in genome_dict[region_type]:
            genome_dict[region_type][chrom] = {}
        region_info = line[8].split('; ')
        if region_type == "gene":
            gene_type = region_info[1].split(" ")[1][1:-1]
            gene_name = region_info[2].split(" ")[1][1:-1]
            if gene_name not in genome_dict[region_type][chrom]:
                genome_dict[region_type][chrom][gene_name] = {"gene_type": gene_type, "loc": (pos1, pos2)}
    return genome_dict


def parse_cosmic_for_gene_dict(cosmic_file):
    """
    parse_cosmic_for_gene_dict()

    This function returns a dictionary of genes and their locations, from the appropriate reference.

    Inputs:
    -- cosmic_file (file) --> file with gene names from Cosmic

    Outputs:
    -- genes_dict (dict) --> {1: {gene1: {"type": type1, "loc": (pos1, pos2)}, gene2: {"type": type2, "loc": (pos1, pos2)}, ...}, ...}
    """
    genes_dict = {}
    with open(cosmic_file, mode='r', encoding='utf-8') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        next(reader)
        for row in reader:
            gene_ID = row[0]
            gene_chr = row[3]
            if gene_chr not in genes_dict:
                genes_dict[gene_chr] = {}
            gene_pos1 = row[4]
            gene_pos2 = row[5]
            gene_classes = row[14]
            gene_classes = gene_classes.split(',')
            genes_dict[gene_chr][gene_ID] = {"type": gene_classes, "loc": (gene_pos1, gene_pos2)}
    return genes_dict


def parse_cycle_data_table_for_cycle_segment_dict(cycle_data_file):
    """
    parse_cycle_data_table_for_cycle_segment_dict()

    This function returns a dictionary of all cycle segments from the full cycle data file.

    Inputs:
    -- cycle_data_file (file) --> path to cycle data file

    Outputs:
    -- cycle_segment_dict (dict) --> {sample: {amp1: {cycle1: {chrom1: [(start1, end1), (start2, end2), ...], ...}, ...}, ...}, ...}
    """
    cycle_segment_dict = {}
    with open(cycle_data_file, 'r') as file:
        file_reader = csv.reader(file, delimiter=',')
        next(file_reader, None)
        for row in file_reader:
            sample = row[0]
            amplicon = row[1]
            cycle = row[2]
            # format (string): [(chrom1, start1, end1), (chrom1, start2, end2), ..., (chrom2, start1, end1), ...]
            # need to parse out data to return to list of lists
            cycle_regions_list = ast.literal_eval(row[7])
            if sample not in cycle_segment_dict:
                cycle_segment_dict[sample] = {}
            if amplicon not in cycle_segment_dict[sample]:
                cycle_segment_dict[sample][amplicon] = {}
            cycle_segment_dict[sample][amplicon][cycle] = {}
            for segment in cycle_regions_list:
                chrom = segment[0]
                start = int(segment[1])
                end = int(segment[2])
                if chrom not in cycle_segment_dict[sample][amplicon][cycle]:
                    cycle_segment_dict[sample][amplicon][cycle][chrom] = []
                cycle_segment_dict[sample][amplicon][cycle][chrom].append((start, end))
    return cycle_segment_dict


def parse_cycle_data_table_for_cycle_genes(cycle_data_file):
    """
    parse_cycle_data_table_for_cycle_genes()

    This function returns a dictionary of all genes in the cycle.

    Inputs:
    -- cycle_data_file (file) --> path to cycle data file

    Outputs:
    -- cycle_gene_dict (dict) --> {sample: {amp1: {cycle1: [gene1, gene2, ...], ...}, ...}, ...}
    """
    cycle_gene_dict = {}
    with open(cycle_data_file, 'r') as file:
        file_reader = csv.reader(file, delimiter=',')
        next(file_reader, None)
        for row in file_reader:
            sample = row[0]
            amplicon = row[1]
            cycle = row[2]
            cycle_genes_list = ast.literal_eval(row[-6])
            if sample not in cycle_gene_dict:
                cycle_gene_dict[sample] = {}
            if amplicon not in cycle_gene_dict[sample]:
                cycle_gene_dict[sample][amplicon] = {}
            cycle_gene_dict[sample][amplicon][cycle] = cycle_genes_list
    return cycle_gene_dict


def parse_cycle_data_table_for_cycle_unique_bes(cycle_data_file):
    """
    parse_cycle_data_table_for_cycle_unique_bes()

    This function returns a dictionary of the unique breakends in each cycle.

    Inputs:
    -- cycle_data_file (file) --> path to cycle data file

    Outputs:
    -- cycle_unique_bes_dict (dict) --> {sample: {amp1: {cycle1: {chrom1: [be1, be2, ...], ...}, ...}, ...}, ...}
    """
    cycle_unique_bes_dict = {}
    with open(cycle_data_file, 'r') as file:
        file_reader = csv.reader(file, delimiter=',')
        next(file_reader, None)
        for row in file_reader:
            sample = row[0]
            amplicon = row[1]
            cycle = row[2]
            cycle_paired_bes = ast.literal_eval(row[8])
            cycle_bes = {}
            for pe in cycle_paired_bes:
                chrom1 = pe[0]
                if chrom1 not in cycle_bes:
                    cycle_bes[chrom1] = []
                be1 = int(pe[1])
                if be1 not in cycle_bes[chrom1]:
                    cycle_bes[chrom1].append(be1)
                chrom2 = pe[2]
                if chrom2 not in cycle_bes:
                    cycle_bes[chrom2] = []
                be2 = int(pe[3])
                if be2 not in cycle_bes[chrom2]:
                    cycle_bes[chrom2].append(be2)
            if sample not in cycle_unique_bes_dict:
                cycle_unique_bes_dict[sample] = {}
            if amplicon not in cycle_unique_bes_dict[sample]:
                cycle_unique_bes_dict[sample][amplicon] = {}
            cycle_unique_bes_dict[sample][amplicon][cycle] = cycle_bes
    return cycle_unique_bes_dict

def normalize_chrom(chrom):
    """
    Return the bare chromosome name. Strips a leading 'chr' only if actually present.
    """
    if chrom is None:
        return chrom
    c = str(chrom).strip()
    if c[:3].lower() == "chr":
        c = c[3:]
    return c

def parse_blacklist_file_for_blacklist_dict(blacklist_file):
    """
    parse_blacklist_file_for_blacklist_dict()

    This function returns a dictionary of blacklist regions.

    Inputs:
    -- blacklist_file (bed file) --> contains information on blacklisted regions (make sure to use the correct assembly!)
    file option 1 retrieved from: http://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/problematic/
    file option 2 retrieved from: https://github.com/Boyle-Lab/Blacklist/tree/master/lists

    Outputs:
    -- blacklist_regions (dict) --> {chrom1: [(start1, end1), (start2, end2), ...], chrom2: [...], ...}
    """
    blacklist_regions_dict = {}
    with open(blacklist_file, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            chrom = normalize_chrom(fields[0])
            region_name = fields[3] if len(fields) > 3 else "blacklist"
            blacklist_regions_dict.setdefault(chrom, []).append(
                (int(fields[1]), int(fields[2]), region_name))
    return blacklist_regions_dict


def parse_SV_file_for_SVs(sv_file):
    """
    parse_SV_file_for_SVs()

    This function returns a dictionary of detected structural variants for the full dataset.

    Inputs:
    -- sv_file (file) --> file to SV information (must be in specific format: see sample SV sample)

    Outputs:
    -- sv_dict (dict) --> {sample: {sv_class1: [{chrom1: chr1, pos1: pos1, orient1: orient1, chrom2: chr2, ...}]}}
    """
    sv_dict = {}
    with open(sv_file, 'r') as file:
        file_reader = csv.reader(file, delimiter='\t')
        next(file_reader, None)
        for row in file_reader:
            sample = row[0]
            chrom_1 = row[1]
            pos_1 = int(row[2])
            orient_1 = row[3]
            chrom_2 = row[4]
            pos_2 = int(row[5])
            orient_2 = row[6]
            sv_class = row[7]
            sv_info_dict = {'chrom1': chrom_1, 'pos1': pos_1, 'orient1': orient_1,
                            'chrom2': chrom_2, 'pos2': pos_2, 'orient2': orient_2}
            if sample not in sv_dict:
                sv_dict[sample] = {}
            if sv_class not in sv_dict[sample]:
                sv_dict[sample][sv_class] = []
            sv_dict[sample][sv_class].append(sv_info_dict)
    return sv_dict


# STRUCTURAL METRIC CALCULATIONS


def calc_cycle_size(cycle_regions_list):
    """
    calc_cycle_size()

    This function returns the total size of a cycle.

    Inputs:
    -- cycle_regions_list (list) --> list of cycle regions in format: [(chr1, start1, end1), (chr1, start2, end2), ...]

    Outputs:
    -- cycle_size (int)
    """
    total_size = 0
    for region_tuple in cycle_regions_list:
        start_bp = region_tuple[1]
        end_bp = region_tuple[2]
        total_size += np.abs(int(end_bp) - int(start_bp) + 1)
    return total_size


def calc_cycle_num_breakpoints(cycle_regions_list):
    """
    calc_cycle_num_breakpoints()

    This function returns the total size of a cycle.

    Inputs:
    -- cycle_regions_list (list) --> list of cycle regions in format: [(chr1, start1, end1), (chr1, start2, end2), ...]

    Outputs:
    -- number of cycle breakpoints (int)
    """
    return len(cycle_regions_list)


def calc_cycle_num_chroms(cycle_regions_list):
    """
    calc_cycle_num_chroms()

    This function returns the total number of unqiue chromosomes in a cycle.

    Inputs:
    -- cycle_regions_list (list) --> list of cycle regions in format: [(chr1, start1, end1), (chr1, start2, end2), ...]

    Outputs:
    -- number of unique chromosomes (int)
    """
    unique_chroms = []
    for region in cycle_regions_list:
        chrom = region[0]
        if chrom not in unique_chroms:
            unique_chroms.append(chrom)
    return len(unique_chroms)


def calc_basepair_intersect(bps1, bps2):
    """
    calc_basepair_intersect()

    This function calculates the number of intersecting base pairs between two lists of base pair regions.

    Inputs:
    -- bps1 (list) --> list of base pair region tuples: [(bp1, bp2), (bp3, bp4), ...]
    -- bps2 (list) --> same as bps1

    Outputs:
    -- intersection (integer)
    """
    sorted_bps1 = [sorted(bp_range) for bp_range in bps1]
    sorted_bps2 = [sorted(bp_range) for bp_range in bps2]

    I = set()
    for start, end in sorted_bps1:
        for i in range(int(start), int(end) + 1):
            I.add(i)
    temp = set()
    for start, end in sorted_bps2:
        for i in range(int(start), int(end) + 1):
            if i in I:
                temp.add(i)
    intersection = len(temp)
    return intersection


def calc_basepair_union(bps1, bps2):
    """
    calc_basepair_union()

    This function calculates the union of base pairs between two lists of base pair regions.

    Inputs:
    -- bps1 (list) --> list of base pair region tuples: [(bp1, bp2), (bp3, bp4), ...]
    -- bps2 (list) --> same as bps1

    Outputs:
    -- union (integer)
    """
    sorted_bps1 = [sorted(bp_range) for bp_range in bps1]
    sorted_bps2 = [sorted(bp_range) for bp_range in bps2]

    U = set()
    for start, end in sorted_bps1:
        for i in range(int(start), int(end) + 1):
            U.add(i)
    for start, end in sorted_bps2:
        for i in range(int(start), int(end) + 1):
            U.add(i)
    union = len(U)
    return union


def calc_gene_intersect(list1, list2):
    """
    calc_gene_intersect()

    This function calculates the number of intersecting genes between two lists of genes.

    Inputs:
    -- list1 (list) --> list of genes: [gene1, gene2, ...]
    -- list2 (list) --> same as list2

    Outputs:
    -- oncogene_intersection (integer)
    -- intersect_list (list) --> list of common genes
    """
    if len(list1) == 0 or len(list2) == 0:
        return 0, []
    intersect_list = []
    for oncogene in list1:
        if oncogene in list2:
            intersect_list.append(oncogene)
    oncogene_intersection = len(intersect_list)
    return oncogene_intersection, intersect_list


def calc_gene_union(list1, list2):
    """
    calc_gene_union()

    This function calculates the length of the union set of genes between two gene lists.

    Inputs:
    -- list1 (list) --> list of genes: [gene1, gene2, ...]
    -- list2 (list) --> same as list2

    Outputs:
    -- oncogene_union (integer)
    """
    if len(list1) == 0:
        return len(list2)
    if len(list2) == 0:
        return len(list1)
    union_list = []
    for oncogene in list1:
        union_list.append(oncogene)
    for oncogene in list2:
        if oncogene not in union_list:
            union_list.append(oncogene)
    oncogene_union = len(union_list)
    return oncogene_union


def calc_breakend_intersection(list1, list2, be_overlap_buffer):
    """
    calc_breakend_intersection()

    This function calculates the number of breakends within range of each other across two lists.

    Inputs:
    -- list1 (list) --> list of breakends: [be1, be2, ...]
    -- list2 (list) --> same as list2
    -- be_overlap_buffer (int) --> maximum distance between breakends for them to be considered overlapping

    Outputs:
    -- intersect_count (integer) --> count of breakend intersections between the two cycles
    -- intersect_list (list) --> [(be1, be2), (be3, be4), ...]
    """
    intersect_count = 0
    intersect_list = []
    for be1 in list1:
        for be2 in list2:
            if np.abs(int(be1) - int(be2)) <= be_overlap_buffer:
                intersect_count += 1
                intersect_list.append((be1, be2))
    return intersect_count, intersect_list


def calc_breakend_union(list1, list2, be_intersection):
    """
    calc_breakend_union()

    This function calculates the length of the union set of all breakends.

    Inputs:
    -- list1 (list) --> list of breakends: [be1, be2, ...]
    -- list2 (list) --> same as list2
    -- be_intersection (int) --> count of breakend intersections between the two cycles

    Outputs:
    -- union (integer) --> count of breakend union between the two cycles
    """
    union = len(list1) + len(list2) - be_intersection
    return union


def calc_mapping_errors(cycle_paired_ends, blacklist_regions_dict, blacklist_buffer):
    """
    calc_mapping_errors()

    This function counts and records any breakends with associated mapping errors.

    Inputs:
    -- cycle_paired_ends (list) --> [(chrom1, be1, chrom2, be2), ...]
    -- blacklist_regions (dict) --> {chrom1: [(start1, end1), (start2, end2), ...], chrom2: [...], ...}
    -- blacklist_buffer (int) --> buffer region around blacklisted region for position to still be considered included

    Outputs:
    -- unmappable_breakends_count (int) --> number of unmappable breakends in the cycle
    -- unmappable_breakends_list (list) --> [(chrom1, bp1), (chrom2, bp2), (chrom2, bp3), ...]
    """
    unmappable_breakends_count = 0
    unmappable_breakends_list = []
    unmappable_reasons = []
    for pe in cycle_paired_ends:
        pe1_unmappable = False
        pe2_unmappable = False
        chrom1 = pe[0]
        bp1 = int(pe[1])
        chrom1_blacklist = blacklist_regions_dict[chrom1]
        for region in chrom1_blacklist:
            if region[0] - blacklist_buffer <= bp1 <= region[1] + blacklist_buffer:
                pe1_unmappable = True
                unmappable_reasons.append(region[2])
        if pe1_unmappable:
            unmappable_breakends_count += 1
            unmappable_breakends_list.append((chrom1, bp1))
        chrom2 = pe[2]
        bp2 = int(pe[3])
        chrom2_blacklist = blacklist_regions_dict[chrom2]
        for region in chrom2_blacklist:
            if region[0] - blacklist_buffer <= bp2 <= region[1] + blacklist_buffer:
                pe2_unmappable = True
                unmappable_reasons.append(region[2])
        if pe2_unmappable:
            unmappable_breakends_count += 1
            unmappable_breakends_list.append((chrom2, bp2))
    return unmappable_breakends_count > 0, unmappable_breakends_list, unmappable_reasons


def remove_small_deletions(cycle_paired_ends, small_del_size):
    """
    remove_small_deletions()

    This function removes any suspected small deletions from the list of cycle paired ends.

    Inputs:
    -- cycle_paired_ends (list) --> [(chrom1, be1, chrom2, be2), ...]
    -- small_del_size (int) --> maximum distance between bps on same chromosome to be considered small deletion

    Outputs:
    -- reduced_cycle_paired_ends (list) --> original cycle_paired_ends list with small deletions removed
    """
    # if there is just one circularized segment, we do not consider this a small deletion (even if the gap is < buffer)
    if len(cycle_paired_ends) == 1:
        return cycle_paired_ends
    reduced_cycle_paired_ends = []
    for i in range(len(cycle_paired_ends)):
        pe = cycle_paired_ends[i]
        if pe[0] == pe[2]:
            if np.abs(int(pe[1]) - int(pe[3])) <= small_del_size:
                continue
            else:
                reduced_cycle_paired_ends.append(pe)
        else:
            reduced_cycle_paired_ends.append(pe)
    return reduced_cycle_paired_ends


def calc_in_buffer_range(buffer, val_1, val_2):
    """
    calc_in_buffer_range()

    This function checks if two values are within a buffer region of each other.

    Inputs:
    -- buffer (int)
    -- val_1 (int)
    -- val_2 (int)

    Outputs:
    -- True/False
    """
    if abs(int(val_1) - int(val_2)) <= buffer:
        return True
    else:
        return False


def calc_bp_jaccard(cycle1_segments, cycle2_segments):
    """
    calc_bp_jaccard()

    This function computes the basepair jaccard between two cycles - the length of their common base pairs divided by
    the union length of the cycles.

    Inputs:
    -- cycle1_segments (dict) --> {chrom1: [(start1, end1), (start2, end2), ...], chrom2: [...], ...}
    -- cycle2_segments (dict) --> {chrom1: [(start1, end1), (start2, end2), ...], chrom2: [...], ...}

    Outputs:
    -- jaccard (float)
    """
    total_intersection = 0
    total_union = 0
    for chrom in cycle1_segments:
        if chrom in cycle2_segments:
            total_intersection += calc_basepair_intersect(cycle1_segments[chrom], cycle2_segments[chrom])
            total_union += calc_basepair_union(cycle1_segments[chrom], cycle2_segments[chrom])
    if total_union == 0:
        jaccard = 0
    else:
        jaccard = total_intersection / total_union
    return jaccard


def calc_gene_jaccard(cycle1_genes, cycle2_genes):
    """
    calc_gene_jaccard()

    This function computes the gene jaccard between two cycles - the number of common genes divided by the union of
    their gene lists.

    Inputs:
    -- cycle1_genes (list) --> [gene1, gene2, ...]
    -- cycle2_genes (list) --> [gene1, gene2, ...]

    Outputs:
    -- jaccard (float)
    -- intersect_list (list) --> list of common genes
    """
    intersection, intersect_list = calc_gene_intersect(cycle1_genes, cycle2_genes)
    union = calc_gene_union(cycle1_genes, cycle2_genes)
    if union == 0:
        jaccard = 0
    else:
        jaccard = intersection / union
    return jaccard, intersect_list


def calc_be_jaccard(cycle1_bes, cycle2_bes, be_overlap_buffer):
    """
    This function computes the breakend jaccard between two cycles - the number of unique breakends within range of each
    other across cycles, divided by the union number of unique breakends between the cycles.

    Inputs:
    -- cycle1_bes (list) --> [be1, be2, ...]
    -- cycle2_bes (list) --> [be1, be2, ...]
    -- be_overlap_buffer (int) --> maximum distance between breakends for them to be considered overlapping

    Outputs:
    -- jaccard (float)
    -- total_intersect_list (list) --> list of common breakends
    """
    total_intersection = 0
    total_union = 0
    total_intersect_list = []
    for chrom in cycle1_bes:
        if chrom in cycle2_bes:
            intersection, intersect_list = calc_breakend_intersection(cycle1_bes[chrom], cycle2_bes[chrom], be_overlap_buffer)
            total_intersection += intersection
            for be_intersect in intersect_list:
                total_intersect_list.append((chrom, be_intersect[0], be_intersect[1]))
            total_union += calc_breakend_union(cycle1_bes[chrom], cycle2_bes[chrom], total_intersection)
    if total_union == 0:
        jaccard = 0
    else:
        jaccard = total_intersection / total_union
    return jaccard, total_intersect_list


def calc_TPR_and_FPR(cycle_paired_ends, SVs_dict, sample_name, in_range_buffer, small_del_size):
    """
    calc_TPR_and_FPR()

    This function calculates the true positive rate and false positive rate of paired ends in a single cycle.

    Inputs:
    -- cycle_paired_ends (list) --> [(chrom1, be1, chrom2, be2), ...]
    -- SVs_dict (dict) --> {sample: {sv_class1: [{chrom1: chr1, pos1: pos1, orient1: orient1, chrom2: chr2, ...}]}}
    -- in_range_buffer (int) --> maximum distance between bps to be considered within sufficient range of each other
    -- small_del_size (int) --> user-defined threshold size for small deletions.

    Outputs:
    -- cycle_TPR (float)
    -- cycle_FPR (float)
    -- TP_list (list) --> list of true positive paired ends --> [(chrom1, be1, chrom2, be2), ...]
    -- FP_list (list) --> list of false positive paired ends
    -- small_del_removed (int) --> number of small dels removed in filtering
    """

    TP_count = 0
    TP_list = []
    FP_count = 0
    FP_list = []
    pre_small_del_removal_PE_count = len(cycle_paired_ends)
    cycle_paired_ends = remove_small_deletions(cycle_paired_ends, small_del_size)
    small_del_removed = pre_small_del_removal_PE_count - len(cycle_paired_ends)
    if sample_name not in SVs_dict:
        print("No SV information for sample.")
        return "NF", "NF", "NF", "NF", small_del_removed
    if len(cycle_paired_ends) == 0:
        print("No paired ends for validation.")
        return "NF", "NF", "NF", "NF", small_del_removed
    for cycle_pe in cycle_paired_ends:
        supported = False
        pe_chrom_1 = cycle_pe[0]
        pe_pos_1 = cycle_pe[1]
        pe_chrom_2 = cycle_pe[2]
        pe_pos_2 = cycle_pe[3]
        SVs_dict_by_type = SVs_dict[sample_name]
        for SV_type in SVs_dict_by_type:
            SVs_list = SVs_dict_by_type[SV_type]
            for SV_pe_dict in SVs_list:
                sv_chrom1 = SV_pe_dict["chrom1"]
                sv_chrom2 = SV_pe_dict["chrom2"]
                sv_pos1 = SV_pe_dict["pos1"]
                sv_pos2 = SV_pe_dict["pos2"]
                if pe_chrom_1 == sv_chrom1 and pe_chrom_2 == sv_chrom2:
                    if calc_in_buffer_range(in_range_buffer, pe_pos_1, sv_pos1) and \
                            calc_in_buffer_range(in_range_buffer, pe_pos_2, sv_pos2):
                        supported = True
                        break
                if pe_chrom_1 == sv_chrom2 and pe_chrom_2 == sv_chrom1:
                    if calc_in_buffer_range(in_range_buffer, pe_pos_1, sv_pos2) and \
                            calc_in_buffer_range(in_range_buffer, pe_pos_2, sv_pos1):
                        supported = True
                        break
        if supported:
            TP_count += 1
            TP_list.append((pe_chrom_1, pe_pos_1, pe_chrom_2, pe_pos_2))
        else:
            FP_count += 1
            FP_list.append((pe_chrom_1, pe_pos_1, pe_chrom_2, pe_pos_2))
    cycle_TPR = TP_count / len(cycle_paired_ends)
    cycle_FPR = FP_count / len(cycle_paired_ends)
    return cycle_TPR, TP_list, cycle_FPR, FP_list, small_del_removed


def calc_PFNR(cycle_paired_ends, small_del_size, bed_files_path, sample_name, SVs_dict, SV_buffer, TP_list,
              amplification_threshold=4):
    """
    calc_PFNR()

    This function calculates the putative false negative rate of breakends in a single cycle. An important distinction is that
    this rate is calculated as a rate of the number of breakends that are in paired ends that are NOT supported by
    the SV calls. Then each breakend is given a boolean of 1 or 0 depending on if it is involved in a different
    qualifying SV.

    Inputs:
    -- cycle_paired_ends (list) --> [(chrom1, be1, chrom2, be2), ...]
    -- small_del_size (int) --> user-defined threshold size for small deletions.
    -- bed_files_path (file path) --> path to standardized bed files for samples
    -- sample_name (string)
    -- SVs_dict (dict) --> {sample: {sv_class1: [{chrom1: chr1, pos1: pos1, orient1: orient1, chrom2: chr2, ...}]}}
    -- cycle (int) --> cycle number
    -- SV_buffer (int) --> buffer region around SVs to be considered overlapping with breakends
    -- amplification_threshold (int) --> minimum amplification for SV to be included in FN calculation

    Outputs:
    -- cycle_PFNR (float)
    -- PFN_list (list) --> list of putative false negative paired ends
    """
    if bed_files_path == "NO_BED_FILES":
        return "NF", "NF"
    if TP_list == "NF":
        return "NF", "NF"
    PFN_count = 0
    PFN_denominator = 0
    PFN_list = []
    cycle_paired_ends = remove_small_deletions(cycle_paired_ends, small_del_size)
    if len(cycle_paired_ends) == 0:
        return "NF", "NF"
    # build dictionary of amplified regions for this sample - only SVs in this region will be counted for PFNs
    amp_regions_dict = {}
    sample_bed_file = bed_files_path + sample_name + ".bed"
    if not os.path.exists(sample_bed_file):
        print("bed file not found")
        return "NF", "NF"
    f = open(sample_bed_file, 'r')
    f.readline()
    for line in f:
        line = line.strip()
        line_info = line.split('\t')
        line_chrom = line_info[0]
        line_pos1 = int(line_info[1])
        line_pos2 = int(line_info[2])
        # check this according to your coverage/bed file
        line_amp = float(line_info[4])
        if line_amp >= amplification_threshold:
            if line_chrom not in amp_regions_dict:
                amp_regions_dict[line_chrom] = []
            amp_regions_dict[line_chrom].append((line_pos1, line_pos2))
    # find PFNs
    for cycle_pe in cycle_paired_ends:
        # we only want to consider the breakpoints that are not part of a TP pair
        if cycle_pe not in TP_list:
            PFN_denominator += 2
            PE_1_PFN = 0
            PE_2_PFN = 0
            # define info on cycle PE
            cycle_chrom1 = cycle_pe[0]
            cycle_chrom2 = cycle_pe[2]
            cycle_pos1 = cycle_pe[1]
            cycle_pos2 = cycle_pe[3]
            # loop over the SVs available
            SVs_dict_by_type = SVs_dict[sample_name]
            for SV_type in SVs_dict_by_type:
                SVs_list = SVs_dict_by_type[SV_type]
                for SV_pe_dict in SVs_list:
                    sv_chrom1 = SV_pe_dict["chrom1"]
                    sv_chrom2 = SV_pe_dict["chrom2"]
                    sv_pos1 = int(SV_pe_dict["pos1"])
                    sv_pos2 = int(SV_pe_dict["pos2"])
                    # need both of the paired ends in SV to be amplified to consider using for PFNR calculation:
                    # one must be amplified if it overlaps with the ecDNA breakpoint, and the other
                    # must be amplified for us to consider it as a denominator-included SV call
                    if sv_chrom1 in amp_regions_dict and sv_chrom2 in amp_regions_dict:
                        sv_chrom_1_regions = amp_regions_dict[sv_chrom1]
                        sv_chrom_2_regions = amp_regions_dict[sv_chrom2]
                        sv1_pos_amped = 0
                        sv2_pos_amped = 0
                        for region in sv_chrom_1_regions:
                            if region[0] - SV_buffer < sv_pos1 < region[1] + SV_buffer:
                                sv1_pos_amped = 1
                        for region in sv_chrom_2_regions:
                            if region[0] - SV_buffer < sv_pos2 < region[1] + SV_buffer:
                                sv2_pos_amped = 1
                        if sv1_pos_amped and sv2_pos_amped:
                            if cycle_chrom1 == sv_chrom1:
                                if cycle_chrom2 == sv_chrom2:
                                    if calc_in_buffer_range(SV_buffer, cycle_pos1,
                                                            sv_pos1) and not calc_in_buffer_range(SV_buffer, cycle_pos2,
                                                                                                  sv_pos2):
                                        PE_1_PFN = 1
                                        PFN_list.append((cycle_chrom1, cycle_pos1, sv_chrom2, sv_pos2))
                            if cycle_chrom1 == sv_chrom2:
                                if cycle_chrom2 == sv_chrom1:
                                    if calc_in_buffer_range(SV_buffer, cycle_pos1,
                                                            sv_pos2) and not calc_in_buffer_range(SV_buffer, cycle_pos2,
                                                                                                  sv_pos1):
                                        PE_1_PFN = 1
                                        PFN_list.append((cycle_chrom1, cycle_pos1, sv_chrom1, sv_pos1))
                            if cycle_chrom1 == sv_chrom1:
                                if cycle_chrom2 == sv_chrom2:
                                    if not calc_in_buffer_range(SV_buffer, cycle_pos1,
                                                                sv_pos1) and calc_in_buffer_range(SV_buffer, cycle_pos2,
                                                                                                  sv_pos2):
                                        PE_2_PFN = 1
                                        PFN_list.append((cycle_chrom2, cycle_pos2, sv_chrom1, sv_pos1))
                            if cycle_chrom1 == sv_chrom2:
                                if cycle_chrom2 == sv_chrom1:
                                    if not calc_in_buffer_range(SV_buffer, cycle_pos1,
                                                                sv_pos2) and calc_in_buffer_range(SV_buffer, cycle_pos2,
                                                                                                  sv_pos1):
                                        PE_2_PFN = 1
                                        PFN_list.append((cycle_chrom2, cycle_pos2, sv_chrom2, sv_pos2))

            if PE_1_PFN:
                PFN_count += 1
            if PE_2_PFN:
                PFN_count += 1

    if PFN_denominator == 0:
        return 0, PFN_list
    else:
        PFNR = PFN_count / PFN_denominator

    return PFNR, PFN_list


# CLUSTERING


def cluster_with_consensus(cycle_level_table_path, cycle_level_data_clustered_table, features, L=2, K=5, H=1000, resample_proportion=0.8, n_init=10):
    """
    cluster_with_consensus()

    This function clusters the cycles based on structural metrics.

    Inputs:
    -- cycle_level_table_path (file path) --> file with all information by cycle
    -- L (int) --> minimum number of clusters to try
    -- K (int) --> maximum number of clusters to try
    -- H (int) --> number of resamples
    -- resample_proportion (float) --> proportion of samples to resample

    Output:
    -- None (adds cluster assignments to the cycle data table)
    """
    cycle_data = pd.read_csv(cycle_level_table_path)
    cycle_data.reset_index(drop=True, inplace=True)
    cycle_data = cycle_data[~cycle_data[features].eq("NF").any(axis=1)]
    clustering_cycle_data = cycle_data[features]
    print("Total cycles eligible for clustering: " + str(len(clustering_cycle_data)))

    consensus_cluster = consensusClustering.ConsensusCluster(
        cluster=partial(KMeans, n_init=n_init),
        L=L,
        K=K,
        H=H,
        resample_proportion=resample_proportion
    )

    consensus_cluster.fit(clustering_cycle_data)

    cluster_labels = consensus_cluster.predict_data(clustering_cycle_data)

    cycle_data['Cycle Cluster Assignment'] = cluster_labels

    cycle_data.to_csv(cycle_level_data_clustered_table, index=False)

    return


def post_clustering_vis_filters(raw_df):
    """
    post_clustering_vis_filters()

    This function applies filters to the cycle dataframe with clustering info, prior to visualizing clusters.

    Inputs:
    -- raw_df (DataFrame) --> unfiltered dataframe with clustering information

    Output:
    -- clustering_df (DataFrame) --> updated clustering df ready for plotting
    """

    # remove any cycles that were unassigned prior to visualization
    clustering_df = raw_df[raw_df["Cycle Cluster Assignment"] != "NF"]

    # remap variables to booleans for visualization
    clustering_df["Cycle Type"] = clustering_df["Cycle Type"].map({'linear': 0, 'circular': 1})
    # clustering_df["Cycle MEB"] = clustering_df["Cycle MEB"].map({False: 0, True: 1})

    # sort the cluster for visualization
    clustering_df = clustering_df.sort_values(by=["Cycle Cluster Assignment", "Cycle Type", "Cycle TPR"],
                                              ascending=[True, False, False])

    # add a name for each cycle to be displayed for visualization
    clustering_df["Sample Name"] = clustering_df[["Sample", "Amplicon", "Cycle"]].astype(str).agg('-'.join, axis=1)

    # confirm all features and metrics are numeric
    clustering_df = clustering_df.apply(pd.to_numeric, errors='coerce')

    return clustering_df


def plot_separated_clustering_heatmap(cohort_df, clustering_features, boolean_metrics, size_metric, structural_metrics,
                                      cluster_feature_labels, boolean_metric_labels, size_metric_labels,
                                      structural_metric_labels, sorting_features, sorting_features_dir,
                                      cluster_plot_path):
    """
    plot_separated_clustering_heatmap()

    Inputs:
    -- cohort_df (DataFrame) --> cohort dataframe with clustering information
    -- clustering_features (list) --> features used in the clustering
    -- boolean_metrics (list) --> Boolean metrics to be displayed
    -- size_metric (string) --> name of size metric to be displayed
    -- structural_metrics (list) --> list of structural metrics to be displayed
    -- cluster_feature_labels (list) --> list of labels for clustering feature metrics (in order with clustering_features list)
    -- boolean_metric_labels (list) --> list of labels for Boolean metrics (in order with boolean_metrics list)
    -- size_metric_labels (list) --> list/string with size metric label(s) for size metrics (in order with size_metric)
    -- structural_metric_labels (list) --> list of labels for structural metrics (in order with structural_metrics list)
    -- sorting_features (list) --> ordered list of features to sort the samples within each cluster by
    -- sorting_features_dir (Boolean) --> if True, samples will be sorted according to features in sorting_features list in ascending order; otherwise, descending
    -- cluster_plot_path (string) --> path to where to save output image

    Output:
    -- None (prints plot)
    """

    # sort the cohort_df by the provided sorting features and order
    cohort_df = cohort_df.sort_values(by=sorting_features, ascending=sorting_features_dir)

    # ensure features used in the clustering are numeric type and construct separated dataframes for each feature type
    cluster_df = cohort_df[clustering_features].apply(pd.to_numeric, errors='coerce')
    boolean_df = cohort_df[boolean_metrics].apply(pd.to_numeric, errors='coerce')
    size_df = cohort_df[size_metric].apply(pd.to_numeric, errors='coerce')
    structural_df = cohort_df[structural_metrics].apply(pd.to_numeric, errors='coerce')

    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, gridspec_kw={'height_ratios': [0.5, 0.5*len(clustering_features), 0.5*len(boolean_metrics), 0.5*len(size_metric), 0.5*len(structural_metrics)]}, figsize=(12, 8))

    cluster_labels = cohort_df["Cycle Cluster Assignment"].astype(str)
    unique_clusters = cluster_labels.unique()
    color_palette = sns.color_palette("husl", len(unique_clusters))
    cluster_colors = {cluster: color_palette[i] for i, cluster in enumerate(unique_clusters)}
    row_colors = cluster_labels.map(cluster_colors)
    color_array = np.array(row_colors.tolist())[np.newaxis, :, :]

    ax1.imshow(color_array, aspect='auto')
    ax1.set_xticks([])
    ax1.set_yticks([])

    cluster_counts = cluster_labels.value_counts().reindex(unique_clusters)

    start = 0
    for i, count in enumerate(cluster_counts):
        box_center_x = start + count / 2
        start += count
        ax1.text(box_center_x, 0.05, f'Cluster {i+1}\n(n = {count})', ha='center', va='center', color='black', fontsize=14)

    sns.heatmap(cluster_df.T, ax=ax2, cmap='Greens', cbar=False, xticklabels=False)
    ax2.set_ylabel('Clustering\nFeatures', fontsize=14)
    ax2.tick_params(axis='y', labelsize=12)
    ax2.set_yticklabels(cluster_feature_labels, ha="right", rotation=0)

    sns.heatmap(boolean_df.T, ax=ax3, cmap='Blues', cbar=False, xticklabels=False)
    ax3.set_ylabel('Boolean \nMetrics', fontsize=14)
    ax3.tick_params(axis='y', labelsize=12)
    ax3.set_yticklabels(boolean_metric_labels, ha="right", rotation=0)

    sns.heatmap(size_df.T, ax=ax4, cmap='Oranges', cbar=False, xticklabels=False)
    ax4.set_ylabel('Cycle \nSize', fontsize=14)
    ax4.tick_params(axis='y', labelsize=12)
    ax4.set_yticklabels(size_metric_labels, ha="right", rotation=0)

    sns.heatmap(structural_df.T, ax=ax5, cmap='Purples', cbar=False, xticklabels=False)
    ax5.set_ylabel('Structural \nMetrics', fontsize=14)
    ax5.tick_params(axis='y', labelsize=12)
    ax5.set_yticklabels(structural_metric_labels, ha="right", rotation=0)

    cbar_ax2 = fig.add_axes([0.2, 0.0, 0.3, 0.02])
    cbar_ax3 = fig.add_axes([0.55, 0.0, 0.3, 0.02])
    cbar_ax4 = fig.add_axes([0.2, -0.15, 0.3, 0.02])
    cbar_ax5 = fig.add_axes([0.55, -0.15, 0.3, 0.02])

    cbar1 = fig.colorbar(plt.cm.ScalarMappable(cmap='Greens', norm=plt.Normalize(vmin=cluster_df.min().min(),
                                                                                 vmax=cluster_df.max().max())),
                         cax=cbar_ax2, orientation='horizontal', fraction=0.2, pad=0.05)
    cbar1.set_label('Clustering Features', labelpad=-55, loc='center', fontsize=14)
    cbar1.ax.tick_params(labelsize=14)

    cbar2 = fig.colorbar(plt.cm.ScalarMappable(cmap='Blues', norm=plt.Normalize(vmin=boolean_df.min().min(),
                                                                                vmax=boolean_df.max().max())),
                         cax=cbar_ax3, orientation='horizontal', fraction=0.2, pad=0.35)
    cbar2.set_ticks([0, 1])
    cbar2.set_ticklabels(['0\n(Linear)', '1\n(Circular)'])
    cbar2.set_label('Boolean Metrics', labelpad=-70, loc='center', fontsize=14)
    cbar2.ax.tick_params(labelsize=14)

    cbar3 = fig.colorbar(plt.cm.ScalarMappable(cmap='Oranges', norm=plt.Normalize(vmin=size_df.min().min(),
                                                                                  vmax=size_df.max().max())),
                         cax=cbar_ax4, orientation='horizontal', fraction=0.2, pad=0.65)
    cbar3.set_label('Cycle Size', labelpad=-55, loc='center', fontsize=14)
    cbar3.ax.tick_params(labelsize=14)

    cbar4 = fig.colorbar(plt.cm.ScalarMappable(cmap='Purples', norm=plt.Normalize(vmin=structural_df.min().min(),
                                                                                  vmax=structural_df.max().max())),
                         cax=cbar_ax5, orientation='horizontal', fraction=0.2, pad=0.95)
    cbar4.set_label('Structural Metrics', labelpad=-55, loc='center', fontsize=14)
    cbar4.ax.tick_params(labelsize=14)

    plt.show()

    plt.savefig(cluster_plot_path, bbox_inches='tight', dpi=300)

    return


def post_clustering_hierarchical_selection(raw_df, selection_values, TPR_threshold, PFNR_threshold,
                                           TPR_scores, PFNR_scores, MEB_scores, ESB_scores,
                                           total_score_thresholds, weights):
    """
    post_clustering_hierarchical_selection()

    Inputs:
    -- raw_df (DataFrame) --> cycle dataframe with clustering results and all selection metrics
    -- selection_values (list) --> list of values to be used in selection process
    -- TPR_threshold (float) --> see ecDNAInspector manuscript figure 2 for details
    -- PFNR_threshold (float) --> see ecDNAInspector manuscript figure 2 for details
    -- TPR_scores (list) --> see ecDNAInspector manuscript figure 2 for details
    -- PFNR_scores (list) --> see ecDNAInspector manuscript figure 2 for details
    -- MEB_scores (list) --> see ecDNAInspector manuscript figure 2 for details
    -- ESB_scores (list) --> see ecDNAInspector manuscript figure 2 for details
    -- total_score_thresholds (list) --> see ecDNAInspector manuscript figure 2 for details
    -- weights (list) --> see ecDNAInspector manuscript figure 2 for details

    Output:
    -- raw_df (DataFrame) --> updated cycle dataframe with confidence assignments based on hierarchical selection
    """

    for value in selection_values:
        # subset to only include samples that have information on all selection values
        raw_df = raw_df[raw_df[value] != "NF"]
        # convert value to float
        raw_df[value] = raw_df[value].astype(float)

    print("Total cycles with complete selection information: " + str(len(raw_df)))

    def convert_TPR(TPR):
        if TPR == 0:
            return TPR_scores[0]
        elif TPR < TPR_threshold:
            return TPR_scores[1]
        elif TPR >= TPR_threshold:
            return TPR_scores[2]
        elif TPR == 1:
            return TPR_scores[3]

    def convert_MEB(boolean):
        if boolean == 0:
            return MEB_scores[0]
        elif boolean == 1:
            return MEB_scores[1]

    def convert_ESB(boolean):
        if boolean == 0:
            return ESB_scores[0]
        elif boolean == 1:
            return ESB_scores[1]

    def convert_PFNR(PFNR):
        if PFNR == 0:
            return PFNR_scores[0]
        elif PFNR < PFNR_threshold:
            return PFNR_scores[1]
        elif PFNR >= PFNR_threshold:
            return PFNR_scores[2]
        elif PFNR == 1:
            return PFNR_scores[3]

    raw_df["Selection Score"] = weights["TPR"] * raw_df["Cycle TPR"].apply(convert_TPR) + \
                                weights["MEB"] * raw_df["Cycle MEB"].apply(convert_MEB) + \
                                weights["ESB"] * raw_df["Cycle ESB"].apply(convert_ESB) + \
                                weights["PFNR"] * raw_df["Cycle PFNR"].apply(convert_PFNR)

    def label_score(score):
        if score < total_score_thresholds[0]: return "high"
        if total_score_thresholds[0] <= score <= total_score_thresholds[1]: return "medium"
        if total_score_thresholds[1] < score: return "low"

    raw_df["Confidence"] = raw_df["Selection Score"].apply(label_score)

    print("Total high confidence cycles: " + str(raw_df["Confidence"].value_counts().get('high', 0)))

    print("Total medium confidence cycles: " + str(raw_df["Confidence"].value_counts().get('medium', 0)))

    print("Total low confidence cycles: " + str(raw_df["Confidence"].value_counts().get('low', 0)))

    return raw_df
