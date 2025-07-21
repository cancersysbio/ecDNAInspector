# ecDNAInspector_cycle_selection
"""
This file includes all steps needed to go from raw cycle prediction files to a csv file
with all calculated metrics and confidence assignments for each prediction.
"""

import ecDNAInspector.core_functions as ecI
from ecDNAInspector.utils import load_config
import pandas as pd
from itertools import combinations
import csv
import random


def cycle_data_calculations(config_file):

    # First unload parameters:

    config = load_config(config_file)

    # File inputs and conversions
    input_cycles_path = config['files']['input_cycles_path']
    amplicon_info = config['files']['amplicon_info']
    SV_path = config['files']['SV_path']
    processed_cycle_files = config['files']['processed_cycle_files']
    bed_files_path = config['files']['bed_files_path']
    gene_file = config['files']['gene_file']
    blacklist_file = config['files']['blacklist_file']

    # Data table inputs
    cycle_level_data_table = config['data']['cycle_level_data_table']

    # Analysis specifications
    cycle_type_to_include = config['analysis']['cycle_type_to_include']
    testing_mode = config['analysis']['testing_mode']
    skip_file_conversion = config['analysis']['skip_file_conversion']

    # Buffers, threshholds, and parameters
    gene_list = config['params']['gene_list']
    SV_in_range_buffer = config['params']['SV_in_range_buffer']
    blacklist_buffer = config['params']['blacklist_buffer']
    small_del_len = config['params']['small_del_len']
    gene_inclusion_prop_buffer = config['params']['gene_inclusion_prop_buffer']
    copy_count_threshold = config['params']['copy_count_threshold']

    # Convert cycle files
    sample_amp_dict = ecI.parse_amp_info_file_for_sample_amp_dict(amplicon_info)

    # Check for testing mode (if true, select random subset of samples - 10%)
    if testing_mode:
        num_samples = max(1, int(len(sample_amp_dict) * .10))
        sampled_keys = random.sample(list(sample_amp_dict.keys()), num_samples)
        subset_dict = {k: sample_amp_dict[k] for k in sampled_keys}
        print(f"[Testing Mode] Using a subset of {num_samples} amplicons out of {len(sample_amp_dict)}.")
        sample_amp_dict = subset_dict

    # Check if processed cycle files are already provided, otherwise proceed with file conversion.
    if not skip_file_conversion:
        ecI.convert_cycle_file_set(sample_amp_dict, input_cycles_path, processed_cycle_files,
                                   copy_count_threshold, cycle_type_to_include)

    # Load blacklist, SV, and gene data
    blacklist_dict = ecI.parse_blacklist_file_for_blacklist_dict(blacklist_file)

    svs_dict = ecI.parse_SV_file_for_SVs(SV_path)

    if gene_list == "Cosmic":
        cosmic_genes_dict = ecI.parse_cosmic_for_gene_dict(gene_file)

    if gene_list == "Reference":
        ref_genes_dict = ecI.parse_reference_for_genome_dict(gene_file)

    output_file_data = [["Sample", "Amplicon", "Cycle", "Cycle Type", "Cycle Size (bp)", "Cycle Breakpoint Count",
                         "Cycle Chromosome Count", "Cycle Regions List", "Cycle Paired Ends List",
                         "Cycle TPR", "Cycle TP BEs", "Cycle FPR", "Cycle FP BEs", "Cycle PFNR", "Cycle PFN BEs",
                         "Cycle MEB", "Cycle Mapping Error Breakends", "Cycle Unmappable Reasons",
                         "Small Deletions Count", "Cycle Genes Info", "Cycle Genes List"]]

    # iterate over every sample, amp, cycle
    print("")
    for sample in sample_amp_dict:
        print("Working on sample: " + str(sample))
        for amp in sample_amp_dict[sample]:
            cycle_file = processed_cycle_files + sample + "_processed_amplicon" + amp + "_cycles.txt"
            cycle_nums = ecI.parse_cycle_file_for_cycle_nums(cycle_file)
            for cycle in cycle_nums:
                cycle_type = ecI.parse_cycle_file_for_cycle_type(cycle_file, cycle)
                cycle_regions_dict = ecI.parse_cycle_file_for_cycle_region_dict(cycle_file, cycle)
                cycle_regions_list = ecI.parse_cycle_file_for_cycle_region_list(cycle_file, cycle)
                cycle_paired_ends = ecI.parse_cycle_file_for_cycle_paired_ends_list(cycle_file, cycle)
                cycle_size = ecI.calc_cycle_size(cycle_regions_list)
                cycle_num_bps = ecI.calc_cycle_num_breakpoints(cycle_paired_ends)
                cycle_num_chroms = ecI.calc_cycle_num_chroms(cycle_regions_list)
                cycle_TPR, TP_list, cycle_FPR, FP_list, small_del_removed = ecI.calc_TPR_and_FPR(cycle_paired_ends,
                                                                                                 svs_dict, sample,
                                                                                                 SV_in_range_buffer,
                                                                                                 small_del_len)
                cycle_PFNR, PFN_list = ecI.calc_PFNR(cycle_paired_ends, small_del_len, bed_files_path, sample,
                                                 svs_dict, SV_in_range_buffer, TP_list, copy_count_threshold)
                cycle_MEB, cycle_mapping_errors, cycle_unmappable_reasons = ecI.calc_mapping_errors(cycle_paired_ends,
                                                                                                    blacklist_dict,
                                                                                                    blacklist_buffer)
                if gene_list == "Cosmic":
                    cycle_genes_dict, cycle_genes_list = ecI.parse_cosmic_dict_for_genes_in_cycle(cycle_regions_dict,
                                                                                                  cosmic_genes_dict,
                                                                                                  gene_inclusion_prop_buffer)
                if gene_list == "Reference":
                    cycle_genes_dict, cycle_genes_list = ecI.parse_genome_dict_for_features_in_cycle(cycle_regions_dict,
                                                                                                     ref_genes_dict,
                                                                                                     gene_inclusion_prop_buffer)
                cycle_info_list = [sample, amp, cycle, cycle_type, cycle_size, cycle_num_bps, cycle_num_chroms,
                                   cycle_regions_list, cycle_paired_ends, cycle_TPR, TP_list, cycle_FPR,
                                   FP_list, cycle_PFNR, PFN_list, cycle_MEB, cycle_mapping_errors,
                                   cycle_unmappable_reasons, small_del_removed,
                                   cycle_genes_dict, cycle_genes_list]
                output_file_data.append(cycle_info_list)

    with open(cycle_level_data_table, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(output_file_data)

    cohort_df = pd.read_csv(cycle_level_data_table)

    # update dataframe with extreme size (ESB) boolean information

    # The cohort percentiles below and above which (respective) cycle size should be labeled as "extreme".
    # We recommend adjusting these cutoffs according to the distribution of your dataset.
    esb_percentiles = (0.1, 0.9)

    lower_size_bound = cohort_df["Cycle Size (bp)"].quantile(esb_percentiles[0])
    upper_size_bound = cohort_df["Cycle Size (bp)"].quantile(esb_percentiles[1])

    cohort_df["Cycle ESB"] = ((cohort_df["Cycle Size (bp)"] < lower_size_bound) |
                              (cohort_df["Cycle Size (bp)"] > upper_size_bound)).astype(int)

    # re-write data to csv file so ESB can be included in clustering visualization

    cohort_df.to_csv(cycle_level_data_table, index=False)


def cycle_clustering(config_file):

    config = load_config(config_file)

    # Data table inputs
    cycle_level_data_table = config['data']['cycle_level_data_table']
    cycle_level_data_clustered_table = config['data']['cycle_level_data_clustered_table']

    # Buffers, threshholds, and parameters
    min_clusters = config['params']['min_clusters']
    max_clusters = config['params']['max_clusters']
    cluster_num_resamples = config['params']['cluster_num_resamples']
    cluster_resample_prop = config['params']['cluster_resample_prop']
    n_init = config['params']['n_init']

    # define clustering features
    features = ["Cycle TPR", "Cycle FPR", "Cycle PFNR", "Cycle MEB"]

    # cluster
    ecI.cluster_with_consensus(cycle_level_data_table, cycle_level_data_clustered_table, features, min_clusters,
                               max_clusters, cluster_num_resamples, cluster_resample_prop, n_init)


def visualize_clustering(clustered_table, cluster_image_output):

    cohort_w_clust_df = pd.read_csv(clustered_table)

    cohort_w_clust_df["Cycle Type Boolean"] = (cohort_w_clust_df["Cycle Type"] == "circular").astype(int)
    cohort_w_clust_df["Cycle MEB"] = cohort_w_clust_df["Cycle MEB"].astype(int)

    cohort_w_clust_df['Cluster'] = cohort_w_clust_df['Cycle Cluster Assignment']
    cohort_w_clust_df = cohort_w_clust_df.drop(columns=['Cycle Cluster Assignment'])
    cohort_w_clust_df = cohort_w_clust_df.rename(columns={"Cluster": "Cycle Cluster Assignment"})

    clustering_features = ["Cycle TPR", "Cycle FPR", "Cycle PFNR", "Cycle MEB"]
    boolean_metrics = ["Cycle Type Boolean", "Cycle ESB"]
    size_metric = ["Cycle Size (bp)"]
    structural_metrics = ["Cycle Breakpoint Count", "Cycle Chromosome Count"]

    cluster_feature_labels = ["TPR", "FPR", "PFNR", "MEB"]
    boolean_metric_labels = ["Cycle Type", "ESB"]
    size_metric_labels = ["Cycle Size (bp)"]
    structural_metric_labels = ["BP Count", "Chrom. Count"]

    sorting_features = ['Cycle Cluster Assignment', 'Cycle Type Boolean', 'Cycle TPR']
    sorting_features_dir = [True, False, False]

    ecI.plot_separated_clustering_heatmap(cohort_w_clust_df,
                                          clustering_features, boolean_metrics, size_metric, structural_metrics,
                                          cluster_feature_labels, boolean_metric_labels, size_metric_labels,
                                          structural_metric_labels,
                                          sorting_features, sorting_features_dir, cluster_image_output)


def assign_cycle_confidence_by_cluster(config_file, hc_clusters, mc_clusters, lc_clusters):

    config = load_config(config_file)

    # Data table inputs
    cycle_level_data_clustered_table = config['data']['cycle_level_data_clustered_table']
    cycle_level_data_w_conf_table = config['data']['cycle_level_data_w_conf_table']

    cohort_w_clust_df = pd.read_csv(cycle_level_data_clustered_table)

    def assign_confidence_by_cluster(cluster):
        if str(cluster) in hc_clusters:
            return "high"
        elif str(cluster) in mc_clusters:
            return "medium"
        elif str(cluster) in lc_clusters:
            return "low"

    cohort_w_clust_df["Confidence"] = cohort_w_clust_df["Cycle Cluster Assignment"].apply(assign_confidence_by_cluster)
    cohort_w_clust_df.to_csv(cycle_level_data_w_conf_table, index=False)

    print("Total high confidence cycles: " + str(cohort_w_clust_df["Confidence"].value_counts().get('high', 0)))

    print("Total medium confidence cycles: " + str(cohort_w_clust_df["Confidence"].value_counts().get('medium', 0)))

    print("Total low confidence cycles: " + str(cohort_w_clust_df["Confidence"].value_counts().get('low', 0)))


def assign_cycle_confidence_by_selection(config_file):

    config = load_config(config_file)

    # Data table inputs
    cycle_level_data_clustered_table = config['data']['cycle_level_data_clustered_table']
    cycle_level_data_w_conf_table = config['data']['cycle_level_data_w_conf_table']

    # Hierarchical selection thresholds and weights
    TPR_threshold = config['params']['TPR_threshold']
    PFNR_threshold = config['params']['PFNR_threshold']
    TPR_scores = config['params']['TPR_scores']
    PFNR_scores = config['params']['PFNR_scores']
    MEB_scores = config['params']['MEB_scores']
    ESB_scores = config['params']['ESB_scores']
    weights = config['params']['weights']
    total_score_thresholds = config['params']['total_score_thresholds']

    selection_values = ["Cycle TPR", "Cycle MEB", "Cycle ESB", "Cycle PFNR"]

    cohort_w_clust_df = pd.read_csv(cycle_level_data_clustered_table)

    cohort_w_conf_df = ecI.post_clustering_hierarchical_selection(cohort_w_clust_df, selection_values,
                                                                  TPR_threshold, PFNR_threshold,
                                                                  TPR_scores, PFNR_scores, MEB_scores, ESB_scores,
                                                                  total_score_thresholds, weights)

    cohort_w_conf_df.to_csv(cycle_level_data_w_conf_table, index=False)


def intrasample_filtering(config_file):
    config = load_config(config_file)

    # File and data table inputs
    input_cycles_path = config['files']['input_cycles_path']
    intrasample_filt_cycle_level_data_w_conf_table = config['data']['intrasample_filt_cycle_level_data_w_conf_table']
    cycle_level_data_table = config['data']['cycle_level_data_table']
    cycle_level_data_w_conf_table = config['data']['cycle_level_data_w_conf_table']

    # Params
    threshold_value = config['params']['intrasample_filtering_threshold_value']

    cohort_df = pd.read_csv(cycle_level_data_table)
    cohort_w_conf_df = pd.read_csv(cycle_level_data_w_conf_table)

    intra_sample_jaccards_results = []

    for sample, group in cohort_df.groupby('Sample'):
        pairs = list(combinations(group.itertuples(index=False), 2))

        for (row1, row2) in pairs:

            row1_cycle_sample_name = row1._asdict()['Sample']
            row1_cycle_amp = row1._asdict()['Amplicon']
            row1_cycle_num = row1._asdict()['Cycle']
            row1_cycle_file = input_cycles_path + row1_cycle_sample_name + "_amplicon" + str(row1_cycle_amp) + "_cycles.txt"
            row1_cycle_regions =ecI.parse_cycle_file_for_cycle_region_dict(row1_cycle_file, row1_cycle_num)

            row2_cycle_sample_name = row2._asdict()['Sample']
            row2_cycle_amp = row2._asdict()['Amplicon']
            row2_cycle_num = row2._asdict()['Cycle']
            row2_cycle_file = input_cycles_path + row2_cycle_sample_name + "_amplicon" + str(row2_cycle_amp) + "_cycles.txt"
            row2_cycle_regions = ecI.parse_cycle_file_for_cycle_region_dict(row2_cycle_file, row2_cycle_num)

            jaccard = ecI.calc_bp_jaccard(row1_cycle_regions, row2_cycle_regions)

            intra_sample_jaccards_results.append({
                'Sample': sample,
                'Amplicon1': row1.Amplicon,
                'Cycle1': row1.Cycle,
                'Amplicon2': row2.Amplicon,
                'Cycle2': row2.Cycle,
                'Jaccard': jaccard
            })

    intrasample_jaccards_df = pd.DataFrame(intra_sample_jaccards_results)

    intrasample_to_filter_df = intrasample_jaccards_df[intrasample_jaccards_df["Jaccard"] > threshold_value]

    def filter_sample(sample_df):
        to_remove = []

        # for each unique sample, pull the pairwise comparisons with a Jaccard > threshold
        sample_name = sample_df['Sample'].iloc[0]
        sample_pairs = intrasample_to_filter_df[intrasample_to_filter_df['Sample'] == sample_name]

        for _, pair in sample_pairs.iterrows():

            amplicon_1, cycle_1 = pair['Amplicon1'], pair['Cycle1']
            amplicon_2, cycle_2 = pair['Amplicon2'], pair['Cycle2']

            pair_rows_1 = sample_df[(sample_df['Amplicon'] == amplicon_1) &
                                    (sample_df['Cycle'] == cycle_1)]
            pair_rows_2 = sample_df[(sample_df['Amplicon'] == amplicon_2) &
                                    (sample_df['Cycle'] == cycle_2)]

            if not pair_rows_1.empty and not pair_rows_2.empty:

                # remove the cycle with the lower TPR
                tpr_1, tpr_2 = pair_rows_1['Cycle TPR'].iloc[0], pair_rows_2['Cycle TPR'].iloc[0]
                if tpr_1 > tpr_2:
                    to_remove.append(pair_rows_2)
                else:
                    to_remove.append(pair_rows_1)

        for rows in to_remove:
            sample_df = sample_df[~sample_df.index.isin(rows.index)]

        return sample_df

    intrasample_filtered_cycle_data_df = pd.concat([filter_sample(sample_df)
                                                    for _, sample_df in cohort_w_conf_df.groupby('Sample')])

    intrasample_filtered_cycle_data_df.to_csv(intrasample_filt_cycle_level_data_w_conf_table, index=False)


