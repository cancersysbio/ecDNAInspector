# ecDNAInspector

ecDNAInspector (ecI) provides a systematic method for ecDNA structural validation and analysis, enabling higher-confidence utilization of sequencing data-based predictions in studies of ecDNA. ecI can be deployed downstream of other ecDNA callers (e.g. [AmpliconArchitect (AA)](https://github.com/virajbdeshpande/AmpliconArchitect), [JaBbA](https://github.com/mskilab-org/JaBbA)) - the current version offers automatic file processing for AA output. ecI can be run through the command line or through a series of Jupyter notebooks (please see "Usage"). Ultimately, ecI takes as input a set of ecDNA cycle predictions and provides a confidence evaluation of each prediction. You can then select a subset of your cycles based on these evaluations for downstream analysis.

## Installation (for command line use)

```bash
# clone the repo
# with https:
git clone https://github.com/cancersysbio/ecDNAInspector.git

# with ssh:
git clone git@github.com:cancersysbio/ecDNAInspector.git

# navigate into the tool directory
cd ecDNAInspector

# create the environment from the YAML file (if desired, or create your own environment)
conda env create -f environment.yml

# activate the environment (note: ecI_env is the name from the environment.yml file; if you want a different name, please change it in the file and in the command below).
conda activate ecI_env

# finally, install the tool!
pip install -e .
```

## Usage

Running ecI requires for each sample: cycle prediction file(s), a .bed coverage file, and a list of SV calls. These files should be provided in a standard format (please see the provided test data folder to inspect the formats and ensure your files conform). If using AA output, you can proceed directly with the raw cycle prediction output files and ecI will automatically convert them to the processed cycle file format. If not using AA output, please manually convert your cycle prediction files to conform with the processed cycle file format. 

### Jupyter notebook version (recommended for initial use)

We recommend using the Jupyter notebook version of ecI for initial runs because you are able to pause and visualize results after steps more easily. As there are numerous parameters that must be carefully considered for ecI to work effectively, we recommend first running through Jupyter notebook and testing a variety of parameters. 

To run this version, please navigate to the "ecDNAInspector_Jupyter_notebooks" folder and download all notebooks and files. Update paths within the Jupyter notebook before running. Start by running "ecDNAInspector_cycle_selection.ipynb", then "ecDNAInspector_analysis.ipynb". Run "ecDNAInspector_Jaccard_calculations.ipynb" before the analysis notebook if you wish to make structural conservation comparisons (note that this step may take a long time to complete, O(n^2)). Please see the manuscript for more details on the different notebooks and their functionalities.

### Command Line Version

Once you are comfortable with your selected parameters, you can run ecI in a more streamlined fashion using the command line interface version. ecI is written as a package so you can import it to any code you like, or run directly from the command line. 

Each command line run requires including the path to a config file. A default file is provided ("default_congif.yaml"), or you may copy and update your own. Please update all file paths within the config according to your directory organization. Please see more information below about the parameters and recommended values, or read our paper (to be linked). 

```bash
# start by running help to view all the commands
ecI --help

# run the cycle metric calculations
ecI --config configs/default_config.yaml --run-metric-calc

# run the cycle clustering
ecI --config configs/default_config.yaml --run-cluster

# visualize your cycle clustering (start with unfiltered)
ecI --config configs/default_config.yaml --visualize_cluster unfiltered

# run the cycle confidence assignments, specifying how you want confidence assignments done
# OPTION 1
ecI --config configs/default_config.yaml --run-conf-assignment -conf_type by_selection
# OPTION 2 (replace numbers following high/med/low_conf_clusters with the cluster numbers from your clustering results. You can type multiple numbers separated by a space, or exclude a category completely by removing the flag)
ecI --config configs/default_config.yaml --run-conf-assignment -conf_type by_cluster --high_conf_clusters 0 --med_conf_clusters 1 --low_conf_clusters 2

# optionally, filter the cycles to remove highly similar cycles in the same sample (note: this step can only be run after confidence has been assigned)
ecI --config configs/default_config.yaml --run-intrasample-filter

# after running an intrasample filtering step, you can re-visualize the same cycle clustering with only the cycles remaining after intrasample filtering
ecI --config configs/default_config.yaml --visualize_cluster intrasample_filtered

# you can also run multiple steps at once!

ecI --config configs/default_config.yaml --run-metric-calc --run-cluster --visualize_cluster unfiltered --run-conf-assignment -conf_type by_selection

ecI --config configs/default_config.yaml --run-metric-calc --run-cluster --run-conf-assignment -conf_type by_selection --run-intrasample-filter --visualize_cluster intrasample_filtered
```

**Important!** Some steps MUST be run before others. Cycle metric calculations must be done before clustering and confidence assignments. You do not need to cluster before confidence assignments if you use the "by_selection" confidence assignment method; if using the "by_cluster" method, you must cluster first (and we recommend visualizing the clusters to help select your high, medium, and low confidence cluster(s)). 

Here is a full table of the flags and their descriptions:

| Flag | Description | Notes | Default/options | 
|----------|----------|----------|----------|
| -h, --help |  Displays all flags and descriptions |   |   |
| -c, --config |  Path to YAML configuration file | **Required.** You can use the default config file in the configs folder, or modify the user config file in the configs folder.  |   |
| --run-metric-calc |  Run cycle metric calculations | Must be run before --run-cluster, --run-conf-assignment, or --run-intrasample-filter. Can be run independently or together with these downstream steps. |   |
| --run-cluster |  Run cycle clustering | Must be run together with or after --run-metric-calc. Must be run before --run-conf-assignment. |   |
| --run-conf-assignment |  Run cycle confidence assignment | Must be run together with or after --run-metric-calc, --run-cluster. Must be run before --run-intrasample-filter. Must be run together with -conf_type. |   |
| --run-intrasample-filter |  Run cycle intrasample filtering | Must be run together with or after --run-metric-calc, --run-cluster, --run-conf-assignment.  |   |
| --visualize-cluster |  Run cycle clustering visualization and specify which cycle cluster data to use. | Must be run together with or after --run-metric-calc, --run-cluster. To use full original cycle set, run with unfiltered mode. To use the intrasample-filtered cycle set, run with intrasample_filtered mode.  |  unfiltered, intrasample_filtered |
| -conf_type, --confidence_assignment_type |  Type of cycle confidence assignment process. | Must be run together with --run-conf-assignment. To assign confidence by cluster, run with by_cluster mode. To assign confidence with the hierarchical selection method, run with by_selection mode. |  by_selection, by_cluster |
| --high_conf_clusters |  List of manually selected high confidence clusters. | Can be included with --run-conf-assignment -conf_type by_cluster. Exclude if you do not want to assign any clusters as high confidence. To assign 1+ clusters as high confidence, type the cluster numbers (0-indexed) after the flag with spaces in-between if multiple are provided.  |   |
| --med_conf_clusters |  List of manually selected medium confidence clusters. | Can be included with --run-conf-assignment -conf_type by_cluster. Exclude if you do not want to assign any clusters as medium confidence. To assign 1+ clusters as medium confidence, type the cluster numbers (0-indexed) after the flag with spaces in-between if multiple are provided.  |   |
| --low_conf_clusters |  List of manually selected low confidence clusters. | Can be included with --run-conf-assignment -conf_type by_cluster. Exclude if you do not want to assign any clusters as low confidence. To assign 1+ clusters as low confidence, type the cluster numbers (0-indexed) after the flag with spaces in-between if multiple are provided.  |   |


## Parameters

ecI uses numerous parameters, which are organized in a config file. A default config file is provided, but the file paths must be updated to reflect your directory organization. A blank user config file is also provided for your convenience. Below is a table of each parameter with descriptions. 

**Important!** **Your choice of input buffers, threshholds, and parameters matters!** Please read the ecDNAInspector manuscript for more information. Below we provide suggestions for each metric, but as each dataset is highly variable, we recommend trying multiple combinations of metrics to assess how decisions affect output and may need to be tailored to your specific dataset. This is most conveniently done through the Jupyter notebook version of the tool. 

| Parameter name | Description/options | Notes | Default | 
|----------|----------|----------|----------|
| **Input files** |   |   |   |
| input_cycles_path | path to directory with raw cycle files (from AA)  | **Optional.** Do not use if you are not starting from AA cycle files.  |   |
| amplicon_info    | path to file with information on sample name & amplicon number  | **Required.** See provided test amplicon file for required format.  |   |
| SV_path    | path to file with SV calls for samples  | **Required.** See provided test SV file for required format. **Please note: the quality of your structural variant calls matters!** The confidence assignment of each ecDNA is extremely depending on your selection of SVs, so the quality of your "ground truth" SV list is very important. We recommend using the consensus of at least two different callers to produce a high quality list (see the ecDNAInspector manuscript for additional details; our original study used the consensus of four callers.) |   |
| processed_cycle_files    | path to directory for processed cycle files  | **Required.** Please make sure this directory is established before running the cycle metric calculations. If not using AA for initial cycle prediction, please store your cycle predictions here in the default format. See provided test processed cycle files for required format. |  |
| bed_files_path    | path to directory with .bed coverage files for samples | **Optional, but required for pFNR calculation.** If not found, pFNR calculations will fail and pFNR must be manually excluded as a metric for any downstream confidence assignments. |   |
| gene_file    | path to file with gene location information  | **Required.** We provide functionality for using the Cosmic cancer gene list and the reference gene list - please replace this value with your desired list from the corresponding assembly. Please confirm the format matches the provided test gene files. Some options include: [Cosmic hg19](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v102/cancergenecensus), [Cosmic hg38](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v102/cancergenecensus), [refgene hg37](https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/)), [refgene hg38](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/) |   |
| blacklist_file    | path to file with blacklist information  | **Required.** Files for different assemblies can be found [here, organized by the Boyle Lab.](https://github.com/Boyle-Lab/Blacklist/tree/master/lists)  |   |
| unfiltered_cluster_image_output    | path to file where unfiltered clustering image will be stored  | **Required.**  |   |
| intrasamp_filtered_cluster_image_output    | path to file where intrasample-filtered clustering image will be stored  | **Required.**  |   |
| **Output data tables**    |   |   |   |
| cycle_level_data_table    | path to file with all cycle metric information (no cluster or confidence information)  | **Required.**  |   |
| cycle_level_data_clustered_table    | path to file with all cycle metric information and cluster assignment (no confidence information)  | **Required.**  |   |
| cycle_level_data_w_conf_table    | path to file with all cycle metric information, cluster assignment, and confidence assignment | **Required.**  |   |
| intrasample_filt_cycle_level_data_w_conf_table    | path to file with intrasample filtered cycle metric information, cluster assignment, and confidence assignment  | **Optional, but required for intrasample filtering step.**  |   |
| **Analysis specifications**    |   |  |  |
| cycle_type_to_include  | "all", "circular", or "linear"  | **Optional, specifically for AA input.** AA distinguishes between cycles with a complete circular contig ("circular") and those missing a read for one paired breakend connection ("linear"). Choose "all" to include all cycles regardless of this distinction, "circular" if you only want to analyze circular cycles, and "linear" if you only want to analyze linear cycles. Note that while circular cycles trend higher in confidence, some linear cycles may be higher confidence than circular ones, so we recommend investigations all cycles to start. Please see the manuscript for further details on this distinction.  | "all" |
| testing_mode    |  True/False (boolean) | Choose True to test the tool on a smaller fraction of your dataset (10%). We recommend initially running in testing mode as some steps are lengthy and testing mode enables quicker debugging of any issues. | False |
| skip_file_conversion    |  True/False (boolean) | Choose True if not using AA input files (i.e., if you need to start directly from your manually processed cycle files). | False |
| **Buffers, threshholds, and parameters**    |   |  |  |
| gene_list    | "Cosmic", "Reference"  | Update the gene_file path accordingly!  | "Cosmic" |
| be_overlap_buffer    | breakend overlap buffer; maximum distance between breakends for them to be considered overlapping  | a larger number will be more lenient in comparing breakends, a smaller number more strict (i.e. at 0, the breakends must be exactly the same)  | 100 |
| SV_in_range_buffer    | SV/BE overlap buffer; maximum distance between an SV breakend and an ecDNA breakend for them to be considered overlapping  |  a larger number will be more lenient (expect higher TPRs), a smaller number more strict (expect lower TPRs) | 100 |
| blacklist_buffer    | Blacklist/BE overlap buffer; maximum distance between blacklist region and an ecDNA breakend for the ecDNA breakend to be considered part of a blacklist region  |  a larger number will be more strict (expect more MEB=1), a smaller number more lenient (expect more MEB=0) | 100 |
| small_del_len    | small deletion length; maximum length of a "deletion"-type paired breakend to be excluded from SV validations; i.e., if an ecDNA cycle contained a paired breakend where the breakends are on the same chromosome within the small_del_len, this paired breakend will not be included in SV validation as SV callers generally perform poorly on small deletions. |  a larger number will call more paired breakends as small deletions and may inflate the TPRs, a smaller number will call fewer paired breakends as small deletions and may deflate the TPRs. | 100 |
| gene_inclusion_prop_buffer    | proportion of total gene length that must be present in ecDNA cycle for full gene to be considered included in the cycle |  e.g., when gene_inclusion_prop_buffer = 1.0, the full gene length must be present in the ecDNA cycle for the full gene to be considered included. when gene_inclusion_prop_buffer = 0.5, only half the gene length must be present. | 1.0 |
| copy_count_threshold    | minimum copy count of cycle to be included in analysis, as well as minimum coverage for regions in .bed file to search for pFN breakpoints. |  Higher numbers are stricter (expect lower ecDNA counts if using AA input and reduced pFNR), lower numbers are more lenient (expect higher ecDNA counts if using AA input and increased pFNR). If using AA input, we recommend using a copy number of at least 4.0 to avoid pulling a high quantity of lower-confidence cycle predictions with low copy number (recall ecDNA are expected to have high copy number).  | 4.0 |
| duplicate_buffer    | maximum distance between matched ends of paired ends for paired ends to be considered duplicates. |  It is possible that a single cycle can contain copies of the same segment, with multiple copies of the same/very similar paired end pair. To avoid inflating cycle quality metrics, these copies are excluded from TPR/FPR/pFNR calculations. Higher numbers are more lenient in calling duplicates (expect more duplicates), lower numbers are more strict (expect fewer duplicates).  | 100 |
| min_clusters    | minimum number of clusters to attempt in clustering |  We recommend using at least 3 clusters (to assign a high, medium, and low confidence group from the clusters), but some smaller datasets may not recognize 3.  | 3 |
| max_clusters    | maximum number of clusters to attempt in clustering |   | 5 |
| cluster_num_resamples    | number of resamples for clustering |  Required by consensusClustering method. [See github page for further details.](https://github.com/ZigaSajovic/Consensus_Clustering) | 1000 |
| cluster_resample_prop    | proportion of samples to resample |  Required by consensusClustering method. [See github page for further details.](https://github.com/ZigaSajovic/Consensus_Clustering) | 0.8 |
| n_init    | number of KMeans runs with different centroid seeds |  A higher number will take more time, but likely lead to a better clustering [See github page for further details.](https://github.com/ZigaSajovic/Consensus_Clustering) | 10 |
| TPR_threshold    | hierarchical selection TPR threshold |  A higher threshold is stricter (expect higher total scores, fewer high confidence cycles). See manuscript for further details on selecting the threshold. | 0.5 |
| PFNR_threshold    | hierarchical selection pFNR threshold |  A higher threshold is more lenient (expect lower total scores, more high confidence cycles). See manuscript for further details on selecting the threshold. | 0.5 |
| TPR_scores    | hierarchical selection score values for TPR |  Higher values are stricter (expect higher total scores, fewer high confidence cycles). Greater distinction between the scores will separate cycles more extremely. Values must be provided as a list of 4 integers/floats. | [1.5, 1.0, 0.5, 0] |
| pFNR_scores    | hierarchical selection score values for pFNR |  Higher values are more lenient (expect lower total scores, more high confidence cycles). Greater distinction between the scores will separate cycles more extremely. Values must be provided as a list of 4 integers/floats. | [1.5, 1.0, 0.5, 0] |
| MEB_scores    | hierarchical selection score values for MEB |  Higher values are more strict (expect higher total scores, fewer high confidence cycles). Value must be provided as a list of 2 integers/floats. | [0, 0.75] |
| ESB_scores    | hierarchical selection score values for ESB |  Higher values are more strict (expect higher total scores, fewer high confidence cycles). Value must be provided as a list of 2 integers/floats. | [0, 0.75] |
| weights    | hierarchical selection weights for metric scores |  Higher values are more strict (expect higher total scores, fewer high confidence cycles). Set value for metric to 0 to exclude a metric from the total score calculation. Value must be provided as a dictionary with the metric as the key and the value as the weight. |  {"TPR": 1, "MEB": 0.75, "ESB": 0.5, "PFNR": 0.25} |
| total_score_thresholds    | hierarchical selection total score thresholds |  Thresholds should be in range of the total scores across the cohort; we recommend centering them around the median (i.e., the range between the thresholds should align with the median of the range of total scores). Value must be provided as a list of 2 integers/floats; cycles with a total score less than the first value will be labeled as high confidence, cycles with a total score greater than the second value will be labeled as low confidence, cycles with a total score between the two values will be labeled as medium confidence.  |  [1, 2] |
| intrasample_filtering_threshold_value    | minimum basepair Jaccard value between two cycles from the same sample for those cycles to be considered "similar" and the lower confidence cycle removed |  **Optional, but required for intrasample filtering step.** We recommend using our analysis Jupyter notebook module to explore the intrasample similarity before selecting the cutoff Jaccard value. Please see our manuscript for more details on choosing intrasample filtering.  |  [1, 2] |
