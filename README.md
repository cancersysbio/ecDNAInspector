# ecDNAInspector

ecDNAInspector (ecI) provides a systematic method for ecDNA structural validation and analysis, enabling higher-confidence utilization of sequencing data-based predictions in studies of ecDNA. ecI can be deployed downstream of other ecDNA callers (e.g. [AmpliconArchitect [AA]]([url](https://github.com/virajbdeshpande/AmpliconArchitect)), [JaBbA]([url](https://github.com/mskilab-org/JaBbA))) - the current version offers automatic file processing for AA output. ecI can be run through the command line or through a series of Jupyter notebooks (please see "Usage").

## Installation (for command line use)

```bash
# clone the repo
git clone gh:cancersysbio/ecDNAInspector
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

To run this version, please navigate to the "ecDNAInspector_Jupyter_notebooks" folder and download all notebooks and files. Update paths within the Jupyter notebook before running. Start by running "ecDNAInspector_cycle_selection.ipynb", then "ecDNAInspector_analysis.ipynb". Run "ecDNAInspector_Jaccard_calculations.ipynb" before the analysis notebook if you wish to make structural conservation comparisons (note that this step may take a long time to complete, O(n^2)). 

### Command Line Version

Once you are comfortable with your selected parameters, you can run ecI in a more streamlined fashion using the command line interface version. ecI is written as a package so you can import it to any code you like, or run directly from the command line. 

Each command line run requires including the path to a config file. A default file is provided ("default_congif.yaml"), or you may copy and update your own. Please update all file paths within the config according to your directory organization. Please see more information below about the parameters and recommended values, or read our paper (to be linked). 

```bash
# start by running help to view all the commands
ecI --help

# run the cycle metric calculations
ecI --config default_config.yaml --run-metric-calc

# run the cycle clustering
ecI --config default_config.yaml --run-cluster

# visualize your cycle clustering
ecI --config default_config.yaml --visualize_cluster

# run the cycle confidence assignments, specifying how you want confidence assignments done
# OPTION 1
ecI --config default_config.yaml --run-conf-assignment -conf_type by_selection
# OPTION 2 (replace numbers following high/med/low_conf_clusters with the cluster numbers from your clustering results. You can type multiple numbers separated by a space, or exclude a category completely by removing the flag)
ecI --config default_config.yaml --run-conf-assignment -conf_type by_cluster --high_conf_clusters 0 --med_conf_clusters 1 --low_conf_clusters 2

# OR, run multiple steps at once!

ecI --config default_config.yaml --run-metric-calc --run-cluster --visualize_cluster --run-conf-assignment -conf_type by_selection
```

**Important!** Some steps MUST be run before others. Cycle metric calculations must be done before clustering and confidence assignments. You do not need to cluster before confidence assignments if you use the "by_selection" confidence assignment method; if using the "by_cluster" method, you must cluster first (and we recommend visualizing the clusters to help select your high, medium, and low confidence cluster(s)). 

