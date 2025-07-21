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
```

## Usage

Running ecI requires for each sample: cycle prediction file(s), a .bed coverage file, and a list of SV calls. These files should be provided in a standard format (please see the provided test data folder to inspect the formats and ensure your files conform). If using AA output, you can proceed directly with the raw cycle prediction output files and ecI will automatically convert them to the processed cycle file format. If not using AA output, please manually convert your cycle prediction files to conform with the processed cycle file format. 

### Jupyter notebook version (recommended for initial use)

We recommend using the Jupyter notebook version of ecI for initial runs because you are able to pause and visualize results after steps more easily. As there are numerous parameters that must be carefully considered for ecI to work effectively, we recommend first running through Jupyter notebook and testing a variety of parameters. 

To run this version, please navigate to the "ecDNAInspector_Jupyter_notebooks" folder and download all notebooks and files. Update paths within the Jupyter notebook before running.

### Command Line Version

Once you are comfortable with your selected parameters, you can run ecI in a more streamlined fashion using the command line interface version. ecI is written as a package so you can import it to any code you like, or run directly from the command line as follows:

```bash



```



