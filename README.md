# Requirements
The pipeline is Linux compatible and it was tested under Ubuntu 18 (or newer). We recommend to run the pipeline under Docker in other Linux distributions.

# Workflow

Before executing, it is necessary to generate a STAR (2.6.1b) index directory from the annotation of interest (for more info visit https://github.com/alexdobin/STAR/blob/2.6.1b/doc/STARmanual.pdf).

The steps for CORALL analysis on demultiplexed FASTQ files are as follows:

1. UMI extraction
2. Adapter trimming (Cutadapt 2.5)
3. Alignment (STAR 2.6.1b)
4. UMI deduplicating (UMI-tools 1.0.0)
5. Transcript quantification (mix-square 1.4.0.2)

# Pipeline set-up

## Download Mix² RNA-Seq Data Analysis Software

Before you do any of the other steps for the pipeline set-up, you need to download Mix² from https://www.lexogen.com/mix-analysis-software/ and extract the binary with the name mix-square to the main folder of this repository.

## Using installation in local miniconda environment

If your miniconda environment is located in ~/miniconda3 you can, for instance, install the pipeline in your local environment with

`./install.sh ~/miniconda3`

and run, e.g. with

`./run_corall_local.sh -a gtf_file -o corall_out -s star_index_dir -m ~/miniconda3 rd1_fastq_gz [rd2_fastq_gz]`

or

`./run_corall_local.sh -a gtf_file -o corall_out -s star_index_dir -m ~/miniconda3 sample_csv`

The first version takes one gzipped fastq file for read 1 and possibly another one for read 2, while the second version takes a csv file as input which contains on each line the path of a gzipped fastq file for read 1 and possibly another one for read 2. If two files are present on a line they are separated by a comma.

## Using docker installation

For building the docker image use, for instance,

`./build_corall_dockerimg.sh corall:v1.0.0-customer`

This requires that environment.yml for conda installation, and the mix-square binary are in the current working directory. The command line above builds an image with tag corall:v1.0.1

For running the image use, e.g.

`./run_corall_docker.sh -a gtf_file -o ${input_dir}/corall_out -s star_index_dir -d corall:v1.0.0-customer ${input_dir}/sample_fastq.csv`

Note that ${input_dir} has to be an absolute path name.


In case you want make a backup of the build image and export it run

`docker save corall:v1.0.0-customer | gzip > corallWF_v1.0.0-customer.tar.gz`


## Using an existing docker image

In case you want use our pre-build docker image just import it with

`docker load < corallWF_v1.0.0-customer.tar.gz`

The pipeline can then be run in the same way as shown above e.g.

`./run_corall_docker.sh -a gtf_file -o ${input_dir}/corall_out -s star_index_dir -d corall:v1.0.0-customer ${input_dir}/sample_fastq.csv`

## Running local tests
In case you need to validate your set-up you can run

`./test_data/test_pipeline.sh`

from the cloned root directory. Conda needs to be accessible in order to run the tests (e.g. by running `conda init` or the conda folder needs to be included in PATH environment variable). 
