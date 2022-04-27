# Pipeline set-up

## Using installation in local miniconda environment

If your miniconda environment is located in ~/miniconda3 you can, for instance, install the pipeline in your local environment with

`./install.sh ~/miniconda3`

and run, e.g. with

`./run_corall_local.sh -o corall_out -m ~/miniconda3 <rd1_fastq_gz> [<rd2_fastq_gz>]`

or

`./run_corall_local.sh -o corall_out -m ~/miniconda3 <sample_csv>`

The first version takes one gzipped fastq file for read 1 and possibly another one for read 2, while the second version takes a csv file as input which contains on each line the path of a gzipped fastq file for read 1 and possibly another one for read 2. If two files are present on a line they are separated by a comma.

## Using docker installation

For building the docker image use, for instance,

`./build_corall_dockerimg.sh corall:v1.0.0-customer`

This requires that environment.yml for conda installation, and the mix-square binary are in the current working directory. The command line above builds an image with tag corall:v1.0.1

For running the image use, e.g.

`./run_corall_docker.sh -o ${input_dir}/corall_out -d corall:v1.0.0-customer ${input_dir}/sample_fastq.csv`

Note that ${input_dir} has to be an absolute path name.


In case you want make a backup of the build image and export it run

`docker save corall:v1.0.0-customer | gzip > corallWF_v1.0.0-customer.tar.gz`


## Using an existing docker image

In case you want use our pre-build docker image just import it with

`docker load < corallWF_v1.0.0-customer.tar.gz`

The pipeline can then be run in the same way as shown above e.g.

`./run_corall_docker.sh -o ${input_dir}/corall_out -d corall:v1.0.0-customer ${input_dir}/sample_fastq.csv`


