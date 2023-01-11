#!/bin/bash -e

usage() { echo "Usage: $0 [-a <gtf_file>] [-s <star_index] [-o <out_dir>] [-m <miniconda_install_dir>] <rd1_fastq_gz|sample_csv> [<rd2_fastq_gz>]" 1>&2; exit 1; }

out_dir=corall_out
miniconda_install_dir=~/miniconda3

while getopts ":a:s:o:m:" o; do
    case "${o}" in
        a)
            gtf_file=${OPTARG}
            ;;
        s)
            star_index=${OPTARG}
            ;;
        o)
            out_dir=${OPTARG}
            ;;
        m)
            miniconda_install_dir=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

[ $# -eq 0 ] && usage

if [[ $1 =~ \.fastq.gz$ ]]
then
	[ $# -gt 2 ] && usage
	if [ $# -eq 1 ]
	then
		corall_opt="--r1 $1"
	else
		corall_opt="--r1 $1 --r2 $2"
	fi
elif [[ $1 =~ \.csv$ ]]
then
	[ $# -gt 1 ] && usage
	corall_opt="--csv $1"
else
	usage
fi

source ${miniconda_install_dir}/bin/activate corall

./corall_batch.sh ${corall_opt} \
  --gtfFile ${gtf_file} \
  --starGenomeDir ${star_index} \
  --baseDir ${out_dir}

conda deactivate
