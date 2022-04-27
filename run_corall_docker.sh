#!/bin/bash

usage() { echo "Usage: $0  [-a <gtf_file>] [-s <star_index] [-o <output_dir>] [-d <docker_image>] <rd1_fastq_gz|sample_csv> [<rd2_fastq_gz>]" 1>&2; exit 1; }


output_dir=corall_out
docker_img=corall:v1.0.1

while getopts ":a:s:o:d:" o; do
    case "${o}" in
        a)
            gtf_file=${OPTARG}
            ;;
        s)
            star_index=${OPTARG}
            ;;
        o)
            output_dir=${OPTARG}
            ;;       
        d)
            docker_img=${OPTARG}
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
	#use dirname twice, don't want to write in "fastq" directory
	input_dir=$(dirname $(dirname $1))
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
	input_dir=$(dirname $1)
else
	usage
fi

[ -z ${output_dir} ] && output_dir=${input_dir}/corall_batch_out
mkdir -p ${output_dir}

docker run -u $(id -u):$(id -g) --rm \
    -v ${input_dir}:${input_dir} \
    -v ${output_dir}:${output_dir} \
    -v ${gtf_file}:${gtf_file} \
    ${docker_img} \
    /corall_batch.sh ${corall_opt} \
    --gtfFile ${gtf_file} \
    --starGenomeDir ${star_index} \
    --baseDir ${output_dir} 