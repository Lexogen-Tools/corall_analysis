#!/bin/bash

usage() { echo "Usage: $0 [-e <conda_environment>] [-m <mix2_binary>] <docker_image_tag>" 1>&2; exit 1; }

conda_env=./environment.yml
mix2_bin=./mix-square

while getopts ":e:m:" o; do
    case "${o}" in
        e)
            conda_env=${OPTARG}
            ;;
        m)
            mix2_bin=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

[ $# -eq 1 ] || usage

docker_image_tag=$1

docker build --tag ${docker_image_tag} --build-arg CONDA_ENV=${conda_env} --build-arg MIX2_BIN=${mix2_bin} . 
