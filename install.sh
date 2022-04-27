#!/bin/bash -e

#requires existing installation of miniconda3 in home directory

usage() { echo "Usage: $0 <path_to_miniconda_installation>" 1>&2; exit 1; }

[ $# -ne 1 ] && usage

miniconda_install_dir=$1

export PATH="${miniconda_install_dir}/bin:$PATH"

conda env create -f ./environment.yml

cp corall_batch.sh ${miniconda_install_dir}/envs/corall/bin

cp mix-square ${miniconda_install_dir}/envs/corall/bin

rm -rf ${extract_dir}

