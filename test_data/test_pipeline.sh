#!/bin/bash

echo "Generating STAR index.."
mkdir ./test_data/star_index

conda_dir=$(which conda | grep -oP "^/.*(?=bin)")

echo "conda dir detected at: ${conda_dir}"

source ${conda_dir}/bin/activate corall
STAR --runMode genomeGenerate --genomeFastaFiles ./test_data/test_genome.fasta --genomeSAindexNbases 7 --sjdbGTFfile ./test_data/test_genome.gtf --sjdbOverhang 99 --genomeDir ./test_data/star_index --outTmpDir ./test_data/star_index/tmp

if [ ! $? -eq 0 ]; then
	echo "an error occured during STAR index creation" >&2
	exit 1
fi

conda deactivate

echo "Running CORALL pipeline.."

./run_corall_local.sh -m ${conda_dir} -a ./test_data/test_genome.gtf -s ./test_data/star_index/ -o ./test_data/test_output/ ./test_data/samples.csv

if [ ! $? -eq 0 ]; then
        echo "an error occured during corall pipeline" >&2
        exit 1
fi

cd test_data/test_output && find . > ../file_list
cd ../..

echo "Comparing file structure of the test output folder.."
d=$(diff <(sort -h ./test_data/file_list) <(sort -h ./test_data/expected_file_list))
if [ ! -z "$d" ]; then
	echo "Unexpected output file structure generated.. differences were found as follows:" >&2
	echo -e $d
	exit 1
fi

echo "Comparing gene quantification with expected output.."
d=$(diff <(awk '{print $3"\t"$9}' ./test_data/test_output/test_R1/genes_summary_umi_deduplicated.dat) ./test_data/gene_table_reference)
if [ ! -z "$d" ]; then
	echo "gene quantification table has changed.. please check if this is intended" >&2
	exit 1
fi

echo "All tests passed successfully"

rm -rf ./test_output
rm -rf ./star_index
