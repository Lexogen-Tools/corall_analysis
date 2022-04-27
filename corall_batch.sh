#!/usr/bin/env bash

######################################
#
# dependencies and versions tested:
#       umi_tools 1.0.0
#       samtools 1.9
#       cutadapt 2.5
#       STAR 2.6.1b
#       GNU Awk 5.0.0
#       FastQC 0.11.8
#       multiQC 1.7
#	    mix-square 1.4.0.1
#	
#	    last change: 15/01/2021
#
#####################################

usage() { echo "Usage: $0 [-1|--r1 <rd1_fastq_gz>]  [-2|--r2 <rd2_fastq_gz>] [-c|--csv <csv_file>] [-g|--gtfFile <gtfFile>] [-s|--starGenomeDir <dir>] [-b|--baseDir <dir>] [-B|--mix2Blocks <nr_blocks>] [-C|--mix2Cores <nr_cores>]" 1>&2; exit 1; }

my_needed_commands="umi_tools samtools cutadapt STAR gawk fastqc multiqc"

missing_counter=0
for needed_command in $my_needed_commands; do
  if ! hash $needed_command >/dev/null 2>&1; then
    printf "Command not found in PATH or environment: %s\n" $needed_command >&2
    ((missing_counter++))
  fi
done

if ((missing_counter > 0)); then
  printf "Minimum %d commands are missing in PATH or environment, aborting\n" $missing_counter >&2
  exit 1
fi

mix2Blocks=3
mix2Cores=1

OPT=$(getopt -o 1:2:c:g:s:b:B:C: --long r1:,r2:,csv:,gtfFile:,starGenomeDir:,baseDir:,mix2Blocks:,mix2Cores: -- "$@")

[ $? != 0 ] && usage

eval set -- "$OPT"

while true; do
  case "$1" in
    -1 | --r1 ) inFastqRead1="$2"; shift 2;;
    -2 | --r2 ) inFastqRead2="$2"; shift 2;;
    -c | --csv ) CSVFILE="$2"; readarray -t INPUTFILES < ${CSVFILE}; shift 2;;
    -g | --gtfFile ) gtfFile="$2"; shift 2;;
    -s | --starGenomeDir ) starGenomeDir="$2"; shift 2;;
    -b | --baseDir ) baseDir="$2"; shift 2;;
    -B | --mix2Blocks ) mix2Blocks="$2"; shift 2;;
    -C | --mix2Cores ) mix2Cores="$2"; shift 2;;
    -- ) shift; break;;
    * ) usage;;
  esac
done

[ $# -gt 0 ] && { echo "no command line arguments outside options!" 1>&2; usage; }

[ -z ${inFastqRead1} ] && [ -z ${CSVFILE} ] && { echo "either --r1 or --csv have to be specified!" 1>&2; usage; }

[ -z ${gtfFile} ] && { echo "--gtfFile has to be specified!" 1>&2; usage; }

[ -z ${starGenomeDir} ] && { echo "--starGenomeDir has to be specified" 1>&2; usage; }

if [ -z ${CSVFILE} ]; then
	INPUTFILES=("${inFastqRead1},${inFastqRead2}")
fi

###################

set -o errexit

# for trimming on nextSeq or NovaSeq data
nextSeqTrim="--nextseq-trim=10"

# def. adapter sequences for corall protocol
rd2Adapter="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
rd1Adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"

for file in ${INPUTFILES[@]}; do

	row=${file// /} # remove whitespaces
	if [ -z $row ]; then
		break
	fi
        IFS="," read -r -a rowArr <<< $row
	inFastqRead1=${rowArr[0]}
	if [ -z $inFastqRead1 ]; then
		echo "no input file(s) entered..."
		exit 10
	fi

	outName=$(basename $inFastqRead1 .fastq.gz)
	outDir=${baseDir}/${outName} 
	mkdir -p ${outDir}

	inFastqRead2=${rowArr[1]}

	if [ -z $inFastqRead2 ]; then # SE

		echo "processing SE data file: $inFastqRead1"

		echo "extracting UMIs..."
		mkdir -p ${outDir}/fastq
		rd1_extract_fq=${outDir}/fastq/rd1_extract.fastq.gz
		awk '
		NR%4==1{ rd_name=$1; rd_info=$2 }
		NR%4==2{ umi=substr($1,1,10); rd_seq=substr($1,13) }
		NR%4==0{ print rd_name"_"umi" "rd_info; print rd_seq; print "+"; print substr($1,13) }' <( zcat $inFastqRead1 ) | gzip > ${rd1_extract_fq}

		echo "adapter trimming..."
		rd1_trim_fq=${outDir}/fastq/rd1_trim.fastq.gz
		cutadapt -m 20 -O 20 -a "QUALITY=G{20}" ${rd1_extract_fq}  | cutadapt  -m 20 $nextSeqTrim  -a $rd1Adapter - | cutadapt  -m 20 -O 3 -a "r1polyA=A{18}"  - | cutadapt  -m 20 -O 20 -g $rd1Adapter  --discard-trimmed  -o ${rd1_trim_fq} -
		echo "creating quality report..."
        	fastqc -o ${outDir}/fastq ${rd1_trim_fq}

		echo "alignment..."
		align_out_dir=${outDir}/align
		mkdir -p ${align_out_dir}
		STAR --runThreadN 12  --readFilesCommand zcat --genomeDir $starGenomeDir --readFilesIn ${rd1_trim_fq} --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
			--outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.6 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD --outSAMtype BAM SortedByCoordinate \
			--outFileNamePrefix ${align_out_dir}/ \
		     --limitBAMsortRAM 25000000000 --limitIObufferSize 200000000 --limitOutSJcollapsed 5000000 --seedPerWindowNmax 10

		# umi_tools dedup needs indexed bam files
		samtools index ${align_out_dir}/Aligned.sortedByCoord.out.bam

		echo "umi deduplicating..."
		umi_tools dedup -I ${align_out_dir}/Aligned.sortedByCoord.out.bam -S ${outDir}/umi_deduplicated.bam --multimapping-detection-method=NH --output-stats=${outDir}/deduplicated --log=${outDir}/deduplication.log 
		samtools index ${outDir}/umi_deduplicated.bam

	else

		echo "processing PE data files: $inFastqRead1 $inFastqRead2 " 
		echo "extracting UMIs..."
		rd1_fq_lines=$(wc -l < <(zcat $inFastqRead1))
        	rd2_fq_lines=$(wc -l < <(zcat $inFastqRead2))
        	if ((${rd1_fq_lines} != ${rd2_fq_lines}))
		then
			echo "ERROR: rd1/rd2 fastq files have different size"
			exit 1
        	fi

		mkdir -p ${outDir}/fastq
		extract_fq_prefix=${outDir}/fastq/rd
		extract_fq_suffix="_extract.fastq.gz"

		# checked window at end of read2 for UMI trimming
		regExWin=24
		matchLen=10

		split -l $rd1_fq_lines --numeric-suffixes=1 --additional-suffix=${extract_fq_suffix} --filter='gzip -c - > $FILE' <(awk -v win=$regExWin -v matchLen=$matchLen '
		function RC_match(match1,seq){
			m1=RC_fct(match1);
				return index(seq,m1);
		}
		function RC_fct(umi){
			split(umi,str,"")
			RC="";
			for(i=0;i<=length(umi);++i) {
				switch (str[i]) {
					case "A":
						RC="T"RC;
						break
					case "T":
						RC="A"RC;
						break
					case "G":
						RC="C"RC;
						break
					case "C":
						RC="G"RC;
						break
					default:
						break
				}
			}
			return RC;
		}

		BEGIN{count=0;winM=win-1;winP=win+1-matchLen}
		{if(NR==FNR){
		if(FNR%4==1){ rd_name=$1; rd_info=$2 }
		if(FNR%4==2){ umi=substr($1,1,10); matchRd=substr($1,13,matchLen);rd_seq=substr($1,13); rd_name2umi[rd_name]=umi;matchPt[rd_name]=matchRd;}
		if(FNR%4==0){ print rd_name"_"umi" "rd_info; print rd_seq; print "+"; print substr($1,13) }
		next }
		{
		if(FNR%4==1){ if(!rd_name2umi[$1]){ print "ERROR: no read1 for read2="$1 > "/dev/stderr"; exit 1 }else{ umi=rd_name2umi[$1];matchRd2=matchPt[$1];print $1"_"umi" "$2} }
		#if(FNR%4==2){val=RC_match(matchRd2,substr($1,length($1)-winM)); if(val!=0){print FNR, "found: "matchRd2 "(RC:"RC_fct(matchRd2)") in ", substr($1,length($1)-winM), "RC UMI:" RC_fct(umi),substr($1,1,length($1)-winP+  val), $1 > "/dev/stderr"}}
		if(FNR%4==2){val=RC_match(matchRd2,substr($1,length($1)-winM)); if(val!=0){print substr($1,1,length($1)-winP+val);count++}else {print}}
		if(FNR%4==3){print}
		if(FNR%4==0){if(val!=0){print substr($1,1,length($1)-winP+val)}else {print}}
		}}END{print "removed UMI seqences in rd2: " count > "/dev/stderr"}' <(zcat $inFastqRead1) <(zcat $inFastqRead2) ) ${extract_fq_prefix}

		echo "adapter trimming..."
		rd1_trim_fq=${outDir}/fastq/rd1_trim.fastq.gz
		rd2_trim_fq=${outDir}/fastq/rd2_trim.fastq.gz
		cutadapt -m 20 -O 20 --interleaved -n 2 -a "QUALITY=G{20}" -A "QUALITY=G{20}" "${extract_fq_prefix}01${extract_fq_suffix}" "${extract_fq_prefix}02${extract_fq_suffix}" \
			| cutadapt -m 20 --interleaved -n 3 $nextSeqTrim -a $rd1Adapter -A N{12}$rd2Adapter";min_overlap=15;max_error_rate=0.045455" -G XT{18} - | cutadapt  -m 20 -O 3 --interleaved -n 1  -a "r1polyA=A{18}"  - \
			| cutadapt -m 20 -O 20 --interleaved -g $rd1Adapter -G $rd2Adapter --discard-trimmed -o ${rd1_trim_fq} -p ${rd2_trim_fq} - 

		echo "creating quality report..."
		fastqc -o ${outDir}/fastq ${rd1_trim_fq}
		fastqc -o ${outDir}/fastq ${rd2_trim_fq}

		echo "alignment..."
		align_out_dir=${outDir}/align
		mkdir -p ${align_out_dir}
		STAR --runThreadN 12  --peOverlapNbasesMin 40 --peOverlapMMp 0.8 --readFilesCommand zcat --genomeDir $starGenomeDir --readFilesIn ${rd1_trim_fq} ${rd2_trim_fq} --outFilterType BySJout \
		     --outFilterMultimapNmax 200 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.6 --alignIntronMin 20 --alignIntronMax 1000000 \
		     --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD --outSAMtype BAM SortedByCoordinate --outFileNamePrefix "${align_out_dir}/" \
		     --limitBAMsortRAM 25000000000 --limitIObufferSize 200000000 --limitOutSJcollapsed 5000000 --seedPerWindowNmax 10

		# umi_tools dedup needs indexed bam files
		samtools index ${align_out_dir}/Aligned.sortedByCoord.out.bam

		echo "umi deduplicating..."
		umi_tools dedup -I ${align_out_dir}/Aligned.sortedByCoord.out.bam -S ${outDir}/umi_deduplicated.bam --multimapping-detection-method=NH --output-stats=${outDir}/deduplicated --paired --log=${outDir}/deduplication.log 
		samtools index ${outDir}/umi_deduplicated.bam
 
	fi

	echo "mixSquare analysis of deduplicated bam file" 
	if (($mix2Cores > 1))
	then
		mix-square -G $gtfFile -p ${mix2Cores} -B ${outDir}/umi_deduplicated.bam -D fwd -o ${outDir} -b ${mix2Blocks} 
	else
		mix-square -G $gtfFile -B ${outDir}/umi_deduplicated.bam -D fwd -o ${outDir} -b ${mix2Blocks} 
	fi

	dat_temp=$(mktemp -p ${outDir})
	awk 'NR==1{ printf("%s\t",$0); print "expected_comp_frags_trans" }NR>1{ printf("%s\t",$0); printf("%d\n",$7*$11) }' <(cut -f1-13 ${outDir}/transcripts_summary_umi_deduplicated.dat) > ${dat_temp}
	mv ${dat_temp} ${outDir}/transcripts_summary_umi_deduplicated.dat

	echo "finished analysis for $inFastqRead1 $inFastqRead2"
	echo ""
done

echo "running multiQC.."
multiqc -d -dd 2 ${baseDir} -o ${baseDir}

echo "done"

exit 0
