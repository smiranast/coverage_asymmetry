#!/bin/bash
 
while getopts ":S:I:O:" flag;
do
	case $flag in
		S) sample_names=$OPTARG ;;
		I) input_dir=$OPTARG ;;
		O) output_dir=$OPTARG ;;
	esac
										         
done

tmp=$output_dir/tmp/
mkdir -p $output_dir
mkdir -p $output_dir/tmp/

while read i
do
    i=$(echo $i | cut -d";" -f1)
        
	echo $i 

	ls $input_dir/$i""Aligned.sortedByCoord.out.bam
	samtools view -b $input_dir/$i""Aligned.sortedByCoord.out.bam --fetch-pair "2:29192774-29921566" > $tmp/$i""Aligned.sortedByCoord.out.alk.bam
	if samtools view -H $tmp/$i""Aligned.sortedByCoord.out.alk.bam | grep '@RG'
	then 
	   pass;
	else
           samtools addreplacerg -r "@RG\tID:RG1\tSM:SampleName\tPL:Illumina\tLB:Library.fa" -o $tmp/$i""Aligned.sortedByCoord.out.alk.2.bam $tmp/$i""Aligned.sortedByCoord.out.alk.bam;
        fi
	mv $tmp/$i""Aligned.sortedByCoord.out.alk.2.bam $tmp/$i""Aligned.sortedByCoord.out.alk.bam

	samtools view $tmp/$i""Aligned.sortedByCoord.out.alk.bam | wc -l 
	picard MarkDuplicates REMOVE_DUPLICATES=true I=$tmp/$i""Aligned.sortedByCoord.out.alk.bam O=$tmp/$i""Aligned.sortedByCoord.out.bam M=$output_dir/metrices.txt VERBOSITY=ERROR QUIET=true
	samtools index $tmp/$i""Aligned.sortedByCoord.out.bam

	python3 $(dirname "$0")/TKR_script_sense_merged.paper.py -i $tmp/$i""Aligned.sortedByCoord.out.bam -o $output_dir/ --names $sample_names --genes ALK --stats $output_dir/coverage_stats.SnS.csv  --plot_all --table $output_dir/exons_stats.paper.csv --original $input_dir/$i""Aligned.sortedByCoord.out.bam

	rm $tmp/$i""Aligned.sortedByCoord.out.bam* $tmp/$i""Aligned.sortedByCoord.out.alk.bam
done < $sample_names

