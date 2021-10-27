#!/bin/bash

method=`sed '/^method=/!d;s/.*=//' conf_splicing.txt`


################
if [ "$method" == "DEXseq" ]; then
gtf_1=`sed '/^gtf_1=/!d;s/.*=//' conf_splicing.txt`
script=../scripts
output_path_1=`sed '/^output_path_1=/!d;s/.*=//' conf_splicing.txt`
sam_file_path_1=`sed '/^sam_file_path_1=/!d;s/.*=//' conf_splicing.txt`
metaData_1=`sed '/^metaData_1=/!d;s/.*=//' conf_splicing.txt`
number_cores_1=`sed '/^number_cores_1=/!d;s/.*=//' conf_splicing.txt`
pvalue_1=`sed '/^pvalue_1=/!d;s/.*=//' conf_splicing.txt`
Gene_to_plot_1=`sed '/^Gene_to_plot_1=/!d;s/.*=//' conf_splicing.txt`



raw_dir=$(pwd)

cd $script
python_script=$(pwd)

python $python_script/dexseq_prepare_annotation.py \
$gtf_1 \
$output_path_1/DEXSeq.chr.gff

cd $sam_file_path_1
for file in $(ls | grep .sam)
        do
pre_name=${file%.sam}
python $python_script/dexseq_count.py \
$output_path_1/DEXSeq.chr.gff \
$file $output_path_1/${pre_name}.txt
        done

Rscript $raw_dir/DEXseq.R -c $output_path_1 -m $metaData_1 -p $pvalue_1 -o $output_path_1 -g $Gene_to_plot_1
fi
##################

##################



if [ "$method" == "StringTie+ballgown" ]; then
sam_path_2=`sed '/^sam_path_2=/!d;s/.*=//' conf_splicing.txt`
#samtools_path_2=`sed '/^samtools_path_2=/!d;s/.*=//' conf_splicing.txt`
out_path_2=`sed '/^out_path_2=/!d;s/.*=//' conf_splicing.txt`
#stringtie_path_2=`sed '/^stringtie_path_2=/!d;s/.*=//' conf_splicing.txt`
gtf_2=`sed '/^gtf_2=/!d;s/.*=//' conf_splicing.txt`
metaData_2=`sed '/^metaData_2=/!d;s/.*=//' conf_splicing.txt`
qvalue_2=`sed '/^qvalue_2=/!d;s/.*=//' conf_splicing.txt`
Gene_to_plot_2=`sed '/^Gene_to_plot_2=/!d;s/.*=//' conf_splicing.txt`

raw_dir=$(pwd)

cd $sam_path_2 || exit 1
for file in $(ls | grep .sam)
        do
pre_name=${file%.sam}
#$samtools_path_2/
samtools sort -@ 8 -o $out_path_2/${pre_name}.bam $file
        done
cd $out_path_2 || exit 1
mkdir assembly
for file in $(ls | grep .bam)
        do
pre_name=${file%.bam}
#$stringtie_path_2/
stringtie $file -l $pre_name -p 8 -G $gtf_2 -o assembly/${pre_name}.gtf
echo "assembly/${pre_name}.gtf" >> mergelist.txt
        done
#$stringtie_path_2/
stringtie --merge -p 8 -G $gtf_2 -o stringtie_merged.gtf mergelist.txt

echo "The number of transcripts"
cat stringtie_merged.gtf | grep -v "^#" | awk '$3=="transcript" {print}' | wc -l

mkdir ballgown
for file in $(ls | grep .bam)
        do
pre_name=${file%.bam}
mkdir ballgown/$pre_name
#$stringtie_path_2/
stringtie -e -B -p 8 -G stringtie_merged.gtf -o ballgown/$pre_name/${pre_name}.gtf $file
        done
        
Rscript $raw_dir/ballgown.R -m $metaData_2 -q $qvalue_2 -o $out_path_2 -g $Gene_to_plot_2
fi
####################

###################

if [ "$method" == "rMATS" ]; then
bam_path_3=`sed '/^bam_path_3=/!d;s/.*=//' conf_splicing.txt`
read_length_3=`sed '/^read_length_3/!d;s/.*=//' conf_splicing.txt`
condition1_files_3=`sed '/^condition1_files_3=/!d;s/.*=//' conf_splicing.txt`
condition2_files_3=`sed '/^condition2_files_3=/!d;s/.*=//' conf_splicing.txt`
gtf_3=`sed '/^gtf_3=/!d;s/.*=//' conf_splicing.txt`
out_path_3=`sed '/^out_path_3=/!d;s/.*=//' conf_splicing.txt`
paired_or_single_3=`sed '/^paired_or_single_3=/!d;s/.*=//' conf_splicing.txt`
tmp_3=`sed '/^tmp_3=/!d;s/.*=//' conf_splicing.txt`
event_type_3=`sed '/^event_type_3=/!d;s/.*=//' conf_splicing.txt`
number_of_top_events_to_plot_3=`sed '/^number_of_top_events_to_plot_3=/!d;s/.*=//' conf_splicing.txt`

RNASeq-MATS.py --b1 $condition1_files_3 --b2 $condition2_files_3 --gtf $gtf_3  --od $out_path_3 --tmp $tmp_3 \
-t $paired_or_single_3 --nthread 3 --readLength $read_length_3

if [ ! -d "rmats2sashimiplot" ]; then
  git clone https://github.com/Xinglab/rmats2sashimiplot
fi

cd rmats2sashimiplot
./2to3.sh
cd ..

file1=$(cat $condition1_files_3)
file2=$(cat $condition2_files_3)

if [ "$event_type_3" == "SE" ]; then
  event_file=$out_path_3/SE.MATS.JC.txt
fi
if [ "$event_type_3" == "A5SS" ]; then
  event_file=$out_path_3/A5SS.MATS.JC.txt
fi
if [ "$event_type_3" == "A3SS" ]; then
  event_file=$out_path_3/A3SS.MATS.JC.txt
fi
if [ "$event_type_3" == "MXE" ]; then
  event_file=$out_path_3/MXE.MATS.JC.txt
fi
if [ "$event_type_3" == "RI" ]; then
  event_file=$out_path_3/RI.MATS.JC.txt
fi


python rMATS_results_filter.py -i $event_file -o $out_path_3 -n $number_of_top_events_to_plot_3

raw_dir=$(pwd)

cd $out_path_3
if [ ! -d "event_plot" ]; then
  mkdir event_plot
fi

cd $raw_dir
python rmats2sashimiplot/src/rmats2sashimiplot/rmats2sashimiplot.py \
--l1 condition1  --l2 condition2 \
--b1 $file1 \
--b2 $file2 \
-t $event_type_3 \
-e $out_path_3/top_events.txt \
--exon_s 1 \
--intron_s 5 \
-o $out_path_3/event_plot


fi


