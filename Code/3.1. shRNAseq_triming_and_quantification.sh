#!/bin/bash
# STAR version=2.7.6a
# Stringtie version=2.1.4
# Trimmomatic version=0.39

## Trim 
trim='/mnt/data132Tp3/public/201212/softwares/Trimmomatic-0.39/trimmomatic-0.39.jar'
adapter='/mnt/data132Tp3/public/201212/softwares/Trimmomatic-0.39/adapters/PE_adapters.fa'
dir_wk='/mnt/data132Tp3/public/201212/K562'
dir_fq='/mnt/data132Tp3/public/201212/K562/K562_fq' 

ls $dir_fq |sed 's/\///' |while read sample
do
ls ${dir_fq}/${sample}|grep .fastq.gz|sed "s/${sample}_//" |sed "s/_.*.fastq.gz$//"|uniq|while read id
do
fq1=${dir_fq}/${sample}/`ls ${dir_fq}/${sample} |grep $id|head -n 1`
fq2=${dir_fq}/${sample}/`ls ${dir_fq}/${sample} |grep $id|tail -n 1`
dir_trim=${dir_wk}'/Trim_results'
mkdir -p ${dir_trim}/${sample}'_'${id}
filt_fn_r1=${sample}'_'${id}'_filtered_R1.fastq.gz'
unp_fn_r1=${sample}'_'${id}'_unpaired_R1.fastq.gz'
filt_fn_r2=${sample}'_'${id}'_filtered_R2.fastq.gz'
unp_fn_r2=${sample}'_'${id}'_unpaired_R2.fastq.gz'

echo "java -jar ${trim} PE ${fq1} ${fq2} ${dir_trim}/${sample}_${id}/${filt_fn_r1} \
        ${dir_trim}/${sample}_${id}/${unp_fn_r1} ${dir_trim}/${sample}_${id}/${filt_fn_r2} ${dir_trim}/${sample}_${id}/${unp_fn_r2} \
        'ILLUMINACLIP:'${adapter}':2:30:10' LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> ${dir_trim}/${sample}_${id}/read_surviving_stat.txt" 
done
done>Trim.sh

nohup bash Trim.sh &

## STAR 2pass RNA-seq read alignment
star='/mnt/data132Tp3/public/201212/softwares/STAR/bin/Linux_x86_64/STAR'
gdir='/mnt/data132Tp3/public/201212/ref/star_gencodev35_grch38.p13_index'
file_gtf='/mnt/data132Tp3/public/201212/ref/gencode.v35.annotation.gtf'
dir_fq='/mnt/data132Tp3/public/201212/K562/Trim_results' 
dir_alignment='/mnt/data132Tp3/public/201212/K562/K562_star_alignments'

ls $dir_fq |sed 's/\///' |while read sample
do

mkdir -p ${dir_alignment}/${sample}
fq1=${dir_fq}/${sample}/`ls ${dir_fq}/${sample} |grep ${sample}'_filtered'|head -n 1`
fq2=${dir_fq}/${sample}/`ls ${dir_fq}/${sample} |grep ${sample}'_filtered'|tail -n 1`

echo "${star} --twopassMode Basic --twopass1readsN -1 --runThreadN 5 --genomeDir ${gdir} --readFilesCommand zcat --readFilesIn ${fq1} ${fq2} \
--sjdbGTFfile ${file_gtf} --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
--outFileNamePrefix ${dir_alignment}/${sample}/${sample}'_alignment' --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 58172905921"
done >star.sh

nohup bash star.sh &

## Stringtie quantification
stringtie='/mnt/data132Tp3/public/201212/softwares/stringtie-2.1.4.Linux_x86_64/stringtie'
dir_align='/mnt/data132Tp3/public/201212/K562/K562_star_alignments'
dir_stringtie='/mnt/data132Tp3/public/201212/K562/stringtie_quantification'

ls $dir_align | while read sample
do
file_bam=${dir_align}'/'${sample}'/'${sample}'_alignmentAligned.sortedByCoord.out.bam'
mkdir -p ${dir_stringtie}'/'${sample}
ref_gtf='/mnt/data132Tp3/public/201212/K562/script_K562/CCLE_1017_merged.gtf'
out_gtf=${dir_stringtie}'/'${sample}'/'${sample}'_stringtie.gtf'
file_gene=${dir_stringtie}/${sample}/${sample}'_gene_abundance.tab'

echo "${stringtie} ${file_bam} -p 8 -o ${out_gtf} -G ${ref_gtf} -A ${file_gene} -e"
done >stringtie.sh

nohup bash stringtie.sh &
