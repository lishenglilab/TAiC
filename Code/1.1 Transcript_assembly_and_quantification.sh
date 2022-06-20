## Transcript assembly and quantification from RNA-seq data.

#!/bin/bash
# STAR version=2.7.6a
# Stringtie version=2.1.4
# genome version gencodev35_grch38.p13
# gtf version gencode.v35.annotation.gtf

sample=$1

## STAR 2pass RNA-seq read alignment

star='/mnt/data132Tp3/public/201212/softwares/STAR/bin/Linux_x86_64/STAR'
gdir='/mnt/data132Tp3/public/201212/ref/star_gencodev35_grch38.p13_index'
file_gtf='/mnt/data132Tp3/public/201212/ref/gencode.v35.annotation.gtf'
dir_fq='/mnt/data132Tp3/public/201212/fastq'
dir_alignment='/mnt/data132Tp3/public/201212/star_alignments'
mkdir -p ${dir_alignment}/${sample}
fq1=${dir_fq}/${sample}/${sample}'_R1.fastq.gz'
fq2=${dir_fq}/${sample}/${sample}'_R2.fastq.gz'
 
${star} --twopassMode Basic --twopass1readsN -1 --runThreadN 5 --genomeDir ${gdir} --readFilesIn ${fq1} ${fq2} --readFilesCommand zcat \
--sjdbGTFfile ${file_gtf} --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
--outFileNamePrefix ${dir_alignment}/${sample}/${sample}'_alignment' --outSAMtype BAM SortedByCoordinate


## Stringtie reference-based transcript assembly

stringtie='/mnt/data132Tp3/public/201212/softwares/stringtie-2.1.4.Linux_x86_64/stringtie'
dir_align='/mnt/data132Tp3/public/201212/star_alignments'
file_bam=${dir_align}'/'${sample}'/'${sample}'_alignmentAligned.sortedByCoord.out.bam'
dir_stringtie='/mnt/data132Tp3/public/201212/stringtie_assembly'
mkdir -p ${dir_stringtie}'/'${sample}
ref_gtf='/mnt/data132Tp3/public/201212/ref/gencode.v35.annotation.gtf'
out_gtf=${dir_stringtie}'/'${sample}'/'${sample}'_stringtie.gtf'
${stringtie} ${file_bam} -p 5 -o ${out_gtf} -G ${ref_gtf}

## merge transcript assembly from individual samples

stringtie='/mnt/data132Tp3/public/201212/softwares/stringtie-2.1.4.Linux_x86_64/stringtie'
gtf_files='/mnt/data132Tp3/public/201212/data/gtf_files.txt' # list of transcript assembly of all samples
ref_gtf='/mnt/data132Tp3/public/201212/ref/gencode.v35.annotation.gtf'
out_gtf='/mnt/data132Tp3/public/201212/data/CCLE_1017_merged.gtf'
${stringtie} --merge -p 10 -o ${out_gtf} -G ${ref_gtf} ${gtf_files}

## Stringtie transcript quantification
## this quantification workflow is for both CCLE and TCGA samples

stringtie='/mnt/data132Tp3/public/201212/softwares/stringtie-2.1.4.Linux_x86_64/stringtie'
dir_align='/mnt/data132Tp3/public/201212/star_alignments'
file_bam=${dir_align}/${sample}/${sample}'_alignmentAligned.sortedByCoord.out.bam'
dir_stringtie='/mnt/data132Tp3/public/201212/stringtie_quantification'
mkdir -p ${dir_stringtie}/${sample}
ref_gtf='/mnt/data132Tp3/public/201212/data/CCLE_1017_merged.gtf'
out_gtf=${dir_stringtie}/${sample}/${sample}'_stringtie.gtf'
file_gene=${dir_stringtie}/${sample}/${sample}'_gene_abundance.tab'
${stringtie} ${file_bam} -p 8 -o ${out_gtf} -G ${ref_gtf} -A ${file_gene} -e
