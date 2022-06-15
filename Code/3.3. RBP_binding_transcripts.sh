#!/bin/bash
# bedtools version=2.29.2

# overlaps of RBP binding sites and exons of transcripts
dir="/home/weihu/RBP_transcript_drug_project/4_RBP_regulation/data/eCLIPseq_bed/HepG2_eCLIPseq_bed/"
ls $dir |while read id
do
sample=`echo $id|sed "s/_.*//"`
bedtools intersect -a /home/weihu/RBP_transcript_drug_project/data/bed/CCLE_exon.bed -b ${dir}${id} -s -wa -wb > '/home/weihu/RBP_transcript_drug_project/4_RBP_regulation/data/HepG2_intersect_bed/'${sample}'_exon_eCLIP_intersect'
done
