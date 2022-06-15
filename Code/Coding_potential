#!/bin/bash
# ORFfinder version=0.4.3
# CPC2 version=1.0.1
# CPAT version=3.0.4
# HMMER version=3.3

# ORFfinder
nohup /home/weihu/softwares/ORFfinder -ml 30 -strand plus -in /home/weihu/RBP_transcript_drug_project/data/fa/CCLE_1017_merged_novel_noncoding.fa -out novel_noncoding_transcript_ORFfinder > ORFfinder.log &

# CPC2
nohup /home/weihu/softwares/CPC2_standalone-1.0.1/bin/CPC2.py -i /home/weihu/RBP_transcript_drug_project/data/fa/CCLE_1017_merged_novel_noncoding.fa -o novel_noncoding_transcript_CPC2 &

# CPAT
nohup /home/weihu/.local/bin/cpat.py -x Human_Hexamer.tsv -d Human_logitModel.RData --top-orf=5 -g /home/weihu/RBP_transcript_drug_project/data/fa/CCLE_1017_merged_novel_noncoding.fa -o novel_noncoding_transcript_CPAT &

# HMMER
nohup hmmsearch --domtblout novel_noncoding_transcript_Pfam /home/Wei.Hu/Task_pfam_scan/PfamScan/Pfam-A.hmm /home/weihu/RBP_transcript_drug_project/data/fa/CCLE_1017_merged_novel_noncoding.fa
