python script/get_DEG_promoter.py fasta/tair10_promoter_2kb.tab DEG/RPKM_table_FDR0.05_FC1.5_DEG.tsv > DEG/RPKM_table_FDR0.05_FC1.5_DEG.tair10_promoter_2kb.fasta

fimo --qv-thresh --thresh 0.01 motif/JASPAR2018_CORE_plants_non-redundant_pfms_meme.txt DEG/RPKM_table_FDR0.05_FC1.5_DEG.tair10_promoter_2kb.fasta

source activate gmatic
echo -e 'TF\ttarget_num\ttarget_list' > table/FIMO_scan_DEG_promoter_2kb_with_JASPAR_motifs.tsv
grep -v '#' fimo_out/fimo.tsv|grep -v 'motif_id' |cut -f2,3|sort -k1,1 -k2,2|groupBy -g 1 -c 2,2 -o count_distinct,distinct|sort -rnk2 >> table/FIMO_scan_DEG_promoter_2kb_with_JASPAR_motifs.tsv
