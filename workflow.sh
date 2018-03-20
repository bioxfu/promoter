# download promoter sequences from Ensembl Plant using Biomart
# http://plants.ensembl.org/biomart/martview?VIRTUALSCHEMANAME=plants_mart&ATTRIBUTES=athaliana_eg_gene.default.sequences.ensembl_gene_id|athaliana_eg_gene.default.sequences.gene_flank|athaliana_eg_gene.default.sequences.upstream_flank."2000"&FILTERS=athaliana_eg_gene.default.filters.with_tair_locus.only&VISIBLEPANEL=resultspanel

# must add "ENDLINE" as the final line in the sequence file
echo 'ENDLINE' >> fasta/tair10_promoter_2kb.fasta

# Searches for the vCGCGb motif
# v:[ACG]
# b:[CGT]
perl script/cisfinder.pl fasta/tair10_promoter_2kb.fasta 2000 ACGCGC ACGCGG ACGCGT CCGCGC CCGCGG CCGCGT GCGCGC GCGCGG GCGCGT
mv cisfinder.v2* cisfinder

Rscript script/enrichment.R
