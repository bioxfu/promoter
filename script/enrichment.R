all_ball <- read.table('fasta/all_TAIR.txt', stringsAsFactors = F)$V1
white_ball <- read.table('cisfinder/cisfinder.v2.2000.positive.list.txt', stringsAsFactors = F)$V1
black_ball <- setdiff(all_ball, white_ball)

enrich_phyper <- function(drawn_ball) {
  q <- length(intersect(white_ball, drawn_ball))
  m <- length(white_ball)
  n <- length(black_ball)
  k <- length(drawn_ball)
  c(percent=q / k * 100, pvalue=phyper(q, m, n, k, lower.tail = FALSE))
}

DEG <- read.table('DEG/expr_table_cpm_DEG.tsv', header = T, sep='\t', quote = '', stringsAsFactors = F)
vs <- grep('_logFC|FDR', grep('_vs_', colnames(DEG), value = T), invert = T, value=T)

result <- NULL
for (i in 1:length(vs)) {
  result <- rbind(result, c(VS=vs[i], regulated='up', enrich_phyper(DEG$Gene[DEG[, vs[i]] == 1])))
  result <- rbind(result, c(VS=vs[i], regulated='down', enrich_phyper(DEG$Gene[DEG[, vs[i]] == -1])))
}

result <- data.frame(result, stringsAsFactors=F)
result$p.adjust <- p.adjust(as.numeric(result$pvalue), method = 'BH')

result$percent <- sprintf("%.2f", as.numeric(result$percent))
result$pvalue <- sprintf("%.2E", as.numeric(result$pvalue))
result$p.adjust <- sprintf("%.2E", as.numeric(result$p.adjust))

write.table(result, 'table/promoter_cis_element_enrichment.tsv', row.names = F, quote = F, sep = '\t')
