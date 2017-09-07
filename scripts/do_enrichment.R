psd <- read.table("psd.tsv", sep="\t", header=TRUE, stringsAsFactors=FALSE, quote = "")

find_overlap <- function(res, psd, alpha=0.05, use_col='cyberT_adj') {

  is_sig = apply(as.matrix(res[,use_col] <= alpha), 1, all)
  up_proteins_psd <- merge(
                           res[is_sig & res$fold_change > 0,], 
                           psd, by.x="accession_number", by.y="Accession") 
  down_proteins_psd <- merge(
                             res[is_sig & res$fold_change < 0,],
                             psd, by.x="accession_number", by.y="Accession")
  is_up <- res$fold_change > 0
  is_down <- res$fold_change < 0

  sig_psd <- sum(!duplicated(rbind(up_proteins_psd, down_proteins_psd)))
  sig <- sum(!duplicated(res[is_sig, "accession_number"]))
  all_seen <- sum(!duplicated(res[,c('accession_number')]))
  all_psd <- length(intersect(psd$Accession, res$accession_number))

  logp <- phyper(sig_psd, all_psd, all_seen - all_psd, sig, log.p=TRUE, lower.tail=FALSE)

  up_proteins <- unique(res$accession_number[is_sig & is_up])
  down_proteins <- unique(res$accession_number[is_sig & is_down])

  return(list(up_proteins_psd=up_proteins_psd,
              down_proteins_psd=down_proteins_psd,
              up_proteins=up_proteins,
              down_proteins=down_proteins,
              log_hyperp=logp,
              num_proteins_up=length(up_proteins),
              num_proteins_up_psd=nrow(up_proteins_psd),
              num_proteins_down=length(down_proteins),
              num_proteins_down_psd=nrow(down_proteins_psd), 
              num_proteins_total=all_seen))
}

plot_venn <- function(overlap, name) {
  # TODO more elegant
  overlap$num_proteins_both_psd <- 0
  p = 1:1556
  u = c(c(1:overlap$num_proteins_up_psd),
        c(1557:(1556 + overlap$num_proteins_up - overlap$num_proteins_up_psd)))
  d = c(c((overlap$num_proteins_up_psd + 1 - overlap$num_proteins_both_psd):
          (overlap$num_proteins_up_psd - overlap$num_proteins_both_psd + overlap$num_proteins_down_psd)),
        c((1558 + overlap$num_proteins_up - overlap$num_proteins_up_psd - overlap$num_proteins_both):
          (1557 + overlap$num_proteins_up - overlap$num_proteins_up_psd - overlap$num_proteins_both + overlap$num_proteins_down - overlap$num_proteins_down_psd)))
  venn.diagram(x=list(
                      "PSD"=p,
                      "Up-regulated"=u,
                      "Down-regulated"=d),
               filename=paste("figures", name, ".tiff", sep=""),
               fill = c("red", "green", "blue"),
               alpha = c(0.5, 0.5, 0.5),
               cex = 2,
               cat.fontface = 4,
               lty =2,
               fontfamily =3,
               euler.d = TRUE)
}

find_reverse <- function(phosphos) {
  return(ngkd_joint_complete[
         (ngkd_joint_complete$adj.P.Val <= 0.05 & ngkd_joint_complete$Average.Log2.Expression > 0 & ngkd_joint_complete$Norm.adj.P.Val <= 0.05 & ngkd_joint_complete$Log2.Norm.Avg < 0)
         | (ngkd_joint_complete$adj.P.Val <= 0.05 & ngkd_joint_complete$Average.Log2.Expression < 0 & ngkd_joint_complete$Norm.adj.P.Val <= 0.05 & ngkd_joint_complete$Log2.Norm.Avg > 0),
         c('Average.Log2.Expression', 'adj.P.Val', 'Gene.Symbol', 'Accession.Number', 'Tot.Log2.Med.Rep1', 'Tot.Log2.Med.Rep2', 'Tot.Log2.Med.Rep3', 'Log2.Norm.Avg', 'Norm.adj.P.Val')])
}

# overlap <- find_overlap(test_res, psd)

setup_enrich <- function() {
  source("https://bioconductor.org/biocLite.R")
  library(clusterProfiler)
  library(org.Mm.eg.db)
}

convert_entrez <- function(gene_list, fromType='ACCNUM', filename=NULL) {
  entrez <- bitr(gene_list, fromType=fromType, toType=c("ENTREZID", "SYMBOL"), OrgDb="org.Mm.eg.db")
  if (!is.null(filename)) {
    write.table(entrez, paste("enrichment/", filename, sep=''), sep='\t', row.names=FALSE, quote=FALSE)
  }
  return(entrez$ENTREZID[!duplicated(entrez$ACCNUM)])
}

go_enrich <- function(gene_list, ont="MF", is_entrez=FALSE) {
  if (!is_entrez) {
    entrez_genes <- convert_entrez(gene_list)
  } else {
    entrez_genes <- gene_list
  }
  return(enrichGO(gene=entrez_genes, OrgDb="org.Mm.eg.db", ont=ont, readable=TRUE))
}

go_group <- function(gene_list, ont="MF", level=3, is_entrez=FALSE) {
  # Provide gene_list as list of 
  if (!is_entrez) {
    entrez_genes <- convert_entrez(gene_list)
  } else {
    entrez_genes <- gene_list
  }
  out <- groupGO(gene=entrez_genes, OrgDb='org.Mm.eg.db', ont=ont, level=level, readable=TRUE)
  out@result <- out@result[order(out@result$Count, decreasing=TRUE),]
  return(out)
}

generate_gsea_gene_list <- function(ng_joint_complete) {
  # TODO: sort by p-values or expression change?
  sorted <- ng_joint_complete[order(ng_joint_complete$Norm.adj.P.Val, decreasing=TRUE),c("Accession.Number.NoIso", "Norm.adj.P.Val")]
  indices <- !duplicated(ng_joint_complete$Accession.Number.NoIso)
  gsea_pvals <- sorted$Norm.adj.P.Val[indices]
  names(gsea_pvals) <- sorted$Accession.Number.NoIso[indices]
  return(gsea_pvals)
}


filter_go_rankings_to <- function(enrich_result, level) {
  # 1 is most general, 5 is most specific
  tmp <- enrich_result
  for (l in 1:level) {
    tmp <- dropGO(tmp, level=l)
  }
  return(tmp)
}

RUN_OVERLAP <- F
RUN_CONVERT <- F
RUN_ENRICH <- F
if (RUN_OVERLAP) {
  psd_overlap <- find_overlap(res_protein, psd, alpha=0.05, use_col='cyberT_adj')
  psd_overlap_strict <- 
    find_overlap(res_protein, psd, alpha=0.01, use_col='modT_2samp_adj')
}
if (RUN_CONVERT) {
  # Convert gene names
  # down_entrez_genes <- convert_entrez(
  #     unique(psd_overlap$down_proteins),
  # )
  # up_entrez_genes <- convert_entrez(
  #     unique(psd_overlap$up_proteins),
  # )

  down_entrez_genes_psd <- convert_entrez(
      unique(psd_overlap$down_proteins_psd$accession_number),
  )
  up_entrez_genes_psd <- convert_entrez(
      unique(psd_overlap$up_proteins_psd$accession_number),
  )
  down_entrez_genes_psd_strict <- convert_entrez(
      unique(psd_overlap_strict$down_proteins_psd$accession_number),
  )
  up_entrez_genes_psd_strict <- convert_entrez(
      unique(psd_overlap_strict$up_proteins_psd$accession_number),
  )
}

if (RUN_ENRICH) {
  # Do enrichment
  # enrich_down <- go_enrich(down_entrez_genes, is_entrez=TRUE)
  # enrich_up <- go_enrich(up_entrez_genes, is_entrez=TRUE)
  enrich_down_psd <- go_enrich(down_entrez_genes_psd, is_entrez=TRUE)
  enrich_up_psd <- go_enrich(up_entrez_genes_psd, is_entrez=TRUE)
  enrich_down_psd_strict <- go_enrich(down_entrez_genes_psd_strict, is_entrez=TRUE)
  enrich_up_psd_strict <- go_enrich(up_entrez_genes_psd_strict, is_entrez=TRUE)
}

# ngkd_down_enrich <- go_enrich(unique(
#   ngkd_joint_complete$Accession.Number.NoIso[
#       ngkd_joint_complete$Norm.adj.P.Val <= 0.05 &
#       ngkd_joint_complete$Log2.Norm.Avg < 0
#   ]), ont="BP")

# ngkd_up_enrich <- go_enrich(unique(
#   ngkd_joint_complete$Accession.Number.NoIso[
#       ngkd_joint_complete$Norm.adj.P.Val <= 0.05 &
#       ngkd_joint_complete$Log2.Norm.Avg > 0
#   ]), ont="BP")

# GSEA FAILS FOR SOME REASON
# TODO figure out why
# ngkd_down_gsea <- gseGO(
#   generate_gsea_gene_list(ngkd_joint_complete),
#   organism="mouse",
#   ont="BP",
#   verbose=TRUE)

# ngkd_down_group <- go_group(unique(
#   ngkd_joint_complete$Accession.Number.NoIso[
#       ngkd_joint_complete$Norm.adj.P.Val <= 0.05 &
#       ngkd_joint_complete$Log2.Norm.Avg < 0
#   ]))
