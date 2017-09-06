library(plyr)

read_psm <- function(filename='psm_raw.csv', med_norm=TRUE) {
  rel_cols = c(
        "TMT_126_131", "TMT_127N_131",
        "TMT_127C_131", "TMT_128N_131",
        "TMT_128C_131", "TMT_129N_131",
        "TMT_129C_131", "TMT_130N_131",
        "TMT_130C_131", "accession_number",
        "sequence", "modifications",
        "accession_numbers", "geneSymbol",
        "filename"
        )
  data <- read.csv(filename)[,rel_cols]
  # Log the value cols
  data[,1:9] <- log2(data[,1:9])

  if (med_norm) {
    # Median normalize each ratio column
    medians = unlist(lapply(data[,1:9], function(x) median(x, na.rm=TRUE)))
    data[,1:9] <- sweep(data[,1:9], 2, medians, '-')
  }

  # Extract the Plex
  data$plex <- unlist(lapply(as.character(data$filename),
                                 function(x) strsplit(x, '_')[[1]][6]))
  data$filename <- NULL 
    
  # Split by plex, then recombine in a more workable form
  plexA <- data[data$plex == 'PlexA',]
  colnames(plexA)[1:9] <- c('GFP_A1', 'GFP_A2',
                            'KO95_A1', 'KO95_A2', 'KO95_A3',
                            'KO93_A1', 'KO93_A2',
                            'DKO_A1', 'DKO_A2')
  plexB <- data[data$plex == 'PlexB',]
  colnames(plexB)[1:9] <- c('GFP_B1', 'GFP_B2',
                            'KO95_B1', 'KO95_B2',
                            'KO93_B1', 'KO93_B2', 'KO93_B3',
                            'DKO_B1', 'DKO_B2')
  data = rbind.fill(plexA, plexB)
  return(data)
}

read_protein <- function(filename='protein_SGS_raw.csv', med_norm=TRUE, group_genes=TRUE) {
  rel_cols <- c(
    # GFP
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepA_TMTRuns_FullPlex.log2_TMT_126_131_median',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepA_TMTRuns_FullPlex.log2_TMT_127N_131_median',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepB_TMTRuns_FullPlex.log2_TMT_126_131_median',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepB_TMTRuns_FullPlex.log2_TMT_127N_131_median',
    # PSD95 KO
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepA_TMTRuns_FullPlex.log2_TMT_127C_131_median',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepA_TMTRuns_FullPlex.log2_TMT_128N_131_median',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepA_TMTRuns_FullPlex.log2_TMT_128C_131_median',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepB_TMTRuns_FullPlex.log2_TMT_127C_131_median',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepB_TMTRuns_FullPlex.log2_TMT_128N_131_median',
    # PSD93KO
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepA_TMTRuns_FullPlex.log2_TMT_129N_131_median',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepA_TMTRuns_FullPlex.log2_TMT_129C_131_median',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepB_TMTRuns_FullPlex.log2_TMT_128C_131_median',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepB_TMTRuns_FullPlex.log2_TMT_129N_131_median',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepB_TMTRuns_FullPlex.log2_TMT_129C_131_median',
    # Double KO
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepA_TMTRuns_FullPlex.log2_TMT_130N_131_median',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepA_TMTRuns_FullPlex.log2_TMT_130C_131_median',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepB_TMTRuns_FullPlex.log2_TMT_130N_131_median',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepB_TMTRuns_FullPlex.log2_TMT_130C_131_median',
    # Quality stats
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepA_TMTRuns_FullPlex.numSpectra',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepB_TMTRuns_FullPlex.numSpectra',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepA_TMTRuns_FullPlex.pct_cov',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepB_TMTRuns_FullPlex.pct_cov',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepA_TMTRuns_FullPlex.unique_peptides',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepB_TMTRuns_FullPlex.unique_peptides',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepA_TMTRuns_FullPlex.TMT_131_total',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepB_TMTRuns_FullPlex.TMT_131_total',
    # 'numSequencesSubgroupSpecificCIDP',
    # 'numSequencesSubgroupSpecificCI',
    'accession_number',
    'accession_numbers',
    'geneSymbol',
    'numPepsUnique',
    'entry_name')
  data <- read.csv(filename)[,rel_cols]
  colnames(data)[1:4] <- c("GFP_A1", "GFP_A2", "GFP_B1", "GFP_B2")
  colnames(data)[5:9] <- c("KO95_A1", "KO95_A2", "KO95_A3", "KO95_B1", "KO95_B2")
  colnames(data)[10:14] <- c("KO93_A1", "KO93_A2", "KO93_B1", "KO93_B2", "KO93_B3")
  colnames(data)[15:18] <- c("DKO_A1", "DKO_A2", "DKO_B1", "DKO_B2")
  colnames(data)[19:26] <- c(
      "numSpectra_A", "numSpectra_B",
      "pctCov_A", "pctCov_B",
      "uniquePeps_A", "uniquePeps_B",
      "totRef_A", "totRef_B")

  data <- data[complete.cases(data),]

  if (group_genes) {
    # Group by gene symbol, keep protein with greatest number of peptides
    # First order dataframe by number of peptides (decreasing order)
    data <- data[order(- data$numSpectra_A - data$numSpectra_B),]
    # Then drop duplicates
    data <- data[!duplicated(data$geneSymbol),]
  }

  if (med_norm) {
    # Median normalize each ratio column
    medians = unlist(lapply(data[,1:18], median))
    data[,1:18] <- sweep(data[,1:18], 2, medians, '-')
  }
  return(data)
}


read_peptide <- function(filename='peptide_raw.csv', med_norm=TRUE) {
  rel_cols <- c(
      # GFP
      'log2_TMT_126_131_mean',
      'log2_TMT_127N_131_mean',
      'log2_TMT_126_131_mean.1',
      'log2_TMT_127N_131_mean.1',
      # PSD 95 KO
      'log2_TMT_127C_131_mean',
      'log2_TMT_128N_131_mean',
      'log2_TMT_128C_131_mean',
      'log2_TMT_127C_131_mean.1',
      'log2_TMT_128N_131_mean.1',
      # PSD 93 KO
      'log2_TMT_129N_131_mean',
      'log2_TMT_129C_131_mean',
      'log2_TMT_128C_131_mean.1',
      'log2_TMT_129N_131_mean.1',
      'log2_TMT_129C_131_mean.1',
      # Double KO
      'log2_TMT_130N_131_mean',
      'log2_TMT_130C_131_mean',
      'log2_TMT_130N_131_mean.1',
      'log2_TMT_130C_131_mean.1',
      # Quality stats
      'num_Spectra',
      'num_Spectra.1',
      'score',
      'score.1',
      # 'percent_scored_peak_intensity',
      # 'percent_scored_peak_intensity.1',
      # Other
      'sequence',
      'modifications',
      'geneSymbol',
      'accession_number',
      'accession_numbers',
      'protein_group_num',
      'protein_score',
      'entry_name')
  data <- read.csv(filename)[,rel_cols]
  colnames(data)[1:4] <- c("GFP_A1", "GFP_A2", "GFP_B1", "GFP_B2")
  colnames(data)[5:9] <- c("KO95_A1", "KO95_A2", "KO95_A3", "KO95_B1", "KO95_B2")
  colnames(data)[10:14] <- c("KO93_A1", "KO93_A2", "KO93_B1", "KO93_B2", "KO93_B3")
  colnames(data)[15:18] <- c("DKO_A1", "DKO_A2", "DKO_B1", "DKO_B2")
  colnames(data)[19:24] <- c(
      "numSpectra_A", "numSpectra_B",
      "score_A", "score_B"
  )

  data <- data[complete.cases(data),]
  if (med_norm) {
    # Median normalize each ratio column
    medians = unlist(lapply(data[,1:18], median))
    data[,1:18] <- sweep(data[,1:18], 2, medians, '-')
  }

  return(data)
}
