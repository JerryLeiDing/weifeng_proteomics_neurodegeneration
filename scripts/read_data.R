library(plyr)
library(ggplot2)

# Plots the distribution of every column in the df
# Calculates mean, std dev, and interquartile range
plot_data_quality <- function(df, filename='raw_quality.png') {
  if (is.null(dim(df))) {
    # Invalid input
    stop("Invalid input to `plot_data_quality`: no dimensions")
  }
  nrows = ceiling(dim(df)[2]/2)

  width = 1080
  height = 480 * nrows
  # Save to file
  png(filename, width=width, height=height)
  # Two columns of subplots
  par(mfrow=c(nrows, 2))

  for (name in colnames(df)) {
    col <- df[,name]

    # Calculate mean, std dev, and interquartile range
    mn <- mean(col, na.rm=T)
    md <- median(col)
    stdev <- sd(col, na.rm=T)
    # 2.25% - 97.75% range (=4x stddev for normal)
    # qt = quantile(col, probs=c(0.02275013, 1-0.02275013))
    # (1x stddev for normal)
    qt = quantile(col, probs=c(0.3085375, 1-0.3085375))
    iq = (qt[2] - qt[1])

    subtitle <- sprintf("Mean=%.3f, Std=%.3f, Median=%.3f, IQ=%.3f", mn, stdev, md, iq)

    # Allow a total of 15 stddevs of x-axis in the graph
    bounds <- c(qt[1] - 7*iq, qt[2] + 7*iq)
    # We want 50 bins between 2% - 98% range
    bins <- seq(bounds[1], bounds[2], (15*iq)/50)
    # Now add extreme breakpoints to bins to capture all data
    bins <- c(min(col), bins, max(col))

    TXT_SCALE <- 2
    h <- hist(col, breaks=bins, include.lowest=TRUE, xlim=bounds, freq=TRUE,
         main=name, sub=subtitle, xlab="",
         cex.lab=TXT_SCALE, cex.axis=TXT_SCALE, cex.main=TXT_SCALE, cex.sub=TXT_SCALE)

    # And plot a normal fit
    scale <- h$counts / h$density
    curve(dnorm(x, md, qt[2] - qt[1])*scale[25],
          add=TRUE,
          col='darkblue',
          lwd=2)
  }
  dev.off()
}


# Encapsulates the median normalization, gene grouping, and plotting functions
# data_cols: indices or colnames corresponding to data columns
# Will filter by complete cases in data cols
clean_data <- function(data, data_cols, med_norm=TRUE, log=TRUE, 
                       group_genes=FALSE, plot_filename=NULL) {
  # Filter to rows where all data cols are complete
  incomplete <- !complete.cases(data[,data_cols])
  if (any(incomplete)) {
    warning(paste('Dropping', sum(incomplete), 'rows with NA values'))
    data <- data[!incomplete,]
  }

  if (log) {
    data[,data_cols] <- log2(data[,data_cols])
  }

  if (group_genes) {
    # Then drop duplicate geneSymbols
    # TODO: should this be median?
    data <- data[!duplicated(data$geneSymbol),]
  }

  if (!is.null(plot_filename)) {
    print(head(data[,data_cols]))
    plot_data_quality(data[,data_cols], plot_filename)
  }

  if (med_norm) {
    # Median normalize each ratio column
    medians = unlist(lapply(data[,data_cols], median))
    data[,data_cols] <- sweep(data[,data_cols], 2, medians, '-')
  }

  return(data)
}

# TODO: fix me to work with clean_data?
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
 
  # TODO: use clean_data instead of existing functions 
  if (FALSE) {
    data <- clean_data(data, c(1:9), med_norm=T, log=T, group_genes=F,
                       plot_filename='figures/raw_quality/Aug_psm.png')
  }

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

read_protein_aug <- function(filename='protein_raw.csv') {
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
  data <- read.csv(paste('data/Aug_2017/',filename,sep=''))[,rel_cols]
  colnames(data)[1:4] <- c("GFP_A1", "GFP_A2", "GFP_B1", "GFP_B2")
  colnames(data)[5:9] <- c("KO95_A1", "KO95_A2", "KO95_A3", "KO95_B1", "KO95_B2")
  colnames(data)[10:14] <- c("KO93_A1", "KO93_A2", "KO93_B1", "KO93_B2", "KO93_B3")
  colnames(data)[15:18] <- c("DKO_A1", "DKO_A2", "DKO_B1", "DKO_B2")
  colnames(data)[19:26] <- c(
      "numSpectra_A", "numSpectra_B",
      "pctCov_A", "pctCov_B",
      "uniquePeps_A", "uniquePeps_B",
      "totRef_A", "totRef_B")

  # First order dataframe by number of peptides (decreasing order)
  # So we keep the more abundant protein when grouping genes
  data <- data[order(- data$numSpectra_A - data$numSpectra_B),]
  data <- clean_data(data, c(1:18), med_norm=T, log=F, group_genes=T,
                     plot_filename='figures/raw_quality/Aug_protein.png')

  return(data)
}


# TODO: Fix data file (extra line above for some reason?)
read_peptide_aug <- function(filename='peptide_raw.csv', med_norm=TRUE) {
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
  data <- read.csv('data/Aug_2017/peptides_raw.csv')[,rel_cols]
  colnames(data)[1:4] <- c("GFP_A1", "GFP_A2", "GFP_B1", "GFP_B2")
  colnames(data)[5:9] <- c("KO95_A1", "KO95_A2", "KO95_A3", "KO95_B1", "KO95_B2")
  colnames(data)[10:14] <- c("KO93_A1", "KO93_A2", "KO93_B1", "KO93_B2", "KO93_B3")
  colnames(data)[15:18] <- c("DKO_A1", "DKO_A2", "DKO_B1", "DKO_B2")
  colnames(data)[19:22] <- c(
      "numSpectra_A", "numSpectra_B",
      "score_A", "score_B"
  )

  data <- clean_data(data, c(1:18), med_norm=T, log=F, group_genes=F,
                     plot_filename='figures/raw_quality/Aug_peptides.png')

  return(data)
}


read_protein_raw_aug <- function() {
  rel_cols <- c(
    # GFP
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepA_TMTRuns_FullPlex.TMT_126_total',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepA_TMTRuns_FullPlex.TMT_127N_total',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepB_TMTRuns_FullPlex.TMT_126_total',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepB_TMTRuns_FullPlex.TMT_127N_total',

    'Melanie.PSDKO_TMT10_Weifeng_052217.RepA_TMTRuns_FullPlex.TMT_127C_total',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepA_TMTRuns_FullPlex.TMT_128N_total',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepA_TMTRuns_FullPlex.TMT_128C_total',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepB_TMTRuns_FullPlex.TMT_127C_total',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepB_TMTRuns_FullPlex.TMT_128N_total',

    'Melanie.PSDKO_TMT10_Weifeng_052217.RepA_TMTRuns_FullPlex.TMT_129N_total',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepA_TMTRuns_FullPlex.TMT_129C_total',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepB_TMTRuns_FullPlex.TMT_128C_total',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepB_TMTRuns_FullPlex.TMT_129N_total',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepB_TMTRuns_FullPlex.TMT_129C_total',
    
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepA_TMTRuns_FullPlex.TMT_130N_total',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepA_TMTRuns_FullPlex.TMT_130C_total',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepB_TMTRuns_FullPlex.TMT_130N_total',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepB_TMTRuns_FullPlex.TMT_130C_total',

    'Melanie.PSDKO_TMT10_Weifeng_052217.RepA_TMTRuns_FullPlex.TMT_131_total',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepB_TMTRuns_FullPlex.TMT_131_total',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepA_TMTRuns_FullPlex.numSpectra',
    'Melanie.PSDKO_TMT10_Weifeng_052217.RepB_TMTRuns_FullPlex.numSpectra',
    'accession_number',
    'accession_numbers',
    'geneSymbol',
    'numPepsUnique',
    'entry_name')
  data <- read.csv('data/Aug_2017/protein_raw.csv')[,rel_cols]
  colnames(data)[1:4] <- c("GFP_A1", "GFP_A2", "GFP_B1", "GFP_B2")
  colnames(data)[5:9] <- c("KO95_A1", "KO95_A2", "KO95_A3", "KO95_B1", "KO95_B2")
  colnames(data)[10:14] <- c("KO93_A1", "KO93_A2", "KO93_B1", "KO93_B2", "KO93_B3")
  colnames(data)[15:18] <- c("DKO_A1", "DKO_A2", "DKO_B1", "DKO_B2")
  colnames(data)[19:20] <- c("REF_A1", "REF_B1")
  colnames(data)[21:22] <- c("numSpectra_A", "numSpectra_B")

  data <- data[order(- data$numSpectra_A - data$numSpectra_B),]
  data <- clean_data(data, c(1:20), med_norm=T, log=F, group_genes=T,
                     plot_filename='figures/raw_quality/Aug_protein_nolog.png')
  return(data)
}


read_protein_march <- function() {
  data <- read.csv('data/March_2017/proteins_raw.csv')

  colnames(data) <- c(
        'reference',
        'accession_number',
        'uniprot',
        'geneSymbol',
        'entry_name',
        'P25F_A1',
        'P25F_A2',
        'P25F_A3',
        'P25F_A4',
        'P25F_A5',
        'P25F_A6',
        'CKF_A1',
        'CKF_A2',
        'CKF_A3',
        'CKF_A4',
        'ratio')

  # TODO: if necessary, dummy columns for data which is not included

  # Reorder data so intensities are at front
  data <- data[,c(6:16,1:5)]
  data <- clean_data(data, c(1:10), med_norm=F, log=T, group_genes=T,
                     plot_filename='figures/raw_quality/Mar_protein_raw.png')

  return(data)
}

read_phosphos_march <- function() {
  data <- read.csv('data/March_2017/phospho_raw.csv')
  colnames(data) <- c(
        'reference',
        'accession_number',
        'uniprot',
        'geneSymbol',
        'entry_name',
        'sites',
        'phosphoresidues',
        'motif',
        'maxscore',
        'sequence',
        'P25F_A1',
        'P25F_A2',
        'P25F_A3',
        'P25F_A4',
        'P25F_A5',
        'P25F_A6',
        'CKF_A1',
        'CKF_A2',
        'CKF_A3',
        'CKF_A4',
        'ratio')

  # Reorder data so intensities are first
  data <- data[,c(11:21,1:7,10)]
  data <- clean_data(data, c(1:10), med_norm=F, log=T, group_genes=F,
                     plot_filename='figures/raw_quality/Mar_phospho_raw.png')

  return(data)
}


read_phosphos_sept <- function() {
  data <- read.csv('data/Sept_2017/phospho_raw.csv')
  colnames(data) <- c(
        'reference',
        'accession_number',
        'uniprot',
        'geneSymbol',
        'entry_name',
        'sites',
        'motif',
        'phosphoresidues',
        'protein_sites',
        'maxscore',
        'P25_EE_A1',
        'P25_EE_A2',
        'P25_EE_A3',
        'CT_EE_A1',
        'CT_EE_A2',
        'CT_EE_A3',
        'P25_HC_A1',
        'P25_HC_A2',
        'CT_HC_A1',
        'CT_HC_A2')

  # Reorder data so intensities are first
  data <- data[,c(11:20,9,2:6,8)]
  data <- clean_data(data, c(1:10), med_norm=F, log=T, group_genes=F,
                     plot_filename='figures/raw_quality/Sept_phospho_raw.png')

  return(data)
}
