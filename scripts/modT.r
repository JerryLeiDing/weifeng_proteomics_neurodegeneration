
# $Id: modT.r 31 2013-01-16 16:49:47Z manidr $
print ("$Id: modT.r 31 2013-01-16 16:49:47Z manidr $", quote=FALSE)


library (limma)
# scatterhist.r or pairs-plot.r needed when plots are required


##
## Example usage
##
#  Data: itraq.csv contains the replicate (log) ratios for a specific comparison
#  Results will be written to itraq-results.csv (table) and itraq-results.pdf (plots)
#  Files needed: modT.r, scatterhist.r or pairs-plot.r (for plots, 2 replicates or more, resp.)
#  To run the moderated T test:
#
#  source ('modT.r')
#  modT.test ('itraq.csv', 'itraq-results',       -- these two arguments are required
#             id.col='id',                        -- to specify id columns in itraq.csv;
#                                                    if not present, one will be created with 1,2,...,N
#             data.col=NULL,                      -- if NULL all columns in itraq.csv after removing id.col is
#                                                    data; else, use these specified columns
#             p.value.alpha=0.05,                 -- p-value threshold below which points will be marked on plot
#             use.adj.pvalue=TRUE,                -- if FALSE nominal p-value will be used for threshold
#             apply.log=FALSE,                    -- if TRUE, data will be log2 transformed
#             na.rm=TRUE,                         -- if TRUE rows with ANY replicate missing will be removed
#             plot=TRUE,                          -- by default no plots are generated; set to TRUE for plotting
#             pairs.plot.2rep=FALSE,              -- if TRUE will create pairs plot even for 2 reps
#             limits=c(-3,3),                     -- plotting axis limits; default uses the range in data
#             xlab='ratio', ylab='ratio',         -- x and y axis labels
#             main='Title',                       -- plot title
#             plot.col='black',                   -- color for data points in scatter plot
#             subset.col='green',                 -- regulated items will be plotted in this color (default: red)
#             hist.breaks=20,                     -- higher number for finer histogram
#             hist.col='grey',                    -- histogram color
#             cex.cor=2,                          -- size for pearson correlation value
#             pch=18, cex=0.7                     -- other graphics options
#  )
# NB: The plotting options above are for using pairs-plot.r (with 3 or more replicates)
#     For 2 replicates, scatterhist.r is used and the options are slightly different;
#     see function definition for scatterhist.ci in scatterhist.r for options.
#
#




##
## run moderated t-test, and plot results
## mainly for iTRAQ, but can be used of other data
##
modT_test <- function (data.file, output.prefix, id.col=NULL, data.col=NULL,
                       p.value.alpha=0.05, use.adj.pvalue=TRUE, apply.log=FALSE,
                       na.rm=FALSE, nastrings=c("NA", "<NA>", "#NUM!", "#DIV/0!", "#NA", "#NAME?"), 
                       plot=FALSE, pairs.plot.2rep=FALSE, limits=NULL, xlab="", ylab="", dframe=FALSE, design=NULL,...) {
  #
  # data.file should contain one peptide in each row.
  # The columns contain the normalized log-ratio from each replicate
  # (technical or biological). The ratio is based on classes of interest
  # that need to be distinguished: i.e., ratio = intensity_A / intensity_B;
  # this test calculates the p-value for determining if peptide p is
  # differentially regulated between classes A and B (i.e., if the log
  # ratio is different from 0).
  # While the standard scatter plot routine can only handle 2 replicates,
  # a pairs plot is created when there are more than 2 replicates.
  # The moderated t-test can be applied to any number of replicates.
  # An id column can be optionally included in the data.file to track
  # peptides (row numbers are used as id if a column is not specified).
  #
  # graphics can be controlled using ...
  #  when using scatterhist (for 2 replicates), this can be arguments to points
  #  when > 2 replicates are present, ... can include arguments to points in 
  #   addition to: plot.col, subset.col, hist.col, hist.breaks,
  #                prefix (for correlation), cex.cor

  # read data file 
  if (!dframe) {
    d <- read.csv (data.file, na.strings=nastrings)
    if (na.rm) {
      # ... and remove rows with any missing values
      no.na <- apply (d, 1, function (x) { all (!is.na (x))} )
      d <- d [no.na,]
    }
  } else {
    d <- data.file
  }

  # extract id column
  # if id column is not specified, use row numbers
  row.num <- 1:nrow (d)
  if (is.null (id.col)) {
    d <- cbind (id=row.num, d)
    id.col <- c ('id')
  }

  id.col <- make.names (id.col)
  id <- d[,id.col]

  # extract data columns
  if (is.null (data.col)) data <- d [, setdiff (colnames (d), id.col)]
  else data <- d [, make.names (data.col)]
  

  # log transform is required
  if (apply.log) data <- log2 (data)

  # moderated t test
  sig <- moderated.t (data, design)
  mod.t.result <- topTable (sig, number=nrow(data), sort.by='none', confint=TRUE)
  if (use.adj.pvalue) mod.sig <- mod.t.result [,'adj.P.Val'] <= p.value.alpha
  else  mod.sig <- mod.t.result [,'P.Value'] <= p.value.alpha
  change <- apply (data, 1,
                   function (x) {
                     x <- x [is.finite (x)]
                     ret.value <- '?'
                     if ( all (x < 0) ) ret.value <- 'down'
                     else if ( all (x > 0)) ret.value <- 'up'
                     return (ret.value)
                   })
  mod.t <- data.frame ( cbind (data.frame (id), data, mod.t.result, change=change, significant=mod.sig) )

  # plot and/or determine prediction CI's, etc.
  # the output/behavior of scatterhist.ci can be controlled using ... for specifying args
  results <- NULL
  if (plot) {
    if (ncol (data)==2 && !pairs.plot.2rep) {
      if (!exists ('libdir')) source ('scatterhist.r')    # source only when not in GenePattern
      # when 2 replicates are present, create a scatterplot with marginal histograms
      pdf (paste (output.prefix, ".pdf", sep=''), width=6, height=6, pointsize=12)
      results <- scatterhist.ci (data[,1], data[,2], id=row.num, ci.level=ci.level,
                                 special.subsets=list (mod.t[,'significant']),
                                 subset.colors=c('red'), limits=limits, xlab=xlab, ylab=ylab, ...)
      if (length (setdiff (colnames (results), c ('x','y','id'))) != 0)
        results <- results [, setdiff (colnames (results), c ('x','y'))]  # remove x and y -- data already has this
      else results <- NULL
      dev.off ()
    } else if (ncol (data) >= 2) {
      # create a pairs plot when more than 2 replicates are present
      if (!exists ('libdir')) source ('pairs-plot.r')     # source only when not in GenePattern
      keep <- !is.na (mod.sig)
      plot.data <- data [keep,]
      significant <<- mod.sig [keep]

      if (is.null (limits)) {
        # use full range of x, y if limits not specified
        r <- apply (plot.data, 2, range, na.rm=TRUE)
        limits <- c ( min(r), max(r) )
      }
      
      size.in <- ncol (data) + 1
      pdf (paste (output.prefix, ".pdf", sep=''), height=size.in, width=size.in, pointsize=8)      
      pairs (plot.data, diag.panel=panel.hist, lower.panel=panel.scatterplot, upper.panel=panel.cor,
             subset=significant, xlim=limits, ylim=limits, ...)
      dev.off()
    } else warning ("No plots generated.")
  }

  # write out / return results
  if (is.null (results)) final.results <- mod.t
  else final.results <- merge (results, mod.t, by='id')
  # write.csv (final.results, paste (output.prefix, ".csv", sep=''), row.names=FALSE)

  return (final.results)
}





##
## Support functions
##


# Moderated t-test for significance testing
moderated.t <- function (data, design=NULL) {
  # data is a table with rows representing peptides/proteins/genes 
  # and columns representing replicates
        
  data.matrix <- data.frame (data)
  if (is.null(design)) {
    # nb: the design matrix for lmFit is the default (unit vector)
    #     which treats each column in data as a replicate
    m <- lmFit (data.matrix)
  } else {
    if (is.vector(design)) {design = model.matrix(~factor(design))}
    m <- lmFit (data.matrix, design=design)
  }
  m <- eBayes (m)
 
  return(m)
}





