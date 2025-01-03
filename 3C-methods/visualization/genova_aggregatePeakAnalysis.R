##########
# Perform APA and view results with GENOVA
# 12/12/23
# Noah Burget
##########
library(GENOVA2)

## Your objects should be stored in a named list 'explist'
## Let's define some parameters for the plots
CHROM <- '8' # Chromosome name to plot
START <-  131e6# Start position to plot
END <-  132e6# End position to plot
SCORE.LIMS <- c(0,20) # The range of scores to plot
DIFF.LIMS <- c(-10,10) # The range of differences to plot when plotting difference matrices
plot.legend <- TRUE
dist.thresh <- c(200e3, Inf) # Range of loop sizes to include in aggregate calculation
outlier.filter <- c(0,1) # Quantiles of data to be used as thresholds. [0,1] performs no outlier filtering
## Quantification parameters
metric <- 'median' # Metric to be used when quantifying aggregates
size <- 3 # The dimensions of the square at the center of the pileups that serve as the 'foreground'
shape <- "circle" # The shape of the central foreground

# Quantification parameters
quant.metric <- 'median' # Metric used when calculating stats, 'mean' or 'median'

## For APA, we need a set of regions to aggregate.
# The dataframe should be in BEDPE (chrom1,start1,end1,chrom2,start2,end2) format
loops <- '../testing/test_dots.tsv'
loops <- read.table(loops, header=T, sep='\t')

##### WARNING: Be sure that the chromosome names in your Contact object matches those in the loops dataframe!
# The below code will replace 'chr' with empty string, making the chromosome names '1','2','3',etc
loops$chrom1 <- gsub('chr', '', loops$chrom1)
loops$chrom2 <- gsub('chr', '', loops$chrom2)

## Define data to be plotted
c1 <- explist[[1]]
c2 <- explist[[2]]
# Or, you can use the explist to plot multiple

########## Perform APA ##########
apa.res <- GENOVA2::APA(explist,
                        dist_thres = dist.thresh,
                        bedpe=loops,
                        outlier_filter = outlier.filter)
### Visualize the results of APA
visualise(apa.res, 
          title = 'Aggregate Peak Analysis',
          metrix = metric)

### Quantify results of APA
# The 'per_loop' feature of the quantification object returns the fold changes, background, and foreground values per loop of c1/c2, i.e. q$per_loop
# 'per_sample' will aggregate these and report fold changes, background, and foreground values, i.e. q$per_sample
q <- quantify(apa.res,
         size = size,
         metric = metric,
         shape = shape)

