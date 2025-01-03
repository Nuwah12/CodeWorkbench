##########
# View contact matrices with GENOVA
# 12/12/23
# Noah Burget
##########
library(GENOVA2)
library(ggplot2)

## Your objects should be stored in a named list 'explist'
## Let's define some parameters for the plots
CHROM <- '17' # Chromosome name to plot
START <-  69e6# Start position to plot
END <-  71e6# End position to plot
LOOPS <- NULL # Dataframe represnting a BEDPE file of loop contacts, will be circled on the matrix
tads <- NULL # Dataframe with TAD boundaries to be plotted
SCORE.LIMS <- c(0,20) # The range of scores to plot
DIFF.LIMS <- c(-10,10) # The range of differences to plot when plotting difference matrices
plot.legend <- TRUE

## Define data to be plotted
c1 <- explist[[1]]
c2 <- explist[[2]]

##### Color palette to be used #####
COLORS <- c()
if (length(COLORS) > 0){
  options("GENOVA.colour.palette" = {COLORS}) 
}

########## Plotting a single matrix ##########
### Here, we can plot 1 or 2 contact matrices and present in different forms
# 1: Plot a single SQUARE matrix 
GENOVA2::hic_matrixplot(exp1 = c1,
               chrom = CHROM, start = START, end = END,
               loops = LOOPS,
               tads = tads,
               colour_lim = SCORE.LIMS,
               colour_bar = plot.legend
               )

# 2: Plot a single PYRAMID-form of a matrix
pdf('figures/MB157_SMC1_WO_HiCmatrix_chr1769000000-71000000.pdf')
GENOVA2::pyramid(exp = c1,
                 chrom = CHROM, start = START, end = END,
                 colour = SCORE.LIMS)+
  ggtitle('MB157 WO, chr17:69000000-71000000 @ 10kb')+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# 3: Plot two matrices about the diagonal (coplot = 'dual'). Each half will be labelled with the objects' name attributes
GENOVA2::hic_matrixplot(exp1 = c1,
                        exp2 = c2,
                        chrom = CHROM, start = START, end = END,
                        loops = LOOPS,
                        tads = tads,
                        colour_lim = SCORE.LIMS,
                        colour_bar = plot.legend,
                        coplot = 'dual')
# 4: Plot the difference of two matrices (coplot = 'diff')
GENOVA2::hic_matrixplot(exp1 = c1,
                        exp2 = c2,
                        chrom = CHROM, start = START, end = END,
                        loops = LOOPS,
                        tads = tads,
                        cut.off = SCORE.CUTOFF,
                        coplot = 'diff',
                        colour_bar = plot.legend,
                        colour_lim = DIFF.LIMS)

# 5: Plot the difference of 2 matrices as a pyramid
GENOVA2::pyramid_difference(c1, c2, chrom = CHROM, start = START, end = END, colour = DIFF.LIMS)








