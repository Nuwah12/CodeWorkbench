##########
# Loading in .cool files for plotting with GENOVA
# 12/12/23
# Noah Burget
##########
library(GENOVA2)

### Paths to all cool files
cools <- list("0hr" = "/mnt/data0/noah/3C_methods/3C-Methods/testing/s38_221203_Granta519EBF1KI_cl27_0hr_MicroC_25U_Nova_5000_unbalanced.cool")

### Define some parameters
RES <- 10000 # Matrix resolution
BALANCE <- TRUE # Apply balancing weights to matrix (VC_SQRT in GENOVA2)
scale_bp <- NULL # 'TRUE' if you want internal coverage-based normalization
colors <- c('red') # Colors the samples will take on when plotting
centromeres <- NULL # Dataframe of centromere positions, if null they are inferred from largest span of empty bins

explist <- list() # This is the list that will hold all of your Contacts objects

### Read in cool files and make an experiment list
for(i in 1:length(cools)) {
  explist[[names(cools[i])]] <- GENOVA2::load_contacts(signal_path = cools[[i]], # Path to .cool file
                                        sample_name=names(cools[i]), # Sample name - this is what this sample will be called in downstream analyses
                                        resolution = RES, # Bin size of matrix
                                        balancing = BALANCE, # Whether or not to apply balancing weights to contact counts 
                                        scale_bp = scale_bp, # NULL = no internal coverage normalization
                                        colour = colors[i], # Color of sample - this is the color the same will appear as im downstream analysis
                                        centromeres = centromeres, # Dataframe of centromere locations. If null, centromeres are inferred from the longest stretch of empty bins                          
                                        
  )
}
