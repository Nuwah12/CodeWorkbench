##########
# Script for storing ggplot variables for styles I like
##########
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
})

#####
# Palettes 
#####


#####
# Themes
#####
min.theme <- function(color="black", xtick=TRUE, )
{
  if (xtick)
  {
    return(theme(axis.ticks.y = element_line(color=color),
                 axis.text.y = element_text(color=color),
                 axis.ticks.x = element_line(color=color),
                 axis.text.x = element_text(color=color),
                 panel.background = element_blank(),
                 axis.line = element_line(color=color)))
  }
  else
  {
    return(theme(axis.ticks.y = element_line(color=color),
                 axis.text.y = element_text(color=color),
                 axis.ticks.x = element_blank(),
                 axis.text.x = element_text(color=color),
                 panel.background = element_blank(),
                 axis.line = element_line(color=color)))
  }
}
