library(foreach)
library(doParallel)
library(ggplot2)
library(RColorBrewer)
library(plotly)

# My colors : 
# -----------

custom.col <- c('#2ECC71', # clear green
                 "#33691E", # olive green
                 '#1E90FF', # dark turquoise
                 "#1A237E", #  dark blue
                 '#6C3483', #  dar purple
                 '#D81B60', # Hot pink
                 "#D16103", # dark orange
                 "#FFD700", # gold
                 '#B22222', # brick red
                 '#626567', # dark gray
                 "#17202A") # dark gray-blue 
custom.col <<- append(custom.col, brewer.pal(n = 12, name = "Paired"))
cols <- function(a) image(1:23, 1, as.matrix(1:23), col=custom.col, axes=FALSE , xlab="", ylab="")
a <- 1:23
cols(23)

# _________________




