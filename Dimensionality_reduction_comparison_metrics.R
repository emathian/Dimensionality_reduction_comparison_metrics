
if(!require(foreach)) install.packages("foreach",repos = "https://CRAN.R-project.org/package=foreach")
if(!require(doParallel)) install.packages("doParallel",repos = "https://CRAN.R-project.org/package=doParallel")
if(!require(ggplot2)) install.packages("ggplot2",repos = "https://CRAN.R-project.org/package=ggplot2")
if(!require(RColorBrewer)) install.packages("RColorBrewer",repos = "https://CRAN.R-project.org/package=RColorBrewer")
if(!require(plotly)) install.packages("plotly",repos = "https://CRAN.R-project.org/package=plotly")
if(!require(spatial)) install.packages("spatial",repos = "https://CRAN.R-project.org/package=plotly")
if(!require(raster)) install.packages("raster",repos = "https://CRAN.R-project.org/package=raster")
if(!require(spdep)) install.packages("spdep",repos = "https://CRAN.R-project.org/package=spdep")
if(!require(FNN)) install.packages("FNN",repos = "https://CRAN.R-project.org/package=FNN")

library(foreach)
library(doParallel)
library(ggplot2)
library(RColorBrewer)
library(plotly)
library(rspatial)
library(spatial)
library(raster)
library(spdep)
library(FNN)
library(rspatial)
library(latex2exp)
library(VennDiagram)
library(venn)

source("../SEQ_DIFF.R") 
source("../CP.R") 
source("../MORAN_I.R")
# My colors : 
# -----------

custom.col <- c( '#1E90FF',"#1A237E", '#6C3483','#D81B60',  '#B22222', "#D16103",  "#FFD700",  '#2ECC71',"#33691E", '#626567',"#17202A") 
custom.col <<- append(custom.col, brewer.pal(n = 12, name = "Paired"))
#cols <- function(a) image(1:23, 1, as.matrix(1:23), col=custom.col, axes=FALSE , xlab="", ylab="")
#a <- 1:23
#cols(23)

# _________________


