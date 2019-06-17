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


################################ SPATIAL AUTO CORRELATION ANALYSIS #############################################if (!require("rspatial")) devtools::install_github('rspatial/rspatial')

spatial_cor_main <-function(l_coords_data , spatial_data, listK, nsim = 99){
  # Merging 
  colnames(spatial_data)[1] <- "Sample_ID"
  L_coords_data <- list()
  Spa_data <- list() # With the same sample IDs than the ones in L_coords_data[[i]]
  for (i in 1:length(l_coords_data)){
    c_data <- l_coords_data[[i]]
    colnames(c_data)[1] <- "Sample_ID"
    if (dim(c_data)[1] != dim(spatial_data)[1]){
      warning(paste("The number of rows between the coordinate data frame [", i,"] and the spatial attribute data frame differs. A inner join is going to be done.", sep=""))
    }
    else if (sum(c_data[, 1] == spatial_data[, 1]) != dim(c_data)[1]){
      warning(paste("Warning : Sample_IDs between the cordinate data frame [", i,"] and the spatial attribute data frame are not the same, or are not in the same order. An inner join is going to be effected.", sep =""))
    }
    data_m <- merge(spatial_data, c_data, by= "Sample_ID")
    c_data <- data_m[, (dim(spatial_data)[2]+1):dim(data_m)[2]]
    c_data <- cbind(data_m[, 1], c_data)
    colnames(c_data)[1] <- "Sample_ID"
    L_coords_data[[i]] <- c_data
    c_spatial_data <- data_m[, 1:dim(spatial_data)[2]]
    Spa_data[[i]] <- c_spatial_data
  }
  l_coords_data <- L_coords_data
  MI_array <- array(rep(NA, length(l_coords_data)* (dim(spatial_data)[2]-1)*length(listK)), dim=c(length(l_coords_data), (dim(spatial_data)[2]-1), length(listK)))
  MS_array <- array(rep(NA, length(l_coords_data)* (dim(spatial_data)[2]-1)*length(listK)), dim=c(length(l_coords_data), (dim(spatial_data)[2]-1), length(listK)))
  for (i in 1:length(l_coords_data)){
    c_data <- l_coords_data[[i]]
    spatial_data <- Spa_data[[i]]
    c_sample_id <- spatial_data[ ,1]
    if (dim(c_data)[2] == 3){
      c_data <- as.matrix(c_data[, 2:dim(c_data)[2]])
      print(head(c_data))
    }
    for (c_k in 1:length(listK)){
      if (dim(c_data)[2] == 2){
        k_neigh <- knn2nb(knearneigh(c_data, k=listK[c_k], RANN=FALSE))
        print(k_neigh)
        ww <- nb2listw(k_neigh, style='B')
      }
      for (j in 2:dim(spatial_data)[2]){
        if (dim(c_data)[2] == 2){
          MI <- moran(spatial_data[, j], ww, n=length(ww$neighbours), S0=Szero(ww))
          MS <- moran.mc(spatial_data[, j], ww, nsim=99)
          MI_array[i,(j-1),c_k] <- MI$I
          MS_array[i,(j-1),c_k] <- MS$p.value
        }
        else{
          c_spatial_data <- data.frame("Sample_ID"= as.character(c_sample_id), "att" = spatial_data[, j] )
          MI <- moran_index_HD(data = c_data, spatial_att = c_spatial_data, K = listK[c_k], merge = NULL)
          MI_array[i,(j-1),c_k] <- MI
          MS <- moran_stat_HD(obs_moran_I = MI, data = c_data, k = listK[c_k], nsim = 500, spatial_att = c_spatial_data)
          MS_array[i,(j-1),c_k] <- MS
          
        }
      } 
    }
  }
  return(list(MI_array, MS_array))
}


######################################################


moran_index_HD <- function(data, spatial_att, K, merge = NULL){
  
   if (is.null(merge) == FALSE){
    colnames(spatial_att)[1] <- 'Sample_ID'
    colnames(data)[1] <- 'Sample_ID'
    data_m <- merge(spatial_att, data, by= "Sample_ID")
    c_data <- data_m[, (dim(spatial_att)[2]+1):dim(data_m)[2]]
    c_data <- cbind(data_m[, 1], c_data)
    c_spatial_att <-  data_m[,1:dim(spatial_att)[2]]
    colnames(c_data)[1] <- "Sample_ID"
    colnames(c_spatial_att)[1] <- "Sample_ID"
    data <- c_data
    spatial_att <- c_spatial_att
   }
  KNN_R  = get.knn(data[, 2:dim(data)[2]], k=K, algorithm=c( "brute"))
  KNN_R = KNN_R$nn.index
  m_neigh <- matrix(0, ncol = dim(KNN_R)[1], nrow =dim(KNN_R)[1])
  for (i in 1:dim(KNN_R)[1]){
    for (j in 1:dim(KNN_R)[2]){
      n_index = KNN_R[i,j]
      m_neigh[i,n_index] = 1
    }
  }
  n <- length(spatial_att[, 2])
  y <- spatial_att[, 2]
  ybar <- mean(y)
  dy <- y - ybar
  g <- expand.grid(dy, dy)
  yiyj <- g[,1] * g[,2]
  pm <- matrix(yiyj, ncol=n)
  pmw <- pm * m_neigh
  spmw <- sum(pmw)
  smw <- sum(wm)
  sw <-spmw/smw
  vr <- n / sum(dy^2)
  MI <- vr * sw
  return(MI)
}


moran_stat_HD <- function(obs_moran_I , data, spatial_att, k , nsim = 99){
  MI_rand <- c()
  for (s in 1:nsim){
    KNN_R <- matrix(0, nrow = dim(data)[1], ncol = k)
    for (I in 1:dim(data)[1]){
      KNN_R[I, ] <- floor(runif(k, min=1, max=(dim(data)[1]+1)))
    }
    m_neigh <- matrix(0, ncol = dim(KNN_R)[1], nrow =dim(KNN_R)[1])
    for (i in 1:dim(KNN_R)[1]){
      for (j in 1:dim(KNN_R)[2]){
        n_index = KNN_R[i,j]
        m_neigh[i,n_index] = 1
      }
    }
    n <- length(spatial_att[, 2])
    y <- spatial_att[, 2]
    ybar <- mean(y)
    dy <- y - ybar
    g <- expand.grid(dy, dy)
    yiyj <- g[,1] * g[,2]
    pm <- matrix(yiyj, ncol=n)
    pmw <- pm * m_neigh
    spmw <- sum(pmw)
    smw <- sum(wm)
    sw <-spmw/smw
    vr <- n / sum(dy^2)
    MI <- vr * sw
    MI_rand<- c(MI_rand,MI)
  }
  MI_rand_df = data.frame('MI' = MI_rand)
  n_greater <- length(MI_rand_df$MI[MI_rand_df$MI > obs_moran_I])
  PV <- n_greater / nsim
  return(PV)
}


#############################################

spa_VISTA <- data.frame("Sample_ID" = as.character(spatial_att$Sample_ID) , "VISTA" = spatial_att$VISTA)
moran_index_HD(data= r_meso_df, spatial_att= spa_VISTA, K = 190, merge =FALSE)
moran_index_HD(data= r_meso_df, spatial_att= spa_VISTA, K =20, merge =TRUE)
moran_index_HD(data= r_meso_df, spatial_att= spa_VISTA, K =20, merge =NULL)



PCA_coords_order <- test_spa_coords[[1]]
PCA_coords_order$Sample_ID <- as.character(PCA_coords_order$Sample_ID)

spatial_vista_order <- test_spa_coords[[2]]
spatial_vista_order$Sample_ID <- as.character(spatial_vista_order$Sample_ID)

moran_index_HD(data= PCA_coords_order, spatial_att= spatial_vista_order, K = 20, merge =T)
moran_index_HD(data= PCA_coords_df, spatial_att= spa_VISTA, K = 20, merge =T)


r_num_meso_df <- apply(R_meso_df[,2:dim(R_meso_df)[2]], 2, as.numeric )
r_meso_df <- cbind(R_meso_df[,1],r_num_meso_df )
r_meso_df <- as.data.frame(r_meso_df)
for (i in 2:dim(r_meso_df)[2]){
  r_meso_df[ ,i] <-as.numeric(r_meso_df[ ,i])
}




PCA_coords_order <- test_spa_coords[[1]]
PCA_coords_order$Sample_ID <- as.character(PCA_coords_order$Sample_ID)

spatial_vista_order <- test_spa_coords[[2]]
spatial_vista_order$Sample_ID <- as.character(spatial_vista_order$Sample_ID)

moran_index_HD(data= PCA_coords_order, spatial_att= spatial_vista_order, K = 20, merge =T)
moran_index_HD(data= PCA_coords_df, spatial_att= spa_VISTA, K = 20, merge =T)


r_num_meso_df <- apply(R_meso_df[,2:dim(R_meso_df)[2]], 2, as.numeric )
r_meso_df <- cbind(R_meso_df[,1],r_num_meso_df )
r_meso_df <- as.data.frame(r_meso_df)
for (i in 2:dim(r_meso_df)[2]){
  r_meso_df[ ,i] <-as.numeric(r_meso_df[ ,i])
}




spatial_att <- read.table("Meso_spatial_cor_att.tsv", header = T)
pca_coords_xy <- as.matrix(PCA_coords_df[, 2:3])
k3 <- knn2nb(knearneigh(pca_coords_xy, k=20, RANN=FALSE))
k3


wm <- nb2mat(k3, style='B') ; wm

ww <- nb2listw(k3, style='B') ; ww
ww$weights
ww$neighbours

MI <- moran(spatial_att$VISTA, ww, n=length(ww$neighbours), S0=Szero(ww))
moran.test(spatial_att$VISTA, ww, randomisation=FALSE)
MS <- moran.mc(spatial_att$VISTA, ww, nsim=99)



spa_VISTA <- data.frame("Sample_ID" = as.character(spatial_att$Sample_ID) , "VISTA" = spatial_att$VISTA)
moran_index_HD(data= r_meso_df, spatial_att= spa_VISTA, K = 190, merge =FALSE)
moran_index_HD(data= r_meso_df, spatial_att= spa_VISTA, K =20, merge =TRUE)
moran_index_HD(data= r_meso_df, spatial_att= spa_VISTA, K =20, merge =NULL)

 
r_meso_df[,1]<- as.character(r_meso_df[,1])
MITEST  <- spatial_cor_main(l_coords_data = list(PCA_coords_df, UMAP_coords_NN230, r_meso_df), spatial_data = spatial_att, listK = c(20, 200), nsim = 101 )
spatial_att_for_sim <- data.frame('Sample_ID' = spatial_att$Sample_ID,  'vista'= spatial_att$VISTA)
moran_stat_HD(obs_moran_I = 0.005, data = PCA_coords_df, k = 20, nsim = 99, spatial_att = spatial_att_for_sim)

