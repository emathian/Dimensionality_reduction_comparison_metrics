library(foreach)
library(doParallel)

CP2 <- function(l_data , list_K , dataRef = NULL , colnames_res_df = NULL , filename = NULL){
  # This function calculs the centrality preservation values, according
  # a list of K values. The typical uses should be the calcul of CP values
  # on the low dimendional space, according data1, and their determination 
  # in the high dimensional space, data2. Results could be written in a 
  # text file.
  # Arguments : data1 and data2 must be arrays of floats and rownames may
  # contain samples' ID. Example:
  # data 1 :  Sample_ID | x | y | ...
  # Furthermore according the previous remark if data2 is defined, it must
  # contain the same number of samples, i.e. the same number of lines.
  
  # __________ Clusters initialization ______
  no_cores <- detectCores() # - 1
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  # _________________________________________
  
  
  # __________ Distance matrices ____________
  l_dist = list()
  for (i in 1:length(l_data)){
    c_data <- l_data[i][[1]]
    dist <- as.matrix(dist(c_data[, 2:dim(c_data)[2]], method = "euclidian", diag = TRUE, upper = TRUE)) 
    rownames(dist) <- as.character(c_data[ ,1])
    colnames(dist) <- as.character(c_data[ ,1]) 
    l_dist[[i]] = dist
  }
  # _________________________________________
  
  # No parallel version
  # __________
  #np.time <- system.time({
  #  list_CP <- list()
  #  for (i in 1:length(l_data)){
  #    c_data <- l_data[[i]]
  #    c_dist <- l_dist[[i]]
  #    CP_c_data <- data.frame()
  #    for (k in list_K){
  #      CP <- data.frame("Sample_ID" = as.character(c_data[, 1]), "CP" = rep(0, length(c_data[, 1])), "K"=rep(k, length(c_data[, 1])))
  #     N <- matrix(ncol = k, nrow = dim(c_data)[1])
  #      rownames(N) <- c_data[, 1] ; colnames(N) <- seq(k)
  #     for (j in 1:dim(c_dist)[1]){
  #        N_dist <- list(c_dist[i, ])[[1]]
  #        names(N_dist) <- c_data[ ,1]
  #       N_dist <- sort(N_dist)
  #        KNeighbors_N <- as.character(names(N_dist)[1:k])
  #        N[j,] <- KNeighbors_N
  #      }
  #     for (l in 1:dim(N)[1]){
  #        c_point <- rownames(N)[l]
  #        CP_j <- 0
  #       for(j in 1:dim(N)[1]){
  #          neighbors_j <- list(N[l, ])[[1]]
  #          if (c_point %in% neighbors_j ){
  #            CP_j = CP_j + (k - match(c_point,neighbors_j))
  #          }
  #        }
  #        CP[l,2] =  CP_j
  #      }
  #      CP_c_data <- rbind(CP_c_data, CP)
  #    }
  #   list_CP[[i]] <- CP_c_data
  #  }
  #})
  #print("nopar")
  #print(np.time)
  
  
  #print(np.time)
  # __________
  
  # Parallel version
  p.time <- system.time({
    list_CP <- list()
    len_list_CP <- list()
    for (I in 1:length(l_data)){
      c_data <- l_data[[I]]
      c_dist <- l_dist[[I]]
      #CP_c_data <- data.frame()
      CP_c_data <- foreach(i=1:length(list_K),.combine=rbind) %dopar% {
        k = list_K[i]
        CP <- data.frame("Sample_ID" = as.character(c_data[, 1]), "CP" = rep(0, length(c_data[, 1])), "K"=rep(k, length(c_data[, 1])))
        N <- matrix(ncol = k, nrow = dim(c_data)[1])
        rownames(N) <- c_data[, 1] ; colnames(N) <- seq(k)
        for (j in 1:dim(c_dist)[1]){
          N_dist <- list(c_dist[i, ])[[1]]
          names(N_dist) <- c_data[ ,1]
          N_dist <- sort(N_dist)
          KNeighbors_N <- as.character(names(N_dist)[1:k])
          N[j,] <- KNeighbors_N
        }
        for (l in 1:dim(N)[1]){
          c_point <- rownames(N)[l]
          CP_j <- 0
          for(j in 1:dim(N)[1]){
            neighbors_j <- list(N[l, ])[[1]]
            if (c_point %in% neighbors_j ){
              CP_j = CP_j + (k - match(c_point,neighbors_j))
            }
          }
          CP[l,2] =  CP_j
        }
        CP
      }
      list_CP[[I]] <- CP_c_data
      len_list_CP[[I]] <- length(CP_c_data[ ,1])
    }
  })
  print("par")
  print(p.time)
  stopCluster(cl)
  
  if (length(unique(len_list_CP))==1){
    df_to_write <- data.frame('Sample_ID' = list_CP[[1]]$Sample_ID, 'K' = list_CP[[1]]$K )
    for (i in 1:length(list_CP)){
      df_to_write <- cbind(df_to_write, list_CP[[i]]$CP)
      colnames(df_to_write)[dim(df_to_write)[2]] <- paste('V', i, sep="")
    }  
  }
  else{
    print("warning :  Input data frames don't the same number of lines a inner join will be done. ")
    df_to_write <- data.frame("Sample_ID"=list_CP[[1]]$Sample_ID, 'K'= list_CP[[1]]$K, 'V1'= list_CP[[1]]$CP)
    for (i in 2:length(list_CP)){
      df_to_write <- merge(df_to_write, list_CP[[i]],  by=c('Sample_ID','K'))
      colnames(df_to_write)[dim(df_to_write)[2]] <- paste('V', i, sep="")
    }  
  }
  
  if (is.null(colnames_res_df) == FALSE) { 
    print('here')
    colnames(df_to_write)[3:length(colnames(df_to_write))] <- colnames_res_df
  }
  
  if (is.null(filename) == FALSE) {
    if (file.exists(as.character(filename))){
      print("Warning : The filename gives as argument exist in the current directory, this name will be 'incremented'.")
      c = 2
      while(file.exists(as.character(filename))){
        filename <- paste(filename, c, sep = "" )
        c = c+1
      }
    }
    write.table(df_to_write, file = filename, sep = "\t")
  }

  
  
  return(df_to_write)
}


L = CP2(list( UMAP_coords_NN230[1:120, ], PCA_coords_df, TM_coords_df) , c(20,50,49), TM_coords_df, c("a","b","c") ,"aname")

LDATA <- list( UMAP_coords_NN230, PCA_coords_df, TM_coords_df)
#def centrality_preservation(dist1 , dist2 , K , filename):
  

PCA_coords_df = read.table("Meso_pca_coords.tab", sep="\t",header = T)
TM_coords_df= read.table("Meso_tm_coords_v2.tab", sep="\t", header = T)
UMAP_coords_NN230 <- read.table( "umap_coords_nn230.tab" ,  sep = "\t", dec="." , header = TRUE,   quote="")

CP(PCA_coords_df , c(3,5), TM_coords_df, FALSE)


Dist <- dist(PCA_coords_df[, 2:3], method = "euclidian", diag = TRUE, upper = TRUE) ; Dist#
m <- as.matrix(Dist) ;  m
rownames(m) <- as.character(PCA_coords_df[ ,1])
colnames(m) <- as.character(PCA_coords_df[ ,1])  
m


CP_R_MOFA_expr<- read.table("CP_pulmo_R_MOFA.txt", sep = "\t", dec="." , header = TRUE,   quote="")
plot(CP_R_MOFA_expr$CP2)
hist(CP_R_MOFA_expr$CP2)

x <- matrix(rnorm(100), nrow = 5)
dist(x)
dist(x, diag = TRUE)
dist(x, upper = TRUE)
m <- as.matrix(dist(x))




d <- as.dist(m)