library(foreach)
library(doParallel)
library(ggplot2)

CP <- function(l_data , list_K , dataRef = NULL , colnames_res_df = NULL , filename = NULL){
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
  if (is.null(dataRef) == FALSE){
    distRef <- as.matrix(dist(c_data[, 2:dim(dataRef)[2]], method = "euclidian", diag = TRUE, upper = TRUE)) 
    rownames(distRef) <- as.character(dataRef[ ,1])
    colnames(distRef) <- as.character(dataRef[ ,1]) 
  }
  # _________________________________________
  
  #_____________ Calculs of CP values _______

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

  if (is.null(dataRef) == FALSE){
    CPRef <- foreach(i=1:length(list_K),.combine=rbind) %dopar% {
      k = list_K[i]
      CPRef_k <- data.frame("Sample_ID" = as.character(dataRef[, 1]), "CP" = rep(0, length(dataRef[, 1])), "K"=rep(k, length(dataRef[, 1])))
      N <- matrix(ncol = k, nrow = dim(dataRef)[1])
      rownames(N) <- dataRef[, 1] ; colnames(N) <- seq(k)
      for (j in 1:dim(distRef)[1]){
        N_dist <- list(distRef[i, ])[[1]]
        names(N_dist) <- dataRef[ ,1]
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
        CPRef_k[l,2] =  CP_j
      }
      CPRef_k
    }
    list_CP[[I+1]] <- CPRef
    len_list_CP[[I+1]] <- length(CPRef[ ,1])
  }
  # _________________________________________
  
  stopCluster(cl)
  
  if (length(unique(len_list_CP))==1){
    df_to_write <- data.frame('Sample_ID' = list_CP[[1]]$Sample_ID, 'K' = list_CP[[1]]$K )
    for (i in 1:length(list_CP)){
      df_to_write <- cbind(df_to_write, list_CP[[i]]$CP)
      colnames(df_to_write)[dim(df_to_write)[2]] <- paste('V', i, sep="")
    }  
  }
  else{
    print("Warning :  Input data frames don't the same number of lines a inner join will be done. ")
    df_to_write <- data.frame("Sample_ID"=list_CP[[1]]$Sample_ID, 'K'= list_CP[[1]]$K, 'V1'= list_CP[[1]]$CP)
    for (i in 2:length(list_CP)){
      df_to_write <- merge(df_to_write, list_CP[[i]],  by=c('Sample_ID','K'))
      colnames(df_to_write)[dim(df_to_write)[2]] <- paste('V', i, sep="")
    }  
  }
  
  if (is.null(colnames_res_df) == FALSE){ 
    print('here')
    if (is.null(dataRef) == FALSE){
      colnames_res_df <- c(colnames_res_df, 'REF_CP')
     }
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





CP_graph_by_k  <-function (data_CP,  ref_CP_data, Names=NULL, list_col=NULL){
  if (dim(data_CP)[1] != dim(ref_CP_data)[1]){
    print("Warning : Input data frames don't have the same number of line")
    break
  }
  else{
    colnames(data_CP)[1] <- "Sample_ID" ; colnames(ref_CP_data)[1] <- "Sample_ID"
    colnames(data_CP)[2] <- "K" ; colnames(ref_CP_data)[2] <- "K"
    data_CP <- cbind(data_CP, ref_CP_data$Real)
    data_diff_mean_k <- data.frame("k" =  unique(data_CP$K))
    if (is.null(Names)==TRUE){
      Names <- colnames(data_CP)[3:length(colnames(data_CP))]
    }
    else if (length(Names) != dim(data_CP)[2]-3){
      print("Warning :  The length of the list of Names is inconsistant in regard of the length of the input data")
      Names <- colnames(data_CP)[3:length(colnames(data_CP))]
      
    }
    for (I in seq(from = 3, to = dim(data_CP)[2]-1, by = 1)) {
      abs_diff <- abs(as.numeric(data_CP[, I]) - as.numeric(data_CP[, dim(data_CP)[2]]))
      c_abs_diff_df <- data.frame("abs_diff" = abs_diff, 'K' = data_CP$K)
      abs_diff_k <- tapply(c_abs_diff_df$abs_diff, c_abs_diff_df$K, mean)
      data_diff_mean_k <- cbind(data_diff_mean_k, abs_diff_k)
    }
    colnames(data_diff_mean_k)[2:length(colnames(data_diff_mean_k))] <- Names
  }
  data_diff_mean_k <- data_diff_mean_k[order(data_diff_mean_k$k), ]
  custom.col <- c( "#F4EDCA","#D16103", "#C4961A",'#B03A2E', "#C3D7A4", "#52854C", "#4E84C4", '#1ABC9C',"#293352", '#6C3483', '#626567', "#17202A") #"#FFDB6D",  
  p1 <- ggplot()
  for (i in 2:dim(data_diff_mean_k)[2]){
    print(i)
    c_df <- data.frame('k'= data_diff_mean_k[,1], 'val'=data_diff_mean_k[ ,i] )
    #colnames(c_df)[2] <- colnames(data_diff_mean_k)[i]
    p1 <- p1 + geom_line(data =  c_df  , aes(y =val  , x = k ), colour=custom.col[i] )
  }
  print(p1)
  #return(data_diff_mean_k)
}  


list_CP_df = list( CP_R_PCA, CP_R_TM ,CP_R_UMAP_MD02)#, CP_R_UMAP_NN150_MD05 , CP_R_UMAP_NN230 , CP_R_UMAP_NN20 , CP_R_UMAP_MD09, CP_PCA_TM ,
Name = c("CP_R_PCA" , "CP_R_TM" ,"CP_R_UMAP_MD02")#, "CP_R_UMAP_NN150_MD05", "CP_R_UMAP_NN230" , "CP_R_UMAP_NN20" , "CP_R_UMAP_MD09", "CP_PCA_TM" 
CP_graph_by_k(list_CP_df,Name)


data_CP <- data.frame("Sample_ID" = CP_R_TM$sample, 'K' = CP_R_TM$K, "PCA" = CP_R_PCA$CP2, 'TM' = CP_R_TM$CP2,
                      'UMAP' = CP_R_UMAP_MD02$CP2, 'UMAP_NN230'= CP_R_UMAP_NN230$CP2, 'CP_R_UMAP_NN20'=CP_R_UMAP_NN20$CP2,
                      'CP_R_UMAP_MD09'= CP_R_UMAP_MD09$CP2, 'CP_R_UMAP_MD02'=CP_R_UMAP_MD02$CP2)
ref_CP_data <- data.frame("Sample_ID" = CP_R_PCA$sample , "K" = CP_R_PCA$K, 'Real' = CP_R_PCA$CPN)
CP_graph_by_k(data_CP,  ref_CP_data, Names=c("PCA_Ref", 'TM_ref', 'Umap_ref'), list_col=NULL)

  #par(mfrow=c(1,2))
  c= 1 
  c_data_frame = as.data.frame(list_data_CP[[1]])
  CP_diff_mean_df = data.frame("k" =unique(c_data_frame$K ))
  for (i in 1:length(list_data_CP)){
    c_data_frame = as.data.frame(list_data_CP[i])
    diff_CP2_CPN = abs(c_data_frame$CP2 - c_data_frame$CPN )
    c_data_frame = cbind(c_data_frame ,diff_CP2_CPN )
    colnames(c_data_frame)[dim(c_data_frame)[2]] <- "Abs_diff"
    Mean_by_k =tapply(c_data_frame$Abs_diff ,c_data_frame$K, mean)
    print( Mean_by_k)
    CP_diff_mean_df <- cbind(CP_diff_mean_df, Mean_by_k)
    colnames(CP_diff_mean_df)[dim(CP_diff_mean_df)[2]] <- Name[i]
  }
  #print(CP_diff_mean_df)

    p1 <- ggplot() 
    for (i in 2:dim(CP_diff_mean_df)[2]){
      print(i)
        c_df <- data.frame('k'= CP_diff_mean_df[,1], 'diff_cp'= CP_diff_mean_df[,i])
        p1 <- p1 + geom_line(data = c_df, aes(y = diff_cp, x = k))
    }
    print(p1)
   # p1 <- p1 + geom_line(aes(y = CP_diff_mean_df[,3], x = CP_diff_mean_df$k), colour = 'blue')
  
#  return(p1)






L = CP(list( UMAP_coords_NN230[1:120, ], PCA_coords_df, TM_coords_df) , c(20,50,49), TM_coords_df, c("a","b","c") ,"aname")

TM_coords  <- read.table("Meso_tm_coords_v2.tab", sep = "\t", dec="." , header = TRUE,   quote="")
colnames(TM_coords)[1] <- "sample"
PCA_coords  <- read.table("Meso_pca_coords.tab", sep = "\t", dec="." , header = TRUE,   quote="")
colnames(PCA_coords)[1] <- "sample"
UMAP_coords_NN230 <- read.table( "umap_coords_nn230.tab" ,  sep = "\t", dec="." , header = TRUE,   quote="")
colnames(UMAP_coords_NN230)[1]<- "sample"
UMAP_coords_NN20 <-  read.table( "umap_coords_nn20.tab" ,  sep = "\t", dec="." , header = TRUE,   quote="")
colnames(UMAP_coords_NN20)[1]<- "sample"
UMAP_coords_MD09 <-  read.table( "umap_coords_md09.tab" ,  sep = "\t", dec="." , header = TRUE,   quote="")
colnames(UMAP_coords_MD09)[1]<- "sample"
UMAP_coords_MD02 <-  read.table( "umap_coord_md_02.tab" ,  sep = "\t", dec="." , header = TRUE,   quote="")
colnames(UMAP_coords_MD02 )[1]<- "sample"
UMAP_coords_NN150_MD05 <- read.table( "umap_coords_nn150_md05.tab" ,  sep = "\t", dec="." , header = TRUE,   quote="")
colnames(UMAP_coords_NN150_MD05 )[1]<- "sample"

CP_PCA_TM <- read.table("CP_PCA_TM_v2.txt", sep = "\t", dec="." , header = TRUE,   quote="")  # HERE
CP_R_PCA <- read.table("CP_PCA_real_v2.txt", sep = "\t", dec="." , header = TRUE,   quote="")
CP_R_TM <- read.table("CP_meso_TM_real_v2.txt", sep = "\t", dec="." , header = TRUE,   quote="")
#CP_R_UMAP_NN150_MD05_nearly_all <- read.table("CP_meso_UMAP_NN150_MD05_real_nearly_all_k.txt", sep = "\t", dec="." , header = TRUE,   quote="")
#CP_R_UMAP_NN150_MD05_Next <- read.table("CP_meso_UMAP_NN150_MD05_real_Next.txt", sep = "\t", dec="." , header = TRUE,   quote="")
#P_R_UMAP_NN150_MD05 <- CP_R_UMAP_NN150_MD05[CP_R_UMAP_NN150_MD05$K %in% unique(CP_R_TM$K) ,  ]
CP_R_UMAP_NN230 <- read.table("CP_meso_UMAP_NN230_real_v2.txt", sep = "\t", dec="." , header = TRUE,   quote="")
CP_R_UMAP_NN20 <- read.table("CP_meso_UMAP_NN20_real_v2.txt", sep = "\t", dec="." , header = TRUE,   quote="")
CP_R_UMAP_MD09 <- read.table("CP_meso_UMAP_MD09_real_v2.txt", sep = "\t", dec="." , header = TRUE,   quote="")
CP_R_UMAP_MD02 <- read.table("CP_meso_UMAP_MD02_real_v2.txt", sep = "\t", dec="." , header = TRUE,   quote="")



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