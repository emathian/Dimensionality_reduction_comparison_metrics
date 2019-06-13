library(foreach)
library(doParallel)
library(ggplot2)
library(RColorBrewer)

# My colors : 
# -----------

custom.col <- c( "#D16103", # dark orange
                 "#C4961A", # gold
                 '#B03A2E', # brick red
                 '#2ECC71', # clear green
                 "#33691E", # olive green
                 '#1ABC9C', # dark turquoise
                 "#1A237E", #  dark blue
                 '#6C3483', #  dar purple
                 '#D81B60', # Hot pink
                 '#626567', # dark gray
                 "#17202A") # dark gray-blue 
custom.col <<- append(brewer.pal(n = 12, name = "Paired"), custom.col)
cols <- function(a) image(1:23, 1, as.matrix(1:23), col=custom.col, axes=FALSE , xlab="", ylab="")
a <- 1:23
cols(23)

# _________________




CP <- function(l_data , list_K , dataRef = NULL , colnames_res_df = NULL , filename = NULL , graphics = FALSE, stats = FALSE){
  
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
    distRef <- as.matrix(dist(dataRef[, 2:dim(dataRef)[2]], method = "euclidian", diag = TRUE, upper = TRUE)) 
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
    if (is.null(dataRef) == FALSE){
      colnames_res_df <- c(colnames_res_df, 'REF')
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
  if (graphics == FALSE & stats == FALSE){
    cp_df <- df_to_write
    return(cp_df)  
  }
  else{
    data_CP <- df_to_write
    data_diff_mean_k <- data.frame("k" =  unique(data_CP$K))
    for (I in seq(from = 3, to = dim(data_CP)[2]-1, by = 1)) {
      abs_diff <- abs(as.numeric(data_CP[, I]) - as.numeric(data_CP[, dim(data_CP)[2]]))
      c_abs_diff_df <- data.frame("abs_diff" = abs_diff, 'K' = data_CP$K)
      abs_diff_k <- tapply(c_abs_diff_df$abs_diff, c_abs_diff_df$K, mean)
      data_diff_mean_k <- cbind(data_diff_mean_k, abs_diff_k)
    }
    colnames(data_diff_mean_k)[2:length(colnames(data_diff_mean_k))] <- colnames(data_CP)[3:(dim(data_CP)[2]-1)]
   
     if (graphics == TRUE){
      data_diff_mean_k_graph <- data.frame('k' = data_diff_mean_k$k , 'diff_cp' = data_diff_mean_k[, 2], 'Method' = rep(as.character(colnames(data_diff_mean_k)[2]), length(data_diff_mean_k$k)))
      for (i in 3:(dim(data_diff_mean_k)[2])){
        c_df <- data.frame('k' = data_diff_mean_k$k , 'diff_cp' = data_diff_mean_k[, i], 'Method' = rep(as.character(colnames(data_diff_mean_k)[i]), length(data_diff_mean_k$k)))
        data_diff_mean_k_graph <- rbind(data_diff_mean_k_graph, c_df)
      }
      theme_set(theme_bw())
      p <- ggplot(data_diff_mean_k_graph, aes(x=k, y=diff_cp,  color=Method)) + geom_line() + geom_point()+
        scale_colour_manual(values=custom.col[1:length(unique(data_diff_mean_k_graph$Method))])
      p <- p +  labs(title="Centrality preservation", caption = "Means of absolulte differences by k, of CP values' between each method and the reference one. ",
                     y="mean(|CPi - CP_ref|)", x="K") +theme(plot.title=element_text(size=18, face="bold", color="#17202A", hjust=0.5,lineheight=1.2),  # title
                                                             plot.subtitle =element_text(size=13, color="#17202A", hjust=0.5),  # caption
                                                             plot.caption =element_text(size=10, color="#17202A", hjust=0.5),  # caption
                                                             axis.title.x=element_text(size=12, face="bold"),  # X axis title
                                                             axis.title.y=element_text(size=12, face="bold"),  # Y axis title
                                                             axis.text.x=element_text(size=12),  # X axis text
                                                             axis.text.y=element_text(size=12))  # Y axis text
      print(p)
     }
  }
  if (graphics == TRUE & stats == FALSE){
    cp_df <- df_to_write
    return (list(cp_df, p))
  }
  else{
    data_diff_mean_k
  }
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
      Names <- colnames(data_CP)[3:(length(colnames(data_CP))-1)]
    }
    else if (length(Names) != dim(data_CP)[2]-3){
      print("Warning :  The length of the list of Names is inconsistant in regard of the length of the input data")
      Names <- colnames(data_CP)[3:(length(colnames(data_CP))-1)]
    }
    for (I in seq(from = 3, to = dim(data_CP)[2]-1, by = 1)) {
      abs_diff <- abs(as.numeric(data_CP[, I]) - as.numeric(data_CP[, dim(data_CP)[2]]))
      c_abs_diff_df <- data.frame("abs_diff" = abs_diff, 'K' = data_CP$K)
      abs_diff_k <- tapply(c_abs_diff_df$abs_diff, c_abs_diff_df$K, mean)
      data_diff_mean_k <- cbind(data_diff_mean_k, abs_diff_k)
    }
    colnames(data_diff_mean_k)[2:length(colnames(data_diff_mean_k))] <- Names
  }
  
  data_diff_mean_k_graph <- data.frame('k' = data_diff_mean_k$k , 'diff_cp' = data_diff_mean_k[, 2], 'Method' = rep(as.character(colnames(data_diff_mean_k)[2]), length(data_diff_mean_k$k)))
  for (i in 3:(dim(data_diff_mean_k)[2]-1)){
    c_df <- data.frame('k' = data_diff_mean_k$k , 'diff_cp' = data_diff_mean_k[, i], 'Method' = rep(as.character(colnames(data_diff_mean_k)[i]), length(data_diff_mean_k$k)))
    data_diff_mean_k_graph <- rbind(data_diff_mean_k_graph, c_df)
  }
  theme_set(theme_bw())
  p <- ggplot(data_diff_mean_k_graph, aes(x=k, y=diff_cp,  color=Method)) + geom_line() + geom_point()+
    scale_colour_manual(values=custom.col[1:length(unique(data_diff_mean_k_graph$Method))])
  p <- p +  labs(title="Centrality preservation", caption = "Means of absolulte differences by k, of CP values' between each method and the reference one. ",
                 y="mean(|CPi - CP_ref|)", x="K") +theme(plot.title=element_text(size=18, face="bold", color="#17202A", hjust=0.5,lineheight=1.2),  # title
                   plot.subtitle =element_text(size=13, color="#17202A", hjust=0.5),  # caption
                   plot.caption =element_text(size=10, color="#17202A", hjust=0.5),  # caption
                   axis.title.x=element_text(size=12, face="bold"),  # X axis title
                   axis.title.y=element_text(size=12, face="bold"),  # Y axis title
                   axis.text.x=element_text(size=12),  # X axis text
                   axis.text.y=element_text(size=12))  # Y axis text
  print(p)
}  

PCA_coords_df = read.table("Meso_pca_coords.tab", sep="\t",header = T)
TM_coords_df= read.table("Meso_tm_coords_v2.tab", sep="\t", header = T)
UMAP_coords_NN230 <- read.table( "umap_coords_nn230.tab" ,  sep = "\t", dec="." , header = TRUE,   quote="")
R_meso_df = read.table("feature_data_with_lv_2.tsv", sep="\t", header = T)
R_meso_df = t(R_meso_df)
colnames(R_meso_df) <- R_meso_df[1, ]
R_meso_df = R_meso_df[-1, ]
R_meso_df = cbind(rownames(R_meso_df), R_meso_df)
colnames(R_meso_df)[1] <- "Sample_ID"
rownames(R_meso_df) <- NULL

L = CP(list( UMAP_coords_NN230[1:100,], PCA_coords_df, TM_coords_df) , c(20,50,49), R_meso_df, c("a","b","c") ,"aname", TRUE , FALSE)













list_CP_df = list( CP_R_PCA, CP_R_TM ,CP_R_UMAP_MD02)#, CP_R_UMAP_NN150_MD05 , CP_R_UMAP_NN230 , CP_R_UMAP_NN20 , CP_R_UMAP_MD09, CP_PCA_TM ,
Name = c("CP_R_PCA" , "CP_R_TM" ,"CP_R_UMAP_MD02")#, "CP_R_UMAP_NN150_MD05", "CP_R_UMAP_NN230" , "CP_R_UMAP_NN20" , "CP_R_UMAP_MD09", "CP_PCA_TM" 
CP_graph_by_k(list_CP_df,Name)


data_CP <- data.frame("Sample_ID" = CP_R_TM$sample, 'K' = CP_R_TM$K, "PCA" = CP_R_PCA$CP2, 'TM' = CP_R_TM$CP2,
                      'UMAP' = CP_R_UMAP_MD02$CP2, 'UMAP_NN230'= CP_R_UMAP_NN230$CP2, 'CP_R_UMAP_NN20'=CP_R_UMAP_NN20$CP2,
                      'CP_R_UMAP_MD09'= CP_R_UMAP_MD09$CP2, 'CP_R_UMAP_MD02'=CP_R_UMAP_MD02$CP2)
ref_CP_data <- data.frame("Sample_ID" = CP_R_PCA$sample , "K" = CP_R_PCA$K, 'Real' = CP_R_PCA$CPN)
CP_graph_by_k(data_CP,  ref_CP_data, Names=c("PCA_Ref", 'TM_ref', 'Umap_ref'), list_col=NULL)








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



cols <- function(a) image(1:11, 1, as.matrix(1:11), col=custom.col, axes=FALSE , xlab="", ylab="")
a <- 1:11
cols(11)


mypalette<-brewer.pal(7,"set3") 
display.brewer.all(type="div")
display.brewer.pal(8, "Dark2")
display.brewer.pal(12, "Paired")
brewer.pal(n = 12, name = "Paired")





for (i in 1:length(custom.col)){
  plot(CP_R_UMAP_NN230$CP2 , col = custom.col[i])
}


LDATA <- list( UMAP_coords_NN230, PCA_coords_df, TM_coords_df)
#def centrality_preservation(dist1 , dist2 , K , filename):
  


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