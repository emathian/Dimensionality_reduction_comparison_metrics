
# -----------------
#  Tools function :
# ----------------
Merging_function <- function(l_data, dataRef){
  colnames(dataRef)[1] <- "Sample_ID"
  res_l_data <- list()
  for (i in 1:length(l_data)){
    c_data <- l_data[[i]]
    colnames(c_data)[1] <- "Sample_ID"
    if (dim(c_data)[1] != dim(dataRef)[1]){
      warning(paste("The data frame[" ,  i ,")] doesn't have the same number of lines than `dataRef`. An inner join will be effecte") )
    }
    else if(sum(c_data[, 1] == dataRef[, 1]) != length(c_data[, 1])){
      warning(paste("Sample_IDs in data frame [", i,"] and dataRef are not the same, or are not in the same order. An inner join will be effected.", sep =""))
    }
    data_m <- merge(dataRef, c_data, by = "Sample_ID")
    dataRef <- dataRef[, 1:dim(dataRef)[2]]
    r_data <- data_m[,(dim(dataRef)[2] + 1):dim(data_m)[2]]
    r_data <- cbind(data_m[, 1], r_data)
    colnames(r_data)[1] <- 'Sample_ID'
    res_l_data[[i]] <- r_data
  }
  return(list(res_l_data, dataRef))
}

# ------------------


###########################################################################################################
CP_main <- function(l_data , list_K , dataRef = NULL , colnames_res_df = NULL , filename = NULL , graphics = FALSE, stats = FALSE){
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
    #for (i in 1:length(list_K)){
      k = list_K[i]
      CP <- data.frame("Sample_ID" = as.character(c_data[, 1]), "CP" = rep(0, length(c_data[, 1])), "K"=rep(k, length(c_data[, 1])))
      N <- matrix(ncol = k, nrow = dim(c_data)[1])
      rownames(N) <- c_data[, 1] ; colnames(N) <- seq(k)
      for (j in 1:dim(c_dist)[1]){
        N_dist <- list(c_dist[j, ])[[1]]
        names(N_dist) <- c_data[ ,1]
        N_dist <- sort(N_dist)
        KNeighbors_N <- as.character(names(N_dist)[1:k])
        N[j,] <- KNeighbors_N
      }
      print("head de N")
      print(head(N))
      for (l in 1:dim(N)[1]){
        c_point <- rownames(N)[l] # select a point J
        CP_j <- 0
        for(j in 1:dim(N)[1]){ # pour chaque i
          neighbors_j <- list(N[j, ])[[1]] # considerer le kème viosinage de i
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
        N_dist <- list(distRef[j, ])[[1]]
        names(N_dist) <- dataRef[ ,1]
        N_dist <- sort(N_dist)
        KNeighbors_N <- as.character(names(N_dist)[1:k])
        N[j,] <- KNeighbors_N
      }
      for (l in 1:dim(N)[1]){
        c_point <- rownames(N)[l]
        CP_j <- 0
        for(j in 1:dim(N)[1]){
          neighbors_j <- list(N[j, ])[[1]]
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
  
  # __________________ Writing option and CP Data Frame __________
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
      df_to_write <- merge(df_to_write, list_CP[[i]],  by=c('Sample_ID'))
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
  # ____________________________________________________
  
  
  # ___________________  NO Stats, NO graphic, NO Ref _____
  if ( is.null(dataRef) == TRUE){
    if (graphics == TRUE | stats == TRUE ){
      print("Warning:  Neither graphics nor sattistics could be computed if any reference data frame is defined ")
    }
    return(list("CP_Data_frame" = cp_df))
  }
  else{ # (is.null(dataRef) == FALSE)
    data_CP <- df_to_write
    data_diff_mean_k <- data.frame("k" =  unique(data_CP$K))
    for (I in seq(from = 3, to = dim(data_CP)[2]-1, by = 1)) {
      abs_diff <- abs(as.numeric(data_CP[, I]) - as.numeric(data_CP[, dim(data_CP)[2]]))
      c_abs_diff_df <- data.frame("abs_diff" = abs_diff, 'K' = data_CP$K)
      abs_diff_k <- tapply(c_abs_diff_df$abs_diff, c_abs_diff_df$K, mean)
      data_diff_mean_k <- cbind(data_diff_mean_k, abs_diff_k)
    }
    colnames(data_diff_mean_k)[2:length(colnames(data_diff_mean_k))] <- colnames(data_CP)[3:(dim(data_CP)[2]-1)]
    print('data_diff_mean_k')
    print(head(data_diff_mean_k))
    if (graphics == FALSE & stats == FALSE){
      return(list("CP_Data_frame" = data_CP,'CP_Diff_mean_by_K' = data_diff_mean_k))
    }
    if (graphics == TRUE){
      data_diff_mean_k_graph <- data.frame('k' = data_diff_mean_k$k , 'diff_cp' = data_diff_mean_k[, 2], 'Method' = rep(as.character(colnames(data_diff_mean_k)[2]), length(data_diff_mean_k$k)))
      if (dim(data_diff_mean_k)[2] > 2){
        for (i in 3:(dim(data_diff_mean_k)[2])){
          c_df <- data.frame('k' = data_diff_mean_k$k , 'diff_cp' = data_diff_mean_k[, i], 'Method' = rep(as.character(colnames(data_diff_mean_k)[i]), length(data_diff_mean_k$k)))
          data_diff_mean_k_graph <- rbind(data_diff_mean_k_graph, c_df)
        }
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
      if (dim(data_diff_mean_k)[2] == 2){ # Only one method was define
        if (stats == TRUE){
          print('Warning : Statistics cannot be compute if only one method was given as input e.g `l_data` contains only one element.')
        }
        return(list("CP_Data_frame" =  data_CP,'CP_Diff_mean_by_K' = data_diff_mean_k, "Graphic" = p))
      }
    }
    if (graphics == TRUE & stats == FALSE){
      return (list("CP_Data_frame" = data_CP, 'CP_Diff_mean_by_K' = data_diff_mean_k ,"Graphic" = p))
    }
    else{
      cp_df <- df_to_write
      # Cheeck if there if more than two ech
      if (dim(data_diff_mean_k)[2] == 3){
        WT =wilcox.test(data_diff_mean_k[,2], data_diff_mean_k[ ,3])
        print(WT)
      }
      # Kruskal test
      else{
        ks_df <- data.frame('mean_diff_cp' = data_diff_mean_k[, 2], 'method'= rep(colnames(data_diff_mean_k)[2], dim(data_diff_mean_k)[1]))
        
        for (i in 3:dim(data_diff_mean_k)[2]){
          c_df <- data.frame('mean_diff_cp' = data_diff_mean_k[, i], 'method'=rep(colnames(data_diff_mean_k)[i], dim(data_diff_mean_k)[1]))
          ks_df <- rbind(ks_df, c_df )
          KST = kruskal.test(mean_diff_cp ~ method, data = ks_df)
          print(KST)
        }
        paired_test_m <- matrix(nrow = (dim(data_diff_mean_k)[2]-1) , ncol = (dim(data_diff_mean_k)[2]-1))
        for (i in 2:dim(data_diff_mean_k)[2]){
          for (j in 2:dim(data_diff_mean_k)[2]){
            if (j < i){
              c_WT <- wilcox.test(data_diff_mean_k[,i], data_diff_mean_k[,j])
              paired_test_m[(i-1),(j-1)] <- c_WT$p.value
            }
          }
        }
        colnames(paired_test_m) <- colnames(data_diff_mean_k)[2:dim(data_diff_mean_k)[2]]
        rownames(paired_test_m) <- colnames(data_diff_mean_k)[2:dim(data_diff_mean_k)[2]]
        paired_test_m[is.na(paired_test_m)]   <- '-'
        print(paired_test_m)
      }
      if (graphics == FALSE & stats == TRUE & dim(data_diff_mean_k)[2] == 3){
        return (list("CP_Data_frame" = cp_df, 'CP_Diff_mean_by_K' = data_diff_mean_k , 'Wilcoxon_test' = WT))
      }
      else if (graphics == TRUE & stats == TRUE & dim(data_diff_mean_k)[2] == 3){
        return (list("CP_Data_frame" = cp_df, 'CP_Diff_mean_by_K' = data_diff_mean_k , 'Wilcoxon_test' = WT, 'Graphic' = p))
      }
      else if (graphics == FALSE & stats == TRUE & dim(data_diff_mean_k)[2] != 3){
        return (list("CP_Data_frame" = cp_df, 'CP_Diff_mean_by_K' = data_diff_mean_k , 'Kruskal_test' = KST, 'Paired_wilocoxon_test' = paired_test_m ))
      }
      else{
        return (list("CP_Data_frame" = cp_df, 'CP_Diff_mean_by_K' = data_diff_mean_k , 'Kruskal_test' = KST, 'Paired_wilocoxon_test' = paired_test_m , 'Graphic' = p))
      }
    }
  }
}
###########################################################################################################
List_projection <- list(data.frame(acp_fig1_li_df), data.frame(umap_md01_res_df))

gene_expr_df_filter <- merge(gene_expr_df, acp_fig1_li_df[,1], by = "Sample_ID")

Main_CP_res <- CP_main(l_data = List_projection , list_K = c(20,100) , dataRef = gene_expr_df_filter , colnames_res_df = c("pca", "umap_md=0.1") , filename = NULL , graphics = TRUE, stats = TRUE)


CP_K = Main_CP_res[[1]][which(Main_CP_res[[1]]$K ==20  | Main_CP_res[[1]]$K ==100),]
CP_PCA = CP_K[,1:3]
CP_map(CP_PCA , acp_fig1_li_df, list(20,100), Title = 'CP map :  resulting from the gene expression and pojected on UMAP layout ')


CP_K_UMAP =CP_K[,1:2]
CP_K_UMAP = cbind(CP_K_UMAP ,"UMAP_CP" =CP_K[,4])
head(CP_K_UMAP)
CP_map(CP_K_UMAP ,  data.frame(umap_md07_res_df), list(20,100), Title = 'CP map :  resulting from the gene expression and pojected on UMAP layout ')

CP_map(CP_K_R, acp_fig1_li_df, list(66,111), Title = 'Centrality preservation map for PCA')


###########################################################################################################
CP_graph_by_k  <-function (data_CP,  ref_CP_data, Names=NULL, list_col=NULL){
  if (dim(data_CP)[1] != dim(ref_CP_data)[1]){
    print("Warning : Input data frames don't have the same number of line")
    break
  }
  else{
    colnames(data_CP)[1] <- "Sample_ID" ; colnames(ref_CP_data)[1] <- "Sample_ID"
    colnames(data_CP)[2] <- "K" ; colnames(ref_CP_data)[2] <- "K"
    data_CP <- cbind(data_CP, ref_CP_data[, 3])
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
  if (dim(data_diff_mean_k)[2] > 3){
    for (i in 3:(dim(data_diff_mean_k)[2]-1)){
      c_df <- data.frame('k' = data_diff_mean_k$k , 'diff_cp' = data_diff_mean_k[, i], 'Method' = rep(as.character(colnames(data_diff_mean_k)[i]), length(data_diff_mean_k$k)))
      data_diff_mean_k_graph <- rbind(data_diff_mean_k_graph, c_df)
    }
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
 
  return(p)
}

###########################################################################

CP_calcul <- function(data, list_K, parallel = TRUE ){
  custom.col <- c( '#1E90FF',"#1A237E", '#6C3483','#D81B60',  '#B22222', "#D16103",  "#FFD700",  '#2ECC71',"#33691E", '#626567',"#17202A") 
  
  no_cores <- detectCores() # - 1
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  
  dist <- as.matrix(dist(data[, 2:dim(data)[2]], method = "euclidian", diag = TRUE, upper = TRUE))
  rownames(dist) <- data[ ,1][[1]]
  colnames(dist) <- data[ ,1][[1]]

  CP_data <- foreach(i=1:length(list_K),.combine=rbind) %dopar% {
    k = list_K[i]
    print(k)
    CP <- data.frame("Sample_ID" = as.character(data[, 1]), "CP" = rep(0, length(data[, 1])), "K"=rep(k, length(data[, 1][[1]])))
    print(dim(CP))
    N <- matrix(ncol = k, nrow = dim(data)[1])
    rownames(N) <- data[, 1] ; colnames(N) <- seq(k)
    for (j in 1:dim(dist)[1]){ 
      N_dist <- data.frame("Sample_ID" = as.character(data[, 1]), "dist" = list(dist[j, ])[[1]] ) 
      print('N_dist after sorting')
      N_dist <- N_dist[order(N_dist$dist),]
      KNeighbors_N <- as.character(N_dist$Sample_ID[1:k])
      N[j,] <- KNeighbors_N
      
    }
    for (l in 1:dim(N)[1]){
      c_point <- rownames(N)[l]
      CP_j <- 0
      for(J in 1:dim(N)[1]){
        neighbors_j <- list(N[J, ])[[1]]
        if (c_point %in% neighbors_j ){
          CP_j = CP_j + (k - match(c_point,neighbors_j))
        }
      }
      CP[l,2] =  CP_j
      print(CP_j)
    }
    CP
  }
  
  stopCluster(cl)
  return(CP_data)
}
###########################################################################################################
#library(doSNOW)



###########################################################################################################


Main_CP_res <- CP_calcul(acp_fig1_li_df , list_K = c(100) )

CP_K = Main_CP_res[which(Main_CP_res$K ==100 ),]

CP_K = data.frame("Sample_ID" = CP_K$Sample_ID, "K" = CP_K$K , "CP" = CP_K$CP)
#CP_PCA = CP_K[,1:3]
CP_map(CP_K , acp_fig1_li_df, list(100), Title = 'CP map :  resulting from the gene expression and pojected on UMAP layout ')


CP_K_UMAP =CP_K[,1:2]
CP_K_UMAP = cbind(CP_K_UMAP ,"UMAP_CP" =CP_K[,4])
head(CP_K_UMAP)
CP_map(CP_K_UMAP ,  data.frame(umap_md07_res_df), list(20), Title = 'CP map :  resulting from the gene expression and pojected on UMAP layout ')

CP_map(CP_K_R, acp_fig1_li_df, list(66,111), Title = 'Centrality preservation map for PCA')



###########################################################################################################
CP_permutation_test <- function(data, data_ref, list_K, n=30, graph = TRUE){
  if (n > 30){
    warning(" The calcul could be long !")
  }
  colnames(data)[1] <- 'Sample_ID' ; colnames(data_ref)[1] <- 'Sample_ID'
  if (dim(data)[1] != dim(data_ref)[1]){
    warning(" Sample IDs don't match between `data` and `data_ref` a merge will be effected.")
    data_m <- merge(data, data_ref, by=c('Sample_ID'))
    data <- data_m[, 1:dim(data)[2]]
    data_ref <- data_m[, (dim(data)[2]+1):dim(data_m)[2]]
    data_ref <- cbind(data_m[, 1], data_ref)
  }
  else if( dim(data)[1] == dim(data_ref)[1] & sum(as.character(data[, 1]) == as.character(data_ref[, 1])) != length(data_ref[, 1])){
    warning("Sample IDs don't match between `data` and `data_ref` a merge will be effected.")
    data_m <- merge(data, data_ref, by=c('Sample_ID'))
    data <- data_m[, 1:dim(data)[2]]
    data_ref <- data_m[, (dim(data)[2]+1):dim(data_m)[2]]
    data_ref <- cbind(data_m[, 1], data_ref)

  }

  CP_data <- CP_calcul(data, list_K)
  CP_ref <- CP_calcul(data_ref, list_K)
  abs_diff <- abs(CP_data$CP - CP_ref$CP)
  abs_diff_df <- data.frame('k'= CP_data$K, "abs_diff" = abs_diff)
  abs_diff_k <- tapply(abs_diff_df$abs_diff, abs_diff_df$k, mean)
  main_diff_df <- data.frame('k' = unique(abs_diff_df$k) , "abs_diff_ref" = abs_diff_k)

  for (i in 1:n){
    print(i)
    data_shuffle <- data[,2:dim(data)[2]]
    data_shuffle <- data_shuffle[,sample(ncol(data_shuffle))]
    data_shuffle <- data_shuffle[sample(nrow(data_shuffle)),]
    data_shuffle <- cbind(data[,1], data_shuffle, row.names = NULL)
    colnames(data_shuffle)[1] <- "Sample_ID"
    CP_data_A <- CP_calcul(data_shuffle, list_K)
    abs_diff <- abs(CP_data_A$CP - CP_ref$CP)
    abs_diff_df <- data.frame('k'= CP_data$K, "abs_diff" = abs_diff)
    abs_diff_k <- tapply(abs_diff_df$abs_diff, abs_diff_df$k, mean)
    main_diff_df <- cbind(main_diff_df , abs_diff_k)
  }
  theme_set(theme_bw())
  p <- ggplot()
  for (i in 3:dim(main_diff_df)[2]){
    print(i)
    c_df <- data.frame('k' = main_diff_df[ ,1] , 'abs_diff' = main_diff_df[ ,i])
    p <- p + geom_line(data = c_df, aes(x=k, y=abs_diff), colour = '#848484')+geom_point(data = c_df, aes(x=k, y=abs_diff), colour = '#848484')

  }
  c_df <- data.frame('k' = main_diff_df[ ,1] , 'abs_diff' = main_diff_df[ ,2])
  p <- p + geom_line(data = c_df, aes(x=k, y=abs_diff), colour = '#B40404')+geom_point(data = c_df, aes(x=k, y=abs_diff), colour = '#B40404')
  p <- p +  labs(title="Centrality preservation signicance test",
                 y="mean(|CPi - CP_ref|)", x="K") +theme(plot.title=element_text(size=18, face="bold", color="#17202A", hjust=0.5,lineheight=1.2),  # title
                                                         plot.subtitle =element_text(size=13, color="#17202A", hjust=0.5),  # caption
                                                         plot.caption =element_text(size=10, color="#17202A", hjust=0.5),  # caption
                                                         axis.title.x=element_text(size=12, face="bold"),  # X axis title
                                                         axis.title.y=element_text(size=12, face="bold"),  # Y axis title
                                                         axis.text.x=element_text(size=12),  # X axis text
                                                         axis.text.y=element_text(size=12))  # Y axis text
  print(p)

  by_k_alea <- main_diff_df[,3:dim(main_diff_df)[2]]
  Means_alea <- rowMeans(by_k_alea)
  WT  = wilcox.test( Means_alea,main_diff_df[ ,1])
  print(WT)
  return(WT)
}
###########################################################################################################

###########################################################################################################
CP_map <- function(data_CP, data_coords, listK, Title = NULL){
  # MERGE Data et Data_coords
  colnames(data_CP) <- c("Sample_ID", "K", "CP")
  colnames(data_coords) <- c("Sample_ID","x", "y")
  if (length(unique(data_CP$Sample_ID)) != dim(data_coords)[1]){
    print("Warning : `data_CP` and `data_coords` don't have the same number of line. An inner join will be effected.")
  }
  else if (sum(data_CP$Sample_ID == data_coords$Sample_ID) != dim(data_CP)[1]){
    print("Warning : Sample_IDs in `data_CP` and `data_coords` differ, or are not in the same order. An inner join will be effected.")
  }

  data_m <- merge(data_CP, data_coords, by = 'Sample_ID')
  data_CP <- data_m[, 1:3]
  data_coords <- data_m[, 4:5]
  data_coords <- cbind(data_m[, 1],data_coords)
  colnames(data_coords)[1] <- "Sample_ID"
  data_CP$CP <- as.numeric(data_CP$CP)
  print(dim(data_CP))
  #print(data_CP$CP[min(which(data_CP$K == 21)):max(which(data_CP$K == 21))] == data_CP$CP[min(which(data_CP$K == 1)):max(which(data_CP$K == 1))])

  L_unique_K <- unique(data_CP$K)

  while (length(listK)!=0){

    #  print(data_CP$CP[min(which(data_CP$K == listK[1])):max(which(data_CP$K == listK[1]))] == data_CP$CP[min(which(data_CP$K == listK[2])):max(which(data_CP$K == listK[2]))]  )
    if (listK[1] %in% L_unique_K){
      CP_k1st = data_CP[min(which(data_CP$K == listK[1])):max(which(data_CP$K == listK[1])),]
      CP_k1st_tm <- merge(CP_k1st, data_coords, by="Sample_ID" )
      Title1 = as.character(paste("k = ", as.character(listK[1])))
      p1_seq<- plot_ly(CP_k1st_tm, x = ~x, y = ~y , type="scatter", mode = "markers",
                       marker=list( size=10 , opacity=1), color = ~CP, text = ~paste('Sample: ', Sample_ID))%>%
        layout(annotations=  list(x = 0.2 , y = 1.05, text = Title1, showarrow = F, xref='paper', yref='paper'),  showlegend =FALSE )

      if( (length(listK)-1) < 1 ){
        if (is.null(Title) == FALSE){
          mytitle <- Title
        }
        else{
          mytitle <- "Centrality preservation"
        }
        p_seq <- subplot(p1_seq)%>%
          layout(title = mytitle,  margin = 0.04)
        print(p_seq)
        break
      }
      else{
        CP_k2nd = data_CP[min(which(data_CP$K == listK[2])):max(which(data_CP$K == listK[2])),]
        CP_k2nd_tm <- merge(CP_k2nd, data_coords ,by="Sample_ID" )

        Title2 = as.character(paste("k = ", as.character(listK[2])))
        p2_seq<- plot_ly(CP_k2nd_tm, x = ~x, y = ~y , type="scatter", mode = "markers",
                         marker=list( size=10 , opacity=1), color = ~CP, text = ~paste('Sample: ', Sample_ID))%>%
          layout(annotations=  list(x = 0.2 , y = 1.05, text = Title2, showarrow = F, xref='paper', yref='paper'),  showlegend =FALSE )
      }
      if( (length(listK)-2) < 1 ){
        if (is.null(Title) == FALSE){
          mytitle <- Title
        }
        else{
          mytitle <- "Centrality preservation"
        }
        p_seq <- subplot(p1_seq, p2_seq)%>%
          layout(title = mytitle,  margin = 0.04)
        print(p_seq)
        break
      }
      else{
        CP_k3rd = data_CP[min(which(data_CP$K == listK[3])):max(which(data_CP$K == listK[3])),]
        CP_k3rd_tm <- merge(CP_k3rd, data_coords, by="Sample_ID" )
        Title3 = as.character(paste("k = ", as.character(listK[3])))
        p3_seq<- plot_ly(CP_k3rd_tm, x = ~x, y = ~y , type="scatter", mode = "markers",
                         marker=list( size=10 , opacity=1), color = ~CP, text = ~paste('Sample: ', Sample_ID))%>%
          layout(annotations=  list(x = 0.2 , y = 1.05, text = Title3, showarrow = F, xref='paper', yref='paper'),  showlegend =FALSE )
      }
      if( (length(listK)-3) < 1 ){
        if (is.null(Title) == FALSE){
          mytitle <- Title
        }
        else{
          mytitle <- "Centrality preservation"
        }
        p_seq <- subplot(p1_seq, p2_seq, p3_seq)%>%
          layout(title = mytitle,  margin = 0.04)
        print(p_seq)
        break
      }
      else{
        CP_k4th = data_CP[min(which(data_CP$K == listK[4])):max(which(data_CP$K == listK[4])),]
        CP_k4th_tm <- merge(CP_k4th, data_coords,by="Sample_ID" )
        Title4 = as.character(paste("k = ", as.character(listK[4])))
        p4_seq<- plot_ly(CP_k4th_tm, x = ~x, y = ~y , type="scatter", mode = "markers",
                         marker=list( size=10 , opacity=1), color = ~CP, text = ~paste('Sample: ', Sample_ID))%>%
          layout(annotations=  list(x = 0.2 , y = 1.05, text = Title4, showarrow = F, xref='paper', yref='paper'),  showlegend =FALSE)
        if (is.null(Title) == FALSE){
          mytitle <- Title
        }
        else{
          mytitle <- "Centrality preservation"
        }
        p_seq <- subplot(p1_seq, p2_seq, p3_seq, p4_seq)%>%
          layout(title = mytitle,  margin = 0.04)
        print(p_seq)
        listK <- listK[-c(1:4)]

      }
    }
    else{
      warning( paste("Warning : CP values for this K level (", listK[1],") was not computed."))
      if (length(listK) == 1){
        listK <- list()
      }
      else{
        listK <- listK[2:length(listK)]
      }
    }
  }
  return(p_seq)
}
