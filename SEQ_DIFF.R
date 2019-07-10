# __________
#  Tools :
# _________
Merging_function <- function(l_data, dataRef){
  colnames(dataRef)[1] <- "Sample_ID"
  res_l_data <- list()
  for (i in 1:length(l_data)){
    c_data <- l_data[[i]]
    colnames(c_data)[1] <- "Sample_ID"
    if (dim(c_data)[1] != dim(dataRef)[1]){
      warning(paste("Warning : the data frame[" ,  i ,")] doesn't have the same number of lines than `dataRef`. An inner join will be effecte") )
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

# ------------------------------

############################################################################################
Seq_calcul <- function(l_data, dataRef, listK){
  # __________ Clusters initialization ______
  no_cores <- detectCores() # - 1
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  # _________________________________________
  global_seq_list <- list()
  for (I in 1:length(l_data)){
    c_data <- l_data[[I]]
    colnames(c_data)[1] <- 'Sample_ID'  ; colnames(dataRef)[1] <- 'Sample_ID'
    if (dim(c_data)[1] != dim(dataRef)[1]){
      warning("The number of lines between `c_data` and `dataRef` differs. A merge will be effected.")
    }
    else if (sum(c_data[, 1] == dataRef[, 1]) != length(c_data[, 1])){
      warning("Sample_IDs in `c_data` and `dataRef` are not the same, or are not in the same order. An inner join will be effected.")
    }
    data_m <- merge(c_data, dataRef, by = 'Sample_ID')
    ncol_c_data <- dim(c_data)[2]
    c_data <- data_m[, 1:dim(c_data)[2]]
    dataRef <- data_m[, (dim(c_data)[2]+1):dim(data_m)[2]]
    dataRef <- cbind(data_m[, 1], dataRef)
    colnames(dataRef)[1] <- 'Sample_ID'
    #________________ Distances matrices __________
    dist1 <- as.matrix(dist(c_data[, 2:dim(c_data)[2]], method = "euclidian", diag = TRUE, upper = TRUE))
    rownames(dist1) <- as.character(c_data[ ,1])
    colnames(dist1) <- as.character(c_data[ ,1])
    
    dist2 <- as.matrix(dist(dataRef[, 2:dim(dataRef)[2]], method = "euclidian", diag = TRUE, upper = TRUE))
    rownames(dist2) <- as.character(dataRef[ ,1])
    colnames(dist2) <- as.character(dataRef[ ,1])
    # ____________________________________________
    seq_c_data <- foreach(i=1:length(listK),.combine=rbind) %dopar% {
    k <- listK[i]
    colnames(c_data)[1] <- 'Sample_ID'  ; colnames(dataRef)[1] <- 'Sample_ID'
    if (dim(c_data)[1] != dim(dataRef)[1]){
      warning("The number of lines between `c_data` and `dataRef` differs. A merge will be effected")
    }
    else if (sum(c_data[, 1] == dataRef[, 1]) != length(c_data[, 1])){
      warning("Sample_IDs in `c_data` and `dataRef` are not the same, or are not in the same order. An inner join will be effected.")
    }
    data_m <- merge(c_data, dataRef, by = 'Sample_ID')
    ncol_c_data <- dim(c_data)[2]
    c_data <- data_m[, 1:dim(c_data)[2]]
    dataRef <- data_m[, (dim(c_data)[2]+1):dim(data_m)[2]]
    dataRef <- cbind(data_m[, 1], dataRef)
    colnames(dataRef)[1] <- 'Sample_ID'
    #________________ Distances matrices __________
    dist1 <- as.matrix(dist(c_data[, 2:dim(c_data)[2]], method = "euclidian", diag = TRUE, upper = TRUE))
    rownames(dist1) <- as.character(c_data[ ,1])
    colnames(dist1) <- as.character(c_data[ ,1])
      
    dist2 <- as.matrix(dist(dataRef[, 2:dim(dataRef)[2]], method = "euclidian", diag = TRUE, upper = TRUE))
    rownames(dist2) <- as.character(dataRef[ ,1])
    colnames(dist2) <- as.character(dataRef[ ,1])
    # ____________________________________________

    seq_diff_l <- c()
    n <- dim(dist1)[1]
    for (i in 1:n){
      c_point <- rownames(dist1)[i]
      N1_dist_l <- list(dist1[i, ])[[1]]
      N2_dist_l <- list(dist2[i, ])[[1]]
        
      names(N1_dist_l) <- rownames(dist1)
      names(N2_dist_l) <- rownames(dist2)
      N1_dist_l <- sort(N1_dist_l)
      N2_dist_l <- sort(N2_dist_l)
        
      N1_rank_l <- seq(length(N1_dist_l))
      N2_rank_l <- seq(length(N2_dist_l))
      names(N1_rank_l) <- names(N1_dist_l)
      names(N2_rank_l) <- names(N2_dist_l)
        
      N1_rank_l <- N1_rank_l[1:k]
      N2_rank_l <- N2_rank_l[1:k]
        
      N1_df <- data.frame("Sample_ID" = names(N1_rank_l) , "Rank1" = N1_rank_l)
      N2_df <- data.frame("Sample_ID" = names(N2_rank_l) , "Rank2" = N2_rank_l)
        
      All_neighbors <- unique(c(as.character(N1_df$Sample_ID),as.character(N2_df$Sample_ID)))
        
      s1 = 0
      s2 = 0
      for (j in 1:length( All_neighbors)){
        if (All_neighbors[j] %in%  N1_df$Sample_ID & All_neighbors[j] %in%  N2_df$Sample_ID ){
          N1_index_j = which(N1_df$Sample_ID  == All_neighbors[j]  )
          N2_index_j = which(N1_df$Sample_ID  == All_neighbors[j]  )
          s1 = s1 + ((k - N1_df$Rank1[N1_index_j]) * abs(N1_df$Rank1[N1_index_j] - N2_df$Rank2[N2_index_j]))
          s2 = s2 + ((k - N2_df$Rank2[N2_index_j]) * abs(N1_df$Rank1[N1_index_j] - N2_df$Rank2[N2_index_j]))
        }
        else if (All_neighbors[j] %in%  N1_df$Sample_ID){ 
          N1_index_j = which(N1_df$Sample_ID  == All_neighbors[j]  )
          s1 = s1 + ((k - N1_df$Rank1[N1_index_j]) * abs(N1_df$Rank1[ N1_index_j]))
          s2 = s2 
        }
        else{ 
          N2_index_j = which(N2_df$Sample_ID  == All_neighbors[j]  )
          s1 = s1 
          s2 = s2 + ((k - N2_df$Rank2[N2_index_j]) * abs( N2_df$Rank2[N2_index_j]))
        }
      }
      S = 0.5 * s1 + 0.5 * s2
      seq_diff_l <- c(seq_diff_l,  S)
      }
      seq_diff_k_df <- data.frame('Sample_ID' = c_data$Sample_ID, 'K' = rep(k, length(c_data$Sample_ID)), 'Seq' = seq_diff_l)
      seq_diff_k_df
    }
    global_seq_list[[I]] <- seq_c_data
  }
  stopCluster(cl)
  return(global_seq_list)
}
############################################################################################

############################################################################################
Seq_main <- function(l_data, dataRef, listK, colnames_res_df = NULL , filename = NULL , graphics = FALSE, stats = FALSE){
  l_data <- Merging_function(l_data, dataRef)
  L_data <- list()
  for (i in 1:length(l_data[[1]])){
    L_data[[i]] <- l_data[[1]][[i]]
  }
  dataRef <- l_data[[2]]
  l_data <- L_data
  
  if (length(l_data) == 1 & stats == TRUE){
    warning("Statistical option are not available if `l_data` length is equal to 1.")
    stats =  FALSE
  }
  global_seq_list <- Seq_calcul(l_data , dataRef , listK )
  for (i in 1:length(global_seq_list)){
    global_seq_list[[i]] <- global_seq_list[[i]][complete.cases(global_seq_list[[i]]), ]
  }
  
  # _______________ Writing _________________
  df_to_write <- data.frame('Sample_ID' = global_seq_list[[1]]$Sample_ID, 'K' = global_seq_list[[1]]$K )
  for (i in 1:length(global_seq_list)){
    df_to_write <- cbind(df_to_write, global_seq_list[[i]]$Seq)
    colnames(df_to_write)[dim(df_to_write)[2]] <- paste('V', i, sep="")
  }
  if (is.null(colnames_res_df) == FALSE){
    colnames(df_to_write)[3:length(colnames(df_to_write))] <- colnames_res_df
  }
  if (is.null(filename) == FALSE) {
    if (file.exists(as.character(filename))){
      warning("The filename gives as argument exist in the current directory, this name will be 'incremented'.")
      c = 2
      while(file.exists(as.character(filename))){
        filename <- paste(filename, c, sep = "" )
        c = c+1
      }
    }
    write.table(df_to_write, file = filename, sep = "\t")
  }
  
  data_Seq <- df_to_write
  data_diff_mean_k <- data.frame("k" =  unique(data_Seq$K))
  for (j in seq(from = 3, to = dim(data_Seq)[2], by = 1)) {
    mean_by_k <- tapply(data_Seq[, j], data_Seq$K, mean)
    data_diff_mean_k <- cbind(data_diff_mean_k, mean_by_k)
  }
  colnames(data_diff_mean_k)[2:length(colnames(data_diff_mean_k))] <- colnames(data_Seq)[3:dim(data_Seq)[2]]
  if (graphics == FALSE & stats == FALSE){
    return(list('Seq_df' = df_to_write, 'Seq_mean_by_k' = data_diff_mean_k))
  }
  if (graphics == TRUE){
    p <- Seq_graph_by_k('nothing', Names=colnames_res_df, list_col=NULL, data_diff_mean_K = data_diff_mean_k)
  }
  else{ # graphics == False
    p <- 0 # Only to respect algorithm structure
  }
  if (graphics == TRUE & stats == FALSE){
    return(list('Seq_df' = df_to_write, 'Seq_mean_by_k' = data_diff_mean_k, 'graphics' = p))
  }
  if (stats == TRUE){
    if(dim(data_diff_mean_k)[2] == 2){
      warning("Statics cannot be computed if length list of `l_data` is smaller than two.")
      if (graphics == TRUE){
        return(list('Seq_df' = df_to_write, 'Seq_mean_by_k' = data_diff_mean_k, 'graphics' = p))
      }
      else{
        return(list('Seq_df' = df_to_write, 'Seq_mean_by_k' = data_diff_mean_k))
      }
    }
    if(dim(data_diff_mean_k)[2] == 3){
      WT = wilcox.test(data_diff_mean_k[,2], data_diff_mean_k[, 3], paired = TRUE)
      print(WT)
      if (graphics == TRUE){
        return(list('Seq_df' = df_to_write, 'Seq_mean_by_k' = data_diff_mean_k, 'graphics' = p, 'stats' = list("WT" = WT)))
      }
      else{
        return(list('Seq_df' = df_to_write, 'Seq_mean_by_k' = data_diff_mean_k, 'stats' = list("WT" = WT)))
      }
    }
    else{
      pwt_df <- data.frame('mean_seq' = data_diff_mean_k[, 2], 'method'= rep(paste(colnames(data_diff_mean_k)[2], 2, sep = ""), dim(data_diff_mean_k)[1]))
      for (i in 3:dim(data_diff_mean_k)[2]){
        c_df <- data.frame('mean_seq' = data_diff_mean_k[, i], 'method'=rep(paste(colnames(data_diff_mean_k)[i], i, sep = ""), dim(data_diff_mean_k)[1]))
        pwt_df <- rbind(pwt_df, c_df )
      }
      paired_test_m  <- pairwise.wilcox.test(pwt_df$mean_seq, pwt_df$method,  p.adj = "holm",  paired = TRUE)$p.value
      
      if (graphics == TRUE){
        return(list('Seq_df' = df_to_write, 'Seq_mean_by_k' = data_diff_mean_k, 'graphics' = p,  'pWT' = paired_test_m ))
      }
      else{
        return(list('Seq_df' = df_to_write, 'Seq_mean_by_k' = data_diff_mean_k,  'pWT' = paired_test_m ))
      }
    }
    warning('Unexpected request ')
  }
}
############################################################################################

############################################################################################
Seq_graph_by_k  <-function (data_Seq, Names=NULL, list_col=NULL, data_diff_mean_K = NULL, log = FALSE){
  if (is.null(data_diff_mean_K) == TRUE) {
    data_diff_mean_k <- data.frame("k" =  unique(data_Seq$K))
    for (j in seq(from = 3, to = dim(data_Seq)[2], by = 1)) {
      mean_by_k <- tapply(data_Seq[, j], data_Seq$K , mean)
      data_diff_mean_k <- cbind(data_diff_mean_k, mean_by_k)
    }
    colnames(data_diff_mean_k)[2:length(colnames(data_diff_mean_k))] <- colnames(data_Seq)[3:dim(data_Seq)[2]]
    if (is.null(Names) == FALSE){
      if (length(Names) != (dim(data_Seq)[2] - 3)){
        warning("The list of names gave as input doesn't match with the number of curve.")
      }
      else{
        colnames(data_diff_mean_k)[2:length(colnames(data_diff_mean_k))] <- Names
      }
    }
  }
  else{
    data_diff_mean_k <- data_diff_mean_K
  }
  
  if (log == FALSE){
    data_diff_mean_k_graph <- data.frame('k' = data_diff_mean_k$k , 'diff_seq' =( data_diff_mean_k[, 2]), 'Method' = rep(as.character(colnames(data_diff_mean_k)[2]), length(data_diff_mean_k$k)))
    if (dim(data_diff_mean_k)[2]>3){
      for (i in 3:(dim(data_diff_mean_k)[2])){
        c_df <- data.frame('k' = data_diff_mean_k$k , 'diff_seq' = (data_diff_mean_k[, i]), 'Method' = rep(as.character(colnames(data_diff_mean_k)[i]), length(data_diff_mean_k$k)))
        data_diff_mean_k_graph <- rbind(data_diff_mean_k_graph, c_df)
      }
    }
  
    theme_set(theme_bw())
    p <- ggplot(data_diff_mean_k_graph, aes(x=k, y=diff_seq,  color=Method)) + geom_line() + geom_point()+
      scale_colour_manual(values=custom.col[1:length(unique(data_diff_mean_k_graph$Method))])
    p <- p +  labs(title="Sequence difference metric", caption = "Means of sequence difference values k  between each method and the reference one. ",
         y="mean(Seq_diff)", x="K") +theme(plot.title=element_text(size=18, face="bold", color="#17202A", hjust=0.5,lineheight=1.2),  # title
                                                   plot.subtitle =element_text(size=13, color="#17202A", hjust=0.5),  # caption
                                                   plot.caption =element_text(size=10, color="#17202A", hjust=0.5),  # caption
                                                   axis.title.x=element_text(size=12, face="bold"),  # X axis title
                                                   axis.title.y=element_text(size=12, face="bold"),  # Y axis title
                                                   axis.text.x=element_text(size=12),  # X axis text
                                                   axis.text.y=element_text(size=12))  # Y axis text
  print(p)
  return(p)
  }
  else {
    data_diff_mean_k_graph <- data.frame('k' = data_diff_mean_k$k , 'diff_seq' = log( data_diff_mean_k[, 2]), 'Method' = rep(as.character(colnames(data_diff_mean_k)[2]), length(data_diff_mean_k$k)))
    if (dim(data_diff_mean_k)[2]>3){
      for (i in 3:(dim(data_diff_mean_k)[2])){
        c_df <- data.frame('k' = data_diff_mean_k$k , 'diff_seq' = log(data_diff_mean_k[, i]), 'Method' = rep(as.character(colnames(data_diff_mean_k)[i]), length(data_diff_mean_k$k)))
        data_diff_mean_k_graph <- rbind(data_diff_mean_k_graph, c_df)
      }
    }
    theme_set(theme_bw())
    p <- ggplot(data_diff_mean_k_graph, aes(x=k, y=diff_seq,  color=Method)) + geom_line() + geom_point()+
      scale_colour_manual(values=custom.col[1:length(unique(data_diff_mean_k_graph$Method))])
    p <- p +  labs(title="Sequence difference metric", caption = "Means of sequence difference values k  between each method and the reference one. ",
                   y="mean(Seq_diff)", x="K") +theme(plot.title=element_text(size=18, face="bold", color="#17202A", hjust=0.5,lineheight=1.2),  # title
                                                     plot.subtitle =element_text(size=13, color="#17202A", hjust=0.5),  # caption
                                                     plot.caption =element_text(size=10, color="#17202A", hjust=0.5),  # caption
                                                     axis.title.x=element_text(size=12, face="bold"),  # X axis title
                                                     axis.title.y=element_text(size=12, face="bold"),  # Y axis title
                                                     axis.text.x=element_text(size=12),  # X axis text
                                                     axis.text.y=element_text(size=12))  # Y axis text
    print(p)
    return(p)
  } 
}
############################################################################################

############################################################################################
seq_permutation_test <- function(data, data_ref, list_K, n=30, graph = TRUE){
  
  if (n > 30){
    warning("the calcul could be long !")
  }
  colnames(data)[1] <- 'Sample_ID' ; colnames(data_ref)[1] <- 'Sample_ID'
  if (dim(data)[1] != dim(data_ref)[1]){
    warning("Sample IDs don't match between `data` and `data_ref` a merge will be effected.")
  }
  else if( dim(data)[1] == dim(data_ref)[1] & sum(as.character(data[, 1]) == as.character(data_ref[, 1])) != length(data_ref[, 1])){
    warning("Sample IDs don't match between `data` and `data_ref` a merge will be effected.")
  }
  data_m <- merge(data, data_ref, by=c('Sample_ID'))
  data <- data_m[, 1:dim(data)[2]]
  data_ref <- data_m[, (dim(data)[2]+1):dim(data_m)[2]]
  data_ref <- cbind(data_m[, 1], data_ref)
  colnames(data_ref)[1] <- "Sample_ID"
  global_seq_df <- Seq_calcul(list(data), dataRef = data_ref , listK = list_K)[[1]]
  mean_k <- tapply(global_seq_df$Seq, global_seq_df$K, mean)
  main_df <- data.frame('k' = unique(global_seq_df$K) , "means_ref" = mean_k)
  for (i in 1:n){
    data_shuffle <- data[,2:dim(data)[2]]
    data_shuffle <- data_shuffle[,sample(ncol(data_shuffle))]
    data_shuffle <- data_shuffle[sample(nrow(data_shuffle)),]
    data_shuffle <- cbind(data[,1], data_shuffle, row.names = NULL)
    colnames(data_shuffle)[1] <- "Sample_ID"
    Seq_data_A <- Seq_calcul(list(data_shuffle), dataRef = data_ref , listK = list_K)[[1]]
    mean_k <- tapply(Seq_data_A$Seq, Seq_data_A$K, mean)
    main_df <- cbind(main_df , mean_k)
  }
  by_k_alea <- main_df[,3:dim(main_df)[2]]
  Means_alea <- rowMeans(by_k_alea)
  WT  = wilcox.test(main_df[ ,1], Means_alea, alternative = "less")
  print(WT)
  
  theme_set(theme_bw())
  p <- ggplot()
  for (i in 3:(dim(main_df)[2])){
    c_df <- data.frame('k' = main_df[ ,1] , 'main_df' = main_df[ ,i])
    p <- p + geom_line(data = c_df, aes(x=k, y=main_df), colour = '#848484')+geom_point(data = c_df, aes(x=k, y= main_df), colour = '#848484')
  }
  c_df <- data.frame('k' = main_df[ ,1] , 'main_df' = main_df[ ,2])
  p <- p + geom_line(data = c_df, aes(x=k, y = main_df), colour = '#B40404')+geom_point(data = c_df, aes(x=k, y=main_df), colour = '#B40404')
  
  c_MA_df <- data.frame('k' = main_df[ ,1] , 'main_df' = Means_alea)
  p <- p + geom_line(data = c_MA_df, aes(x=k, y = main_df), colour = '#388E3C')+geom_point(data = c_MA_df, aes(x=k, y=main_df), colour ='#388E3C')
  
  p <- p +  labs(title="Significance test of the Sequence difference metric",
                 y=TeX("$\\bar{SD_k} = \\frac{1}{N} \\sum_{i=1}^N SD_{i,k}$"), x="K") +theme(plot.title=element_text(size=18, face="bold", color="#17202A", hjust=0.5,lineheight=1.2),  # title
                                                   plot.subtitle =element_text(size=13, color="#17202A", hjust=0.5),  # caption
                                                   plot.caption =element_text(size=10, color="#17202A", hjust=0.5),  # caption
                                                   axis.title.x=element_text(size=12, face="bold"),  # X axis title
                                                   axis.title.y=element_text(size=12, face="bold"),  # Y axis title
                                                   axis.text.x=element_text(size=12),  # X axis text
                                                   axis.text.y=element_text(size=12))  # Y axis text
  print(p)
  

  return(list(WT,p))
}


