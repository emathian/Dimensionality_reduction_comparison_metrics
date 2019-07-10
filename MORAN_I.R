
moran_I_main <-function(l_coords_data , spatial_data, listK, nsim = 500, Stat=FALSE, Graph = FALSE, methods_name = NULL){
  colnames(spatial_data)[1] <- "Sample_ID"
  L_coords_data <- list()
  Spa_data <- list() 
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
    }
    for (c_k in 1:length(listK)){
      if (dim(c_data)[2] == 2){
        k_neigh <- knn2nb(knearneigh(c_data, k=listK[c_k], RANN=FALSE))
        ww <- nb2listw(k_neigh, style='B')
      }
      for (j in 2:dim(spatial_data)[2]){
        if (dim(c_data)[2] == 2){
          MI <- moran(spatial_data[, j], ww, n=length(ww$neighbours), S0=Szero(ww))
          MI_array[i,(j-1),c_k] <- MI$I
          if (Stat != FALSE){
            MS <- moran.mc(spatial_data[, j], ww, nsim=99)
            MS_array[i,(j-1),c_k] <- MS$p.value
          }
        }
        else{
          c_spatial_data <- data.frame("Sample_ID"= as.character(c_sample_id), "att" = spatial_data[, j] )
          MI <- moran_index_HD(data = c_data, spatial_att = c_spatial_data, K = listK[c_k], merge = FALSE)
          MI_array[i,(j-1),c_k] <- MI 
          if (Stat != FALSE){
            MS <- moran_stat_HD(data = c_data, K = listK[c_k], spatial_att = c_spatial_data, obs_moran_I = MI, nsim = nsim)$p.value
            MS_array[i,(j-1),c_k] <- MS
          }
        }
      } 
    }
  }
  for (k in 1:length(listK)){
    if(is.null(methods_name)==F & length(methods_name) == dim(MI_array)[1]){
       rownames(MI_array) <- methods_name
    }
    colnames(MI_array) <- colnames(spatial_data)[2:dim(spatial_data)[2]]
    dimnames(MI_array)[[3]]  <- listK
  }
  if (Graph != FALSE){
    moran_I_scatter_plot(MI_array)
  }
  
  if (Stat != FALSE){
    for (k in 1:length(listK)){
      if(is.null(methods_name)==F & length(methods_name) == dim(MS_array)[1]){
        rownames(MS_array) <- methods_name
      }
      colnames(MS_array) <- colnames(spatial_data)[2:dim(spatial_data)[2]]
      dimnames(MS_array)[[3]]  <- listK
    }
    return(list('MoranIndex' = MI_array,'MoranStat'= MS_array))
  }
  else{
    return(list('MoranIndex' = MI_array))
  }
}

####################################################################
################################################################
moran_index_HD <- function(data, spatial_att, K, merge = TRUE){
  if (merge == TRUE){
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
  smw <- sum(m_neigh)
  sw <-spmw/smw
  vr <- n / sum(dy^2)
  MI <- vr * sw
  return(MI)
}
#########################################################################################
moran_stat_HD <- function(data, K, spatial_att, obs_moran_I, nsim = 99){
  MI_rand <- numeric(length=nsim+1)
  for (s in 1:nsim){
    spatial_att_shuffle <- spatial_att[sample(length(spatial_att[,2])),2]
    spatial_att[ ,2] <- spatial_att_shuffle
    MI <- moran_index_HD(data, spatial_att, K, merge = FALSE)
    
    MI_rand[s]<- MI
  }
  MI_rand[nsim +1]<- obs_moran_I
  rankres <- rank(MI_rand)

  xrank <- rankres[length(MI_rand)]
  diff <- nsim - xrank
  diff <- ifelse(diff > 0, diff, 0)
  pval <- punif((diff + 1)/(nsim + 1))
  if (!is.finite(pval) || pval < 0 || pval > 1) 
    warning("Out-of-range p-value: reconsider test arguments")
  statistic <- MI_rand[nsim+1]
  names(statistic) <- "statistic"
  parameter <- xrank
  names(parameter) <- "observed rank"
  method <- "Monte-Carlo simulation of Moran I"
  lres <- list(statistic=statistic, parameter=parameter,
               p.value=pval, alternative="greater", method=method)
  print(lres)
  return(lres)
}

#####################################################################
moran_I_scatter_plot <- function(data, Xlab = NULL, Ylab=NULL, Title= NULL){
  if (dim(data)[3] == 1){
    data = data[,,1]
    vect_metod <- c()
    moranI <- c()
    if (is.null(rownames(data))){
      rownames(data) <- seq(dim(data)[1])
    }
    if (is.null(colnames(data))){
      colnames(data) <- seq(dim(data)[2])
    }
    for (i in 1:dim(data)[1]){
      cm <- rep(rownames(data)[i], dim(data)[2])
      vect_metod <- c(vect_metod, cm)
      moranI <- c(moranI, data[i, ])
    }
    att <- rep(colnames(data), dim(data)[1])
    
    df_graph <- data.frame("Methods" = vect_metod, "Attributes" =att, "moranI" = moranI )
    theme_set(theme_bw())
    if (is.null(Title)){
      Title <- "Moran indexes by method"
    }
    if(is.null(Xlab)){
      Xlab <- "variable"  
    }
    if(is.null(Ylab)){
      Ylab <-"Moran.Index"
    }
    
    
    p <- ggplot(df_graph, aes(x=Attributes, y=moranI,  group=Methods, color = Methods)) +  geom_point(size = 4)+
      scale_colour_manual(values=custom.col[1:length(unique(df_graph$Methods))])
    p <- p +  labs(title=Title, caption = "Moran indexes for each variable and for each method. ",
                   y=Ylab, x= Xlab) +theme(plot.title=element_text(size=18, face="bold", color="#17202A", hjust=0.5,lineheight=1.2),  # title
                                           plot.subtitle =element_text(size=13, color="#17202A", hjust=0.5),  # caption
                                           plot.caption =element_text(size=10, color="#17202A", hjust=0.5),  # caption
                                           axis.title.x=element_text(size=12, face="bold"),  # X axis title
                                           axis.title.y=element_text(size=12, face="bold"),  # Y axis title
                                           axis.text.x=element_text(size=12),  # X axis text
                                           axis.text.y=element_text(size=12))  # Y axis text
    print(p) 
    return(list(p,df_graph))
  }
  else{
    vect_metod_k <- c() 
    moranI_k <- c()
    att_k <-c()
    for (k in 1:dim(data)[3]){
      vect_metod <- c()
      moranI <- c()
      att_I <- c()
      if (is.null(rownames(data))){
        rownames(data) <- seq(dim(data)[1])
      }
      if (is.null(colnames(data))){
        colnames(data) <- as.character(seq(dim(data)[2]))
      }
      for (i in 1:dim(data)[1]){
        cm <- rownames(data)[i]
        vect_metod <- c(vect_metod, cm)
        moranI <- c(moranI, data[i, , k])
        att <- colnames(data)
        att_I <- c(att_I, att)
      }
      vect_metod_k <- c(vect_metod_k, vect_metod)
      moranI_k <- c(moranI_k, moranI)
      att_k <- c(att_k, att_I)
    }
    
    df_graph <- data.frame("Methods" = vect_metod_k, "Attributes" =att_k, "moranI" = moranI_k )
    if (is.null(Title)){
      Title <- "Moran indexes by method"
    }
    if(is.null(Xlab)){
      Xlab <- "variable"  
    }
    if(is.null(Ylab)){
      Ylab <-"Moran.Index"
    }
    p <- ggplot(df_graph, aes(x=Attributes, y=moranI,   fill = Methods)) + geom_boxplot(notch=F)+
      scale_fill_manual(values=custom.col[1:length(unique(df_graph$Methods))])
    p <- p +  labs(title=Title, caption = "Moran indexes distribution by k level for each variable and for each method. ",
                   y=Ylab, x= Xlab) +theme(plot.title=element_text(size=18, face="bold", color="#17202A", hjust=0.5,lineheight=1.2),  # title
                                           plot.subtitle =element_text(size=13, color="#17202A", hjust=0.5),  # caption
                                           plot.caption =element_text(size=10, color="#17202A", hjust=0.5),  # caption
                                           axis.title.x=element_text(size=12, face="bold"),  # X axis title
                                           axis.title.y=element_text(size=12, face="bold"),  # Y axis title
                                           axis.text.x=element_text(size=12),  # X axis text
                                           axis.text.y=element_text(size=12))  # Y axis text
    print(p) 
    return(list(p,df_graph))
    
  }
}

########################################################################################
moran_I_scatter_plot_by_k <- function(data, Xlab = NULL, Ylab=NULL, Title= NULL){
  custom.col <- c( '#1E90FF', '#6C3483','#D81B60',  '#B22222', "#D16103",  "#FFD700",  '#2ECC71',"#33691E", '#626567',"#17202A") 
  
  if (dim(data)[3] == 1){
    data = data[,,1]
    vect_metod <- c()
    moranI <- c()
    if (is.null(rownames(data))){
      rownames(data) <- seq(dim(data)[1])
    }
    if (is.null(colnames(data))){
      colnames(data) <- seq(dim(data)[2])
    }
    for (i in 1:dim(data)[1]){
      cm <- rep(rownames(data)[i], dim(data)[2])
      vect_metod <- c(vect_metod, cm)
      moranI <- c(moranI, data[i, ])
    }
    att <- rep(colnames(data), dim(data)[1])
    df_graph <- data.frame("Methods" = vect_metod, "Attributes" =att, "moranI" = moranI )
    
    theme_set(theme_bw())
    if (is.null(Title)){
      Title <- "Moran indexes by method"
    }
    if(is.null(Xlab)){
      Xlab <- "K"  
    }
    if(is.null(Ylab)){
      Ylab <-"Moran.Index"
    }
    
    
    p <- ggplot(df_graph, aes(x=Attributes, y=moranI,  group=Methods, color = Methods)) +  geom_point(size = 4)+
      scale_colour_manual(values=custom.col[1:length(unique(df_graph$Methods))])
    p <- p +  labs(title=Title, caption = "Moran indexes for each variable and for each method. ",
                   y=Ylab, x= Xlab) +theme(plot.title=element_text(size=18, face="bold", color="#17202A", hjust=0.5,lineheight=1.2),  # title
                                           plot.subtitle =element_text(size=13, color="#17202A", hjust=0.5),  # caption
                                           plot.caption =element_text(size=10, color="#17202A", hjust=0.5),  # caption
                                           axis.title.x=element_text(size=12, face="bold"),  # X axis title
                                           axis.title.y=element_text(size=12, face="bold"),  # Y axis title
                                           axis.text.x=element_text(size=12),  # X axis text
                                           axis.text.y=element_text(size=12))  # Y axis text
    print(p) 
    return(p)
  }
  else{
    vect_metod_k <- c() 
    moranI_k <- c()
    att_k <-c()
    k_vect <- c()
    for (k in 1:dim(data)[3]){
      vect_metod <- c()
      moranI <- c()
      att_I <- c()
      if (is.null(rownames(data))){
        rownames(data) <- seq(dim(data)[1])
      }
      if (is.null(colnames(data))){
        colnames(data) <- as.character(seq(dim(data)[2]))
      }
      for (i in 1:dim(data)[1]){
        cm <- rownames(data)[i]
        vect_metod <- c(vect_metod, cm)
        moranI <- c(moranI, data[i, , k])
        att <- colnames(data)
        att_I <- c(att_I, att)
      }
      rep_k <- rep(dimnames(data)[[3]][k], dim(data)[1] * dim(data)[2])
      k_vect <- c(k_vect, rep_k)
      vect_metod_k <- c(vect_metod_k, vect_metod)
      moranI_k <- c(moranI_k, moranI)
      att_k <- c(att_k, att_I)
    }

    df_graph <- data.frame("Methods" = as.character(vect_metod_k), "Attributes" =att_k, "moranI" = moranI_k, "K_level" =  as.numeric(k_vect))
    if (is.null(Title)){
      Title <- "Moran indexes by attribute according "
    }
    if(is.null(Xlab)){
      Xlab <- "variable"  
    }
    if(is.null(Ylab)){
      Ylab <-"Moran.Index"
    }
    list_p <- list()
    for ( i in 1:length(rownames(data))){
    df_graph_c <- df_graph[ df_graph$Methods == rownames(data)[i],]
    p <- ggplot(df_graph_c, aes(x=as.numeric(K_level), y=moranI, color = Attributes)) + geom_point()+
      scale_colour_manual(values=custom.col[1:length(unique(df_graph$Attributes))])
    p <- p +  labs(title= paste(Title,  rownames(data)[i] ) , caption = "Moran indexes distribution by k level for each variable and for each method. ",
                   y=Ylab, x= Xlab) +theme(plot.title=element_text(size=18, face="bold", color="#17202A", hjust=0.5,lineheight=1.2),  # title
                                           plot.subtitle =element_text(size=13, color="#17202A", hjust=0.5),  # caption
                                           plot.caption =element_text(size=10, color="#17202A", hjust=0.5),  # caption
                                           axis.title.x=element_text(size=12, face="bold"),  # X axis title
                                           axis.title.y=element_text(size=12, face="bold"),  # Y axis title
                                           axis.text.x=element_text(size=12),  # X axis text
                                           axis.text.y=element_text(size=12))  # Y axis text
    print(p) 
    list_p <- c(list_p, p)
    }
    return(list(list_p, df_graph))
    
  }
}
##########################################################################################################
moran_ranking <-function(l_coords_data , spatial_data, K_value, N  ,ref=NULL,methods_name = NULL){
  methods_names <- names(l_coords_data)
  colnames(spatial_data)[1] <- "Sample_ID"
  L_coords_data <- list()
  Spa_data <- list() 
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
  MI_array <- matrix(NA, length(l_coords_data), (dim(spatial_data)[2]) -1)
  for (i in 1:length(l_coords_data)){
    c_data <- l_coords_data[[i]]
    spatial_data <- Spa_data[[i]]
    c_sample_id <- spatial_data[ ,1]
    if (dim(c_data)[2] == 3){
      c_data <- as.matrix(c_data[, 2:dim(c_data)[2]])
    }
    if (dim(c_data)[2] == 2){
      k_neigh <- knn2nb(knearneigh(c_data, k=K_value, RANN=FALSE))
      ww <- nb2listw(k_neigh, style='B')
    }
    for (j in 2:dim(spatial_data)[2]){
      if (dim(c_data)[2] == 2){
        MI <- moran(spatial_data[, j], ww, n=length(ww$neighbours), S0=Szero(ww))
        MI_array[i,(j-1)] <- MI$I
      }
      else{
        c_spatial_data <- data.frame("Sample_ID"= as.character(c_sample_id), "att" = spatial_data[, j] )
        MI <- moran_index_HD(data = c_data, spatial_att = c_spatial_data, K = K_value, merge = FALSE)
        MI_array[i,(j-1)] <- MI 
      }
    } 
  }
 
  if(is.null(methods_name)==F & length(methods_name) == dim(MI_array)[1]){
    rownames(MI_array) <- methods_name
  }
  else{
    rownames(MI_array) <- methods_names 
  }
  colnames(MI_array) <- colnames(spatial_data)[2:dim(spatial_data)[2]]
  MI_array <- t(MI_array)  
  MI_LNEN_DF <- MI_array
  first_cor_list <- c()
  for (i in 1:dim(MI_LNEN_DF)[2]){
    sort_df <- MI_LNEN_DF[order(MI_LNEN_DF[,i]),]
    sortN_l <- rownames(sort_df)[1:N]
    first_cor_list[[i]] <- sortN_l
  }
  names(first_cor_list) <- names(List_coords)
    
  v.table<- venn(first_cor_list)
  return(list('MoranIndex' = MI_array, 'Venn.Diagram' =v.table ))
}
####################################################################




