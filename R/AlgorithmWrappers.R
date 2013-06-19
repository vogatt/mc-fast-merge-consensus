require(cluster);
require(parallel);
require(doMC);
require(foreach);
require(plyr);

.agnes_cut <- function(x, params, clnums) {
  params$x <- x
  result_tree <- do.call(agnes, params)
  clmem <- cutree(result_tree, unlist(clnums))
  rownames(clmem) <- row.names(params$x)
  return(t(clmem))
}

.diana_cut <- function(x, params, clnums) {
  params$x <- x
  result_tree <- do.call(diana, params)
  clmem <- cutree(result_tree, unlist(clnums))
  rownames(clmem) <- row.names(params$x)
  return(t(clmem))
}

.hclust_cut <- function(x, params, clnums) {
  params$diss=NULL;
  
  if(class(x)!='dist'){
    dm <- dist(x);
    params$d=dm;
  } else{
    dm <- x;
    params$d=dm;
  }
  result_tree <- do.call(hclust, params)
  clmem <- cutree(result_tree, unlist(clnums))
  rownames(clmem) <- row.names(x)
  return(t(clmem))
}

.hierarchical_algorithm <- function(x, hierarchical_fun, clnums=list(), params=list(), reps, diss, prop, data_size, numofcores) {
  clmems <- list()
  
  res<-foreach(iter=1:reps,.combine=rbind.fill.matrix) %dopar% {
    if(diss) {
      samp_row_names <- sample(row.names(x), as.integer(prop*data_size));
      samp_x <- x[samp_row_names, samp_row_names];
    } else {
      samp_x <- x[sample(row.names(x),as.integer(prop*data_size)),];
    }
    hierarchical_fun(x=samp_x, params=params,clnums=clnums)
  }

  rownames(res) <- rep(clnums, reps)

  foreach(iter=clnums) %do% {
    clmems[[paste('k',iter,sep='')]] <- list(data=t(res[rownames(res)==iter,]) , k_number = iter)
  }
  rm(res)
  gc()
  return(clmems)
}

.agnes_clmem_fc <- function(x, clnums=list(),params=list(), reps, diss, prop, data_size, numofcores){
  return(.hierarchical_algorithm(x=x, hierarchical_fun=.agnes_cut, clnums=clnums, params=params, prop=prop, reps=reps, diss=diss, data_size=data_size, numofcores=numofcores))
  
}

.hclust_clmem_fc <- function(x, clnums=list(),params=list(), reps, diss, prop, data_size, numofcores){
  
  return(.hierarchical_algorithm(x=x, hierarchical_fun=.hclust_cut, clnums=clnums, params=params, prop=prop, reps=reps, diss=diss, data_size=data_size, numofcores=numofcores))

}

.diana_clmem_fc <- function(x, clnums=list(),params=list(), reps, diss, prop, data_size, numofcores){
  
  return(.hierarchical_algorithm(x=x, hierarchical_fun=.diana_cut, clnums=clnums, params=params, prop=prop, reps=reps, diss=diss, data_size=data_size, numofcores=numofcores))
}

.kmeans_clmem_fc <- function(x, diss, data_size, prop, numofcores, clnums=list(), reps, params=list() ) {
  clmems <- list()
  params$diss = NULL
  params$algorithm = 'Mac'
  res <- foreach(iter=1:reps,.combine=rbind.fill.matrix) %dopar% {
    
    if(diss) {
      samp_row_names <- sample(row.names(x), as.integer(prop*data_size));
      samp_x <- x[samp_row_names, samp_row_names];
    } else {
      samp_x <- x[sample(row.names(x),as.integer(prop*data_size)),];
    }

    res2 <- foreach(iter=clnums,.combine=rbind) %do% {
      params$centers = iter
      params$x<-as.matrix(samp_x);
    
      r<-do.call(kmeans, params)$cluster
      e <- r[rownames(x)]
      names(e) <- rownames(x)
      return(e)
    }
    return(res2);
  }
  
  rownames(res) <- rep(clnums, reps)
  
  foreach(iter=clnums) %do% {
    clmems[[paste('k',iter,sep='')]] <- list(data=t(res[rownames(res)==iter,]) , k_number = iter)
  }
  rm(res);
  gc();
  return(clmems)
}



.pam_clmem_fc <- function(x, diss, data_size, prop, numofcores, clnums=list(), reps, params=list() )  {
  clmems <- list()
  params$diss = NULL

  res <- foreach(iter=1:reps,.combine=rbind.fill.matrix) %dopar% {
    
    if(diss) {
      samp_row_names <- sample(row.names(x), as.integer(prop*data_size));
      samp_x <- x[samp_row_names, samp_row_names];
    } else {
      samp_x <- x[sample(row.names(x),as.integer(prop*data_size)),];
    }
    
    res2 <- foreach(iter=clnums,.combine=rbind) %do% {
      if(length(params)==0){
        params <- list(metric='euclidean',k=iter)
      }
      else{
        if(sum(regexpr('^metric$',names(params))==1)<1){
          params$metric='euclidean'
        }
        params$k <- iter
      }
      params$x<-as.matrix(samp_x);
      r<-do.call(pam, params)$clustering
      e <- r[rownames(x)]
      names(e) <- rownames(x)
      return(e)
    }
    return(res2)
    
  }
  rownames(res) <- rep(clnums, reps)
  
  foreach(iter=clnums) %do% {
    singlek <- res[rownames(res)==iter,]
    clmems[[paste('k',iter,sep='')]] <- list(data=t(singlek) , k_number = iter)
  }
  rm(res);
  gc();
  return(clmems)
}


