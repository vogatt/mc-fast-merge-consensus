require(clusterCons);

mcmergecons <- function(x,diss=FALSE,algorithms=list('agnes'),alparams=list(),alweights=list(),clmin=2,clmax=10,prop=0.8,reps=50,merge=0, numofcores=1){
  #CHECK INPUTS
  #first check if the input data is an ExpressionSet if so extract the data
  if(class(x)=='ExpressionSet'){x <- expSetProcess(x);};
  #check integrity of data (what about NAs)
  if(data_check(x)!=TRUE){stop('the provided data fails the data integrity check')}
  #check that the algorithms are passed as a list
  if(class(algorithms)!='list'){stop('You must pass the algorithm names in a list')};
  #check that all of the algorithms specified are available
  for(i in 1:length(algorithms)){
    if(existsFunction(algorithms[[i]])!=1){stop('One or more of the specified algorithms does not exist')} else {}
  }
  #alweight checking
  #first check for alweights if NULL then create a list of 1's as long as the number of cms
  if(length(alweights)==0) {
    alweights = as.list(rep('1',length(algorithms)))
  }
  #otherwise check that number of weighting elements equals the number of algorithms specified
  else{
    if(length(alweights)!= length(algorithms)){stop('If you sepcify algorithm weightings their number must equal the number of algorithms being used')}
  }
  #clmin,clmax,reps integers
  if(clmin < 2){stop('cannot have less than two clusters')}
  #prop between 0 and 1
  if(prop <= 0 | prop >1){stop('resampling proportion must be greater than 0 and less than 1')}
  #check merge
  if(merge != 0 & merge !=1){stop('merge can only be 0 or 1')}
  
  #if we get past all this now call the re-sampling
  sample_number = dim(x)[1];
  
  #list to hold all of the consmatrix objects
  cmlist <- list();
  
  registerDoMC(cores=numofcores)
  #cmlistret <- apply_function(clmin:clmax, clst=cmlist);
  result <- foreach(iter=clmin:clmax) %dopar% {
    clst <- list()
    #normalise algorithm weightings
    total = sum(as.numeric(alweights));
    norm = as.numeric(alweights)/total;
    if(merge>0){
      mm <- matrix(0,nrow=dim(x)[1],ncol=dim(x)[1],dimnames=list(row.names(x),row.names(x)))
    }
    
    #for each algorithm
    for(a in 1:length(algorithms)){
      final_conn <- matrix(0,nrow=dim(x)[1],ncol=dim(x)[1],dimnames=list(row.names(x),row.names(x)));
      final_id <- final_conn;
      
      #get the algorithm to call for cm
      algo_clmem <- paste(algorithms[[a]],'_clmem',sep='');
      current_algo <- get(algo_clmem);
      #check if params have been specified
      if(length(alparams)!=0){
        current_params <- alparams[[a]];
        #check if the diss parameter is set and true, in case user forgets diss flag
        if(!is.null(current_params$diss)){
          diss=TRUE;
        }
      }
      else{
        #here someone has specified that x is a distance matrix, but not included it in the params
        if(diss==TRUE){
          current_params = list(diss=TRUE);
        }
        #otherwise x not a distance matrix params stay empty
        else{
          current_params = list();
        }
      }
      for(i in 1:reps){
        #the generic clustering algorithm calling method
        #checking if data is in fact a distance matrix (data.frame)
        #either by flag diss being set or diss being set in parameters
        if (diss==TRUE){
          #we have already put the object through the data_check
          #we need to check that the row and column names are the same in the distance data.frame
          if(unique(names(x)==row.names(x))==FALSE){stop('the row and column names of the distance matrix are not the same or the distance matrix is not square !')}
          #if OK perform the sample
          else{
            samp_row_names <- sample(row.names(x),as.integer(prop*sample_number));
            samp_x <- x[samp_row_names,samp_row_names];
            #convert to object of class 'dissimilarity'
            #samp_x <- as.dist(samp_x);
          }
        }
        else{
          samp_x <- x[sample(row.names(x),as.integer(prop*sample_number)),]; #ORIGINAL
        }
        clmem <- current_algo(samp_x,iter,params=current_params);
        conn <- matrix(0,nrow=dim(x)[1],ncol=dim(x)[1],dimnames=list(c(row.names(x)),c(row.names(x))))
        
        conn[] <- 0
        for(cluster in unique(clmem$cm)) {
          conn[rownames(subset(clmem,cm==cluster)),rownames(subset(clmem,cm==cluster))] <- 1
        }
        
        id = conn;
        id[] <- 0;
        id[row.names(clmem),row.names(clmem)] <- 1;
        final_conn = final_conn+conn;
        final_id = final_id+id;
      }
      #this is the consensus matrix
      consensus = final_conn/final_id;
      #check for NAs (this is safer than pseudo counting, replace NAs with 0 i.e. they are never drawn together and/or never connected)
      ind <- is.na(consensus)
      consensus[ind] <- 0
      #perform the full clustering as reference
      rm <- current_algo(x,iter,params=current_params)
      
      
      #now create the S4 class
      current_consmatrix <- new('consmatrix',cm=as.matrix(consensus),rm=rm,a=algorithms[[a]],k=iter);
      #add the consmatrix object to the list
      clst[paste('e',a,'_',algorithms[[a]],'_k',iter,sep='')] <- current_consmatrix
      #now add to the running matrix for merge weighting by the correct algorithm weighting if merge !=0
      if(merge!=0){
        weighted <- current_consmatrix@cm*norm[a]
        mm <- mm+weighted
      }
    }
    if(merge!=0){
      mm <- new('mergematrix',cm=as.matrix(mm),k=iter,a='merge')
      clst[[paste('merge_k',iter,sep='')]] <- mm
    }
    
    
    return (clst)
  }
  
  return(unlist(result));
}
mcfastmergecons <- function(x, diss=FALSE, algorithms=list('agnes'), alparams=list(), alweights=list(), clmin=2, clmax=10, prop=0.8, reps=50, merge=0, numofcores=2) {

  if(data_check(x)!=TRUE){stop('the provided data fails the data integrity check')}
  if(class(algorithms)!='list'){stop('You must pass the algorithm names in a list')};
  for(i in 1:length(algorithms)){
    if(existsFunction(algorithms[[i]])!=1){stop('One or more of the specified algorithms does not exist')} else {}
  }  
  if(length(alweights)==0) {
    alweights = as.list(rep('1',length(algorithms)))
  } else {
    if(length(alweights)!= length(algorithms)){stop('If you sepcify algorithm weightings their number must equal the number of algorithms being used')}
  }

  if(clmin < 2){stop('cannot have less than two clusters')}
  if(prop <= 0 | prop >1){stop('resampling proportion must be greater than 0 and less than 1')}
  if(merge != 0 & merge !=1){stop('merge can only be 0 or 1')}
  

  total = sum(as.numeric(alweights));
  norm = as.numeric(alweights)/total;
  

  sample_number = dim(x)[1];
  
  clst = list();
  registerDoMC(cores=numofcores)

  for(a in 1:length(algorithms)){
    
    algo_clmem <- paste('.', algorithms[[a]],'_clmem_fc',sep='');
    current_algo <- get(algo_clmem);
    
   
    if(length(alparams)!=0){
      current_params <- alparams[[a]];

      if(!is.null(current_params$diss)){
        diss=TRUE;
      }
    }
    else{

      if(diss==TRUE){
        current_params = list(diss=TRUE);
      }
      
      else{
        current_params = list();
      }
    }
    
    
    clmems <- current_algo(x, clnums=as.list(clmin:clmax), params=current_params, diss=diss, data_size=sample_number, reps=reps, prop=prop, numofcores=numofcores)
    
    consensus_matrices <- foreach(iter=clmin:clmax) %dopar% {
      final_conn <- matrix(0,nrow=length(row.names(x)),ncol=length(row.names(x)),dimnames=list(row.names(x),row.names(x)));
      final_id <- final_conn;
    
      for(m in 1:reps) {
        rn <- row.names(x);
        rl <- length(rn);
        result <- list(connectivity= matrix(0,nrow=rl,ncol=rl,dimnames=list(c(rn),c(rn))), id=matrix(0,nrow=rl,ncol=rl,dimnames=list(c(rn),c(rn))))
        index <- paste('k',iter,sep='')
        
        colnames(clmems[[index]]$data)<-c(rep('cm', dim(clmems[[index]]$data)[2]))
        
        for(clustnum in unique(na.omit(clmems[[index]]$data[,m])) ) {
          result$connectivity[names(subset(clmems[[index]]$data[,m],clmems[[index]]$data[,m]==clustnum)), names(subset(clmems[[index]]$data[,m],clmems[[index]]$data[,m]==clustnum))] <- 1
        }
        
        result$id[names(na.omit(clmems[[index]]$data[,m])),names(na.omit(clmems[[index]]$data[,m]))] <- 1;
        
        
        final_conn <- final_conn + result$connectivity;
        final_id <- final_id + result$id;
        rm(result)
        gc();
      }
      
      consensus = final_conn/final_id;
      ind <- is.na(consensus)
      consensus[ind] <- 0
      k_num <- clmems[[paste('k',iter,sep='')]]$k_number
      algo_clmem <- paste(algorithms[[a]],'_clmem',sep='');
      algorithm_for_reference <- get(algo_clmem);
      reference_clustering <- algorithm_for_reference(x,k_num,params=current_params)
      current_consmatrix <- new('consmatrix',cm=as.matrix(consensus),rm=reference_clustering,a=algorithms[[a]],k=k_num);
      
      if(merge>0){
        mm <- matrix(0,nrow=dim(x)[1],ncol=dim(x)[1],dimnames=list(row.names(x),row.names(x)))
      }
      if(merge!=0){
        mm <- new('mergematrix',cm=mm,k=k_num,a='merge')
        ret <- list(mm=mm,cm=current_consmatrix)
      } else {
        ret <- list(cm=current_consmatrix)
      } 
      rm(consensus);
      rm(reference_clustering);
      return(ret)
    }
    
    for(cmat in consensus_matrices) {
      clst[paste('e',a,'_',algorithms[[a]],'_',paste('k',cmat$cm@k, sep=''),sep='')] <- cmat$cm
      if(merge!=0) {
        clst[[paste('merge_k',cmat$mm@k,sep='')]] <- cmat$mm
      }
    }
    rm(clmems);
    rm(consensus_matrices);
    gc();
  }
  
  if(merge!=0) {
    for(a in 1:length(algorithms)){
      for(clnum in as.list(clmin:clmax)) {
        weighted <- clst[[paste('e',a,'_',algorithms[[a]],'_',paste('k',clnum, sep=''),sep='')]]@cm*norm[a]
        clst[[paste('merge_k',clnum,sep='')]]@cm <- as.matrix(clst[[paste('merge_k',clnum,sep='')]]@cm + weighted)
      }
    }
  }
  
  return(clst);
}