#' MEREDITH combines a list of input data (data.frame or matrix) and perform tSNE
#' @param data list of matrix or data frame (Note that even when input is one matrix/dataframe, use list())
#' @param scale option for PCA (see ?prcomp). Default is FALSE
#' @param transfer logical. TRUE if the input matrix/data frame have to be transfered. Default is FALSE.
#' @param takepc the number of PC to keep. Default is 50.
#' @param explainedvar FALSE or percentage of explained variance to keep. Default is FALSE. If this argument is value between 0 and 100, takepc will be ignored.
#' @param nTSNE The number of run for tSNE. Default is 100.
#' @param dims,initial_dims,preplexity,theta,max_iter option of Rtsne (see ?Rtse).
#' @return list same as Rtsne output (see ?Rtsne for details)
#' @export
meredith <- function(data, scale=FALSE, transfer=FALSE, takepc=50, explainedVar=FALSE, nTSNE=100, dims=2, initial_dims=50, perplexity = 30, theta=0.5, max_iter=1000){
  #require(Rtsne)
  
  if(transfer){
    for(i in 1:length(data)){
      data[[i]]<-t(data[[i]])
    }
  }
  
  mincol = ncol(data[[1]])
  if(length(data)>1){
    for(i in 2:length(data)){
      if(nrow(data[[i-1]])!=nrow(data[[i]])){
        stop("Input sample size have to match\n")
      }
      if(mincol > ncol(data[[i]])){mincol = ncol(data[[i]])}
    }
  }
  
  if(explainedVar){
    message("Use explainedVar rather than takepc\nPlease set explainedVar as FALSE to use takepc")
    if(explainedVar <= 0 | explainedVar > 100){
      stop("explainedVar has to be between 0 and 100%")
    }
  }
  
  if(explainedVar==FALSE & takepc > mincol){
    stop("takepc has to be smaller than the number of column of input data")
  }
  
  #PCA and combind data
  if(length(data)>1){
    for(i in 1:length(data)){
      tmpPC = prcomp(data[[i]], scale.=scale)
      if(explainedVar){
        percentVar= cumsum((tmpPC$sdev)^2)/sum((tmpPC$sdev)^2)
        n = length(which(percentVar<=explainedVar))
        pcScore = tmpPC$x[,1:n]
      }else{
        pcScore = tmpPC$x[,1:takepc]
      }
      
      tmpvar = apply(pcScore, 2, var)
      pcScore = pcScore/sqrt(sum(tmpvar))
      if(i==1){
        PCAcombined = pcScore
      }else{
        PCAcombined = cbind(PCAcombined, pcScore)
      }
    }
    cat("PCA is done!!\n",length(data)," data sets were combined\n")
    
    cat("tSNE\n")
    pb = txtProgressBar(min=0, max=100, initial=0, style=3)
    for(i in 1:nTSNE){
      tmpTSNE = Rtsne(PCAcombined, pca=FALSE, dims=dims, initial_dims=initial_dims, perplexity=perplexity, theta=theta, max_iter=max_iter)
      if(i==1){
        tSNE.out = tmpTSNE
        cost = tmpTSNE$itercost[length(tmpTSNE$itercost)]
      }else{
        if(tmpTSNE$itercost[length(tmpTSNE$itercost)] < cost){
          tSNE.out = tmpTSNE
          cost = tmpTSNE$itercost[length(tmpTSNE$itercost)]
        }
      }
      setTxtProgressBar(pb, i*100/nTSNE)
    }
    cat("\nCompleted!!\n")
  }else{
    data = data[[1]]
    tmpPC = prcomp(data, scale. = scale)
    if(explainedVar){
      percentVar= cumsum((tmpPC$sdev)^2)/sum((tmpPC$sdev)^2)
      n = length(which(percentVar<=explainedVar))
      data.pc = tmpPC$x[,1:n]
    }else{
      data.pc = tmpPC$x[,1:takepc]
    }
    #tSNE
    cat("tSNE\n")
    pb = txtProgressBar(min=0, max=100, initia=0, style=3)
    for(i in 1:nTSNE){
      tmpTSNE = Rtsne(data.pc, pca=FALSE, dims=dims, initial_dims=initial_dims, perplexity=perplexity, theta=theta, max_iter=max_iter)
      if(i==1){
        tSNE.out = tmpTSNE
        cost = tmpTSNE$itercost[length(tmpTSNE$itercost)]
      }else{
        if(tmpTSNE$itercost[length(tmpTSNE$itercost)] < cost){
          tSNE.out = tmpTSNE
          cost = tmpTSNE$itercost[length(tmpTSNE$itercost)]
        }
      }
      setTxtProgressBar(pb, i*100/nTSNE)
    }
    cat("\nCompleted!!\n")
  }
  return(tSNE.out)
}

#` DBSCAN with silhouette score optimization
#' @param 
#' @return
#' @export
dbscan_SH<-function(data,eps=NULL,showplot=F,prop_outliers=.1,eps_res=500,eps_range=NULL){
  #require(fpc)
  #require(cluster)
  data<-data.frame(data)
  if(is.null(eps)){
    cat("Optimising eps: ")
    if(is.null(eps_range)){
      eps_scale<-mean(apply(data,2,sd)) # makes the search scale independent
      epsvec<-seq(0,4,length.out=eps_res)*eps_scale # space to search for eps parameter
    } else epsvec<-seq(eps_range[1],eps_range[2],length.out=eps_res)
    silvec<-numeric(length(epsvec))
    for(i in 1:length(epsvec)){
      eps<-epsvec[i]
      DBcl<-dbscan(data,eps) # quite fast
      cl<-DBcl$cluster
      cat(".")
      if(all(cl==1)) break else if(max(cl)==1) silvec[i]<-0 else
        if(all(cl==0)) silvec[i]<-0 else
          if(mean(cl==0)>prop_outliers) silvec[i]<-0 else{
            S<-silhouette(x=cl[cl!=0],dist=dist(data[cl!=0,])) # exclude the 0's
            silvec[i]<-summary(S)$avg.width
          }
    }
    cat("\n")
    if(showplot){
      if(ncol(data)==2)par(mfrow=c(1,2))
      end<-length(silvec)-which(cumsum(rev(silvec))>0)[1]+10
      plot(epsvec[1:end],silvec[1:end],xlab="eps value",ylab="silhouette score")
    }
    eps<-epsvec[which.max(silvec)]
    if(showplot){
      abline(h=max(silvec),lty=2)
      abline(v=eps,lty=2)
    }
  }
  DBcl<-dbscan(data,eps)$cluster
  if(showplot){
    plot(data,col=c("lightgrey",sample(rainbow(max(DBcl),v=.8)))[DBcl+1],pch=16)
    par(mfrow=c(1,1))
  }
  cat("Used eps: ",eps,"\n")
  return(DBcl)
}