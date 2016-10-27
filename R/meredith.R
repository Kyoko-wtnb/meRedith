#' @title MEREDITH (multi omic data integration approach)
#' @description MEREDITH combines a list of input data (data.frame or matrix) and perform tSNE
#' @param data List of matrix or data frame. (Note that even when input is one matrix/dataframe, use list())
#' @param Logical. The scale option for PCA (see ?prcomp). Default is FALSE
#' @param transfer Logical. TRUE if the input matrix/data frame have to be transfered. Default is FALSE.
#' @param takepc Integer. The number of PC to keep. Default is 50.
#' @param explainedvar FALSE or a value. The percentage of explained variance to keep. Default is FALSE. If this argument is value between 0 and 100, takepc will be ignored.
#' @param nTSNE Integer. The number of run for tSNE. Default is 100.
#' @param dims,initial_dims,preplexity,theta,max_iter option of Rtsne (see ?Rtsne).
#' @return list of data same as the Rtsne output (see ?Rtsne for details)
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