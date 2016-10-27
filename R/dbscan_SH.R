#' @title dbscan_SH
#' @description DBSCAN with silhouette score optimization
#' @param data $Y of Rtsne output
#' @param eps
#' @param showplot Logical. TRUE to plot silhouette score and eps. Default is FALSE.
#' @param prop_outliers A value between 0 to 1. The proportion of data points allowed to be unclustered. Default is 0.1.
#' @param eps_res Integer. The resolution of simulation. Default is 500.
#' @param eps_range NULL or vectore of two values. The space of search for eps parameter. If a vectore is given, simulation is performed between the given range. Otherwise, te serch space is defined based on input data size.
#' @return A vector of cluster label. 0 is for outlyers.
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