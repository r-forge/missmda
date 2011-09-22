imputePCA <- function (X, ncp = 2, scale=TRUE, method="Regularized",threshold = 1e-6,seed = NULL,nb.init=1,maxiter=1000,row.w=NULL,...){

impute <- function (X, ncp = 4, scale=TRUE, method=NULL,threshold = 1e-6,seed = NULL,init=1,maxiter=1000,row.w=NULL,...){

    moy.p <- function(V, poids) {
        res <- sum(V * poids,na.rm=TRUE)/sum(poids[!is.na(V)])
    }
    ec <- function(V, poids) {
        res <- sqrt(sum(V^2 * poids,na.rm=TRUE)/sum(poids[!is.na(V)]))
    }

   if (!is.null(seed)) set.seed(seed)
   X <- as.matrix(X)
   ncp <- min(ncp,ncol(X),nrow(X)-1)
   missing <- which(is.na(X))
   mean.p <- apply(X, 2, moy.p,row.w)
   Xhat <- sweep(X, 2,mean.p,FUN="-")
   et <- apply(Xhat, 2, ec,row.w)
   if (scale) Xhat <- sweep(Xhat, 2,et,FUN="/")
   if (any(is.na(X))) Xhat[missing] <- 0
   if (init>1) Xhat[missing] <- rnorm(length(missing)) ## random initialization
   recon <- Xhat

   nb.iter <- 1
   old <- Inf
   while (nb.iter > 0) {
       Xhat[missing] <- recon[missing]
       if (scale) Xhat=sweep(Xhat,2,et, FUN="*")
       Xhat <- sweep(Xhat,2,mean.p, FUN="+")
       mean.p <- apply(Xhat, 2, moy.p,row.w)
       Xhat <- sweep(Xhat,2,mean.p, FUN="-")
       et <- apply(Xhat, 2, ec,row.w)
       if (scale) Xhat <- sweep(Xhat,2,et, FUN="/")

       svd.res <- svd.triplet(Xhat,row.w=row.w)
       sigma2 <- mean(svd.res$vs[-(1:ncp)]^2)
       if (method=="em") sigma2 <-0

       if (ncp==1) {
         lambda.shrinked=(svd.res$vs[1]^2-sigma2)/svd.res$vs[1]
##         recon=(svd.res$U[,1]%*%diag(lambda.shrinked[1],1)%*%(t(svd.res$V[,1])))
         recon=tcrossprod(sweep(svd.res$U[,1],1,row.w,FUN="*")*lambda.shrinked[1],svd.res$V[,1])
       } else {
         lambda.shrinked=(svd.res$vs[1:ncp]^2-sigma2)/svd.res$vs[1:ncp]
##         recon=(svd.res$U[,1:ncp]%*%lambda.shrinked%*% (t(svd.res$V[,1:ncp])))
         recon = tcrossprod(sweep(sweep(svd.res$U[,1:ncp],1,row.w,FUN="*"),2,lambda.shrinked,FUN="*"),svd.res$V[,1:ncp])
       }
       recon <- sweep(recon,1,row.w,FUN="/")
	   diff <- Xhat-recon
	   diff[missing] <- 0
       objective <- mean(sweep(diff^2,1,row.w,FUN="*"))
#       objective <- mean((Xhat[-missing]-recon[-missing])^2)
       criterion <- abs(1 - objective/old)
       old <- objective
       nb.iter <- nb.iter + 1
       if (!is.nan(criterion)) {
         if ((criterion < threshold) && (nb.iter > 5))  nb.iter <- 0
       }
       if (nb.iter>maxiter) {
         nb.iter <- 0
         warning(paste("Stopped after ",maxiter," iterations"))
       }
   }
   if (scale) Xhat <- sweep(Xhat,2,et, FUN="*")
   Xhat <- sweep(Xhat,2,mean.p, FUN="+")
   completeObs <- X
   completeObs[missing] <- Xhat[missing]
   if (scale) recon <- sweep(recon,2,et, FUN="*")
   recon <- sweep(recon,2,mean.p, FUN="+")

   result <- list()
   result$completeObs <- completeObs
   result$objective <- objective
   result$recon <- recon
   return(result) 
}

#### Main program
 obj=Inf
 method <- tolower(method)
 if (ncp>=min(nrow(X)-2,ncol(X)-1)) stop("ncp is too large")
 if (is.null(row.w)) row.w = rep(1,nrow(X))/nrow(X)
 for (i in 1:nb.init){
  if (!any(is.na(X))) return(X)
  res.impute=impute(X, ncp=ncp, scale=scale, method=method, threshold = threshold,seed=seed,init=i,maxiter=maxiter,row.w=row.w)
  if (mean((res.impute$recon[!is.na(X)]-X[!is.na(X)])^2) < obj){
    res <- res.impute
    obj <- mean((res.impute$recon[!is.na(X)]-X[!is.na(X)])^2)
  }
 }
return(res)
}

