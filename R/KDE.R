
library(ks)
KDE_setup<-function(orig_resid,orig_simflow,warmup){
  orig_resid_nowarm<-orig_resid[-(1:(warmup))]
  orig_simflow_nowarm<-orig_simflow[-(1:(warmup))]
  orig_resid_nowarm_tminus1<-c(NA,orig_resid_nowarm[-length(orig_resid_nowarm)])

  X <- cbind(orig_resid_nowarm_tminus1, orig_simflow_nowarm, orig_resid_nowarm)
  sums<-rowSums(X)
  if(any(is.na(sums))){
    X<-X[-which(is.na(sums)),]
  }
  H <- Hpi.diag(X) # plug-in bandwidth matrix
  fhat <- kde(X, H=H)

  return(list(KDE_data=X,fhat=fhat))
}

KDE_sim<-function(fhat,sim,initial_resid=0){
  grid_y <- seq(min(fhat$x[,3])-10*sd(fhat$x[,3]), max(fhat$x[,3])+10*sd(fhat$x[,3]), length=200)
  # This uses conditional=1:2 to condition on the first and second variable, and gives the conditional density of resid error
  orig_inflow<-new_inflow<-sim
  resid_sim<-rep(NA,length(sim))
  x0<-initial_resid
  for(dd in 1:length(orig_inflow)){
    # x0 <- 20
    x1 <- orig_inflow[dd]
    cond <- predict(fhat, x=cbind(rep(x0, length(grid_y)),rep(x1, length(grid_y)), grid_y), conditional=c(1,2))
    if(!all(cond==0)){
      # normalize density
      cond <- cond / sum(cond) / (grid_y[2] - grid_y[1])
      # plot(y=cond,x=grid_y,type="l")

      # sample by inverse transform
      cdf <- cumsum(cond)
      cdf <- cdf / max(cdf)
      # samples <- approx(cdf, grid_y, runif(5000))$y
      samp<-approx(cdf, grid_y, runif(1))$y

      # adjust samp in case cause negative flow
      samp<-min(orig_inflow[dd],samp)

    } else {
      samp<-0
    }
    resid_sim[dd]<-samp
    new_inflow[dd]<-orig_inflow[dd]-samp
    # set up prev resid error for next time step
    x0<-samp

  }
  return(resid_sim)
}
