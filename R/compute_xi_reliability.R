#' Calculates the xi metric from Renard et al. (2010)
#'
#' @param sim matrix; simulated data in different columns
#' @param obs vector; observed data, length should be the same as ncol(sim)
#'
#' @return numeric; the xi value
#' @export
#'
#' @examples
#' library(CDNP)
#' set.seed(42)
#' obs<-rnorm(1000,sd=1)
#'
#' sd<-1 #  # try with different sd values: 0.5, 1, 2
#' sim<-NULL
#' for(i in 1:100){
#'   sim<-cbind(sim,rnorm(1000,sd=sd))
#' }
#' compute_xi_reliability(sim,obs)
compute_xi_reliability<-function(sim,obs){

  all_obs_pvalues<-rep(NA,nrow(sim))
  all_obs_pvalues_theo<-rep(NA,nrow(sim))
  sim_ecdf_theo<-ecdf(obs)
  for(i in 1:nrow(sim)){
    sim_ecdf<-ecdf(sim[i,])
    # plot(sim_ecdf)
    all_obs_pvalues[i]<-sim_ecdf(obs[i])
    all_obs_pvalues_theo[i]<-sim_ecdf_theo(obs[i])
  }

  quant_theoret<-ppoints(length(all_obs_pvalues)) #seq(0,1,0.01)
  quant_obs<-quantile(all_obs_pvalues,quant_theoret,type=1,na.rm=T)

  indicator_fun<-function(x){
    out<-rep(0,length(x))
    out[which(x==0 | x==1)]<-1
    return(out)
  }

  xi_prime<-sum(indicator_fun(quant_obs))/length(quant_obs)
  xi<-1-xi_prime

  return(xi)
}
