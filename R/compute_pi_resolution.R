#' Calculates the pi abs metric from Renard et al. (2010)
#'
#' @param sim matrix; simulated data in different columns
#'
#' @return numeric; the pi abs value
#' @export
#'
#' @examples
#' library(CDNP)
#' set.seed(42)
#'
#' sd<-1 #  # try with different sd values: 0.5, 1, 2
#' sim<-NULL
#' for(i in 1:100){
#'   sim<-cbind(sim,rnorm(1000,sd=sd))
#' }
#' compute_pi_abs_resolution(sim)
compute_pi_abs_resolution<-function(sim){
  sd_sim<-apply(sim,1,sd)
  pi_abs<-1/length(sd_sim[sd_sim>0])*sum(1/(sd_sim[sd_sim>0]))
  return(pi_abs)
}

#' Calculates the pi rel metric from Renard et al. (2010)
#'
#' @param sim matrix; simulated data in different columns
#'
#' @return numeric; the pi rel value
#' @export
#'
#' @examples
#' library(CDNP)
#' set.seed(42)
#'
#' sd<-1 #  # try with different sd values: 0.5, 1, 2
#' mean<-1 # try with different mean values: 0.5, 1, 2
#' sim<-NULL
#' for(i in 1:100){
#'   sim<-cbind(sim,rnorm(1000,sd=sd,mean=mean))
#' }
#' compute_pi_rel_resolution(sim)
compute_pi_rel_resolution<-function(sim){
  sd_sim<-apply(sim,1,sd)
  mean_sim<-apply(sim,1,mean)
  # mean_sim<-pmax(mean_sim,0)
  mean_sim<-abs(mean_sim)
  pi_rel<-1/length(sd_sim[sd_sim>0])*sum(mean_sim[sd_sim>0]/sd_sim[sd_sim>0])
  return(pi_rel)
}
