# https://journals.ametsoc.org/view/journals/hydr/15/4/jhm-d-14-0008_1.xml
# equation 9 and 15
# s^2
unbiased_variance_estimator<-function(x){
  # R: num_ens_members
  num_ens_members<-length(x)
  s_squared<-1/(num_ens_members-1)*sum((mean(x)-x)^2)
  return(s_squared)
}
ensemble_spread<-function(xx){
  # each ensemble member is in different column
  # each time step is in different row
  s_squared<-apply(xx,1,unbiased_variance_estimator)
  es<-sqrt(((ncol(xx)+1)/ncol(xx))*(1/nrow(xx))*sum(s_squared)) # equation 15
  return(es)
}

mean_ensemble_rmse<-function(sim_ens,obs){
  mean_ens<-apply(sim_ens,1,mean)
  mean_ens_rmse<-hydroGOF::rmse(sim=mean_ens,obs=obs)
  return(mean_ens_rmse)
}

spread_skill_ratio<-function(sim,obs){
  ens_spread<-ensemble_spread(sim)
  mean_ens_rmse<-mean_ensemble_rmse(sim,obs)
  return(ens_spread/mean_ens_rmse)
}



