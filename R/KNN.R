
library(FNN)

setup_KNN_data<-function(orig_resid,orig_simflow,warmup){
  orig_resid_nowarm<-orig_resid[-(1:(warmup))]
  orig_simflow_nowarm<-orig_simflow[-(1:(warmup))]
  orig_resid_nowarm_tminus1<-c(NA,orig_resid_nowarm[-length(orig_resid_nowarm)])

  X <- cbind(orig_resid_nowarm_tminus1, orig_simflow_nowarm, orig_resid_nowarm)
  sums<-rowSums(X)
  if(any(is.na(sums))){
    X<-X[-which(is.na(sums)),]
  }


  # min-max scaling
  col.names<-names(X)
  # only for prev resid and current simflow
  prevresid_norm_dat<-normalise(X[,1])
  simflow_norm_dat<-normalise(X[,2])
  X<-cbind(prevresid_norm_dat$normalised_data,
           simflow_norm_dat$normalised_data,
           X[,3])
  colnames(X)<-col.names

  return(list(KNN_data=X,prevresid_norm_dat=prevresid_norm_dat,simflow_norm_dat=simflow_norm_dat))
}


kNN_sim<-function(data,sim,k=10,prevresid_norm_dat,simflow_norm_dat,initial_resid=0,zero_handling_option=2){
  # set.seed(45)

  resid_sim<-rep(NA,length(sim))
  prev_error<-initial_resid
  for(dd in 1:length(sim)){
    if(dd>1){
      prev_error<-resid_sim[dd-1]
    }
    cur_simflow<-sim[dd]

    new_data<-matrix(c(normalise(prev_error,min_vals = prevresid_norm_dat$min_vals,max_vals = prevresid_norm_dat$max_vals)$normalised_data,
                       normalise(cur_simflow,min_vals = simflow_norm_dat$min_vals,max_vals = simflow_norm_dat$max_vals)$normalised_data),ncol=2)
    knx<-get.knnx(data=data[,1:2], query=new_data, k=k)
    # knx$nn.dist

    err_to_sample<-data[c(knx$nn.index),3]
    if(zero_handling_option==1){
      if(any(err_to_sample>cur_simflow)){
        err_to_sample<-err_to_sample[-which(err_to_sample>cur_simflow)]
        err_to_sample<-c(err_to_sample,cur_simflow)
      }

      if(length(err_to_sample)==0){
        resid_sim[dd]<-0
      } else {
        resid_sim[dd]<-sample(err_to_sample,size=1)
      }
    } else if(zero_handling_option==2){
      resid_sim[dd]<-sample(err_to_sample,size=1)
      if(resid_sim[dd]>cur_simflow) resid_sim[dd]<-cur_simflow
    } else {
      stop("zero_handling_option not available")
    }


    # samp_indices<-apply(knx$nn.index,1,function(x) sample(x,size=1))
    # resid_sim[dd]<-data[samp_indices,3]
  }

  return(resid_sim)
}


get_optimal_kNN<-function(resid,simflow,warmup=1095,nrep=10,itermax=20,
                          lower=5,upper=100){

  if(length(resid)!=length(simflow)) stop("length of residual is different to length of streamflow")
  if(warmup>length(simflow)) stop("warmup is longer than the length of time series")
  if(!check_data_ok(resid)) stop("not enough data to split the time series")
  resid_nowarm<-resid[-(1:(warmup))]
  simflow_nowarm<-simflow[-(1:(warmup))]

  split_index<-ceiling(length(resid_nowarm)/2)

  resid_nowarm_spl1<-resid_nowarm[1:split_index]
  resid_nowarm_spl2<-resid_nowarm[(split_index+1):length(resid_nowarm)]

  simflow_nowarm_spl1<-simflow_nowarm[1:split_index]
  simflow_nowarm_spl2<-simflow_nowarm[(split_index+1):length(simflow_nowarm)]

  simflow_withwarm_spl1<-simflow[1:(split_index+warmup)]
  head(simflow_withwarm_spl1[-(1:warmup)])
  head(simflow_nowarm_spl1)
  tail(simflow_withwarm_spl1[-(1:warmup)])
  tail(simflow_nowarm_spl1)

  simflow_withwarm_spl2<-simflow_nowarm[((split_index+1)-warmup):length(simflow_nowarm)]
  length(((split_index+1)-warmup):(split_index))
  head(simflow_withwarm_spl2[-(1:warmup)])
  head(simflow_nowarm_spl2)
  tail(simflow_withwarm_spl2[-(1:warmup)])
  tail(simflow_nowarm_spl2)

  compute_alpha_for_optim<-function(par,nrep=nrep){

    kk<-ceiling(par[1])

    resid_nowarm_spl1
    simflow_nowarm_spl1
    ts_data_resid_tminus1<-c(NA,resid_nowarm_spl1[-length(resid_nowarm_spl1)])
    all_dat<-cbind(ts_data_resid_tminus1,simflow_nowarm_spl1,resid_nowarm_spl1)
    all_dat<-as.data.frame(all_dat)
    to_remove<-which(is.na(all_dat[,1]) | is.na(all_dat[,2]) | is.na(all_dat[,3]))
    if(length(to_remove)>0) all_dat<-all_dat[-to_remove,]

    # scale
    col.names<-names(all_dat)
    # only for prev resid and current simflow
    prevresid_norm_dat<-normalise(all_dat[,1])
    simflow_norm_dat<-normalise(all_dat[,2])
    all_dat<-cbind(prevresid_norm_dat$normalised_data,
                   simflow_norm_dat$normalised_data,
                   all_dat[,3])
    colnames(all_dat)<-col.names
    all_dat<-as.data.frame(all_dat)

    # kNN_sim(X,X[,1:2],kk)

    # cluster_out_spl<-get_CDNP_clusters(ncluster=ncluster,nbin=nbin,ts_data_resid=resid_nowarm_spl1,
    #                                    ts_data_simflow=simflow_nowarm_spl1,use_quantile_spacing=T,seed=45)
    # all_posterior_lookup_spl<-get_CDNP_posterior_lookup(cluster_out_spl)


    all_knn_sim_out<-matrix(NA,nrow=length(simflow_withwarm_spl2)-warmup,ncol=nrep)
    for(rrr in 1:nrep){
      knn_sim_out<-kNN_sim(all_dat,simflow_withwarm_spl2,kk,prevresid_norm_dat,simflow_norm_dat)
      # CDNP_sim_out<-CDNP_sim(all_posterior_lookup_spl,simflow_withwarm_spl2,seed=rrr)
      # lines(CDNP_sim_out,col=2,lty=2)
      knn_sim_out_nowarm<-knn_sim_out[-(1:warmup)]
      all_knn_sim_out[,rrr]<-knn_sim_out_nowarm
    }

    # compute_alpha_reliability(all_CDNP_sim_out,orig_resid_nowarm_spl1)
    alpha<-compute_alpha_reliability(all_knn_sim_out,resid_nowarm_spl2)
    cat("Sampled pars:",par,"Alpha:",alpha,"\n")
    return(-alpha)
  }

  # opt<-optim(par<-as.integer(c(10,10)),compute_alpha_for_optim,method="SANN",control=list(fnscale=-1,trace=1))
  # opt<-optim(par<-as.integer(c(10,10)),compute_alpha_for_optim,control=list(fnscale=-1,trace=1))

  # library(DEoptim)
  # Mapping function to ensure integers
  # map_fun <- function(x) {
  #   return(ceiling(x))
  # }
  # deopt<-DEoptim::DEoptim(compute_alpha_for_optim,lower=lower,upper=upper,
  #                         control=list(trace=1,itermax=itermax),nrep=nrep) #fnMap = map_fun,
  # library(GenSA)
  # opt<-GenSA(par = 10, fn = compute_alpha_for_optim, lower = lower, upper = upper,nrep=nrep,
  #            control=list(verbose=T))
  all_alpha<-c()
  trial_k<-round(seq(lower,upper,length.out=itermax))
  # trial_k<-c(150,200,300,400,500)
  for(kkk in trial_k){
    alph<-compute_alpha_for_optim(kkk,nrep=nrep)
    all_alpha<-c(all_alpha,alph)
  }
  # all_pars<-as.matrix(expand.grid(10:30,10:30))
  # all_alpha<-rep(NA,nrow(all_pars))
  # for(pp in 1:nrow(all_pars)){
  #   cat(pp,"/",nrow(all_pars),"\n")
  #   pars<-all_pars[pp,]
  #   all_alpha[pp]<-compute_alpha_for_optim(pars)
  # }
  best_k<-trial_k[which.min(all_alpha)]
  best_alpha<-all_alpha[which.min(all_alpha)]
  return(list(best_k=best_k,best_alpha=best_alpha))
}
