#' Normalise values to fit between 0 and 1
#'
#' @param vals numeric values to be normalised
#' @param min_vals numeric value representing the minimum value of the population
#' @param max_vals numeric value representing the maximum value of the population
#'
#' @return list containing normalised values, min_vals and max_vals
#' @export
#'
#' @examples
#' library(CDNP)
#' norm_dat<-normalise(eg_data$resid)
#' norm_dat
normalise<-function(vals,min_vals=NA,max_vals=NA){
  # Normalise values to fit between 0 and 1

  if(is.na(min_vals)){
    min_vals<-min(vals)
    max_vals<-max(vals)
  }

  if(min_vals==max_vals) max_vals<-max_vals+1e-15

  out<-(vals - min_vals) / (max_vals - min_vals)

  return(list(normalised_data=out,min_vals=min_vals,max_vals=max_vals))
}



#' Get the k-means clusters based on error and streamflow data
#'
#' @param nclusters number of clusters to use
#' @param nbin number of bins for the error residual. Determines the intervals if not NA (default). Equal to the number of intervals minus 1.
#' @param ts_data_resid vector of residual errors. Should be the same length as ts_data_simflow.
#' @param ts_data_simflow vector of streamflow. Should be the same length as ts_data_resid.
#' @param use_quantile_spacing logical; if true uses quantiles to determine the error intervals. If false, uses equal size spacing of intervals (default=TRUE).
#' @param normalise_data logical; if true uses min-max normalisation for residual (at t-1) and streamflow (default=TRUE).
#' @param seed numeric; optional seed that can be set if not NA (default).
#' @param bootstrap logical; if true performs a bootstrap of the input data (default=FALSE).
#'
#' @return returns a list of objects which is used by get_CDNP_posterior_lookup:\cr\cr
#' \verb{    }resid_intervals: vector of error residual intervals\cr\cr
#' \verb{    }kmeans_model: k-mean model object\cr\cr
#' \verb{    }orig_resid: vector; error residual (same as input)\cr\cr
#' \verb{    }orig_simflow: vector; streamflow (same as input)\cr\cr
#' \verb{    }normalise_data: logical; if true uses min-max normalisation for residual (at t-1) and streamflow (same as input)\cr\cr
#' \verb{    }prevresid_norm_dat: list; results from normalisation of residual (at t-1)\cr\cr
#' \verb{    }simflow_norm_dat: list; results from normalisation of streamflow\cr\cr
#' \verb{    }bootstrap_indices: vector; resample indices from bootstrapping
#' @export
#'
#' @examples
#' library(CDNP)
#' CDNP_clusters_out<-get_CDNP_clusters(nclusters=10,
#'                                      nbin=8,
#'                                      ts_data_resid=eg_data$resid,
#'                                      ts_data_simflow=eg_data$sim)
#' CDNP_clusters_out
get_CDNP_clusters<-function(nclusters,nbin=NA,ts_data_resid,ts_data_simflow,ts_data_USresid=NA,use_quantile_spacing=T,normalise_data=T,seed=NA,bootstrap=F){

  # nclusters=10
  # ts_data_resid=eg_data$resid
  # ts_data_simflow=eg_data$sim
  # ts_data_USresid=eg_data$resid+rnorm(length(eg_data$resid))
  # use_quantile_spacing=T
  # normalise_data=T
  # seed=NA
  # bootstrap=F

  if(!is.na(seed)) set.seed(seed)

  ts_data_resid_tminus1<-c(NA,ts_data_resid[-length(ts_data_resid)])
  if(length(ts_data_USresid)==length(ts_data_simflow)){
    all_dat<-cbind(ts_data_resid_tminus1,ts_data_simflow,ts_data_USresid,ts_data_resid)
    all_dat<-as.data.frame(all_dat)
    to_remove<-which(is.na(all_dat[,1]) | is.na(all_dat[,2]) | is.na(all_dat[,3]) | is.na(all_dat[,4]))
  } else {
    all_dat<-cbind(ts_data_resid_tminus1,ts_data_simflow,ts_data_resid)
    all_dat<-as.data.frame(all_dat)
    to_remove<-which(is.na(all_dat[,1]) | is.na(all_dat[,2]) | is.na(all_dat[,3]))
  }

  obsflow_dat<-ts_data_simflow-ts_data_resid

  if(length(to_remove)>0){
    all_dat<-all_dat[-to_remove,]
    obsflow_dat<-obsflow_dat[-to_remove]
  }
  if(any(obsflow_dat< -1e-10)) stop("Observed flow values should not be less than zero")


  # resample
  if(bootstrap){
    bootstrap_indices<-sample.int(n=nrow(all_dat),size=nrow(all_dat),replace = T)
    all_dat<-all_dat[bootstrap_indices,]
    obsflow_dat<-obsflow_dat[bootstrap_indices]
  } else {
    bootstrap_indices<-NA
  }

  if(normalise_data){
    if(length(ts_data_USresid)==length(ts_data_simflow)){
      col.names<-names(all_dat)
      # only for prev resid, current simflow, US resid
      prevresid_norm_dat<-normalise(all_dat[,1])
      simflow_norm_dat<-normalise(all_dat[,2])
      USresid_norm_dat<-normalise(all_dat[,3])
      all_dat<-cbind(prevresid_norm_dat$normalised_data,
                     simflow_norm_dat$normalised_data,
                     USresid_norm_dat$normalised_data,
                     all_dat[,4])
      colnames(all_dat)<-col.names
      all_dat<-as.data.frame(all_dat)
    } else {
      USresid_norm_dat<-NA
      col.names<-names(all_dat)
      # only for prev resid and current simflow
      prevresid_norm_dat<-normalise(all_dat[,1])
      simflow_norm_dat<-normalise(all_dat[,2])
      all_dat<-cbind(prevresid_norm_dat$normalised_data,
                     simflow_norm_dat$normalised_data,
                     all_dat[,3])
      colnames(all_dat)<-col.names
      all_dat<-as.data.frame(all_dat)
    }


  } else {
    prevresid_norm_dat<-NA
    simflow_norm_dat<-NA
    USresid_norm_dat<-NA
  }

  if(length(ts_data_USresid)==length(ts_data_simflow)){
    num_unique_points<-length(unique(paste0(all_dat[,1],"_",all_dat[,2],"_",all_dat[,3])))
    if(num_unique_points<nclusters){
      cat("Number of unique data points is less than specified nclusters. Will continue with",num_unique_points,"clusters\n")
    }
  } else {
    num_unique_points<-length(unique(paste0(all_dat[,1],"_",all_dat[,2])))
    if(num_unique_points<nclusters){
      cat("Number of unique data points is less than specified nclusters. Will continue with",num_unique_points,"clusters\n")
    }
  }

  if(length(ts_data_USresid)==length(ts_data_simflow)){
    km<-kmeans(all_dat[,1:3],min(num_unique_points,nclusters),nstart = 10,iter.max = 100)
  } else {
    km<-kmeans(all_dat[,1:2],min(num_unique_points,nclusters),nstart = 10,iter.max = 100)
  }

  # km$cluster
  # km$centers
  # plot(y=all_dat[,2],x=all_dat[,1],log="")
  # points(km$centers, col = 2, pch = 8, cex = 2)

  if(!is.na(nbin)){
    # residual intervals (small buffer for this at top end)
    if(!use_quantile_spacing){
      ts_intervals_curresid<-seq(min(all_dat$ts_data_resid),max(all_dat$ts_data_resid)+1,length.out=nbin+1)
    } else {
      ts_intervals_curresid<-quantile(all_dat$ts_data_resid,probs = seq(0,1,length.out=nbin+1))
      if(length(unique(ts_intervals_curresid))<length(ts_intervals_curresid)){
        cat("Number of error bins changed from",nbin,"to",length(unique(ts_intervals_curresid))-1,"to avoid repeat intervals\n")
        ts_intervals_curresid<-unique(ts_intervals_curresid)
      }
      ts_intervals_curresid[length(ts_intervals_curresid)]<-ts_intervals_curresid[length(ts_intervals_curresid)]+1
    }

    # perform check to see if there is data inside each interval
    num_data_points_within_each_bin<-numeric(length(ts_intervals_curresid)-1)
    for(interv in 1:(length(ts_intervals_curresid)-1)){
      indices_within_bin<-which(all_dat$ts_data_resid>=ts_intervals_curresid[interv] & all_dat$ts_data_resid<ts_intervals_curresid[interv+1])
      num_data_points_within_each_bin[interv]<-length(indices_within_bin)
    }
    if(any(num_data_points_within_each_bin==0)){
      to_remove<-which(num_data_points_within_each_bin==0)
      ts_intervals_curresid<-ts_intervals_curresid[-to_remove]
      num_data_points_within_each_bin<-num_data_points_within_each_bin[-to_remove]
      cat("Number of error bins changed to",length(ts_intervals_curresid)-1,"to avoid empty bins\n")
    }
  } else {
    ts_intervals_curresid<-NA
  }



  return(list(resid_intervals=ts_intervals_curresid,
              kmeans_model=km,
              orig_resid=ts_data_resid,
              orig_simflow=ts_data_simflow,
              orig_USresid=ts_data_USresid,
              normalise_data=normalise_data,
              prevresid_norm_dat=prevresid_norm_dat,
              simflow_norm_dat=simflow_norm_dat,
              USresid_norm_dat=USresid_norm_dat,
              bootstrap_indices=bootstrap_indices,
              all_dat=all_dat,
              obsflow_dat=obsflow_dat))
}


# get_CDNP_clusters_output<-get_CDNP_clusters(nclusters=10,
#                                      nbin=NA,
#                                      ts_data_resid=eg_data$resid,
#                                      ts_data_simflow=eg_data$sim)



CDNP_clusters_sim<-function(get_CDNP_clusters_output,simflow,USresid=NA,initial_resid=0,seed=NA,recompute_all_dat=F,prevent_neg_flow_after_sample=T,
                            simflow_trans=NULL,resid_trans=NULL,simflow_invtrans=NULL,resid_invtrans=NULL,sample_option=1,...){

  # below is just to ensure old scripts work with inaccurate argument name - should eventually be removed
  if(exists("prevent_zeroflow_after_sample")){
    prevent_neg_flow_after_sample<-prevent_zeroflow_after_sample
  }

  orig_simflow<-get_CDNP_clusters_output$orig_simflow
  orig_resid<-get_CDNP_clusters_output$orig_resid
  orig_USresid<-get_CDNP_clusters_output$orig_USresid
  resid_intervals<-get_CDNP_clusters_output$resid_intervals
  kmeans_model<-get_CDNP_clusters_output$kmeans_model
  normalise_data<-get_CDNP_clusters_output$normalise_data
  prevresid_norm_dat<-get_CDNP_clusters_output$prevresid_norm_dat
  simflow_norm_dat<-get_CDNP_clusters_output$simflow_norm_dat
  USresid_norm_dat<-get_CDNP_clusters_output$USresid_norm_dat
  bootstrap_indices<-get_CDNP_clusters_output$bootstrap_indices
  obsflow_dat<-get_CDNP_clusters_output$obsflow_dat

  if(recompute_all_dat){
    ts_data_resid_tminus1<-c(NA,orig_resid[-length(orig_resid)])
    if(!all(is.na(orig_USresid))){
      all_dat<-cbind(ts_data_resid_tminus1,orig_simflow,orig_USresid,orig_resid)
      all_dat<-as.data.frame(all_dat)
      to_remove<-which(is.na(all_dat[,1]) | is.na(all_dat[,2]) | is.na(all_dat[,3]) | is.na(all_dat[,4]))
    } else {
      all_dat<-cbind(ts_data_resid_tminus1,orig_simflow,orig_resid)
      all_dat<-as.data.frame(all_dat)
      to_remove<-which(is.na(all_dat[,1]) | is.na(all_dat[,2]) | is.na(all_dat[,3]))
    }

    obsflow_dat<-orig_simflow-orig_resid

    if(length(to_remove)>0){
      all_dat<-all_dat[-to_remove,]
      obsflow_dat<-obsflow_dat[-to_remove]
    }
    if(any(obsflow_dat< -1e-10)) stop("Observed flow values should not be less than zero")

    if(!is.na(bootstrap_indices[1])){
      all_dat<-all_dat[bootstrap_indices,]
      obsflow_dat<-obsflow_dat[bootstrap_indices]
    }

    if(normalise_data){
      col.names<-names(all_dat)
      if(!all(is.na(orig_USresid))){
        # only for prev resid, current simflow, US resid
        prevresid_norm_dat<-normalise(all_dat[,1])
        simflow_norm_dat<-normalise(all_dat[,2])
        USresid_norm_dat<-normalise(all_dat[,3])
        all_dat<-cbind(prevresid_norm_dat$normalised_data,
                       simflow_norm_dat$normalised_data,
                       USresid_norm_dat$normalised_data,
                       all_dat[,4])
      } else {
        # only for prev resid and current simflow
        prevresid_norm_dat<-normalise(all_dat[,1])
        simflow_norm_dat<-normalise(all_dat[,2])
        all_dat<-cbind(prevresid_norm_dat$normalised_data,
                       simflow_norm_dat$normalised_data,
                       all_dat[,3])
      }

      colnames(all_dat)<-col.names
      all_dat<-as.data.frame(all_dat)
    }
  } else {
    all_dat<-get_CDNP_clusters_output$all_dat
  }


  if(!is.na(seed)) set.seed(seed)
  resid_sim<-rep(NA,length(simflow))
  prev_error<-initial_resid
  for(dd in 1:length(simflow)){
    # if(dd==40) browser()
    if(dd>1){
      prev_error<-resid_sim[dd-1]
    }
    cur_simflow<-simflow[dd]
    if(!all(is.na(orig_USresid))){
      cur_USerror<-USresid[dd]
    }
    # find the cluster condition
    if(normalise_data){
      # removed data.frame to speed up processing
      # new_data<-data.frame(ts_data_resid_tminus1=normalise(prev_error,
      #                                                      min_vals = prevresid_norm_dat$min_vals,
      #                                                      max_vals = prevresid_norm_dat$max_vals)$normalised_data,
      #                      ts_data_simflow=normalise(cur_simflow,
      #                                                min_vals = simflow_norm_dat$min_vals,
      #                                                max_vals = simflow_norm_dat$max_vals)$normalised_data)
      if(!all(is.na(orig_USresid))){
        new_data<-matrix(c(normalise(prev_error,min_vals = prevresid_norm_dat$min_vals,max_vals = prevresid_norm_dat$max_vals)$normalised_data,
                           normalise(cur_simflow,min_vals = simflow_norm_dat$min_vals,max_vals = simflow_norm_dat$max_vals)$normalised_data,
                           normalise(cur_USerror,min_vals = USresid_norm_dat$min_vals,max_vals = USresid_norm_dat$max_vals)$normalised_data),ncol=3)
      } else {
        new_data<-matrix(c(normalise(prev_error,min_vals = prevresid_norm_dat$min_vals,max_vals = prevresid_norm_dat$max_vals)$normalised_data,
                           normalise(cur_simflow,min_vals = simflow_norm_dat$min_vals,max_vals = simflow_norm_dat$max_vals)$normalised_data),ncol=2)
      }

    } else {
      # removed data.frame to speed up processing
      # new_data<-data.frame(ts_data_resid_tminus1=prev_error,ts_data_simflow=cur_simflow)
      if(!all(is.na(orig_USresid))){
        new_data<-matrix(c(prev_error,cur_simflow,cur_USerror),ncol=3)
      } else {
        new_data<-matrix(c(prev_error,cur_simflow),ncol=2)
      }
    }

    predicted_cluster<-predict_kmeans(new_data,kmeans_model)

    cluster_indices<-which(kmeans_model$cluster==predicted_cluster)
    if(length(cluster_indices)==0) stop("Couldn't find the cluster condition - perhaps run get_CDNP_clusters again")

    if(sample_option==3){
      # sample_option==3 samples both the residual error and the observed flow assuming they are both as likely
      # The samples are treated differently depending whether the residual error or observed flow is sampled
      resid_or_obs<-sample(2,size = 1)
    } else {
      resid_or_obs<-0
    }
    if(sample_option==1 | (sample_option==3 & resid_or_obs==1)){
      # sample_option==1 samples the residual error and subtracts it from the simulated flow
      if(!all(is.na(orig_USresid))){
        errors_to_sample<-all_dat[cluster_indices,4]
      } else {
        errors_to_sample<-all_dat[cluster_indices,3]
      }
      if(!prevent_neg_flow_after_sample){
        if(any(errors_to_sample>cur_simflow)){
          errors_to_sample<-errors_to_sample[-which(errors_to_sample>cur_simflow)]
          errors_to_sample<-c(errors_to_sample,cur_simflow)
        }
      }
      error_sample<-sample(x=errors_to_sample,size=1)

      if(prevent_neg_flow_after_sample){
        # if necessary inverse transforms the data and then see if there are negatives then retransforms
        error_sample_tmp<-error_sample
        cur_simflow_tmp<-cur_simflow
        if(!is.null(simflow_invtrans)){
          if(is.null(simflow_trans)) stop("Missing simflow transformation")
          cur_simflow_tmp<-simflow_invtrans(cur_simflow_tmp)
        }
        if(!is.null(resid_invtrans)){
          if(is.null(resid_trans)) stop("Missing resid transformation")
          error_sample_tmp<-resid_invtrans(error_sample_tmp)
        }
        if(error_sample_tmp>cur_simflow_tmp){
          error_sample_tmp<-cur_simflow_tmp
        }
        if(!is.null(simflow_trans)){
          if(is.null(simflow_invtrans)) stop("Missing simflow inv transformation")
          cur_simflow_tmp<-simflow_trans(cur_simflow_tmp)
        }
        if(!is.null(resid_trans)){
          if(is.null(resid_invtrans)) stop("Missing resid inv transformation")
          error_sample_tmp<-resid_trans(error_sample_tmp)
        }
        error_sample<-error_sample_tmp
        cur_simflow<-cur_simflow_tmp

      }
      resid_sim[dd]<-error_sample

    } else if(sample_option==2 | (sample_option==3 & resid_or_obs==2)){
      # sample_option==2 samples the observed flow and determines what the residual error should be based on that
      obs_to_sample<-obsflow_dat[cluster_indices]
      obs_sample<-sample(x=obs_to_sample,size=1)
      resid_sim[dd]<-cur_simflow-obs_sample

    } else {
      stop("invalid sample_option")
    }

  }
  return(resid_sim)
}

# CDNP_clusters_sim(get_CDNP_clusters_output,simflow=runif(15,0,3),initial_resid=0,seed=1,recompute_all_dat=F)


#' Create lookup table that determines the conditional probabilities
#'
#' @param get_CDNP_clusters_output list; output from get_CDNP_clusters function
#'
#' @return list containing:\cr\cr
#' \verb{    } posterior_lookup: data.frame; table of conditions and probabilities\cr\cr
#' \verb{    } bin_residuals: list; contains the residuals for each bin (for resampling)\cr\cr
#' \verb{    } normalise_data: logical; whether the residual and streamflow are normalised
#' \verb{    } prevresid_norm_dat: normalised previous day residual data
#' \verb{    } simflow_norm_dat: normalised streamflow data
#' @export
#'
#' @examples
#' library(CDNP)
#' # first run get_CDNP_clusters function
#' CDNP_clusters_out<-get_CDNP_clusters(nclusters=10,
#'                                      nbin=8,
#'                                      ts_data_resid=eg_data$resid,
#'                                      ts_data_simflow=eg_data$sim)
#' # get CNDP lookup
#' CDNP_posterior_lookup<-get_CDNP_posterior_lookup(CDNP_clusters_out)
#' CDNP_posterior_lookup
get_CDNP_posterior_lookup<-function(get_CDNP_clusters_output){
  orig_simflow<-get_CDNP_clusters_output$orig_simflow
  orig_resid<-get_CDNP_clusters_output$orig_resid
  resid_intervals<-get_CDNP_clusters_output$resid_intervals
  kmeans_model<-get_CDNP_clusters_output$kmeans_model
  normalise_data<-get_CDNP_clusters_output$normalise_data
  prevresid_norm_dat<-get_CDNP_clusters_output$prevresid_norm_dat
  simflow_norm_dat<-get_CDNP_clusters_output$simflow_norm_dat
  bootstrap_indices<-get_CDNP_clusters_output$bootstrap_indices

  ts_data_resid_tminus1<-c(NA,orig_resid[-length(orig_resid)])
  all_dat<-cbind(ts_data_resid_tminus1,orig_simflow,orig_resid)
  all_dat<-as.data.frame(all_dat)
  to_remove<-which(is.na(all_dat[,1]) | is.na(all_dat[,2]) | is.na(all_dat[,3]))
  if(length(to_remove)>0) all_dat<-all_dat[-to_remove,]

  if(!is.na(bootstrap_indices[1])){
    all_dat<-all_dat[bootstrap_indices,]
  }

  if(normalise_data){
    col.names<-names(all_dat)
    # only for prev resid and current simflow
    prevresid_norm_dat<-normalise(all_dat[,1])
    simflow_norm_dat<-normalise(all_dat[,2])
    all_dat<-cbind(prevresid_norm_dat$normalised_data,
                   simflow_norm_dat$normalised_data,
                   all_dat[,3])
    colnames(all_dat)<-col.names
    all_dat<-as.data.frame(all_dat)
  }

  # compute all priors, likelihoods and normalising constants
  all_posterior_lookup<-matrix(NA,nrow=(length(resid_intervals)-1)*nrow(kmeans_model$centers),ncol=11)

  all_bin_residuals<-list()
  # cat("Computing",nrow(kmeans_model$centers),"clusters. Remaining:")
  counter<-0
  for(rr in 1:nrow(kmeans_model$centers)){
    # cat(".",nrow(kmeans_model$centers)-rr,".")
    for(res_cur_indx in 1:(length(resid_intervals)-1)){
      # res_cur_indx<-1
      # res_prev_indx<-2
      # sf_indx<-17

      # probability of current error is in resid bin - prior
      indices_prior<-which(all_dat$orig_resid>=resid_intervals[res_cur_indx] & all_dat$orig_resid<resid_intervals[res_cur_indx+1])

      num_prior<-length(indices_prior)
      prob_prior<-num_prior/length(all_dat$orig_simflow)

      ind_prior_4like<-indices_prior
      # if(any(ind_prior_4like==1)){
      #   ind_prior_4like<-ind_prior_4like[-which(ind_prior_4like==1)]
      # }

      # probability of prev error is in resid bin A and current sim is in flow bin GIVEN current error is in resid bin B - likelihood

      indices_likelihood<-which(kmeans_model$cluster[ind_prior_4like]==rr)

      num_likelihood<-length(indices_likelihood)
      prob_likelihood<-num_likelihood/length(indices_prior)

      # probability of prev error is in resid bin A and current sim is in flow bin - normalising constant
      indices_normconst<-which(kmeans_model$cluster==rr)
      num_normconst<-length(indices_normconst)
      prob_normconst<-num_normconst/length(all_dat$orig_simflow)

      # posterior probability
      prob_posterior<-prob_prior*prob_likelihood/prob_normconst

      # calculate posterior directly
      indices_posterior<-which(all_dat$orig_resid[indices_normconst]>=resid_intervals[res_cur_indx] & all_dat$orig_resid[indices_normconst]<resid_intervals[res_cur_indx+1])

      # collect the residuals in each bin for resampling later
      bin_residuals<-all_dat$orig_resid[indices_normconst][indices_posterior]

      num_posterior<-length(indices_posterior)
      prob_posterior_direct<-num_posterior/length(indices_normconst)

      # if prob_normconst is zero then the posterior in undefined (not zero)
      # if(is.na(prob_posterior_direct)) prob_posterior_direct<-0
      # if(prob_normconst==0) prob_posterior<-0
      # allow undefined posteriors
      # if(is.na(prob_posterior_direct)) stop("NA posterior")
      # cat(prob_posterior,prob_posterior_direct,"\n")
      # browser()

      counter<-counter+1
      all_posterior_lookup[counter,]<-c(res_cur_indx,
                                        rr,
                                        resid_intervals[res_cur_indx],resid_intervals[res_cur_indx+1],
                                        kmeans_model$centers[rr,1],kmeans_model$centers[rr,2],
                                        prob_prior,
                                        prob_likelihood,
                                        prob_normconst,
                                        prob_posterior,
                                        prob_posterior_direct)
      # all_posterior_direct[counter]<-prob_posterior_direct

      all_bin_residuals[[counter]]<-bin_residuals
    }
  }
  # cat("\n")
  colnames(all_posterior_lookup)<-c("current_error_bin","cluster_row",
                                    "current_error_lower","current_error_upper",
                                    "previous_error_center","simflow_center",
                                    "prior_probability",
                                    "likelihood_probability",
                                    "normalising_constant_probability",
                                    "bayes_posterior_probability",
                                    "direct_posterior_probability")

  # check to see if all probabilities add to 1
  for(cl in 1:nrow(kmeans_model$centers)){
    index<-which(all_posterior_lookup[,2]==cl)
    # all_posterior_lookup[index,]
    if(round(sum(all_posterior_lookup[index,11]),6)!=1) stop("Posteriors for cluster",cl,"do not add to 1 (",sum(all_posterior_lookup[index,11]),")\n")
  }

  return(list(posterior_lookup=as.data.frame(all_posterior_lookup),
              bin_residuals=all_bin_residuals,
              kmeans_model=kmeans_model,
              normalise_data=normalise_data,
              prevresid_norm_dat=prevresid_norm_dat,
              simflow_norm_dat=simflow_norm_dat))
}

#' Predict k-means cluster for a new point based on smallest Euclidean distance to cluster centroids
#'
#' @param new_data data.frame; first column: new error at t-1, second column: new streamflow at t
#' @param kmeans_model output from kmeans function, should be given from get_CDNP_clusters function
#'
#' @return predicted cluster numbers that the new data belongs to
#' @export
#'
#' @examples
#' library(CDNP)
#' # first run get_CDNP_clusters function
#' CDNP_clusters_out<-get_CDNP_clusters(nclusters=10,
#'                                      nbin=8,
#'                                      ts_data_resid=eg_data$resid,
#'                                      ts_data_simflow=eg_data$sim)
#'
#' kmeans_model<-CDNP_clusters_out$kmeans_model
#' kmeans_model$centers
#' # by default the data is normalised.
#' new_data<-data.frame(ts_data_resid_tminus1=0.7,ts_data_simflow=0.11)
#'
#' predict_kmeans_out<-predict_kmeans(new_data,kmeans_model)
#' # predicted cluster should be close to (0.7,0.11)
#' kmeans_model$centers[predict_kmeans_out,]
#'
predict_kmeans <- function(new_data, kmeans_model) {
  centers <- kmeans_model$centers
  n_new_data <- nrow(new_data)
  n_clusters <- nrow(centers)

  # Initialize a vector to store predicted cluster assignments
  predicted_clusters <- numeric(n_new_data)

  for (j in 1:n_new_data) {
    distances <- apply(centers, 1, function(x) {
      dist(rbind(new_data[j, ], x))
    })
    predicted_clusters[j] <- which.min(distances) # index of the FIRST minimum
  }
  return(predicted_clusters)
}


#' Simulate a daily time series of errors based on CDNP lookup table
#'
#' @param get_CDNP_posterior_lookup_output output from get_CDNP_posterior_lookup function
#' @param simflow vector; streamflow
#' @param initial_resid numeric; initial residual error value (default=0)
#' @param seed numeric; optional seed
#' @param sampling_method sampling_method=1 assumes uniform distribution within bin (with truncation to avoid negative flows),
#' sampling_method=2 (default) resamples from bin_residuals (excluding those that cause negative flows), sampling_method=3 picks mid value of bin (avoids negatives)
#'
#' @return vector; predicted errors the same length as simflow
#' @export
#'
#' @examples
#' library(CDNP)
#' # first run get_CDNP_clusters function
#' CDNP_clusters_out<-get_CDNP_clusters(nclusters=10,
#'                                      nbin=8,
#'                                      ts_data_resid=eg_data$resid,
#'                                      ts_data_simflow=eg_data$sim)
#' # get CNDP lookup
#' CDNP_posterior_lookup<-get_CDNP_posterior_lookup(CDNP_clusters_out)
#'
#' # perform error simulation
#' CDNP_sim_out<-CDNP_sim(CDNP_posterior_lookup,eg_data$sim,seed=42)
#' plot(eg_data$resid,type="l") # original residual
#' lines(CDNP_sim_out,col=2,lty=2) # a new simulated residual
CDNP_sim<-function(get_CDNP_posterior_lookup_output,simflow,initial_resid=0,seed=NA,sampling_method=2){
  # simflow<-orig_simflow
  # initial_resid<-0
  posterior_lookup<-get_CDNP_posterior_lookup_output$posterior_lookup
  bin_residuals<-get_CDNP_posterior_lookup_output$bin_residuals
  kmeans_model<-get_CDNP_posterior_lookup_output$kmeans_model
  normalise_data<-get_CDNP_posterior_lookup_output$normalise_data
  prevresid_norm_dat<-get_CDNP_posterior_lookup_output$prevresid_norm_dat
  simflow_norm_dat<-get_CDNP_posterior_lookup_output$simflow_norm_dat

  if(!is.na(seed)) set.seed(seed)
  resid_sim<-rep(NA,length(simflow))
  prev_error<-initial_resid
  for(dd in 1:length(simflow)){
    # if(dd==40) browser()
    if(dd>1){
      prev_error<-resid_sim[dd-1]
    }
    cur_simflow<-simflow[dd]
    # find the cluster condition
    if(normalise_data){
      # removed data.frame to speed up processing
      # new_data<-data.frame(ts_data_resid_tminus1=normalise(prev_error,
      #                                                      min_vals = prevresid_norm_dat$min_vals,
      #                                                      max_vals = prevresid_norm_dat$max_vals)$normalised_data,
      #                      ts_data_simflow=normalise(cur_simflow,
      #                                                min_vals = simflow_norm_dat$min_vals,
      #                                                max_vals = simflow_norm_dat$max_vals)$normalised_data)
      new_data<-matrix(c(normalise(prev_error,min_vals = prevresid_norm_dat$min_vals,max_vals = prevresid_norm_dat$max_vals)$normalised_data,
                         normalise(cur_simflow,min_vals = simflow_norm_dat$min_vals,max_vals = simflow_norm_dat$max_vals)$normalised_data),ncol=2)
    } else {
      # removed data.frame to speed up processing
      # new_data<-data.frame(ts_data_resid_tminus1=prev_error,ts_data_simflow=cur_simflow)
      new_data<-matrix(c(prev_error,cur_simflow),ncol=2)
    }

    predicted_cluster<-predict_kmeans(new_data,kmeans_model)

    lookup_row_indices<-which(posterior_lookup$cluster_row==predicted_cluster)
    if(length(lookup_row_indices)==0) stop("Couldn't find the cluster condition - perhaps rerun get_CDNP_clusters again")

    # posterior_lookup[posterior_lookup$cluster_row==predicted_cluster,]
    # sum(posterior_lookup$direct_posterior_probability[lookup_row_indices])

    sampled_lookup_row_index<-sample(lookup_row_indices,1,prob = posterior_lookup$direct_posterior_probability[lookup_row_indices])
    # posterior_lookup[sampled_lookup_row_index,]

    # ensure the bounds can't go negative
    # obs = sim - error
    if(sampling_method==1){
      # assumes uniform distribution
      min_bound<-min(posterior_lookup$current_error_lower[sampled_lookup_row_index],cur_simflow)
      max_bound<-min(posterior_lookup$current_error_upper[sampled_lookup_row_index],cur_simflow)
      error_sample<-runif(1,min_bound,max_bound)
    } else if(sampling_method==2){
      # resamples from the bin
      # ensure no negative flows
      bin_to_sample<-bin_residuals[[sampled_lookup_row_index]]
      if(any(bin_to_sample>cur_simflow)){
        bin_to_sample<-bin_to_sample[-which(bin_to_sample>cur_simflow)]
        bin_to_sample<-c(bin_to_sample,cur_simflow)
      }
      error_sample<-sample(bin_to_sample,1)
    } else if(sampling_method==3){
      # middle value
      min_bound<-min(posterior_lookup$current_error_lower[sampled_lookup_row_index],cur_simflow)
      max_bound<-min(posterior_lookup$current_error_upper[sampled_lookup_row_index],cur_simflow)
      error_sample<-mean(c(min_bound,max_bound))
    } else if(sampling_method==4){
      # median
      bin_to_sample<-bin_residuals[[sampled_lookup_row_index]]
      error_sample<-min(median(bin_to_sample),cur_simflow)
    } else if(sampling_method==5){
      # mean
      bin_to_sample<-bin_residuals[[sampled_lookup_row_index]]
      error_sample<-min(mean(bin_to_sample),cur_simflow)
    } else if(sampling_method==6){
      # resamples from the bin
      # ensure no negative flows by setting to zero
      bin_to_sample<-bin_residuals[[sampled_lookup_row_index]]
      error_sample<-sample(bin_to_sample,1)
      if(error_sample>cur_simflow) error_sample<-cur_simflow
    } else {
      stop("Not legitimate sampling_method")
    }

    resid_sim[dd]<-error_sample
    # if(error_sample < (-20)) browser()
  }
  return(resid_sim)
}

# just checks to see if period is split (accounting for warmup) if there's valid data in each period
check_data_ok<-function(obs_or_resid,warmup=1095,number_required=1){
  if(warmup>0){
    obs_or_resid_nowarm<-obs_or_resid[-(1:(warmup))]
  } else {
    obs_or_resid_nowarm<-obs_or_resid
  }

  split_index<-ceiling(length(obs_or_resid_nowarm)/2)

  obs_or_resid_spl1<-obs_or_resid_nowarm[1:split_index]
  obs_or_resid_spl2<-obs_or_resid_nowarm[(split_index+1):length(obs_or_resid_nowarm)]
  if(number_required==1){
    data_ok<-any(!is.na(obs_or_resid_spl1)) & any(!is.na(obs_or_resid_spl2))
  } else {
    data_ok<-min(length(which(!is.na(obs_or_resid_spl1))),length(which(!is.na(obs_or_resid_spl2))))>=number_required
  }

  return(data_ok)
}

#' Uses Differential Evolution optimiser (DEoptim) to determine the number of clusters and bins that give the most reliable results
#'
#' @description
#' This splits the input series into two halves. The first half is used to compute the clusters
#' and lookup tables. These are used to simulate residuals for the second half and these are
#' checked against the actual residuals of the second using the  reliability metric (alpha) from
#' Renard et al. (2010).
#'
#' @param resid vector; residual series
#' @param simflow vector; streamflow series
#' @param warmup integer; length of data that will be exluded at the beginning of each split period (default=1095)
#' @param nrep integer; number of residual simulation replicates (default=10)
#' @param itermax integer; maximum number of iterations for DEoptim (default=50)
#' @param lower vector of integers; lower bounds for ncluster and nbin (default=c(2,2))
#' @param upper vector of integers; upper bounds for ncluster and nbin (default=c(100,100))
#'
#' @return DEoptim object; optimal ncluster and nbin are in $optim$bestmem
#' @export
#'
#' @references Renard, B., Kavetski, D., Kuczera, G., Thyer, M., & Franks, S. W. (2010). Understanding predictive uncertainty in hydrologic modeling: The challenge of identifying input and structural errors. Water Resources Research, 46(5).
#'
#' @examples
#' library(CDNP)
#' optimal_ncluster_and_nbins_out<-get_optimal_ncluster_and_nbins(resid=eg_data$resid,
#'                                      simflow=eg_data$sim,warmup=30,nrep=2,itermax=3)
#' # Best member (ncluster and nbin)
#' optimal_ncluster_and_nbins_out$optim$bestmem
#'
get_optimal_ncluster_and_nbins<-function(resid,simflow,warmup=1095,nrep=10,itermax=50,
                                         lower=c(2,2),upper=c(100,100),sampling_method=2,
                                         optim_method=4){
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

  compute_alpha_for_optim<-function(par,nrep=nrep,use_round=F){
    if(use_round){
      ncluster<-round(par[1])
      nbin<-round(par[2])
    } else {
      ncluster<-ceiling(par[1])
      nbin<-ceiling(par[2])
    }


    cluster_out_spl<-get_CDNP_clusters(ncluster=ncluster,nbin=nbin,ts_data_resid=resid_nowarm_spl1,
                                       ts_data_simflow=simflow_nowarm_spl1,use_quantile_spacing=T,seed=45)
    all_posterior_lookup_spl<-get_CDNP_posterior_lookup(cluster_out_spl)

    all_CDNP_sim_out<-matrix(NA,nrow=length(simflow_withwarm_spl2)-warmup,ncol=nrep)
    for(rrr in 1:nrep){
      CDNP_sim_out<-CDNP_sim(all_posterior_lookup_spl,simflow_withwarm_spl2,seed=rrr,sampling_method=sampling_method)
      # lines(CDNP_sim_out,col=2,lty=2)
      CDNP_sim_out_nowarm<-CDNP_sim_out[-(1:warmup)]
      all_CDNP_sim_out[,rrr]<-CDNP_sim_out_nowarm
    }

    # compute_alpha_reliability(all_CDNP_sim_out,orig_resid_nowarm_spl1)
    alpha<-compute_alpha_reliability(all_CDNP_sim_out,resid_nowarm_spl2)
    cat("Sampled pars:",par,"Alpha:",alpha,"\n")
    return(-alpha)
  }

  # opt<-optim(par<-as.integer(c(10,10)),compute_alpha_for_optim,method="SANN",control=list(fnscale=-1,trace=1))
  # opt<-optim(par<-as.integer(c(10,10)),compute_alpha_for_optim,control=list(fnscale=-1,trace=1))


  if(optim_method==1){
    # library(DEoptim)
    # Mapping function to ensure integers
    map_fun <- function(x) {
      return(ceiling(x))
    }
    deopt<-DEoptim::DEoptim(compute_alpha_for_optim,lower=lower,upper=upper,
                            fnMap = map_fun,control=list(trace=1,itermax=itermax),nrep=nrep)
  } else if(optim_method==2){

    var_length<-round(sqrt(itermax))

    trial_cluster<-round(seq(lower[1],upper[1],length.out=var_length))
    trial_bin<-round(seq(lower[2],upper[2],length.out=var_length))
    exgrid<-as.matrix(expand.grid(trial_cluster,trial_bin))
    all_alpha<-c()
    for(rrr in 1:nrow(exgrid)){
      alph<-compute_alpha_for_optim(c(exgrid[rrr,]),nrep=nrep)
      all_alpha<-c(all_alpha,alph)
    }
    best_nclust_nbin<-c(exgrid[which.min(all_alpha),])
    best_alpha<-all_alpha[which.min(all_alpha)]

    deopt<-list(optim=list(bestmem=best_nclust_nbin,bestval=best_alpha))
  } else if(optim_method==3){
    map_fun <- function(x) {
      return(round(x))
    }
    deopt<-DEoptim::DEoptim(compute_alpha_for_optim_round,lower=lower,upper=upper,
                            fnMap = map_fun,control=list(trace=1,itermax=itermax),nrep=nrep)
  } else if(optim_method==4){
    map_fun <- function(par) {
      par_j <- par + runif(length(par), -0.49, 0.49)  # jitter before rounding


      random_int_vector <- function(lower, upper) {
        mapply(function(lo, up) sample(lo:up, 1), lo = lower, up = upper)
      }
      # occasionally re-randomize (say 5% chance)
      if (runif(1) < 0.05) {
        par_j <- random_int_vector(lower, upper)
      }

      return(pmin(pmax(round(par_j), lower), upper))          # round + clip to bounds
    }
    deopt<-DEoptim::DEoptim(compute_alpha_for_optim,lower=lower,upper=upper,
                            fnMap = map_fun,control=list(trace=1,itermax=itermax,CR=0.9,NP=30,F=0.6),nrep=nrep,
                            use_round=T)

  } else if(optim_method==5){
    random_int_vector <- function(lower, upper) {
      mapply(function(lo, up) sample(lo:up, 1), lo = lower, up = upper)
    }
    map_fun <- function(par) {
      par_j <- par + runif(length(par), -0.49, 0.49)  # jitter before rounding

      # occasionally re-randomize (say 5% chance)
      # if (runif(1) < 0.05) {
      #   par_j <- random_int_vector(lower, upper)
      # }

      return(pmin(pmax(round(par_j), lower), upper))          # round + clip to bounds
    }
    deopt<-DEoptim::DEoptim(compute_alpha_for_optim,lower=lower,upper=upper,
                            fnMap = map_fun,control=list(trace=1,itermax=4,CR=0.9,NP=30,F=0.6),nrep=nrep,
                            use_round=T)
    for(iter in 1:round(itermax/4)){
      cat(iter,"/",round(itermax/4),"\n")
      new_pop<-deopt$member$pop
      indices_to_replace<-which(runif(nrow(new_pop))<0.05)
      if(length(indices_to_replace)>0){
        for(rrr in 1:length(indices_to_replace)){
          new_pop[indices_to_replace[rrr],]<-random_int_vector(lower, upper)
        }

      }
      deopt<-DEoptim::DEoptim(compute_alpha_for_optim,lower=lower,upper=upper,
                              fnMap = map_fun,control=list(trace=1,itermax=4,CR=0.9,NP=30,F=0.6,
                                                           initialpop=new_pop),nrep=nrep,
                              use_round=T)
    }


  } else if(optim_method==6){
    # reseed with totally new polulations
    random_int_vector <- function(lower, upper) {
      mapply(function(lo, up) sample(lo:up, 1), lo = lower, up = upper)
    }
    map_fun <- function(par) {
      par_j <- par + runif(length(par), -0.49, 0.49)  # jitter before rounding

      # occasionally re-randomize (say 5% chance)
      if (runif(1) < 0.05) {
        par_j <- random_int_vector(lower, upper)
      }

      return(pmin(pmax(round(par_j), lower), upper))          # round + clip to bounds
    }
    itermax_part<-max(round(itermax/5),1)
    deopt<-DEoptim::DEoptim(compute_alpha_for_optim,lower=lower,upper=upper,
                            fnMap = map_fun,control=list(trace=1,itermax=itermax_part,CR=0.9,NP=30,F=0.6),nrep=nrep,
                            use_round=T)

    for(iter in 1:max(round(itermax/itermax_part),1)){
      cat(iter,"/",max(round(itermax/itermax_part),1),"\n")
      new_pop<-deopt$member$pop
      best_so_far <- deopt$optim$bestmem
      new_pop[1,]<-best_so_far
      for(rrr in 2:nrow(new_pop)){
        new_pop[rrr,]<-random_int_vector(lower, upper)
      }

      deopt<-DEoptim::DEoptim(compute_alpha_for_optim,lower=lower,upper=upper,
                              fnMap = map_fun,control=list(trace=1,itermax=itermax_part,CR=0.9,NP=30,F=0.6,
                                                           initialpop=new_pop),nrep=nrep,
                              use_round=T)
    }

  } else {
    stop("optim_method not valid")
  }


  # all_pars<-as.matrix(expand.grid(10:30,10:30))
  # all_alpha<-rep(NA,nrow(all_pars))
  # for(pp in 1:nrow(all_pars)){
  #   cat(pp,"/",nrow(all_pars),"\n")
  #   pars<-all_pars[pp,]
  #   all_alpha[pp]<-compute_alpha_for_optim(pars)
  # }
  return(deopt)
}


get_optimal_ncluster<-function(resid,simflow,USresid=NA,warmup=1095,nrep=10,itermax=20,
                                         lower=5,upper=100,prevent_neg_flow_after_sample=T,
                               simflow_trans=NULL,resid_trans=NULL,simflow_invtrans=NULL,resid_invtrans=NULL,
                               sample_option=1,fold=F,...){
  # below is just to ensure old scripts work with inaccurate argument name - should eventually be removed
  if(exists("prevent_zeroflow_after_sample")){
    prevent_neg_flow_after_sample<-prevent_zeroflow_after_sample
  }

  if(length(resid)!=length(simflow)) stop("length of residual is different to length of streamflow")
  if(warmup>length(simflow)) stop("warmup is longer than the length of time series")
  if(!check_data_ok(resid)) stop("not enough data to split the time series")
  resid_nowarm<-resid[-(1:(warmup))]
  simflow_nowarm<-simflow[-(1:(warmup))]
  if(!all(is.na(USresid)) & length(simflow)==length(USresid)){
    USresid_nowarm<-USresid[-(1:(warmup))]
  }

  split_index<-ceiling(length(resid_nowarm)/2)

  resid_nowarm_spl1<-resid_nowarm[1:split_index]
  resid_nowarm_spl2<-resid_nowarm[(split_index+1):length(resid_nowarm)]

  simflow_nowarm_spl1<-simflow_nowarm[1:split_index]
  simflow_nowarm_spl2<-simflow_nowarm[(split_index+1):length(simflow_nowarm)]

  if(!all(is.na(USresid)) & length(simflow)==length(USresid)){
    USresid_nowarm_spl1<-USresid_nowarm[1:split_index]
    USresid_nowarm_spl2<-USresid_nowarm[(split_index+1):length(USresid_nowarm)]
  } else {
    USresid_nowarm_spl1<-NA
    USresid_nowarm_spl2<-NA
  }

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

  if(!all(is.na(USresid)) & length(simflow)==length(USresid)){
    USresid_withwarm_spl1<-USresid[1:(split_index+warmup)]
    USresid_withwarm_spl2<-USresid[((split_index+1)-warmup):length(USresid_nowarm)]
  } else {
    USresid_withwarm_spl1<-NA
    USresid_withwarm_spl2<-NA
  }

  compute_alpha_for_optim<-function(par,nrep=nrep,setup_first_half=T){

    ncluster<-ceiling(par[1])
    # nbin<-ceiling(par[2])
    if(setup_first_half){
      resid_nowarm_setup<-resid_nowarm_spl1
      simflow_nowarm_setup<-simflow_nowarm_spl1
      USresid_nowarm_setup<-USresid_nowarm_spl1

      simflow_withwarm_test<-simflow_withwarm_spl2
      USresid_withwarm_test<-USresid_withwarm_spl2
      resid_nowarm_test<-resid_nowarm_spl2
    } else {
      resid_nowarm_setup<-resid_nowarm_spl2
      simflow_nowarm_setup<-simflow_nowarm_spl2
      USresid_nowarm_setup<-USresid_nowarm_spl2

      simflow_withwarm_test<-simflow_withwarm_spl1
      USresid_withwarm_test<-USresid_withwarm_spl1
      resid_nowarm_test<-resid_nowarm_spl1
    }

    cluster_out_spl<-get_CDNP_clusters(ncluster=ncluster,nbin=NA,ts_data_resid=resid_nowarm_setup,
                                       ts_data_simflow=simflow_nowarm_setup,ts_data_USresid=USresid_nowarm_setup,use_quantile_spacing=T,seed=45)
    # all_posterior_lookup_spl<-get_CDNP_posterior_lookup(cluster_out_spl)

    all_CDNP_sim_out<-matrix(NA,nrow=length(simflow_withwarm_test)-warmup,ncol=nrep)
    for(rrr in 1:nrep){
      CDNP_sim_out<-CDNP_clusters_sim(cluster_out_spl,simflow_withwarm_test,USresid=USresid_withwarm_test,seed=rrr,recompute_all_dat=F,
                                      prevent_neg_flow_after_sample=prevent_neg_flow_after_sample,
                                      simflow_trans=simflow_trans,
                                      resid_trans=resid_trans,
                                      simflow_invtrans=simflow_invtrans,
                                      resid_invtrans=resid_invtrans,
                                      sample_option=sample_option)
      # lines(CDNP_sim_out,col=2,lty=2)
      CDNP_sim_out_nowarm<-CDNP_sim_out[-(1:warmup)]
      all_CDNP_sim_out[,rrr]<-CDNP_sim_out_nowarm
    }

    # compute_alpha_reliability(all_CDNP_sim_out,orig_resid_nowarm_spl1)
    alpha<-compute_alpha_reliability(all_CDNP_sim_out,resid_nowarm_test)
    cat("Sampled pars:",par,"Alpha:",alpha,"setup_first_half:",setup_first_half,"\n")
    return(-alpha)
  }

  # opt<-optim(par<-as.integer(c(10,10)),compute_alpha_for_optim,method="SANN",control=list(fnscale=-1,trace=1))
  # opt<-optim(par<-as.integer(c(10,10)),compute_alpha_for_optim,control=list(fnscale=-1,trace=1))

  all_alpha<-c()
  trial_k<-round(seq(lower,upper,length.out=itermax))
  # trial_k<-c(150,200,300,400,500)
  for(kkk in trial_k){
    alph<-compute_alpha_for_optim(kkk,nrep=nrep)
    if(fold){
      alph2<-compute_alpha_for_optim(kkk,nrep=nrep,setup_first_half=F)
      alph<-(alph+alph2)/2
    }
    all_alpha<-c(all_alpha,alph)
  }

  best_k<-trial_k[which.min(all_alpha)]
  best_alpha<-all_alpha[which.min(all_alpha)]
  return(list(best_k=best_k,best_alpha=best_alpha))

  # optimize(compute_alpha_for_optim, interval = c(lower, upper), tol = 1,nrep=nrep)
  # opt<-optim(50,fn=compute_alpha_for_optim,method="Brent",lower=lower,upper=upper,control=list(reltol=0.01),nrep=2)

  # library(DEoptim)
  # Mapping function to ensure integers
  # map_fun <- function(x) {
  #   return(ceiling(x))
  # }
  # deopt<-DEoptim::DEoptim(compute_alpha_for_optim,lower=lower,upper=upper,
  #                         control=list(trace=1,itermax=itermax),nrep=nrep) # fnMap = map_fun


  # all_pars<-as.matrix(expand.grid(10:30,10:30))
  # all_alpha<-rep(NA,nrow(all_pars))
  # for(pp in 1:nrow(all_pars)){
  #   cat(pp,"/",nrow(all_pars),"\n")
  #   pars<-all_pars[pp,]
  #   all_alpha[pp]<-compute_alpha_for_optim(pars)
  # }
  # return(deopt)

}
