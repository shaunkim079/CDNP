
# transform for streamflow
logplconst<-function(x,e){
  trans_x<-log(x+e)
  return(trans_x)
}
inv_logplconst<-function(trans_x,e){
  x<-exp(trans_x)-e
  return(x)
}

# transform for residual
trans_asinh <- function(x, k = 0.275) asinh(x / k)
inv_trans_asinh <- function(x, k = 0.275) sinh(x)*k
dble_trans_asinh<-function(x, k = 0.275) trans_asinh(trans_asinh(x,k),k)
inv_dble_trans_asinh<-function(x, k = 0.275) inv_trans_asinh(inv_trans_asinh(x,k),k)
