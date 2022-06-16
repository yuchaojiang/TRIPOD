transformBeta<-function(betas){
  sign(betas)*log(1+abs(betas))
}
