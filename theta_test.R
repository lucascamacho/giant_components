library(ggplot2)
library(ggpubr)

theta_variance<-function(n, min, max, n_sim){
  
  r<-lapply(1:n_sim, function(x, n, min, max){
    
    theta<-runif(n=n, min=min, max=max);
    return(data.frame(n=n, var_theta=var(theta)))}, 
    n=n, min=min, max=max)
  
  r<-do.call(rbind, r)
  
  return(r)
  
}

results<-lapply(2:100, theta_variance, n_sim=1000, min=0, max=10)
results<-do.call(rbind, results)

p<-ggline(data=results, x="n", y="var_theta", add="mean_sd")
ggpar(p, xlab="Number of nodes", ylab="Observed variance in theta values")
