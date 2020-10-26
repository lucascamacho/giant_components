n_sp <- 4
M <- matrix(1, nrow = n_sp, ncol = n_sp)
diag(M) = 0
M

m <- 0.9
phi <- 0.2
alpha <- 0.2
theta <- runif(n_sp, 0, 10)
init <- runif(n_sp, 0, 10)
eps <- 0.0001
t_max <- 1000

CoevoMutNet <- function(n_sp, M, m, phi, alpha, theta, init, eps, t_max) {
  z_mat <- matrix(NA, ncol = n_sp, nrow = t_max)
  z_mat[1,] <- init
  
  for(i in 1:(t_max - 1)) {
    z <- z_mat[i, ]
    z_dif <- t(M * z) - (M * z)
    Q <- M * exp((-alpha * (z_dif^2)))
    Q_n <- Q / apply(Q, 1, sum)
    
    env <- phi * (1 - m) * (theta - z)
    
    mut <- Q_n * m * z_dif
    mut <- phi * apply(mut, 1, sum)
    
    z_mat[i+1, ] <- z + env + mut
    
    dif <- mean(abs(z - z_mat[i+1,]))
    if(dif < eps){break}
  }
  return(z_mat[1:i+1, ])
}

CoevoMutNet(n_sp, M, m, phi, alpha, theta, init, eps, t_max)
