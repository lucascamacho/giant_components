coevo_net <- function(x, net, netname, matrix_type="adjacency", mi, phi, alpha, t_max, eps) {
  

  #Obs: x variable in the function definition is a dummy variable to be used in the lapply. Because its not used in the function code, we pass it to lapply n times to run n simulations
  
  netn<-netname #Getting the name of the network to use in the data frame of results
  
  net[net > 0] <- 1 #Changing matrix to binary
  
  #-----------------------------------#
  #-Setting up matrix for simulations-#
  #-----------------------------------#
  
  #We need an adjacency matrix, thus, if the user enters an incidence matrix, we first need to set it up as an adjacency matrix
  
  if(matrix_type=="adjacency"){ A <- net} else if(matrix_type=="incidence"){
    
    n_p <- nrow(net) #Getting the number of species of the row set of the incidence matrix
    n_a <- ncol(net) #Getting the number of species of the column set of the incidence matrix
    
    # Creating adjacency matrix from incidence matrix
    A <- rbind(cbind(matrix(0, n_p, n_p), net),
             cbind(t(net), matrix(0, n_a, n_a)))
  }
  
  n_sp<-nrow(A) #Getting the total number of species
  rownames(A)<-colnames(A)<-paste0("SP", 1:n_sp) #Setting row and column names
  
  m<-rep(mi, n_sp) #Setting up m values for each species - this allows different values for each species
  phi<-rep(phi, n_sp) #Setting up phi values for each species
  
  
  #randomly sampling trait and theta values
  
  theta <- runif(n_sp, min = 0, max = 10) #Sampling theta from a uniform distribution U[0, 10]
  init <- runif(n_sp, min = 0, max = 10) #Sampling initial z values from a uniform distribution U[0, 10]

  # creating matrix to store z values
  z_mat <- matrix(NA, nrow = t_max, ncol = n_sp)
  # adding initial trait values
  z_mat[1, ] <- init  
  
  for (t in 1:(t_max - 1)) {
    # defining current z values
    z <- z_mat[t, ]
    # creating matrix with all trait differences
    z_dif <- (A %*% diag(z)) - (diag(z) %*% A) 
    # calculating matrix Q
    Q <- A * (exp( - alpha * (z_dif^2)))
    # normalizing matrix Q
    Q <- Q / apply(Q, 1, sum)
    Q[is.na(Q)]<-0
    # multiplying each row i of matrix Q by m[i]
    Q <- Q * m
    qij_sum <- rowSums(Q) #Sum of qij=m, thus, for species that do not interact sum_qij=0.
    # multiplying each entry in Q by corresponding trait difference
    Q_dif <- Q * z_dif
    # calculating mutualistic selection for each species i
    mut_sel <- apply(Q_dif, 1, sum)
    # calculating environmental selection for each species i
    env_sel <- (1 - qij_sum) * (theta - z)
    # storing z values of next timestep
    
    z_mat[t+1, ] = z + phi *( mut_sel + env_sel)
    
    dif <- mean(abs(z - z_mat[t + 1, ])) 
    if (dif < eps){break}
    
  }
  
  t_eq<-t #Time until equilibrium
  
  #Indirect effects and equilibrium trait values
  
  degree<-rowSums(A)
  
  I <- diag(1, sum(!!degree)) #Expression sum(!!degree) counts the number of species with at least one interaction
  Psi<-diag((1-mi), sum(!!degree)) #Species completely disconnected will be excluded from indirect effects calculation
  
  z_eq<-z_mat[t, ] #Equilibrium trait values
  
  #---------------------------------------#
  #-Q_matrix at the end of the simulation-#
  #---------------------------------------#
  
  zeq_dif<-(A %*% diag(z_eq)) - (diag(z_eq) %*% A) #Matrix of trait differences
  
  Q_eq_end<-A * (exp( -alpha * (zeq_dif^2))) # Matrix of evolutionary effects
  Q_eq_end<-Q_eq_end / apply(Q_eq_end, 1, sum) #Normalizing qij
  Q_eq_end<-Q_eq_end*m #Multiplying qij by m
  
  Q_eq_end<-Q_eq_end[which(degree!=0), which(degree!=0)] #Subsetting species with at least one interaction
  
  T_mat <- solve((I - Q_eq_end)) %*% Psi #Getting matrix of direct+indirect effects
  diag(T_mat)<-0 #Removing "self" indirect effects
  
  T_mat_ind<-T_mat #Setting up matrix of only indirect effects
  T_mat_ind[Q_eq_end>0]<-0 #Setting to 0 all direct effects in T-matrix
  
  ind_effects<-rowSums(T_mat_ind)/rowSums(T_mat) #Getting proportional contribution of indirect effects to trait evolution of each species
  sp_indirect_effects<-rep(0, nrow(A)) #Creating vector to store indirect effects for all species
  aux<-match(rownames(T_mat), rownames(A)) #Ordering the vector of indirect effects
  sp_indirect_effects[aux]<-ind_effects #Setting indirect effects, species that receive 0 indirect effects remain as 0
  
  graph<-graph_from_adjacency_matrix(A, mode=c("undirected"))
  component_id<-components(graph)$membership
  
  component_tm<-comp_tm(component_id, z_eq, alpha=alpha, sp_id=rownames(A))
  network_tm<-global_tm(z_eq, alpha=alpha)
  partners_tm<-species_tm(A=A, z=z_eq, alpha=alpha) #Average trait matching of one species with its mutualistic partners
  
  output<-data.frame(network_id=netn, component_id, sp_id=rownames(A), 
                     z_init=z_mat[1,], z_eq, theta, degree, sp_indirect_effects, 
                    partners_tm, component_tm, network_tm, m, phi, alpha, t_eq, nsim=x)
  
  return(output)
  
}

tmatch<-function(x,y, alpha){
  tm<-exp(-alpha*(y-x)^2)
}

global_tm<-function(z, alpha){
  tm<-outer(z,z, FUN=tmatch, alpha=alpha)
  diag(tm)<-NA
  return(rowMeans(tm, na.rm=TRUE))
}

species_tm<-function(A,z,alpha){
  tm<-outer(z,z, FUN=tmatch, alpha=alpha)
  tm<-A*tm
  tm[tm==0]<-NA
  return(rowMeans(tm, na.rm=TRUE))
}

comp_tm<-function(comp_id, z, alpha, sp_id){
  
  names(z)<-paste0("SP", 1:length(z))
  tm_list<-list()
  
  for(i in 1:length(unique(comp_id))){
    
    z_i<-z[which(comp_id==unique(comp_id)[i])]
    tm_i<-outer(z_i, z_i, FUN=tmatch, alpha=alpha)
    diag(tm_i)<-NA
    
    mean_tm<-rowMeans(tm_i, na.rm=TRUE)
    tm_list[[i]]<-mean_tm
    
  }
  tm_values<-do.call(c, tm_list)
  aux<-match(sp_id, names(tm_values))
  
  return(tm_values[aux])
  
}

comp_thetavar<-function(comp_id, theta){
  
  var_list<-list()
  
  for(i in 1:length(unique(comp_id))){
    
    theta_i<-theta[which(comp_id==unique(comp_id)[i])]
    var_list[[i]]<-data.frame(comp_id=comp_id[i], var_theta=var(theta_i))
  
  }
  
  return(do.call(rbind, var_list))
  
}

coevo_sim <- function(nsim, net_list, matrix_type="adjacency", mi, phi, alpha, t_max, eps){
  
  results_list<-vector("list", length(net_list))
  
  for(l in 1:length(net_list)){
    
    network<-net_list[[l]]
    net_name<-names(net_list)[l]
    
    results<-sfLapply(1:nsim, fun = coevo_net, net=network, matrix_type=matrix_type, netname=net_name, 
                     mi=mi, phi=phi, alpha=alpha, t_max=t_max, eps=eps)
    
    sim_results<-do.call(rbind, results)
    
    results_list[[l]]<-sim_results
    
    print(l)
    
  }
  
  r_final<-do.call(rbind, results_list)
  
  return(r_final)
  
}