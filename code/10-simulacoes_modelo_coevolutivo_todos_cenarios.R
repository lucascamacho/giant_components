#-------------------#
#-Loading libraries-#
#-------------------#

library(snowfall) #Parallel simulations
library(parallel) #Detect number of CPU cores available for parallel simulations
library(here) #Get directories within the same folder of the .Rproj file
library(igraph) #Tools for network analysis - will be used to detect network components

#-------------------#
#-Loading functions-#
#-------------------#

source(here("code/functions", "coevo_functions.R")) #Coevolution in networks and function that summarize results

#------------------#
#-Loading networks-#
#------------------#

uni_files <- list.files(here("data/uni_matrices"), pattern="*.txt") #Network names
uni_path <- list.files(here("data/uni_matrices"), pattern="*.txt", full.names=TRUE) #Networks file paths
uni_list <- lapply(uni_path, read.table, header=FALSE) #Reading networks
uni_list <- lapply(uni_list, as.matrix) #Converting all networks to matrices
names(uni_list)<-uni_files #Setting networks names

uni_verif<-do.call(c, lapply(uni_list, function(x){return(sum(rowSums(x)))})) #Verifying for networks without any connections
uni_ok<-uni_list[which(uni_verif > 0)]

bip1_files <- list.files(here("data/bip_matrices"), pattern="*.txt") #Network names
bip1_path <- list.files(here("data/bip_matrices"), pattern="*.txt", full.names=TRUE) #Networks file paths
bip1_list <- lapply(bip1_path, read.table, header=FALSE) #Reading networks
bip1_list <- lapply(bip1_list, as.matrix) #Converting all networks to matrices
names(bip1_list)<-bip1_files #Setting networks names

bip1_verif<-do.call(c, lapply(bip1_list, function(x){return(sum(rowSums(x)))})) #Verifying networks - excluding network 74 (no connection)
bip1_ok<-bip1_list[which(bip1_verif>0)]

bip2_files <- list.files(here("data/bip_matrices_asymmetric"), pattern="*.txt") #Network names
bip2_path <- list.files(here("data/bip_matrices_asymmetric"), pattern="*.txt", full.names=TRUE) #Networks file paths
bip2_list <- lapply(bip2_path, read.table, header=FALSE) #Reading networks
bip2_list <- lapply(bip2_list, as.matrix) #Converting all networks to matrices
names(bip2_list)<-bip2_files #Setting networks names

bip2_verif<-do.call(c, lapply(bip2_list, function(x){return(sum(rowSums(x)))})) #Verifying networks - excluding network 74 (no connection)
bip2_ok<-bip2_list[which(bip2_verif>0)]

emp_files <- list.files(here("data/empirical_bip"), pattern="*.txt") #Network names
emp_path <- list.files(here("data/empirical_bip"), pattern="*.txt", full.names=TRUE) #Networks file paths
emp_list <- lapply(emp_path, read.table, header=TRUE) #Reading networks
emp_list <- lapply(emp_list, as.matrix) #Converting all networks to matrices
names(emp_list)<-emp_files #Setting networks names

emp_verif<-do.call(c, lapply(emp_list, function(x){return(sum(rowSums(x)))})) #Verifying networks - excluding network 74 (no connection)
emp_ok<-emp_list[which(emp_verif>0)]

#------------------------------#
#-Setting parallel simulations-#
#------------------------------#

sfInit(parallel = TRUE, cpus = 12, type = "SOCK") #Initializing parallel clusters
sfExportAll() #Export all objects in Global Environment - include user-defined functions
sfLibrary(igraph) #Loading libraries inside clusters

m_values<-seq(from=0.1, to=0.9, by=0.1) #Setting up m values to use in simulations, from 0.1 to 0.9

uni_results<-do.call(rbind, lapply(m_values, FUN=coevo_sim, net_list=uni_ok, matrix_type="adjacency", phi=0.2, alpha=0.2, t_max=1000, eps=0.0001, nsim=100))
bip1_results<-do.call(rbind, lapply(m_values, FUN=coevo_sim, net_list=bip1_ok, matrix_type="incidence", phi=0.2, alpha=0.2, t_max=1000, eps=0.0001, nsim=100))
bip2_results<-do.call(rbind, lapply(m_values, FUN=coevo_sim, net_list=bip2_ok, matrix_type="incidence", phi=0.2, alpha=0.2, t_max=1000, eps=0.0001, nsim=100))
emp_results<-do.call(rbind, lapply(m_values, FUN=coevo_sim, net_list=emp_ok, matrix_type="incidence", phi=0.2, alpha=0.2, t_max=1000, eps=0.0001, nsim=100))


save(file="uni_results.RData", uni_results)
save(file="bip1_results.RData", bip1_results)
save(file="bip2_results.RData", bip2_results)
save(file="emp_results.RData", emp_results)