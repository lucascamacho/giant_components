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

networks_files <- list.files(here("data/uni_matrices"), pattern="*.txt") #Network names
networks_path <- list.files(here("data/uni_matrices"), pattern="*.txt", full.names=TRUE) #Networks file paths
networks_list <- lapply(networks_path, read.table, header=FALSE) #Reading networks
networks_list <- lapply(networks_list, as.matrix) #Converting all networks to matrices
names(networks_list)<-networks_files #Setting networks names

verif<-lapply(networks_list, function(x){return(sum(rowSums(x)))}) #Verifying networks - excluding network 74 (no connection)

networks_list_ok<-networks_list[c(1:73, 75:100)]

#------------------------------#
#-Setting parallel simulations-#
#------------------------------#

sfInit(parallel = TRUE, cpus = detectCores(), type = "SOCK") #Initializing parallel clusters
sfExportAll() #Export all objects in Global Environment - include user-defined functions
sfLibrary(igraph) #Loading libraries inside clusters

m_values<-seq(from=0.1, to=0.9, by=0.1) #Setting up m values to use in simulations, from 0.1 to 0.9

results<-lapply(m_values, FUN=coevo_sim, net_list=networks_list_ok, matrix_type="adjacency", phi=0.2, alpha=0.2, t_max=1000, eps=0.0001, nsim=100)

results_df<-do.call(rbind, results) #Getting full data frame of results
