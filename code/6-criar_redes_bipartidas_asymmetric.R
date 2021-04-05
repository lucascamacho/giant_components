#------------------------------------------------------------------------#
#---- Creates bipartite networks with non-uniform interaction probability
#------------------------------------------------------------------------#

# Author: Kate
# Date: 02/2021
# Interaction probability at two different levels: between networks, and 
# across species in networks.
# Two versions attempted: log-normal (did not generate the structural 
# variability needed), and power-law: the final version).
# Power-law - L30
# Lognormal - L75

#------------------------------------------------------------------------#

library(igraph)
source("./code/functions/sample_coin.R")

#definir numero de especies e numero de redes aleatorias que queremos
n_net <- 100
n_plant <- 100
n_pol <- 100

#criar vetores e listas vazias
nets <- list()
probs <- vector()
tam <- vector()
num <- vector()
var <- vector()

#------------------------------- POWER-LAW ------------------------------#

for (i in 1:n_net) {
  
  x <- round(runif(100, 1, 100), digits = 0)
  e <- runif(n = 1, min = -2.5, max = -0.2)
  
  k_plant <- k_pol <- x^e # creates a degree dist for both sp sets with exponent e
  
  p_plant <- p_pol <- k_plant/max(k_plant) # turns degree dist into prob dist
  
  P <- outer(p_plant, p_pol, "*") # probability matrix
  
  links <- sample.coin(P) # samples interactions using P
  
  M <- matrix(links, nrow = n_plant, ncol = n_pol) # interaction incidence matrix
  
  nets[[i]] <- M
  probs[i] <- e
  var[i] <- var(c(rowSums(M), colSums(M))) # degree variance
  
  tam[i] <- max(igraph::components(graph_from_incidence_matrix(M))$csize) # tamanho do maior componente
  num[i] <- igraph::components(graph_from_incidence_matrix(M))$no # numero de componentes
  
  # write.table(M, paste("./Data/bip_matrices_asymmetric/BIP-A-100-", i, ".txt", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE)
}

result_df <- data.frame(id = paste("BIP-A-100-", 1:100, sep = ""), exponent = probs, ncomp = num, sizelargcomp = tam)
# write.table(result_df, "./Data/Bip-A-100_NetStructure.txt", sep = "\t")

range(num)
range(tam)

plot(probs, num, xlab = "Power-law exponent", ylab = "Number of network's largest component")

plot(probs, tam, xlab = "Power-law exponent", ylab = "Size of network's largest component")

plot(probs, var, xlab = "Power-law exponent", ylab = "Variance in species' degree")

do.call(rbind, lapply(nets, function(x){c(sum(rowSums(x) == 0), sum(colSums(x) == 0))}))

#------------------------------ LOG-NORMAL ------------------------------#
# Attempt of using a lognormal degree distribution which did not result in the require
# variation in number and size of network components.

for (i in 1:n_net) {
  
  mean <- runif(n = 1, min = 0, max = 10)
  sd <- 1
  k_plant <- round(rlnorm(n_plant, meanlog = mean, sdlog = sd))
  k_pol <- round(rlnorm(n_pol, meanlog = mean, sdlog = sd))
  
  p_plant <- k_plant/max(k_plant)
  p_pol <- k_pol/max(k_pol)
  
  P <- outer(p_plant, p_pol, "*")
  
  links <- sample.coin(P)
  
  M <- matrix(links, nrow = n_plant, ncol = n_pol) # crindo matriz M
  
  nets[[i]] <- M # guardar a rede e a probabilidade de conexao em vetores e listas
  probs[i] <- mean
  
  tam[i] <- max(igraph::components(graph_from_incidence_matrix(M))$csize) # tamanho do maior componente
  num[i] <- igraph::components(graph_from_incidence_matrix(M))$no # numero de componentes
}

#plotar alguns resultados
plot(probs, tam)
abline(v = 1/(n_plant + n_pol), col = "red")
abline(v = log(n_plant + n_pol)/(n_plant + n_pol), col = "green")
range(num); range(tam)

# k/max(k)
# mean 0:1; sd = 0.1: 1 component with 200 species
# mean 0:1; sd = 1: 27:166 components with 15:174 species
# mean 0:1; sd = 10: 192:199 components with 2:88 species
# mean 0:1; sd = 100: 198:199 components with 2:3 species
# mean 0:10; sd = 0.5: 1:16 components with 183:200 species
# mean 0:10; sd = 1: 24:164 components with 25:177 species # BEST but still no emergence
# mean 0:10; sd = 5: 191:199 components with 2:10 species
# mean 0:10; sd = 10: 194:199 components with 2:7 species
# mean 0:10; sd = 100: 197:199 components with 2:4 species
# mean 0:100; sd = 0.5: 1:13 components with 188:200 species
# mean 0:100; sd = 5: 188:199 components with 2:11 species
# mean 0:100; sd = 100: 197:199 components with 2:4 species

# k/sum(k)
# mean 0:1; sd = 0.1: all fragmented 196:200 components with 1:3 species
# mean 0:1; sd = 1: all fragmented 195:200 components with 1:3 species
# mean 0:1; sd = 10: all fragmented 196:200 components with 1:5 species
# mean 0:1; sd = 100: all fragmented 198:200 components with 1:3 species
# mean 0:10; sd = 0.5: all fragmented 195:200 components with 1:3 species
# mean 0:10; sd = 1: all fragmented 194:200 components with 1:3 species
# mean 0:10; sd = 5: all fragmented 196:200 components with 1:5 species
# mean 0:10; sd = 10: all fragmented 197:200 components with 1:4 species
# mean 0:10; sd = 100: all fragmented 198:200 components with 1:3 species
# mean 0:100; sd = 5: all fragmented 197:200 components with 1:4 species
# mean 0:100; sd = 0.5: all fragmented 195:200 components with 1:3 species
