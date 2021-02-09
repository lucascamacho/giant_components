library(igraph)

#definir numero de especies e numero de redes aleatorias que queremos
n_net <- 200
n_plant <- 100
n_pol <- 100


#criar vetores e listas vazias
nets <- list()
probs <- vector()
tam <- vector()
num <- vector()

for(i in 1:n_net){
 
  #prob <- runif(1, 0, (log(n_plant + n_pol)/(n_plant + n_pol))) #probabilidade de conexão
  prob <- runif(1, 0, 0.05) #probabilidade de conexão
  
  M <- matrix(prob, nrow = n_plant, ncol = n_pol) # crindo matriz M
  
  P <- matrix(runif(n_plant*n_pol), nrow = n_plant, ncol = n_pol) # criando matriz de probabilidades
  
  M[M >= P] <- 1
  M[M < P] <- 0
  
  nets[[i]] <- M # guardar a rede e a probabilidade de conexao em vetores e listas
  probs[i] <- prob
  
  tam[i] <- max(components(graph_from_incidence_matrix(M))$csize) # obter tamanho do maior componente
  num[i] <- components(graph_from_incidence_matrix(M))$no # obter numero de componentes
  
}

#plotar alguns resultados
plot(probs, tam)
abline(v = 1/(n_plant + n_pol), col = "red")
abline(v = log(n_plant + n_pol)/(n_plant + n_pol), col = "green")
  