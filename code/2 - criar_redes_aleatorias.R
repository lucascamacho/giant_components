# Primeiro codigo para gerar redes aleatorias
# carregando pacotes
library(igraph)

#definir numero de especies e numero de redes aleatorias que queremos
n_net <- 100
n_sp <- 100

#criar vetores e listas vazias
nets <- list()
probs <- vector()
tam <- vector()
num <- vector()


for(i in 1:n_net){
  prob <- runif(1, 0, (log(n_sp)/n_sp)) #probabilidade de conexão

  M <- matrix(prob, nrow = n_sp, ncol = n_sp) # criando matriz com diagonal zerada
  diag(M) <- 0

  P <- matrix(runif(n_sp*n_sp, 0, 1), nrow = n_sp, ncol = n_sp) # criar matriz de probabilidades
  M[upper.tri(M)] <- 0 # usaremos somente o triangulo inferior da matriz

  M[M > P] <- 1 #criando interacoes baseado na matriz de probabilidades 
  M[M < P] <- 0

  M <- M + t(M) # trque para deixar a matriz simetrica e a rede nao direcionada
  M[M > 1] <- 1

  nets[[i]] <- M # guardar a rede e a probabilidade de conexao em vetores e listas
  probs[i] <- prob  
  
  tam[i] <- max(components(graph_from_adjacency_matrix(M))$csize) # obter tamanho do maior componente
  num[i] <- components(graph_from_adjacency_matrix(M))$no # obter numero de componentes
  
}

#plotar alguns resultados
plot(probs, tam)
abline(v = 1/n_sp, col = "red")
abline(v = log(n_sp)/n_sp, col = "green")




