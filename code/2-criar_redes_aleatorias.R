
# Primeiro codigo para gerar redes aleatorias
# carregando pacotes
library(igraph)

#---------- Code by Camacho, Leandro e Dani: for a function version see L45 ----------

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
plot(num, tam)
abline(v = 1/n_sp, col = "red")
abline(v = log(n_sp)/n_sp, col = "green")

#------------------------------ Function by Kate e Dani ------------------------------

# function that returns one random unipartite binary network & its probability of connection
uni_net <- function(n_sp){
  
  prob <- runif(1, 0, (log(n_sp)/n_sp)) #probabilidade de conexao
  
  M <- matrix(prob, nrow = n_sp, ncol = n_sp) #cria matriz com diagonal zerada
  diag(M) <- 0
  
  P <- matrix(runif(n_sp*n_sp, 0, 1), nrow = n_sp, ncol = n_sp) #cria matriz de probabilidades
  
  M[upper.tri(M)] <- 0 #usamos somente o triangulo inferior da matriz
  M[M >= P] <- 1 #cria interacoes baseado na matriz de probabilidades 
  M[M < P] <- 0
  
  M <- M + t(M) #truque para deixar a matriz simetrica e a rede nao direcionada
  if (any(M > 1)) {print("ERROR: M higher than 1")}
  
  return <- list(prob, M) # retorna a rede e a probabilidade de conexao
  return(return)  
}

# running function uni_net to create matrices and a network dataframe
n_net <- 100 # number of networks
net_df <- data.frame(id = character(0), prob = numeric(0), ncomp = numeric(0), 
                     sizelargcomp = numeric(0))

for (i in 1:n_net){
  
  n_sp <- 100
  M <- uni_net(n_sp) # creates network
  id <- paste("UNI-100_", i, sep = "")
  comp <- components(graph_from_adjacency_matrix(M[[2]]))
  
  df <- data.frame(id = id, prob = M[[1]], ncomp = comp$no, sizelargcomp = max(comp$csize))
  net_df <- rbind(net_df, df)
  
  #write.table(M[[2]], paste("./data/uni_matrices/", id, ".txt", sep = ""), sep = "\t", row.names = F, col.names = F)
}

#write.table(net_df, "./data/Uni-100_NetStructure.txt", sep = "\t")
