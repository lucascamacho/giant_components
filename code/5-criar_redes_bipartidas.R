library(igraph)
library(here)

#definir numero de especies e numero de redes aleatorias que queremos
n_net <- 100
n_plant <- 50
n_pol <- 50


#criar vetores e listas vazias
nets <- list()
probs <- vector()
tam <- vector()
num <- vector()

for(i in 1:n_net){
 
  #prob <- runif(1, 0, (log(n_plant + n_pol)/(n_plant + n_pol))) #probabilidade de conexão
  prob <- runif(1, 0, 0.06) #probabilidade de conexão
  
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
plot(probs, num)
abline(v = 1/(n_plant + n_pol), col = "red")
abline(v = log(n_plant + n_pol)/(n_plant + n_pol), col = "green")

#salvando atributos das redes
bip_names<-paste("BIP", paste(100, 1:100, sep="_"), sep="-")
net_str<-data.frame(id=bip_names, prob=probs, ncomp=num, sizelargcomp=tam)

write.table(net_str, here("data", "Bip-100_NetStructure.txt"), sep="\t")

#salvando redes

files_names=paste0(bip_names, ".txt")

for (i in 1:length(nets)){
  
  write.table(nets[[i]], here("data","bip_matrices", files_names[i]), sep="\t", col.names=FALSE, row.names=FALSE)
  
}
