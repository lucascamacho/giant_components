
#------------------------------------------------------------------------#
#------------------------------ criteria() ------------------------------
#------------------------------------------------------------------------#

# Author: Kate
# Date: 06/2020
# Function that calculates "criteria", i.e. the connectivity parameter for 
# bipartite which marks the abrupt transition between fragmented and highly 
# connected networks (presence of GC).
# input: M (interaction INCIDENCE matrix)
# output: c (criteria) 

#------------------------------------------------------------------------

# Os subscritos i/j nao sao sp, mas sim elementos dos vetores que vao do 
# menor ao maior grau dos sets I/J (I = rows, J = columns). 

criteria <- function(M) {
  
  M[M > 0] <- 1 # turns potentially quantitative into binary matrix 
  
  # vector with degrees per species set
  VkI <- unique(rowSums(M)[order(rowSums(M))]) # grau das linhas
  VkJ <- unique(colSums(M)[order(colSums(M))]) # grau das colunas
  
  # frequency of species with degrees in VkI/VkJ
  pkI <- table(rowSums(M))/sum(table(rowSums(M))) # frequencia relativa dos graus das linhas
  pkJ <- table(colSums(M))/sum(table(colSums(M))) # frequencia relativa dos graus das colunas
  
  length(VkI) == length(pkI)
  length(VkJ) == length(pkJ)
  
  # multiplies pkI/pkJ vectors to create PIJ matrix: freq of int between sp with ki and kj 
  pkI <- matrix(pkI, ncol = 1)
  pkJ <- matrix(pkJ, nrow = 1)
  
  PIJ <- pkI %*% pkJ
  sum(PIJ) == 1
  
  # turns degree vectors into degree matrices per set
  KI <- matrix(rep(VkI, length(VkJ)), nrow = length(VkI), byrow = FALSE)
  KJ <- matrix(rep(VkJ, length(VkI)), ncol = length(VkJ), byrow = TRUE)
  dim(KI) == dim(KJ)
  
  c <- sum(KI*KJ*(KI*KJ - KI - KJ)*PIJ)
  return(c)
}
