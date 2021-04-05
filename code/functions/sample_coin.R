
#------------------------------------------------------------------------#
#----------------------------- sample_coin() -----------------------------
#------------------------------------------------------------------------#

# Author: Kate
# Date: 02/2021
# Function which returns 1 or 0 given a probability p.
# Input: a vector of probabilities p (p = probability of getting 1)
# Output: a resulting vector of 1s and 0s 

#------------------------------------------------------------------------

sample.coin <- function(p)
{
  res <- rep(NA, length(p))
  for (i in 1:length(p)) {res[i] <- sample(c(1,0), 1, prob = c(p[i], 1 - p[i]))}
  return(res)
}