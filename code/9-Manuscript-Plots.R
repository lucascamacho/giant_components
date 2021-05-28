
library(data.table)
library(here)
library(ggpubr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(stringr)
library(tidyverse)

source("./code/functions/criteria.R")
source("./code/functions/clean_simdata.R")

#-------------------------------------------------------#
#-Loading simulation results and network structure data-#
#-------------------------------------------------------#

# Loading network structure data
uni_structure <- data.table(read.table(file = here("data", "Uni-100_NetStructure.txt"), header = TRUE))
bip_structure <- data.table(read.table(file = here("data", "Bip-100_NetStructure.txt"), header = TRUE))
bipA_structure <- data.table(read.table(file = here("data", "Bip-A-100_NetStructure.txt"), header = TRUE))
dim(uni_structure); dim(bip_structure); dim(bipA_structure)

# Loading simulation results data
load(here("data", "simdata_Uni.RData")) # unipartite with fix degree
uni_data <- results_df
load(here("data", "simdata_Bip.RData")) # bipartite with fix degree
bip_data <- results; bip_data <- do.call(rbind, bip_data)
load(here("data", "simdata_BipA.RData")) # bipartite with asymmetric degree
bipA_data <- results_df
dim(uni_data); dim(bip_data); dim(bipA_data)

# Cleaning simulation results data
uni_data <- clean_simdata(uni_data, uni_structure)
bip_data <- clean_simdata(bip_data, bip_structure)
bipA_data <- clean_simdata(bipA_data, bipA_structure)
dim(uni_data); dim(bip_data); dim(bipA_data)

# rbinding simdata for all network types
all(colnames(uni_data) == colnames(bip_data) & colnames(bip_data) == colnames(bipA_data))
colnames(bipA_data)[15] <- "prob" # changing from exponent to prop
all(colnames(uni_data) == colnames(bip_data) & colnames(bip_data) == colnames(bipA_data))
all_df <- rbind(uni_data, bip_data, bipA_data)

# Adding network type and size
all_df$type <- c(rep("Uni", nrow(uni_data)), rep("Bip", nrow(bip_data)), rep("BipA", nrow(bipA_data)))
all_df$size <- 100

# Calculating connectivity parameter for all networks
uni_conn_par <- unique(uni_data[, c("id", "degree")])
uni_conn_par <- data.frame(id = uni_conn_par$id, conn_par = uni_conn_par$degree)
range(uni_conn_par$conn_par)

net_names <- dir(here("data/bip_matrices"))
net_list <- lapply(here("data/bip_matrices", net_names), read.table, header = FALSE, sep = "\t")
conn_par <- unlist(lapply(net_list, criteria))
bip_conn_par <- data.frame(id = str_sub(net_names, end = -5), conn_par = conn_par)
range(bip_conn_par$conn_par)

net_names <- dir(here("data/bip_matrices_asymmetric"))
net_list <- lapply(here("data/bip_matrices_asymmetric", net_names), read.table, header = FALSE, sep = "\t")
conn_par <- unlist(lapply(net_list, criteria))
bipA_conn_par <- data.frame(id = str_sub(net_names, end = -5), conn_par = conn_par)
range(bipA_conn_par$conn_par)

all_conn_par <- rbind(uni_conn_par, bip_conn_par, bipA_conn_par)

# Adding connectivity parameter for all networks
all_df <- merge(x = all_df, y = all_conn_par, by = "id", all.x = TRUE)

head(all_df)
dim(all_df)

# ----- Unipartite --------------------------------------------------

pal2 <- brewer.pal(9, "YlOrRd")
uni_df <- all_df %>% filter(type == "Uni")

p1 <- ggplot(uni_df) +
  geom_point(aes(x = degree, y = sp_indirect_effects, fill = m),
             color = "white", shape = 21, size = 3, alpha = 0.9) +
  scale_fill_manual(values = pal2) +
  xlab("<k>") + ylab("Contribution of indirect effects") +
  geom_vline(xintercept = 0, linetype="dashed") +
  theme_pubr()

p2 <- ggplot(uni_df) +
  geom_point(aes(x = degree, y = network_tm, fill = m), 
             color = "white", shape = 21, size = 3, alpha = 0.9) +
  scale_fill_manual(values = pal2) +
  xlab("<k>") + ylab("Network trait matching") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_pubr()

p3 <- ggplot(uni_df) +
  geom_point(aes(x = sp_indirect_effects, y = network_tm, fill = m), 
             color = "white", shape = 21, size = 3, alpha = 0.9) +
  scale_fill_manual(values = pal2) +
  xlab("Indirect effects") + 
  ylab("Network trait matching") +
  theme_pubr()

p4 <- ggplot(uni_df) +
  geom_point(aes(x = degree, y = component_tm, fill = m), 
             color = "white", shape = 21, size = 3, alpha = 0.9) +
  scale_fill_manual(values = pal2) +
  xlab("<k>") + ylab("Component trait matching") +
  theme_pubr()

p1 + p2 + p3 + p4 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

# ----- Bipartite --------------------------------------------------

pal2 <- brewer.pal(9, "YlGnBu")
bip_df <- all_df %>% filter(type == "Bip")

p1 <- ggplot(bip_df) +
  geom_point(aes(x = degree, y = sp_indirect_effects, fill = m),
             color = "white", shape = 21, size = 3, alpha = 0.9) +
  scale_fill_manual(values = pal2) +
  xlab("<k>") + ylab("Contribution of indirect effects") +
  geom_vline(xintercept = 1, linetype="dashed") +
  theme_pubr()

p2 <- ggplot(bip_df) +
  geom_point(aes(x = degree, y = network_tm, fill = m), 
             color = "white", shape = 21, size = 3, alpha = 0.9) +
  scale_fill_manual(values = pal2) +
  xlab("<k>") + 
  ylab("Network trait matching") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  theme_pubr()

p3 <- ggplot(bip_df) +
  geom_point(aes(x = sp_indirect_effects, y = network_tm, fill = m), 
             color = "white", shape = 21, size = 3, alpha = 0.9) +
  scale_fill_manual(values = pal2) +
  xlab("Indirect effects") + 
  ylab("Network trait matching") +
  theme_pubr()

p4 <- ggplot(bip_df) +
  geom_point(aes(x = degree, y = component_tm, fill = m), 
             color = "white", shape = 21, size = 3, alpha = 0.9) +
  scale_fill_manual(values = pal2) +
  xlab("<k>") + 
  ylab("Component trait matching") +
  theme_pubr()

p1 + p2 + p3 + p4 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

# ----- Bipartite_A ------------------------------------------------

pal2 <- brewer.pal(9, "YlGn")
bipA_df <- all_df %>% filter(type == "BipA")

p1 <- ggplot(bipA_df) +
  geom_point(aes(x = degree, y = sp_indirect_effects, fill = m),
             color = "white", shape = 21, size = 3, alpha = 0.9) +
  scale_fill_manual(values = pal2) +
  xlab("<k>") + ylab("Contribution of indirect effects") +
  geom_vline(xintercept = 1, linetype="dashed") +
  theme_pubr()

p2 <- ggplot(bipA_df) +
  geom_point(aes(x = degree, y = network_tm, fill = m), 
             color = "white", shape = 21, size = 3, alpha = 0.9) +
  scale_fill_manual(values = pal2) +
  xlab("<k>") + 
  ylab("Network trait matching") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  theme_pubr()

p3 <- ggplot(bipA_df) +
  geom_point(aes(x = sp_indirect_effects, y = network_tm, fill = m), 
             color = "white", shape = 21, size = 3, alpha = 0.9) +
  scale_fill_manual(values = pal2) +
  xlab("Indirect effects") + 
  ylab("Network trait matching") +
  theme_pubr()

p4 <- ggplot(bipA_df) +
  geom_point(aes(x = degree, y = component_tm, fill = m), 
             color = "white", shape = 21, size = 3, alpha = 0.9) +
  scale_fill_manual(values = pal2) +
  xlab("<k>") + 
  ylab("Component trait matching") +
  theme_pubr()

p1 + p2 + p3 + p4 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
