
library(data.table)
library(here)
library(ggpubr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(stringr)
library(tidyverse)
source("./code/functions/criteria.R")

#-------------------------------------------------------#
#-Loading simulation results and network structure data-#
#-------------------------------------------------------#

load(here("data", "dadosBipA.RData"))

net_structure <- data.table(read.table(file = here("data", "Bip-A-100_NetStructure.txt"), header=TRUE))
results_df <- data.table(results_df)
component_df <- results_df[, setdiff(names(results_df), "sp_id"), with = FALSE]
component_df <- component_df[, lapply(.SD, mean, na.rm = TRUE), by = .(network_id, m)]
names(component_df)[1] <- "id"
component_df$id <- gsub(".txt.*", "", component_df$id)

all_df <- merge(x = component_df, y = net_structure, by = "id", all.x = TRUE)
all_df$m <- as.factor(all_df$m)

# adding connectivity parameter for bipartite networks
net_names <- dir(here("data/bip_matrices_asymmetric"))
net_list <- lapply(here("data/bip_matrices_asymmetric", net_names), read.table, header = FALSE, sep = "\t")
criteria <- unlist(lapply(net_list, criteria))
criteria_df <- data.frame(id = str_sub(net_names, end = -5), criteria = criteria)
all_df <- merge(x = all_df, y = criteria_df, by = "id", all.x = TRUE)

criteria_df[order(criteria_df$criteria, decreasing = TRUE),]

unique(all_df[all_df$criteria > 10000, "id"])

all_df <- all_df %>% filter(criteria < 10000)

pal <- brewer.pal(9, "PuRd")
pal2 <- c("palegreen2", "palegreen3", "palegreen4", "coral1", "coral3", "indianred4", pal[7:9])

all_df %>% mutate(prop_large = sizelargcomp/200) %>% 
ggplot() +
  geom_point(aes(x = criteria, y = prop_large),
             color = "white", shape = 21, size = 3, alpha = 0.9) +
  scale_fill_manual(values = pal2) +
  xlab("c") + ylab("Proportion sp in largest component") +
  geom_vline(xintercept=1, linetype="dashed") +
  theme_pubr()

p1 <- ggplot(all_df) +
  geom_point(aes(x = criteria, y = sp_indirect_effects, fill = m),
             color = "white", shape = 21, size = 3, alpha = 0.9) +
  scale_fill_manual(values = pal2) +
  xlab("c") + ylab("Contribution of indirect effects") +
  geom_vline(xintercept=1, linetype="dashed") +
  theme_pubr()

p2 <- ggplot(all_df) +
  geom_point(aes(x = criteria, y = network_tm, fill = m), 
             color = "white", shape = 21, size = 3, alpha = 0.9) +
  scale_fill_manual(values = pal2) +
  xlab("c") + 
  ylab("Network trait matching") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  theme_pubr()

p3 <- ggplot(all_df) +
  geom_point(aes(x = sp_indirect_effects, y = network_tm, fill = m), 
             color = "white", shape = 21, size = 3, alpha = 0.9) +
  scale_fill_manual(values = pal2) +
  xlab("Indirect effects") + 
  ylab("Network trait matching") +
  theme_pubr()

p4 <- ggplot(all_df) +
  geom_point(aes(x = criteria, y = component_tm, fill = m), 
             color = "white", shape = 21, size = 3, alpha = 0.9) +
  scale_fill_manual(values = pal2) +
  xlab("c") + 
  ylab("Component trait matching") +
  theme_pubr()

p1 + p2 + p3 + p4 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
