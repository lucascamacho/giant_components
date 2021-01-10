library(data.table)
library(here)
library(ggpubr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

#-------------------------------------------------------#
#-Loading simulation results and network structure data-#
#-------------------------------------------------------#

net_structure<-data.table(read.table(file=here("data", "Uni-100_NetStructure.txt"), header=TRUE))
results_df<-data.table(results_df)
component_df<-results_df[,setdiff(names(results_df), "sp_id"), with=FALSE]
component_df<-component_df[, lapply(.SD, mean, na.rm=TRUE), 
                           by=.(network_id, m)]
names(component_df)[1]<-"id"
component_df$id<-gsub(".txt.*", "", component_df$id)

all_df<-merge(x = component_df, y = net_structure, by = "id", all.x = TRUE)
all_df$m<-as.factor(all_df$m)

pal<-brewer.pal(15, "PuRd")
pal2<-c("palegreen2", "palegreen3", "palegreen4", "coral1", "coral3", "indianred4", pal[7:9])

p1<-ggplot(all_df)+
  geom_point(aes(x=degree, y=network_tm, fill=m),color="white", shape=21, size=3, alpha=0.9) +
  scale_fill_manual(values=pal2)+
  xlab("<k>")+ylab("Network trait matching")+
  geom_vline(xintercept=1, linetype="dashed")+
  theme_pubr()

p2<-ggplot(all_df)+
  geom_point(aes(x=degree, y=sp_indirect_effects, fill=m),color="white", shape=21, size=3, alpha=0.9) +
  scale_fill_manual(values=pal2)+
  xlab("<k>")+ylab("Contribution of indirect effects")+
  geom_vline(xintercept=1, linetype="dashed")+
  theme_pubr()

p1+p2+plot_layout(guides="collect") & theme(legend.position="bottom")

p3<-ggscatter(all_df, x="sp_indirect_effects", y="network_tm", fill="m", color="white", shape=21, size=2.5)

p3<-ggpar(p3, xlab="Contribution of indirect effects", ylab="Network trait matching")+
  scale_fill_manual(values=pal2)

p3
