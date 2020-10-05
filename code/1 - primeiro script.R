# primeiro script para teste

# # salvar dados quaisquer
vec = c(1,2,3,4)

save(vec, file = "dados.RData")

pdf("1 - primeiro_plot.pdf")
plot(vec)
dev.off()
#teste commit leandro - teste 2 sem pull request