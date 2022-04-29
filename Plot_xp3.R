Nbit_comp <- 50
#bon model

cbp1 <- c("#999999","#56B4E9","#009999","#0072B2","#CC79A7","#990033","#FFDB6D", "#C4961A", "#F4EDCA", 
          "#D16103","#660033","#00CC00")


load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARz_n100.RData")
ARI_MCAR_probit <- ARI_MCAR
ARI_MNARz_probit <- ARI_MNARz


load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp3a/xp3_logit.RData")
ARI_MCAR_logit <- ARI_MCAR
ARI_MNARz_logit <- ARI_MNARz


load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp3a/xp3_Laplace.RData")
ARI_MCAR_Laplace <- ARI_MCAR
ARI_MNARz_Laplace <- ARI_MNARz



df_plot_main = data.frame(ARI=c(ARI_MCAR_probit,ARI_MCAR_logit,ARI_MCAR_Laplace,ARI_MNARz_probit,ARI_MNARz_logit,ARI_MNARz_Laplace))
df_plot_main['Mechanism']=c(rep('MCAR',Nbit_comp*3),rep('MNARz',Nbit_comp*3))
df_plot_main['Link function']=c(rep(c(rep("probit",Nbit_comp),rep("logit",Nbit_comp),rep("Laplace",Nbit_comp)),2))

plot_main <- ggplot(df_plot_main, aes(x=`Mechanism`, y=ARI)) + geom_hline(yintercept=0.9,color="red",linetype="dashed")+ geom_boxplot(aes(fill=Mechanism))
plot_main <- plot_main + facet_grid(. ~ `Link function` ) + scale_fill_manual(values=cbp1[c(1,2)]) + xlab('')
