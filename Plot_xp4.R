Nbit_comp <- 50
#bon model

cbp1 <- c("#999999","#56B4E9","#009999","#0072B2","#CC79A7","#990033","#FFDB6D", "#C4961A", "#F4EDCA", 
          "#D16103","#660033","#00CC00")


load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp4b/xp4b_30perc.RData")
#load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp4a/xp4_30perc_ARI15_bis.RData")
ARI_MCAR_30perc <- ARI_MCAR[,3]
ARI_MNARz_30perc <- ARI_MNARz[,3]


load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp4a/xp4_10perc.RData")
#load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp4a/xp4_10perc_ARI15.RData")
ARI_MCAR_10perc <- ARI_MCAR
ARI_MNARz_10perc <- ARI_MNARz


load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp4b/xp4b_50perc.RData")
#load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp4a/xp4_50perc_ARI15_bis.RData")
ARI_MCAR_50perc <- ARI_MCAR[,3]
ARI_MNARz_50perc <- ARI_MNARz[,3]



df_plot_main = data.frame(ARI=c(ARI_MCAR_10perc,ARI_MCAR_30perc,ARI_MCAR_50perc,ARI_MNARz_10perc,ARI_MNARz_30perc,ARI_MNARz_50perc))
df_plot_main['Mechanism']=c(rep('MCAR',Nbit_comp*3),rep('MNARz',Nbit_comp*3))
df_plot_main['Percentage NA']=c(rep(c(rep("10% NA",Nbit_comp),rep("30% NA",Nbit_comp),rep("50% NA",Nbit_comp)),2))

plot_main <- ggplot(df_plot_main, aes(x=`Mechanism`, y=ARI)) + geom_hline(yintercept=0.9,color="red",linetype="dashed") + geom_boxplot(aes(fill=Mechanism))
plot_main <- plot_main + facet_grid(. ~ `Percentage NA` ) + scale_fill_manual(values=cbp1[c(1,2)]) + xlab("")


#sum(apply(ICL_MCAR,1,which.max)==3)/nb_it
#sum(apply(ICL_MNARz,1,which.max)==3)/nb_it