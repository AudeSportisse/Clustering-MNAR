Nbit_comp <- 50
#bon model

cbp1 <- c("#999999","#56B4E9","#009999","#0072B2","#CC79A7","#990033","#FFDB6D", "#C4961A", "#F4EDCA", 
          "#D16103","#660033","#00CC00")


load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARz_n100.RData")
ARI_MCAR_rho0 <- ARI_MCAR
ARI_MNARz_rho0 <- ARI_MNARz


load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp3b/xp3b_rho01.RData")
ARI_MCAR_rho01 <- ARI_MCAR
ARI_MNARz_rho01 <- ARI_MNARz


load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp3b/xp3b_rho025.RData")
ARI_MCAR_rho025 <- ARI_MCAR
ARI_MNARz_rho025 <- ARI_MNARz

load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp3b/xp3b_rho05.RData")
ARI_MCAR_rho05 <- ARI_MCAR
ARI_MNARz_rho05 <- ARI_MNARz

df_plot_main = data.frame(ARI=c(ARI_MCAR_rho0,ARI_MCAR_rho01,ARI_MCAR_rho025,ARI_MCAR_rho05,ARI_MNARz_rho0,ARI_MNARz_rho01,ARI_MNARz_rho025,ARI_MNARz_rho05))
df_plot_main['Mechanism']=c(rep('MCAR',Nbit_comp*4),rep('MNARz',Nbit_comp*4))
df_plot_main['Rho']=c(rep(c(rep("l=0",Nbit_comp),rep("l=0.1",Nbit_comp),rep("l=0.25",Nbit_comp),rep("l=0.5",Nbit_comp)),2))

plot_main <- ggplot(df_plot_main, aes(x=`Mechanism`, y=ARI)) + geom_hline(yintercept=0.9,color="red",linetype="dashed")+ geom_boxplot(aes(fill=Mechanism))
plot_main <- plot_main + facet_grid(. ~ `Rho` ) + scale_fill_manual(values=cbp1[c(1,2)])+xlab("")
