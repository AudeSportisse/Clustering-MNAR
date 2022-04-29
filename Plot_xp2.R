library(ggplot2)

Nbit_comp <- 50
#bon model

cbp1 <- c("#999999","#56B4E9","#009999","#0072B2","#CC79A7","#990033","#FFDB6D", "#C4961A", "#F4EDCA", 
          "#D16103","#660033","#00CC00")




###############
##MNARzj
###############

load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp2/xp2_MNARzj_3.RData")

ARI_MNARzj_3 <- ARI_MNAR_true
ARI_MCAR_3_zj <- ARI_MCAR
ARI_MNARz_3_zj <- ARI_MNARz
ARI_Mice_3_zj <- ARI_Mice

load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp2/xp2_MNARzj_6.RData")

ARI_MNARzj_6 <- ARI_MNAR_true
ARI_MCAR_6_zj <- ARI_MCAR
ARI_MNARz_6_zj <- ARI_MNARz
ARI_Mice_6_zj <- ARI_Mice

load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp2/xp2_MNARzj_9.RData")

ARI_MNARzj_9 <- ARI_MNAR_true
ARI_MCAR_9_zj <- ARI_MCAR
ARI_MNARz_9_zj <- ARI_MNARz
ARI_Mice_9_zj <- ARI_Mice


df_plot1 = data.frame(ARI=c(ARI_MCAR_3_zj,ARI_MCAR_6_zj,ARI_MCAR_9_zj,ARI_MNARzj_3,ARI_MNARzj_6,ARI_MNARzj_9,ARI_MNARz_3_zj,ARI_MNARz_6_zj,ARI_MNARz_9_zj,ARI_Mice_3_zj,ARI_Mice_6_zj,ARI_Mice_9_zj))
df_plot1['Mechanism']=c(rep('MCAR',Nbit_comp*3),rep('MNARzj',Nbit_comp*3),rep('MNARz',Nbit_comp*3),rep('Mice',Nbit_comp*3))
df_plot1['Number of variables']=c(rep(c(rep('3',Nbit_comp),rep('6',Nbit_comp),rep('9',Nbit_comp)),4))
plotzj <- ggplot(data=df_plot1,aes(x=`Number of variables`, y=ARI,fill=Mechanism))+geom_boxplot()
plotzj = plotzj #+ theme(legend.position = "none")
plotzj = plotzj + ggtitle('True model is MNARzj') + scale_fill_manual(values=cbp1[c(1,4,2,11)])
plotzj




###############
##MNARy
###############

load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp2/xp2_MNARy_3.RData")

ARI_MNARy_3 <- ARI_MNAR_true
ARI_MCAR_3_y <- ARI_MCAR
ARI_MNARz_3_y <- ARI_MNARz
ARI_Mice_3_y <- ARI_Mice

load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp2/xp2_MNARy_6.RData")

ARI_MNARy_6 <- ARI_MNAR_true
ARI_MCAR_6_y <- ARI_MCAR
ARI_MNARz_6_y <- ARI_MNARz
ARI_Mice_6_y <- ARI_Mice

load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp2/xp2_MNARy_9.RData")

ARI_MNARy_9 <- ARI_MNAR_true
ARI_MCAR_9_y <- ARI_MCAR
ARI_MNARz_9_y <- ARI_MNARz
ARI_Mice_9_y <- ARI_Mice


df_plot1 = data.frame(ARI=c(ARI_MCAR_3_y,ARI_MCAR_6_y,ARI_MCAR_9_y,ARI_MNARy_3,ARI_MNARy_6,ARI_MNARy_9,ARI_MNARz_3_y,ARI_MNARz_6_y,ARI_MNARz_9_y,ARI_Mice_3_y,ARI_Mice_6_y,ARI_Mice_9_y))
df_plot1['Mechanism']=c(rep('MCAR',Nbit_comp*3),rep('MNARy',Nbit_comp*3),rep('MNARz',Nbit_comp*3),rep('Mice',Nbit_comp*3))
df_plot1['Number of variables']=c(rep(c(rep('3',Nbit_comp),rep('6',Nbit_comp),rep('9',Nbit_comp)),4))
ploty <- ggplot(data=df_plot1,aes(x=`Number of variables`, y=ARI,fill=Mechanism))+geom_boxplot()
ploty = ploty #+ theme(legend.position = "none")
ploty = ploty + ggtitle('True model is MNARy') + scale_fill_manual(values=cbp1[c(1,5,2,11)])
ploty




###############
##MNARykzj
###############

load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp2/xp2_MNARykzj_3.RData")

ARI_MNARykzj_3 <- ARI_MNAR_true
ARI_MCAR_3_ykzj <- ARI_MCAR
ARI_MNARz_3_ykzj <- ARI_MNARz
ARI_Mice_3_ykzj <- ARI_Mice

load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp2/xp2_MNARykzj_6.RData")

ARI_MNARykzj_6 <- ARI_MNAR_true
ARI_MCAR_6_ykzj <- ARI_MCAR
ARI_MNARz_6_ykzj <- ARI_MNARz
ARI_Mice_6_ykzj <- ARI_Mice

load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp2/xp2_MNARykzj_9.RData")

ARI_MNARykzj_9 <- ARI_MNAR_true
ARI_MCAR_9_ykzj <- ARI_MCAR
ARI_MNARz_9_ykzj <- ARI_MNARz
ARI_Mice_9_ykzj <- ARI_Mice


df_plot1 = data.frame(ARI=c(ARI_MCAR_3_ykzj,ARI_MCAR_6_ykzj,ARI_MCAR_9_ykzj,ARI_MNARykzj_3,ARI_MNARykzj_6,ARI_MNARykzj_9,ARI_MNARz_3_ykzj,ARI_MNARz_6_ykzj,ARI_MNARz_9_ykzj,ARI_Mice_3_ykzj,ARI_Mice_6_ykzj,ARI_Mice_9_ykzj))
df_plot1['Mechanism']=c(rep('MCAR',Nbit_comp*3),rep('MNARykzj',Nbit_comp*3),rep('MNARz',Nbit_comp*3),rep('Mice',Nbit_comp*3))
df_plot1['Number of variables']=c(rep(c(rep('3',Nbit_comp),rep('6',Nbit_comp),rep('9',Nbit_comp)),4))
plotykzj <- ggplot(data=df_plot1,aes(x=`Number of variables`, y=ARI,fill=Mechanism))+geom_boxplot()
plotykzj = plotykzj #+ theme(legend.position = "none")
plotykzj = plotykzj + ggtitle('True model is MNARykzj') + scale_fill_manual(values=cbp1[c(1,8,2,11)])
plotykzj





###############
##MNARyz
###############

load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp2/xp2_MNARyz_3.RData")

ARI_MNARyz_3 <- ARI_MNAR_true
ARI_MCAR_3_yz <- ARI_MCAR
ARI_MNARz_3_yz <- ARI_MNARz
ARI_Mice_3_yz <- ARI_Mice

load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp2/xp2_MNARyz_6.RData")

ARI_MNARyz_6 <- ARI_MNAR_true
ARI_MCAR_6_yz <- ARI_MCAR
ARI_MNARz_6_yz <- ARI_MNARz
ARI_Mice_6_yz <- ARI_Mice

load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp2/xp2_MNARyz_9.RData")

ARI_MNARyz_9 <- ARI_MNAR_true
ARI_MCAR_9_yz <- ARI_MCAR
ARI_MNARz_9_yz <- ARI_MNARz
ARI_Mice_9_yz <- ARI_Mice


df_plot1 = data.frame(ARI=c(ARI_MCAR_3_yz,ARI_MCAR_6_yz,ARI_MCAR_9_yz,ARI_MNARyz_3,ARI_MNARyz_6,ARI_MNARyz_9,ARI_MNARz_3_yz,ARI_MNARz_6_yz,ARI_MNARz_9_yz,ARI_Mice_3_yz,ARI_Mice_6_yz,ARI_Mice_9_yz))
df_plot1['Mechanism']=c(rep('MCAR',Nbit_comp*3),rep('MNARyz',Nbit_comp*3),rep('MNARz',Nbit_comp*3),rep('Mice',Nbit_comp*3))
df_plot1['Number of variables']=c(rep(c(rep('3',Nbit_comp),rep('6',Nbit_comp),rep('9',Nbit_comp)),4))
plotyz <- ggplot(data=df_plot1,aes(x=`Number of variables`, y=ARI,fill=Mechanism))+geom_boxplot()
plotyz = plotyz #+ theme(legend.position = "none")
plotyz = plotyz + ggtitle('True model is MNARyz') + scale_fill_manual(values=cbp1[c(1,7,2,11)])
plotyz



vec_MNARzj <- c(ARI_MCAR_3_zj,ARI_MCAR_6_zj,ARI_MCAR_9_zj,ARI_MNARzj_3,ARI_MNARzj_6,ARI_MNARzj_9,ARI_MNARz_3_zj,ARI_MNARz_6_zj,ARI_MNARz_9_zj,ARI_Mice_3_zj,ARI_Mice_6_zj,ARI_Mice_9_zj)
vec_MNARy <- c(ARI_MCAR_3_y,ARI_MCAR_6_y,ARI_MCAR_9_y,ARI_MNARy_3,ARI_MNARy_6,ARI_MNARy_9,ARI_MNARz_3_y,ARI_MNARz_6_y,ARI_MNARz_9_y,ARI_Mice_3_y,ARI_Mice_6_y,ARI_Mice_9_y)
vec_MNARyz <- c(ARI_MCAR_3_yz,ARI_MCAR_6_yz,ARI_MCAR_9_yz,ARI_MNARyz_3,ARI_MNARyz_6,ARI_MNARyz_9,ARI_MNARz_3_yz,ARI_MNARz_6_yz,ARI_MNARz_9_yz,ARI_Mice_3_yz,ARI_Mice_6_yz,ARI_Mice_9_yz)
vec_MNARykzj <- c(ARI_MCAR_3_ykzj,ARI_MCAR_6_ykzj,ARI_MCAR_9_ykzj,ARI_MNARykzj_3,ARI_MNARykzj_6,ARI_MNARykzj_9,ARI_MNARz_3_ykzj,ARI_MNARz_6_ykzj,ARI_MNARz_9_ykzj,ARI_Mice_3_ykzj,ARI_Mice_6_ykzj,ARI_Mice_9_ykzj)

df_plot_main = data.frame(ARI=c(vec_MNARzj,vec_MNARy,vec_MNARyz,vec_MNARykzj))
df_plot_main['Mechanism']=rep(c(rep('MCAR',Nbit_comp*3),rep('True',Nbit_comp*3),rep('MNARz',Nbit_comp*3),rep('Mice',Nbit_comp*3)),4) #,c(rep('MCAR',Nbit_comp*3),rep('MNARy',Nbit_comp*3),rep('MNARz',Nbit_comp*3),rep('Mice',Nbit_comp*3)),c(rep('MCAR',Nbit_comp*3),rep('MNARyz',Nbit_comp*3),rep('MNARz',Nbit_comp*3),rep('Mice',Nbit_comp*3)),c(rep('MCAR',Nbit_comp*3),rep('MNARykzj',Nbit_comp*3),rep('MNARz',Nbit_comp*3),rep('Mice',Nbit_comp*3)))
df_plot_main['Number of variables']=rep(rep(c(rep('3',Nbit_comp),rep('6',Nbit_comp),rep('9',Nbit_comp)),4),4)
df_plot_main['True mechanism']=c(rep('MNARzj',Nbit_comp*3*4),rep('MNARy',Nbit_comp*3*4),rep('MNARyz',Nbit_comp*3*4),rep('MNARykzj',Nbit_comp*3*4))

plot_main <- ggplot(df_plot_main, aes(x=`Number of variables`, y=ARI)) + geom_boxplot(aes(fill=Mechanism))
plot_main <- plot_main + facet_grid(. ~ `True mechanism` ) + scale_fill_manual(values=cbp1[c(1,11,2,12)])



plot_main2 <- ggplot(df_plot_main, aes(x=`Mechanism`, y=ARI)) + geom_hline(yintercept=0.9,color="red",linetype="dashed")+geom_boxplot(aes(fill=Mechanism))
plot_main2 <- plot_main2 + facet_grid(`Number of variables` ~ `True mechanism` ) + scale_fill_manual(values=cbp1[c(1,11,2,12)])+xlab('')
