#MNARz, MNARz conc, MNARzj, MNARy, MNARyk, MNARyz, MNARykzj, MNARykz, MNARyzj

Nbit_comp <-50
#bon model

cbp1 <- c("#999999","#56B4E9","#000033","#0072B2","#CC79A7","#990033","#FFDB6D", "#C4961A", "#F4EDCA", 
          "#D16103","#660033","#00CC00")



###############
##MNARz
###############

load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARz_n100.RData")
ARI_MNARz_n100 <- ARI_MNARz
ARI_MNARz_n100_conc <- ARI_MNARz_conc
ARI_MCAR_n100_z <- ARI_MCAR
Tdiff_MNARz_n100 <- Tdiff_MNARz
Tdiff_MNARz_n100_conc <- Tdiff_MNARz_comp
Tdiff_MCAR_n100 <- Tdiff_MCAR

load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARz_250.RData")

ARI_MNARz_n250 <- ARI_MNARz
ARI_MCAR_n250_z <- ARI_MCAR
ARI_MNARz_n250_conc <- ARI_MNARz_conc
Tdiff_MNARz_n250 <- Tdiff_MNARz
Tdiff_MNARz_n250[Tdiff_MNARz_n250<5]=Tdiff_MNARz_n250[Tdiff_MNARz_n250<5]*60
Tdiff_MNARz_n250_conc <- Tdiff_MNARz_comp
Tdiff_MCAR_n250 <- Tdiff_MCAR
Tdiff_MCAR_n250[Tdiff_MCAR_n250<5]=Tdiff_MCAR_n250[Tdiff_MCAR_n250<5]*60


load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARz_500.RData")

ARI_MNARz_n500 <- ARI_MNARz
ARI_MCAR_n500_z <- ARI_MCAR
ARI_MNARz_n500_conc <- ARI_MNARz_conc
Tdiff_MNARz_n500 <- Tdiff_MNARz
Tdiff_MNARz_n500 <- Tdiff_MNARz_n500*60
Tdiff_MNARz_n500_conc <- Tdiff_MNARz_comp
Tdiff_MCAR_n500 <- Tdiff_MCAR*60

df_plot1 = data.frame(ARI=c(ARI_MCAR_n100_z,ARI_MCAR_n250_z,ARI_MCAR_n500_z,ARI_MNARz_n100,ARI_MNARz_n250,ARI_MNARz_n500,ARI_MNARz_n100_conc,ARI_MNARz_n250_conc,ARI_MNARz_n500_conc))
df_plot1['Mechanism']=c(rep('1.MCAR',Nbit_comp*3),rep('2.MNARz',Nbit_comp*3),rep('3.MAR',Nbit_comp*3))
df_plot1['Sample size']=c(rep(c(rep('100',Nbit_comp),rep('250',Nbit_comp),rep('500',Nbit_comp)),3))
plotz <- ggplot(data=df_plot1,aes(x=`Sample size`, y=ARI,fill=Mechanism))+geom_hline(yintercept=0.9,color="red",linetype="dashed") +geom_boxplot()
plotz = plotz #+ theme(legend.position = "none")
plotz = plotz + ggtitle('MNARz') +    scale_fill_manual(values=cbp1[c(1,2,3)],labels=c("MCAR","MNARz","RMixtComp"))
plotz

###############
##MNARzj
###############

load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARzj_100.RData")
ARI_MNARzj_n100 <- ARI_MNARz
ARI_MCAR_n100_zj <- ARI_MCAR
Tdiff_MNARzj_n100 <- Tdiff_MNARz
Tdiff_MNARzj_n100[Tdiff_MNARzj_n100<5] <- Tdiff_MNARzj_n100[Tdiff_MNARzj_n100<5]*60 


load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARzj_250.RData")

ARI_MNARzj_n250 <- ARI_MNARz
ARI_MCAR_n250_zj <- ARI_MCAR
Tdiff_MNARzj_n250 <- Tdiff_MNARz
Tdiff_MNARzj_n250[Tdiff_MNARzj_n250<5]=Tdiff_MNARzj_n250[Tdiff_MNARzj_n250<5]*60


load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARzj_500.RData")

ARI_MNARzj_n500 <- ARI_MNARz
ARI_MCAR_n500_zj <- ARI_MCAR
Tdiff_MNARzj_n500 <- Tdiff_MNARz
Tdiff_MNARzj_n500 <- Tdiff_MNARzj_n500*60


df_plot1 = data.frame(ARI=c(ARI_MCAR_n100_zj,ARI_MCAR_n250_zj,ARI_MCAR_n500_zj,ARI_MNARzj_n100,ARI_MNARzj_n250,ARI_MNARzj_n500))
df_plot1['Mechanism']=c(rep('MCAR',Nbit_comp*3),rep('MNARzj',Nbit_comp*3))
df_plot1['Size']=c(rep(c(rep('100',Nbit_comp),rep('250',Nbit_comp),rep('500',Nbit_comp)),2))
plotzj <- ggplot(data=df_plot1,aes(x=Size, y=ARI,fill=Mechanism))+geom_boxplot()
plotzj = plotzj #+ theme(legend.position = "none")
plotzj = plotzj + ggtitle('True model is MNARzj') + scale_fill_manual(values=cbp1[c(1,4)])
plotzj


###############
##MNARyk
###############

load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARyk_100.RData")
ARI_MNARyk_n100 <- ARI_MNARy
ARI_MCAR_n100_yk <- ARI_MCAR
Tdiff_MNARyk_n100 <- Tdiff_MNARy
Tdiff_MNARyk_n100 <- Tdiff_MNARy*60

load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARyk_250.RData")

ARI_MNARyk_n250 <- ARI_MNARy
ARI_MCAR_n250_yk <- ARI_MCAR
Tdiff_MNARyk_n250 <- Tdiff_MNARy
Tdiff_MNARyk_n250=Tdiff_MNARyk_n250*60


load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARyk_500.RData")

ARI_MNARyk_n500 <- ARI_MNARy
ARI_MCAR_n500_yk <- ARI_MCAR
Tdiff_MNARyk_n500 <- Tdiff_MNARy
Tdiff_MNARyk_n500 <- Tdiff_MNARyk_n500*60

df_plot1 = data.frame(ARI=c(ARI_MCAR_n100_yk,ARI_MCAR_n250_yk,ARI_MCAR_n500_yk,ARI_MNARyk_n100,ARI_MNARyk_n250,ARI_MNARyk_n500))
df_plot1['Mechanism']=c(rep('MCAR',Nbit_comp*3),rep('MNARyk',Nbit_comp*3))
df_plot1['Size']=c(rep(c(rep('100',Nbit_comp),rep('250',Nbit_comp),rep('500',Nbit_comp)),2))
plotyk <- ggplot(data=df_plot1,aes(x=Size, y=ARI,fill=Mechanism))+geom_boxplot()
plotyk = plotyk #+ theme(legend.position = "none")
plotyk = plotyk + ggtitle('True model is MNARyk') + scale_fill_manual(values=cbp1[c(1,6)])
plotyk


###############
##MNARy
###############

load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARy_100.RData")
ARI_MNARy_n100 <- ARI_MNARy
ARI_MCAR_n100_y <- ARI_MCAR
Tdiff_MNARy_n100 <- Tdiff_MNARy

load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARy_250.RData")

ARI_MNARy_n250 <- ARI_MNARy
ARI_MCAR_n250_y <- ARI_MCAR
Tdiff_MNARy_n250 <- Tdiff_MNARy
Tdiff_MNARy_n250=Tdiff_MNARy_n250*60


load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARy_500.RData")

ARI_MNARy_n500 <- ARI_MNARy
ARI_MCAR_n500_y <- ARI_MCAR
Tdiff_MNARy_n500 <- Tdiff_MNARy
Tdiff_MNARy_n500 <- Tdiff_MNARy_n500*60

df_plot1 = data.frame(ARI=c(ARI_MCAR_n100_y,ARI_MCAR_n250_y,ARI_MCAR_n500_y,ARI_MNARy_n100,ARI_MNARy_n250,ARI_MNARy_n500))
df_plot1['Mechanism']=c(rep('MCAR',Nbit_comp*3),rep('MNARy',Nbit_comp*3))
df_plot1['Size']=c(rep(c(rep('100',Nbit_comp),rep('250',Nbit_comp),rep('500',Nbit_comp)),2))
ploty <- ggplot(data=df_plot1,aes(x=Size, y=ARI,fill=Mechanism))+geom_boxplot()
ploty = ploty #+ theme(legend.position = "none")
ploty = ploty + ggtitle('True model is MNARy') + scale_fill_manual(values=cbp1[c(1,5)])
ploty


###############
##MNARyz
###############

load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARyz_100.RData")
ARI_MNARyz_n100 <- ARI_MNARy
ARI_MCAR_n100_yz <- ARI_MCAR
Tdiff_MNARyz_n100 <- Tdiff_MNARy
Tdiff_MNARyz_n100[Tdiff_MNARyz_n100<5] <- Tdiff_MNARyz_n100[Tdiff_MNARyz_n100<5]*60 


load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARyz_250.RData")

ARI_MNARyz_n250 <- ARI_MNARy
ARI_MCAR_n250_yz <- ARI_MCAR
Tdiff_MNARyz_n250 <- Tdiff_MNARy
Tdiff_MNARyz_n250 <- Tdiff_MNARyz_n250*5*60


load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARyz_500.RData")

ARI_MNARyz_n500 <- ARI_MNARy
ARI_MCAR_n500_yz <- ARI_MCAR
Tdiff_MNARyz_n500 <- Tdiff_MNARy
Tdiff_MNARyz_n500 <- Tdiff_MNARyz_n500*5*60

df_plot1 = data.frame(ARI=c(ARI_MCAR_n100_yz,ARI_MCAR_n250_yz,ARI_MCAR_n500_yz,ARI_MNARyz_n100,ARI_MNARyz_n250,ARI_MNARyz_n500))
df_plot1['Mechanism']=c(rep('MCAR',Nbit_comp*3),rep('MNARyz',Nbit_comp*3))
df_plot1['Size']=c(rep(c(rep('100',Nbit_comp),rep('250',Nbit_comp),rep('500',Nbit_comp)),2))
plotyz <- ggplot(data=df_plot1,aes(x=Size, y=ARI,fill=Mechanism))+geom_boxplot()
plotyz = plotyz #+ theme(legend.position = "none")
plotyz = plotyz + ggtitle('True model is MNARyz') + scale_fill_manual(values=cbp1[c(1,7)])
plotyz


###############
##MNARykzj
###############

load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARykzj_100.RData")
ARI_MNARykzj_n100 <- ARI_MNARy
ARI_MCAR_n100_ykzj <- ARI_MCAR
Tdiff_MNARykzj_n100 <- Tdiff_MNARy
Tdiff_MNARykzj_n100 <- (Tdiff_MNARykzj_n100/Nbit_run)*5 
Tdiff_MNARykzj_n100[Tdiff_MNARykzj_n100<30] <- Tdiff_MNARykzj_n100[Tdiff_MNARykzj_n100<30]*60


load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARykzj_250.RData")

ARI_MNARykzj_n250 <- ARI_MNARy
ARI_MCAR_n250_ykzj <- ARI_MCAR
Tdiff_MNARykzj_n250 <- Tdiff_MNARy
Tdiff_MNARykzj_n250 <- (Tdiff_MNARykzj_n250/Nbit_run)*5 
Tdiff_MNARykzj_n250[Tdiff_MNARykzj_n250<30] <- Tdiff_MNARykzj_n250[Tdiff_MNARykzj_n250<30]*60


load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARykzj_500.RData")

ARI_MNARykzj_n500 <- ARI_MNARy
ARI_MCAR_n500_ykzj <- ARI_MCAR
Tdiff_MNARykzj_n500 <- Tdiff_MNARy
Tdiff_MNARykzj_n500 <- (Tdiff_MNARykzj_n500/Nbit_run)*5 
Tdiff_MNARykzj_n500 <- Tdiff_MNARykzj_n500*60

df_plot1 = data.frame(ARI=c(ARI_MCAR_n100_ykzj,ARI_MCAR_n250_ykzj,ARI_MCAR_n500_ykzj,ARI_MNARykzj_n100,ARI_MNARykzj_n250,ARI_MNARykzj_n500))
df_plot1['Mechanism']=c(rep('MCAR',Nbit_comp*3),rep('MNARykzj',Nbit_comp*3))
df_plot1['Size']=c(rep(c(rep('100',Nbit_comp),rep('250',Nbit_comp),rep('500',Nbit_comp)),2))
plotykzj <- ggplot(data=df_plot1,aes(x=Size, y=ARI,fill=Mechanism))+geom_boxplot()
plotykzj = plotykzj #+ theme(legend.position = "none")
plotykzj = plotykzj + ggtitle('True model is MNARykzj') + scale_fill_manual(values=cbp1[c(1,8)])
plotykzj


###############
##MNARykz
###############

load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARykz_100.RData")
ARI_MNARykz_n100 <- ARI_MNARy
ARI_MCAR_n100_ykz <- ARI_MCAR
Tdiff_MNARykz_n100 <- Tdiff_MNARy
Tdiff_MNARykz_n100 <- (Tdiff_MNARykz_n100/Nbit_run)*5 
Tdiff_MNARykz_n100 <- Tdiff_MNARykz_n100*60


load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARykz_250.RData")

ARI_MNARykz_n250 <- ARI_MNARy
ARI_MCAR_n250_ykz <- ARI_MCAR
Tdiff_MNARykz_n250 <- Tdiff_MNARy
Tdiff_MNARykz_n250 <- (Tdiff_MNARykz_n250/Nbit_run)*5 
Tdiff_MNARykz_n250<- Tdiff_MNARykz_n250*60


load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARykz_500.RData")

ARI_MNARykz_n500 <- ARI_MNARy
ARI_MCAR_n500_ykz <- ARI_MCAR
Tdiff_MNARykz_n500 <- Tdiff_MNARy
Tdiff_MNARykz_n500 <- (Tdiff_MNARykz_n500/Nbit_run)*5 
Tdiff_MNARykz_n500 <- Tdiff_MNARykz_n500*60

df_plot1 = data.frame(ARI=c(ARI_MCAR_n100_ykz,ARI_MCAR_n250_ykz,ARI_MCAR_n500_ykz,ARI_MNARykz_n100,ARI_MNARykz_n250,ARI_MNARykz_n500))
df_plot1['Mechanism']=c(rep('MCAR',Nbit_comp*3),rep('MNARykz',Nbit_comp*3))
df_plot1['Size']=c(rep(c(rep('100',Nbit_comp),rep('250',Nbit_comp),rep('500',Nbit_comp)),2))
plotykz <- ggplot(data=df_plot1,aes(x=Size, y=ARI,fill=Mechanism))+geom_boxplot()
plotykz = plotykz #+ theme(legend.position = "none")
plotykz = plotykz + ggtitle('True model is MNARykz') + scale_fill_manual(values=cbp1[c(1,9)])
plotykz


###############
##MNARyzj
###############

load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARyzj_100.RData")
ARI_MNARyzj_n100 <- ARI_MNARy
ARI_MCAR_n100_yzj <- ARI_MCAR
Tdiff_MNARyzj_n100 <- Tdiff_MNARy
Tdiff_MNARyzj_n100 <- (Tdiff_MNARyzj_n100/Nbit_run)*5 
Tdiff_MNARyzj_n100 <- Tdiff_MNARyzj_n100*60


load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARyzj_250.RData")

ARI_MNARyzj_n250 <- ARI_MNARy
ARI_MCAR_n250_yzj <- ARI_MCAR
Tdiff_MNARyzj_n250 <- Tdiff_MNARy
Tdiff_MNARyzj_n250 <- (Tdiff_MNARyzj_n250/Nbit_run)*5 
Tdiff_MNARyzj_n250<- Tdiff_MNARyzj_n250*60


load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp1/xp1_MNARyzj_500.RData")

ARI_MNARyzj_n500 <- ARI_MNARy
ARI_MCAR_n500_yzj <- ARI_MCAR
Tdiff_MNARyzj_n500 <- Tdiff_MNARy
Tdiff_MNARyzj_n500 <- (Tdiff_MNARyzj_n500/Nbit_run)*5 
Tdiff_MNARyzj_n500 <- Tdiff_MNARyzj_n500*60

df_plot1 = data.frame(ARI=c(ARI_MCAR_n100_yzj,ARI_MCAR_n250_yzj,ARI_MCAR_n500_yzj,ARI_MNARyzj_n100,ARI_MNARyzj_n250,ARI_MNARyzj_n500))
df_plot1['Mechanism']=c(rep('MCAR',Nbit_comp*3),rep('MNARyzj',Nbit_comp*3))
df_plot1['Size']=c(rep(c(rep('100',Nbit_comp),rep('250',Nbit_comp),rep('500',Nbit_comp)),2))
plotyzj <- ggplot(data=df_plot1,aes(x=Size, y=ARI,fill=Mechanism))+geom_boxplot()
plotyzj = plotyzj #+ theme(legend.position = "none")
plotyzj = plotyzj + ggtitle('True model is MNARyzj') + scale_fill_manual(values=cbp1[c(1,9)])
plotyzj


###########################
###########################
###########################

#vec_MNARz <- c(ARI_MCAR_n100_z,ARI_MCAR_n250_z,ARI_MCAR_n500_z,ARI_MNARz_n100,ARI_MNARz_n250,ARI_MNARz_n500,ARI_MNARz_n100_conc,ARI_MNARz_n250_conc,ARI_MNARz_n500_conc)
vec_MNARz <- c(ARI_MCAR_n100_z,ARI_MCAR_n250_z,ARI_MCAR_n500_z,ARI_MNARz_n100,ARI_MNARz_n250,ARI_MNARz_n500)
vec_MNARzj <- c(ARI_MCAR_n100_zj,ARI_MCAR_n250_zj,ARI_MCAR_n500_zj,ARI_MNARzj_n100,ARI_MNARzj_n250,ARI_MNARzj_n500)
vec_MNARy <- c(ARI_MCAR_n100_y,ARI_MCAR_n250_y,ARI_MCAR_n500_y,ARI_MNARy_n100,ARI_MNARy_n250,ARI_MNARy_n500)
vec_MNARyk <- c(ARI_MCAR_n100_yk,ARI_MCAR_n250_yk,ARI_MCAR_n500_yk,ARI_MNARyk_n100,ARI_MNARyk_n250,ARI_MNARyk_n500)
vec_MNARyz <- c(ARI_MCAR_n100_yz,ARI_MCAR_n250_yz,ARI_MCAR_n500_yz,ARI_MNARyz_n100,ARI_MNARyz_n250,ARI_MNARyz_n500)
vec_MNARykzj <- c(ARI_MCAR_n100_ykzj,ARI_MCAR_n250_ykzj,ARI_MCAR_n500_ykzj,ARI_MNARykzj_n100,ARI_MNARykzj_n250,ARI_MNARykzj_n500)
vec_MNARykz <- c(ARI_MCAR_n100_ykz,ARI_MCAR_n250_ykz,ARI_MCAR_n500_ykz,ARI_MNARykz_n100,ARI_MNARykz_n250,ARI_MNARykz_n500)
vec_MNARyzj <- c(ARI_MCAR_n100_yzj,ARI_MCAR_n250_yzj,ARI_MCAR_n500_yzj,ARI_MNARyzj_n100,ARI_MNARyzj_n250,ARI_MNARyzj_n500)


df_plot_main = data.frame(ARI=c(vec_MNARz,vec_MNARzj,vec_MNARy,vec_MNARyk,vec_MNARyz,vec_MNARykzj,vec_MNARykz,vec_MNARyzj))
#df_plot_main['Mechanism']=c(c(rep('MCAR',Nbit_comp*3),rep('True',Nbit_comp*3),rep('with mask MAR',Nbit_comp*3)),rep(c(rep('MCAR',Nbit_comp*3),rep('True',Nbit_comp*3)),7))
#df_plot_main['Size']=c(rep(c(rep('100',Nbit_comp),rep('250',Nbit_comp),rep('500',Nbit_comp)),3),rep(rep(c(rep('100',Nbit_comp),rep('250',Nbit_comp),rep('500',Nbit_comp)),2),7))
df_plot_main['Mechanism']=rep(c(rep('MCAR',Nbit_comp*3),rep('True',Nbit_comp*3)),8)
df_plot_main['Sample size']=rep(rep(c(rep('100',Nbit_comp),rep('250',Nbit_comp),rep('500',Nbit_comp)),2),8)
#df_plot_main['True mechanism']=c(rep('MNARz',Nbit_comp*3*3),rep('MNARzj',Nbit_comp*3*2),rep('MNARy',Nbit_comp*3*2),rep('MNARyk',Nbit_comp*3*2),rep('MNARyz',Nbit_comp*3*2),rep('MNARykzj',Nbit_comp*3*2),rep('MNARykz',Nbit_comp*3*2),rep('MNARyzj',Nbit_comp*3*2))
df_plot_main['True mechanism']=c(rep('1.MNARz',Nbit_comp*3*2),rep('2.MNARzj',Nbit_comp*3*2),rep('3.MNARy',Nbit_comp*3*2),rep('4.MNARyk',Nbit_comp*3*2),rep('5.MNARyz',Nbit_comp*3*2),rep('6.MNARykzj',Nbit_comp*3*2),rep('7.MNARykz',Nbit_comp*3*2),rep('8.MNARyzj',Nbit_comp*3*2))



plot_main <- ggplot(df_plot_main, aes(x=`Sample size`, y=ARI)) + geom_boxplot(aes(fill=Mechanism))
plot_main <- plot_main + facet_grid(. ~ `True mechanism`) + scale_fill_manual(values=cbp1[c(1,12,3)])

plot_main2 <- ggplot(df_plot_main, aes(x=`Mechanism`, y=ARI)) +geom_boxplot(aes(fill=Mechanism))
plot_main2 <- plot_main2 + facet_grid(`Sample size` ~ `True mechanism` ) #+ scale_fill_manual(values=cbp1[c(1,11,2,12)])


library(gridExtra)
mecha.labs.1 <- c("MNARz","MNARzj","MNARy","MNARyz")
names(mecha.labs.1) <- c('1.MNARz','2.MNARzj','3.MNARy','5.MNARyz')
mecha.labs.2 <- c("MNARyk","MNARykzj","MNARykz","MNARyzj")
names(mecha.labs.2) <-c('4.MNARyk','6.MNARykzj','7.MNARykz','8.MNARyzj')

grid.arrange(
  ggplot(data = df_plot_main[df_plot_main$`True mechanism` %in% c('1.MNARz','2.MNARzj','3.MNARy','5.MNARyz'),], aes(x = `Sample size`, y = ARI)) + geom_hline(yintercept=0.9,color="red",linetype="dashed") +
    geom_boxplot(aes(fill=Mechanism)) + facet_grid(. ~ `True mechanism`,labeller=labeller(`True mechanism`=mecha.labs.1)) + scale_fill_manual(values=cbp1[c(1,12)]) ,
  ggplot(data = df_plot_main[df_plot_main$`True mechanism` %in% c('4.MNARyk','6.MNARykzj','7.MNARykz','8.MNARyzj'),], aes(x = `Sample size`, y = ARI)) + geom_hline(yintercept=0.9,color="red",linetype="dashed") +
    geom_boxplot(aes(fill=Mechanism)) + facet_grid(. ~ `True mechanism`,labeller=labeller(`True mechanism`=mecha.labs.2)) + scale_fill_manual(values=cbp1[c(1,12)]) ,
  nrow=2)



###########################
###########################
###########################


time_cal <- c(Tdiff_MCAR_n100,Tdiff_MCAR_n250,Tdiff_MCAR_n500,Tdiff_MNARz_n100,Tdiff_MNARz_n250,Tdiff_MNARz_n500,Tdiff_MNARz_n100_conc,Tdiff_MNARz_n250_conc,Tdiff_MNARz_n500_conc,Tdiff_MNARzj_n100,Tdiff_MNARzj_n250,Tdiff_MNARzj_n500,Tdiff_MNARy_n100,Tdiff_MNARy_n250,Tdiff_MNARy_n500,Tdiff_MNARyk_n100,Tdiff_MNARyk_n250,Tdiff_MNARyk_n500,Tdiff_MNARyz_n100,Tdiff_MNARyz_n250,Tdiff_MNARyz_n500,Tdiff_MNARykzj_n100,Tdiff_MNARykzj_n250,Tdiff_MNARykzj_n500,Tdiff_MNARykz_n100,Tdiff_MNARykz_n250,Tdiff_MNARykz_n500,Tdiff_MNARyzj_n100,Tdiff_MNARyzj_n250,Tdiff_MNARyzj_n500)

mecha <- c(rep('1.MCAR',Nbit_comp*3),rep('2.MNARz',Nbit_comp*3), rep('3.with mask MAR',Nbit_comp*3), rep('4.MNARzj',Nbit_comp*3), rep('5.MNARy',Nbit_comp*3), rep('6.MNARyk',Nbit_comp*3), rep('7.MNARyz',Nbit_comp*3), rep('8.MNARykzj',Nbit_comp*3), rep('9.MNARykz',Nbit_comp*3), rep('99.MNARyzj',Nbit_comp*3))

size <- c(rep(c(rep('100',Nbit_comp),rep('250',Nbit_comp),rep('500',Nbit_comp)),10))


df_plot2 = data.frame(Time=time_cal)
df_plot2['Mechanism']=mecha
df_plot2['Sample size']=size
plot2 <- ggplot(data=df_plot2,aes(x=`Sample size`, y=Time,fill=Mechanism))+geom_boxplot()
plot2 = plot2 + ylab("Time (logscale in seconds)") #+ theme(legend.position = "none")
plot2 = plot2 +scale_fill_manual(values=cbp1,labels = c("MCAR","MNARz", "RMixtComp", "MNARzj", "MNARy", "MNARyk", "MNARyz", "MNARykzj", "MNARykz", "MNARyzj"))
plot2+scale_y_continuous(trans='log2')#+ylim(c(0,2000))

##échelle log?



df_plot2 = data.frame(Time=time_cal[1:450])
mecha <- c(rep('1.MCAR',Nbit_comp*3),rep('2.MNARz',Nbit_comp*3), rep('3.with mask MAR',Nbit_comp*3))
size <- c(rep(c(rep('100',Nbit_comp),rep('250',Nbit_comp),rep('500',Nbit_comp)),3))
df_plot2['Mechanism']=mecha
df_plot2['Sample size']=size
plot2 <- ggplot(data=df_plot2,aes(x=`Sample size`, y=Time,fill=Mechanism))+geom_boxplot() + ylab("Time (in seconds)") +scale_fill_manual(labels=c("MCAR",'MNARz',"RMixtComp"),values=cbp1)
legend = get_legend(plot2)
plot2 = plot2 + theme(legend.position = "none")
plot2


grid.arrange(
  plotz + theme(legend.position = "none") ,
  plot2 ,
  legend,
  ncol=3)

get_legend <- function(p) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}