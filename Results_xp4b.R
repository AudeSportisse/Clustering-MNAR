
load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp4b/xp4b_10perc.RData")

MCAR_choiceK_10per <- sum(apply(ICL_MCAR,1,which.max)==3)/nb_it
MNARz_choiceK_10per <- sum(apply(ICL_MNARz,1,which.max)==3)/nb_it


load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp4b/xp4b_30perc.RData")

MCAR_choiceK_30per <- sum(apply(ICL_MCAR,1,which.max)==3)/nb_it
MNARz_choiceK_30per <- sum(apply(ICL_MNARz,1,which.max)==3)/nb_it


load("~/Dropbox/clusteringMNAR/codes/Code pour Github/xp4b/xp4b_50perc.RData")

MCAR_choiceK_50per <- sum(apply(ICL_MCAR,1,which.max)==3)/nb_it
MNARz_choiceK_50per <- sum(apply(ICL_MNARz,1,which.max)==3)/nb_it


