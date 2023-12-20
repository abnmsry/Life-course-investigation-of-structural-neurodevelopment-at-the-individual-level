if (!require("dplyr")){
  install.packages('dplyr')
  library(dplyr)
}


if (!require("lme4")){
  install.packages('lme4')
  library(lme4)
}

if (!require("nlme")){
  install.packages('nlme')
  library(nlme)
}

if (!require("splines")){
  install.packages('splines')
  library(splines)
}

if (!require("mgcv")){
  install.packages('mgcv')
  library(mgcv)
}


if (!require("ggplot2")){
  install.packages('ggplot2')
  library(ggplot2)
}


################################################################################
# 0. Trajectory calculation
################################################################################
data=read.table("demo.txt",sep="\t",header=TRUE)

trajectory_vol_slope=data%>%select(SubID)%>%unique()

for (i in 6:54){ # all ROIs
  fit=lmer(data[,i]~1+age+(1+age|SubID)+ImagenCentreID+handedness+sex+EstimatedTotalIntraCranialVol,
           data=data,control = lmerControl(calc.derivs = FALSE),REML=TRUE)
  fixed=summary(fit)$coefficient
  rand_int=ranef(fit)$SubID[,2]
  trajectory_vol_slope[,(i-4)]=rand_int+fixed[2,1]
}
colnames(trajectory_vol_slope)[2:50]=colnames(data)[6:54]

write.table(trajectory_vol_slope,"trajectory_summary_vol.txt",row.names=F,sep=" ",quote=F)

################################################################################
# 1. Group clustering
################################################################################
trajectory_vol_slope=read.table("trajectory_summary_vol.txt", header = TRUE, sep=" ")
trajectory_vol_ind=trajectory_vol_slope[, c(1, 2:45)]
trajectory_vol_grouped=trajectory_vol_slope[, c(46:50)]

PC_trajectory=prcomp(trajectory_vol_ind[,-1],
                        center = TRUE,
                        scale. = TRUE)
summary(PC_trajectory)
screeplot(PC_trajectory,type="line",lwd=2,npcs=45)

derived_trajectory=cbind(trajectory_vol_ind[, 1], 
                            as.matrix(trajectory_vol_ind[, 2:45]) %*% as.matrix(PC_trajectory$rotation[,1:15]))
derived_trajectory=as.data.frame(derived_trajectory)
colnames(derived_trajectory)=c("SubID", paste0(rep("PC_", 15),1:15))

set.seed(026)

# total_ss=numeric()
# within_ss=numeric()
# between_ss=numeric()
# sil=numeric()
# 
# for(k in 1:10){
#   clustering_k=kmeans(as.matrix(derived_trajectory[, 2:16]), k, iter.max = 100)
#   print(k)
#   total_ss[k]=clustering_k$totss
#   within_ss[k]=clustering_k$tot.withinss
#   between_ss[k]=clustering_k$betweenss
# }

clustering_3 = kmeans(as.matrix(derived_trajectory[, 2:16]), 3, iter.max = 500)
table(clustering_3$cluster)

trajectory_vol_ind$membership = clustering_3$cluster

write.table(trajectory_vol_ind%>%select(SubID,membership),"group_membership.txt",row.names=F,quote=F,na="")



################################################################################
# 2.0 Estimation of GMV development among groups - data loading
################################################################################
IMAGEN_vol=read.table("IMAGEN_full_v4_all.txt", header = T, sep=" ")
ABCD_vol=read.table("ABCD_stats.txt", header = T, sep=",")
HCP_vol=read.table("HCP_Development_stats.txt", header = T, sep=",")
HCP_YA_vol=read.table("HCP_YA_stats.txt", header = T, sep=",")
PNC_vol=read.table("PNC_stats.txt", header = T, sep=",")

IMAGEN_subgroup=read.table("group_membership.txt", header = T)

ctrl = lmeControl(opt='optim')

################################################################################
# 2.1 Estimation of age and region-specific GMV development among groups
################################################################################
IMAGEN_vol$Age_center = IMAGEN_vol$Age
IMAGEN_vol$Age_center_sq = (IMAGEN_vol$Age_center)^2
IMAGEN_vol$Age_center_cub = (IMAGEN_vol$Age_center)^3
IMAGEN_vol_new=IMAGEN_vol

HCP_vol$Age_center = HCP_vol$Age
HCP_vol$Age_center_sq = (HCP_vol$Age_center)^2
HCP_vol$Age_center_cub = (HCP_vol$Age_center)^3

PNC_vol$Age_center = PNC_vol$Age
PNC_vol$Age_center_sq = (PNC_vol$Age_center)^2
PNC_vol$Age_center_cub = (PNC_vol$Age_center)^3

ref_vol=rbind(HCP_vol,PNC_vol)

SUBID_G1 = IMAGEN_subgroup[IMAGEN_subgroup$membership == 2, "SubID"]
SUBID_G2 = IMAGEN_subgroup[IMAGEN_subgroup$membership == 3, "SubID"]
SUBID_G3 = IMAGEN_subgroup[IMAGEN_subgroup$membership == 1, "SubID"]

traj_G1=matrix(NA,nrow=42,ncol=2)
traj_G2=matrix(NA,nrow=42,ncol=2)
traj_G3=matrix(NA,nrow=42,ncol=2)

rownames(traj_G1)=colnames(HCP_vol)[3:44]
colnames(traj_G1)=c("Age_center","Age_center_sq")
rownames(traj_G2)=colnames(HCP_vol)[3:44]
colnames(traj_G2)=c("Age_center","Age_center_sq")
rownames(traj_G3)=colnames(HCP_vol)[3:44]
colnames(traj_G3)=c("Age_center","Age_center_sq")

# build model
for (i in 3:44){ # all cortical regions
  fit_model=glm(ref_vol[,i] ~ 1 + Age_center + Age_center_sq + sex + EstimatedTotalIntraCranialVol, data=ref_vol)
  fix=summary(fit_model)$coefficient
  pre=predict(fit_model,IMAGEN_vol)
  IMAGEN_vol_new[,i]=IMAGEN_vol_new[,i]-pre
  
  tmp=na.omit(IMAGEN_vol_new[,c("SubID","sex","handedness","EstimatedTotalIntraCranialVol",i)])
  colnames(tmp)[5]="ROI"
  fit_model=lme(ROI ~ 1 + Age_center, method="REML", random = ~1+Age_center|SubID, 
                data=tmp,control=ctrl)
  
  rand_int = ranef(fit_model)
  rand_int$SubID = rownames(rand_int)
  rand_int$int = rand_int$`(Intercept)`
  rand_int$`(Intercept)` = NULL
  
  beta_G1 = c(mean(rand_int[rand_int$SubID %in% SUBID_G1, "Age_center"])+fixef(fit_model)[2]+fix[2,1],
              fix[3,1])
  traj_G1[i-2,1:2]=beta_G1
  
  beta_G2 = c(mean(rand_int[rand_int$SubID %in% SUBID_G2, "Age_center"])+fixef(fit_model)[2]+fix[2,1],
              fix[3,1])
  traj_G2[i-2,1:2]=beta_G2
  
  beta_G3 = c(mean(rand_int[rand_int$SubID %in% SUBID_G3, "Age_center"])+fixef(fit_model)[2]+fix[2,1],
              fix[3,1])
  traj_G3[i-2,1:2]=beta_G3
}

traj_G1=data.frame(traj_G1)%>%mutate(traj_5=Age_center+2*Age_center_sq*5,
                                     traj_10=Age_center+2*Age_center_sq*10,
                                     traj_15=Age_center+2*Age_center_sq*15,
                                     traj_20=Age_center+2*Age_center_sq*20,
                                     traj_25=Age_center+2*Age_center_sq*25,
                                     group="G1",
                                     roi=paste0("Left_",rownames(traj_G1)))
traj_G2=data.frame(traj_G2)%>%mutate(traj_5=Age_center+2*Age_center_sq*5,
                                     traj_10=Age_center+2*Age_center_sq*10,
                                     traj_15=Age_center+2*Age_center_sq*15,
                                     traj_20=Age_center+2*Age_center_sq*20,
                                     traj_25=Age_center+2*Age_center_sq*25,
                                     group="G2",
                                     roi=paste0("Left_",rownames(traj_G2)))
traj_G3=data.frame(traj_G3)%>%mutate(traj_5=Age_center+2*Age_center_sq*5,
                                     traj_10=Age_center+2*Age_center_sq*10,
                                     traj_15=Age_center+2*Age_center_sq*15,
                                     traj_20=Age_center+2*Age_center_sq*20,
                                     traj_25=Age_center+2*Age_center_sq*25,
                                     group="G3",
                                     roi=paste0("Left_",rownames(traj_G3)))
traj=rbind(traj_G1,traj_G2,traj_G3)

write.csv(traj,"trajectory_roi_subgroup.csv",row.names=F)

################################################################################
# 2.2 Estimation of group-specific developmental curve of total GMV in IMAGEN
################################################################################
all_vol=rbind(
  ABCD_vol[, c("studyID", "SubID", "Age", "sex", "EstimatedTotalIntraCranialVol",
               "totalGM", "CortexVol", "SubCortGrayVol")],
  HCP_vol[,  c("studyID", "SubID", "Age", "sex", "EstimatedTotalIntraCranialVol",
               "totalGM", "CortexVol", "SubCortGrayVol")],
  PNC_vol[,  c("studyID", "SubID", "Age", "sex", "EstimatedTotalIntraCranialVol",
               "totalGM", "CortexVol", "SubCortGrayVol")])

longitudinal_vol=rbind(IMAGEN_vol[, c("studyID", "SubID", "Age", "sex", "EstimatedTotalIntraCranialVol",
                                      "totalGM", "CortexVol", "SubCortGrayVol")],
                       ABCD_vol[, c("studyID", "SubID", "Age", "sex", "EstimatedTotalIntraCranialVol",
                                    "totalGM", "CortexVol", "SubCortGrayVol")])


SUBID_G1 = IMAGEN_subgroup[IMAGEN_subgroup$membership == 1, "SubID"]
SUBID_G2 = IMAGEN_subgroup[IMAGEN_subgroup$membership == 2, "SubID"]
SUBID_G3 = IMAGEN_subgroup[IMAGEN_subgroup$membership == 3, "SubID"]


all_vol$Age_center = (all_vol$Age - 14)/5
all_vol$Age_center_sq = (all_vol$Age_center)^2
all_vol$Age_center_cub = (all_vol$Age_center)^3

longitudinal_vol$Age_center = (longitudinal_vol$Age - 14)/5
longitudinal_vol$Age_center_sq = (longitudinal_vol$Age_center)^2
longitudinal_vol$Age_center_cub = (longitudinal_vol$Age_center)^3


# 1-step: fit curve for estimated intracranial volume

# civ_model_1 = lme(EstimatedTotalIntraCranialVol ~ 1 + Age_center + sex, method="ML", random = list(~1|SubID,~1|studyID), data=all_vol, control=ctrl)
# civ_model_2 = lme(EstimatedTotalIntraCranialVol ~ 1 + Age_center + Age_center_sq + sex, method="ML", random = list(~1|SubID,~1|studyID), data=all_vol, control=ctrl)
# civ_model_3 = lme(EstimatedTotalIntraCranialVol ~ 1 + Age_center + Age_center_sq + Age_center_cub + sex, method="ML", random = list(~1|SubID,~1|studyID), data=all_vol, control=ctrl)
# anova(civ_model_1,civ_model_2,civ_model_3)

plot_civ = data.frame(Age_center=(5:23-14)/5,
                       Age_center_sq=((5:23-14)/5)^2,
                       Age_center_cub=((5:23-14)/5)^3)
plot_civ$Age=(plot_civ$Age_center)*5+14

civ_model = lme(EstimatedTotalIntraCranialVol ~ 1 + Age_center + Age_center_sq + sex, method="ML", random = list(~1|SubID,~1|studyID), data=all_vol, control=ctrl)
plot_civ$predvalue2=predict(civ_model, plot_civ, level=0) # predict the output for each point based on the model
designmat_civ = model.matrix(eval(eval(civ_model$call$fixed)[-2]), plot_civ) # make design matrix
plot_civ$SDvalue2=sqrt(diag(designmat_civ %*% civ_model$varFix %*% t(designmat_civ))) # calculate standard deviation for each point for each model
plot_civ$lowerCIvalue2=(plot_civ$predvalue2-(1.96*plot_civ$SDvalue2)) # calculate confidence intervals - lower
plot_civ$upperCIvalue2=(plot_civ$predvalue2+(1.96*plot_civ$SDvalue2)) # calculate confidence intervals - upper

plot(x = plot_civ$Age, y = plot_civ$predvalue2, type = "l", col="deeppink",
     xlim=c(5, 23), ylim = c(min(plot_civ$lowerCIvalue2)*0.95, max(plot_civ$lowerCIvalue2)*1.05), lwd=2,pch=18, xlab="Age", ylab="Intra Cranial Volume") 
lines(x = plot_civ$Age, y = plot_civ$lowerCIvalue2, type = "l", lty=2, pch=2, col=adjustcolor("deeppink",alpha=1)) # lower CI
lines(x = plot_civ$Age, y = plot_civ$upperCIvalue2, type = "l", lty=2, pch=2, col=adjustcolor("deeppink",alpha=1)) # upper CI


# 2-step: fit curve for totalGMV

longitudinal_vol$studyID=factor(longitudinal_vol$studyID)
totalGM_model = lme(totalGM ~ 1 + Age_center + Age_center_sq+Age_center_cub + EstimatedTotalIntraCranialVol + sex + studyID, 
                     method="ML", random = ~1+Age_center|SubID, data=longitudinal_vol, control=ctrl)
# extract their random effects ( individual specific coefficients)
rand_int = ranef(totalGM_model)
rand_int$SubID = rownames(rand_int)
rand_int$int = rand_int$`(Intercept)`
rand_int$`(Intercept)` = NULL

plot_totalGM = data.frame(Age_center=(5:23-14)/5,
                           Age_center_sq=((5:23-14)/5)^2,
                           Age_center_cub=((5:23-14)/5)^3,
                           EstimatedTotalIntraCranialVol=plot_civ$predvalue2)
plot_totalGM$Age = (plot_totalGM$Age_center)*5+14

# for subjects in IMAGEN subgroup 1
designmat_totalGM = model.matrix(~ Age_center + Age_center_sq +Age_center_cub+ EstimatedTotalIntraCranialVol, data=plot_totalGM) # make design matrix
beta_totalGM_G1 = c(median(rand_int[rand_int$SubID %in% SUBID_G1, "int"])+fixef(totalGM_model)[1], 
                     median(rand_int[rand_int$SubID %in% SUBID_G1, "Age_center"])+fixef(totalGM_model)[2],
                     fixef(totalGM_model)[3:5])
plot_totalGM$pred_G1 = designmat_totalGM %*% beta_totalGM_G1
beta_totalGM_G1_sub=cbind(rand_int[rand_int$SubID %in% SUBID_G1, "int"]+fixef(totalGM_model)[1],
                          rand_int[rand_int$SubID %in% SUBID_G1, "Age_center"]+fixef(totalGM_model)[2],
                          t(replicate(nrow(rand_int[rand_int$SubID %in% SUBID_G1,]),fixef(totalGM_model)[3:5])))
tmp=designmat_totalGM %*% t(beta_totalGM_G1_sub)
tmp=apply(tmp,1,function(x) quantile(x,probs=c(0.025,0.975)))
plot_totalGM$lowerCI_G1=smooth(tmp[1,])
plot_totalGM$upperCI_G1=smooth(tmp[2,])

# for subjects in IMAGEN subgroup 2
beta_totalGM_G2 = c(median(rand_int[rand_int$SubID %in% SUBID_G2, "int"])+fixef(totalGM_model)[1], 
                     median(rand_int[rand_int$SubID %in% SUBID_G2, "Age_center"])+fixef(totalGM_model)[2],
                     fixef(totalGM_model)[3:5])
plot_totalGM$pred_G2 = designmat_totalGM %*% beta_totalGM_G2
beta_totalGM_G2_sub=cbind(rand_int[rand_int$SubID %in% SUBID_G2, "int"]+fixef(totalGM_model)[1],
                          rand_int[rand_int$SubID %in% SUBID_G2, "Age_center"]+fixef(totalGM_model)[2],
                          t(replicate(nrow(rand_int[rand_int$SubID %in% SUBID_G2,]),fixef(totalGM_model)[3:5])))
tmp=designmat_totalGM %*% t(beta_totalGM_G2_sub)
tmp=apply(tmp,1,function(x) quantile(x,probs=c(0.025,0.975)))
plot_totalGM$lowerCI_G2=smooth(tmp[1,])
plot_totalGM$upperCI_G2=smooth(tmp[2,])

# for subjects in IMAGEN subgroup 3
beta_totalGM_G3 = c(median(rand_int[rand_int$SubID %in% SUBID_G3, "int"])+fixef(totalGM_model)[1], 
                     median(rand_int[rand_int$SubID %in% SUBID_G3, "Age_center"])+fixef(totalGM_model)[2],
                     fixef(totalGM_model)[3:5])
plot_totalGM$pred_G3 = designmat_totalGM %*% beta_totalGM_G3
beta_totalGM_G3_sub=cbind(rand_int[rand_int$SubID %in% SUBID_G3, "int"]+fixef(totalGM_model)[1],
                          rand_int[rand_int$SubID %in% SUBID_G3, "Age_center"]+fixef(totalGM_model)[2],
                          t(replicate(nrow(rand_int[rand_int$SubID %in% SUBID_G3,]),fixef(totalGM_model)[3:5])))
tmp=designmat_totalGM %*% t(beta_totalGM_G3_sub)
tmp=apply(tmp,1,function(x) quantile(x,probs=c(0.025,0.975)))
plot_totalGM$lowerCI_G3=smooth(tmp[1,])
plot_totalGM$upperCI_G3=smooth(tmp[2,])

plot_totalGM$type=as.factor(if_else(plot_totalGM$Age>=14,0,1))


ggplot(data=plot_totalGM)+
  
  geom_ribbon(aes(x=Age,ymin=lowerCI_G1,ymax=upperCI_G1,fill="G1"),alpha=0.11)+
  geom_line(aes(x=Age,y=pred_G1,colour="G1",linetype=type),lwd=1,alpha=1)+
  
  geom_ribbon(aes(x=Age,ymin=lowerCI_G2,ymax=upperCI_G2,fill="G2"),alpha=0.11)+
  geom_line(aes(x=Age,y=pred_G2,colour="G2",linetype=type),lwd=1,alpha=1)+
  
  geom_ribbon(aes(x=Age,ymin=lowerCI_G3,ymax=upperCI_G3,fill="G3"),alpha=0.11)+
  geom_line(aes(x=Age,y=pred_G3,colour="G3",linetype=type),lwd=1,alpha=1)+
  
  scale_colour_manual(values=c(G1="#bf3eff",G2="#008b00", G3="#ff7f00"),
                      labels = c("Group 1 (46.07%)", "Group 2 (49.58%)","Group 3 (4.34%)"))+
  
  scale_fill_manual(values=c(G1="#bf3eff",G2="#008b00", G3="#ff7f00"),
                    labels = c("Group 1", "Group 2","Group 3"))+
  
  coord_cartesian(xlim=c(5, 23), ylim = c(450000, 750000))+
  guides(colour=guide_legend(title=NULL,override.aes=list(size=1.5,width=1)),
         linetype="none",size="none",fill="none")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        legend.position = c(0.95, 0.95),#plot??λ??
        legend.justification = c("right", "top"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.height=unit(0.5,"cm"),
        legend.key.width =unit(1.5,"cm"),
        legend.spacing.y = unit(-0.12, 'cm'),
        legend.text = element_text(size=10),
        axis.title.y = element_text(margin=margin(r=6)),
        axis.text=element_text(size=10),
        axis.title=element_text(size=11))+
  labs(x="Age", y="Total GMV")


################################################################################
# 3. Estimation of peak total GMV in IMAGEN
################################################################################
localMaxima = function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y = diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y = cumsum(rle(y)$lengths)
  y = y[seq.int(1L, length(y), 2L)]
  if (x[[1]] >= x[[2]]) {
    y = y[-1]
  }
  y
}


peak_local=function(x,type){
  tmp=which(diff(sign(diff(x)))==-2)
  tmp=tmp[length(1)]
  if (length(tmp)!=0){
    if (type=="peak_value"){
      as.numeric(x)[tmp+1]
    }else if (type=="peak_time"){
      names(x)[tmp+1]
    }
  }else{
    NA
  }
}

all_vol = rbind(IMAGEN_vol[, c("studyID", "SubID", "Age", "sex", "EstimatedTotalIntraCranialVol",
                                "totalGM", "CortexVol", "SubCortGrayVol")],
                 HCP_vol[,  c("studyID", "SubID", "Age", "sex", "EstimatedTotalIntraCranialVol",
                              "totalGM", "CortexVol", "SubCortGrayVol")],
                 HCP_YA_vol[,  c("studyID", "SubID", "Age", "sex", "EstimatedTotalIntraCranialVol",
                              "totalGM", "CortexVol", "SubCortGrayVol")],
                 PNC_vol[,  c("studyID", "SubID", "Age", "sex", "EstimatedTotalIntraCranialVol",
                              "totalGM", "CortexVol", "SubCortGrayVol")],
                 ABCD_vol[,  c("studyID", "SubID", "Age", "sex", "EstimatedTotalIntraCranialVol",
                               "totalGM", "CortexVol", "SubCortGrayVol")])
all_vol$studyID=as.factor(all_vol$studyID)

longitudinal_vol = rbind(IMAGEN_vol[, c("studyID", "SubID", "Age", "sex","EstimatedTotalIntraCranialVol",
                                         "totalGM", "CortexVol", "SubCortGrayVol")],
                          ABCD_vol[, c("studyID", "SubID", "Age", "sex","EstimatedTotalIntraCranialVol",
                                       "totalGM", "CortexVol", "SubCortGrayVol")])
all_vol$Age_center = all_vol$Age-14
all_vol$Age_center_sq = (all_vol$Age_center)^2
all_vol$Age_center_cub = (all_vol$Age_center)^3

longitudinal_vol$Age_center = (longitudinal_vol$Age-14)
longitudinal_vol$Age_center_sq = (longitudinal_vol$Age_center)^2
longitudinal_vol$Age_center_cub = (longitudinal_vol$Age_center)^3


# 1-step

civ_model=lme(EstimatedTotalIntraCranialVol ~ 1+ns(Age_center,df=5)+Age_center:sex+sex+studyID,
              data=all_vol,method="REML",control=ctrl,random=~1+Age_center|SubID)

plot_civ_trial = data.frame(Age_center=(5:37-14),
                             Age_center_sq=((5:37)-14)^2,
                             Age_center_cub=((5:37)-14)^3,
                             gender=0,
                             studyID=as.factor("IMAGEN"))
plot_civ_trial$Age=(plot_civ_trial$Age_center)+14
plot_civ_trial$predvalue2=predict(civ_model, plot_civ_trial, level=0)


IMAGEN_fit=IMAGEN_vol%>%filter(SubID %in% IMAGEN_subgroup$SubID)%>%select(SubID,gender)%>%unique()
plot_civ=data.frame(Age_center=rep(seq(5,37,by=0.05)-14,nrow(IMAGEN_fit)),
                    SubID=rep(IMAGEN_fit$SubID,each=641),
                    gender=rep(IMAGEN_fit$gender,each=641),
                    studyID=rep(as.factor("IMAGEN"),each=641))
plot_civ$predvalue_civ=predict(civ_model, plot_civ)
plot_civ=plot_civ%>%left_join(IMAGEN_subgroup%>%select(SubID,membership))



# 2-step
longitudinal_vol$studyID=as.factor(longitudinal_vol$studyID)

totalGM_model=lme(totalGM~1+EstimatedTotalIntraCranialVol+bs(Age_center,degree=2,df=12)+Age_center:sex+sex+studyID,random=list(~1+Age_center|SubID),
                  data=longitudinal_vol,method="REML")


# plot
# library(RColorBrewer)
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# set.seed(007)
# col_group1=sample(col_vector, 765,replace=TRUE)
# col_group2=sample(col_vector, 80,replace=TRUE)
# col_group3=sample(col_vector, 711,replace=TRUE)
# colnames(plot_civ)[5]="EstimatedTotalIntraCranialVol"
# plot_civ_longi=plot_civ%>%filter(SubID %in% longitudinal_vol$SubID)
# plot_civ_longi$predvalue_gmv=predict(totalGM_model, plot_civ_longi)
# plot_civ_longi$Age=plot_civ_longi$Age+14
# 
# p1=ggplot(plot_civ_longi%>%filter(membership==1))+
#   geom_line(aes(x=Age,y=predvalue_gmv,col=as.factor(SubID),group=SubID),stat="smooth",
#             method="gam",alpha=0.7,size=0.2,formula=y~ns(x,df=3),se=FALSE)+
#   theme_classic()+
#   coord_cartesian(ylim=c(100000,1000000))+
#   scale_color_manual(values=col_group3)+
#   guides(col="none")+
#   labs(title="Group 1")
# 
# p2=ggplot(plot_civ_longi%>%filter(membership==2))+
#   geom_line(aes(x=Age,y=predvalue_gmv,col=as.factor(SubID),group=SubID),stat="smooth",
#             method="gam",alpha=0.7,size=0.2,formula=y~ns(x,df=3),se=FALSE)+
#   theme_classic()+
#   coord_cartesian(ylim=c(100000,1000000))+
#   scale_color_manual(values=col_group2)+
#   guides(col="none")+
#   labs(title="Group 3")
# 
# 
# p3=ggplot(plot_civ_longi%>%filter(membership==3))+
#   geom_line(aes(x=Age,y=predvalue_gmv,col=as.factor(SubID),group=SubID),stat="smooth",
#             method="gam",alpha=0.5,size=0.2,formula=y~ns(x,df=3),se=FALSE)+
#   theme_classic()+
#   coord_cartesian(ylim=c(100000,1000000))+
#   scale_color_manual(values=col_group1)+
#   guides(col="none")+
#   labs(title="Group 2")
# 
# p1+p2+p3


# peak

peak=plot_civ_longi%>%select(SubID,Age,predvalue_gmv,membership)%>%
  filter(!is.na(membership))%>%mutate(membership=as.factor(membership))
peak=tidyr::spread(peak,Age,predvalue_gmv)
tmp=peak[,3:643]

subgroup=peak%>%select(SubID)
subgroup$peak_value=apply(peak[,3:643],1,function(x) peak_local(x,"peak_value"))
subgroup$peak_time=apply(peak[,3:643],1,function(x) peak_local(x,"peak_time"))

peak=peak%>%left_join(subgroup)

peak$peak_time=as.numeric(peak$peak_time)
peak$peak_value=as.numeric(peak$peak_value)
peak=peak%>%select(SubID,peak_value,peak_time,membership)
peak$membership=as.factor(peak$membership)
peak=na.omit(peak)
peak%>%group_by(membership)%>%summarise(value_sd=quantile(peak_value,probs=c(0.025,0.975)),
                                        peak_value=median(peak_value),
                                        time_sd=quantile(peak_time,probs=c(0.025,0.975)),
                                        peak_time=median(peak_time))

write.table(peak,"peak.txt",sep=" ",row.names=F,quote=F,na="")
