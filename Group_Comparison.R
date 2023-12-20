if (!require("dplyr")){
  install.packages('dplyr')
  library(dplyr)
}

if (!require("lme4")){
  install.packages('lme4')
  library(lme4)
}



group_member=read.table("group_membership_3group.txt", header = TRUE)
colnames(group_member)[2]="membership"

tot_var=read.csv("./Integ_var.csv", header = TRUE)
invar=read.csv("./Integ_invar.csv", header = TRUE)
gray_vol=read.table("./IMAGEN_full_v4_Age.txt", header = TRUE, sep=" ")%>%left_join(invar%>%select(SubID,sex))

var=tot_var%>%select(SubID,Age,PRM_Percent_correct,
                     AGN_Total_omissions_positive,AGN_Total_omissions_negative,
                     SWM_Between_errors,SWM_Strategy,starts_with("CGT"),
                     RVP_A,starts_with("Kest"),GO_RT,GO_ACC,Stop_ACC,SSRT,
                     starts_with("NEO"),starts_with("TCI"),starts_with("LEQ"),
                     sebdtot,semotion,sconduct,shyper,speer,sprosoc,simpact,
                     Lifetime,Frequency_last30,starts_with("Drink"),starts_with("Row"),
                     dep,dcmadep,adhd,dcadhd,padhd)
var=var%>%left_join(group_member)%>%filter(!is.na(membership))%>%
  left_join(invar%>%select(SubID,sex,mo_edu,fa_edu,ImagingCentreID,Handedness,mo_eth,fa_eth,WISCIV))%>%
  left_join(gray_vol%>%select(SubID,Time,EstimatedTotalIntraCranialVol),by=c("Age"="Time","SubID"='SubID'))%>%
  mutate(Age_stan=(Age-14)/5)%>%
  left_join(invar%>%select(SubID,WISCIV))

var=var%>%mutate(across(c(membership,ImagingCentreID,sex),as.factor))
var=as.data.frame(var)

coef_p=matrix(NA,59*3,6)
v=colnames(var)[3:61]
rownames(coef_p)=rep(v,each=3)
colnames(coef_p)=c("Glow_normal_d","Glow_normal_p","Glate_normal_d","Glate_normal_p","Glate_low_d","Glate_low_p")


for (i in 3:60){
  
  ################################low/late vs normal
  fit=lmer(var[,i]~Age_stan*relevel(membership,ref=1)+sex+(1|SubID)+ImagingCentreID+handedness+EstimatedTotalIntraCranialVol, 
           data=var,REML=TRUE,control = lmerControl(calc.derivs = FALSE))
  a=summary(fit)$coefficients
  # trajectory
  coef_p[(i-2)*3-1,2]=a[15,5]
  coef_p[(i-2)*3-1,4]=a[15,6]
  
  
  tmp=var%>%filter(Age==14)
  fit=glm(tmp[,i]~relevel(membership,ref=1)+sex+ImagingCentreID+handedness+EstimatedTotalIntraCranialVol,data=tmp)
  a=summary(fit)$coefficients
  #baseline
  coef_p[(i-2)*3-2,2]=a[2,4]
  coef_p[(i-2)*3-2,4]=a[3,4]
  
  
  tmp=var%>%filter(Age==23)
  if (all(is.na(tmp[,i]))){
    tmp=var%>%filter(Age==19)
    fit=glm(tmp[,i]~relevel(membership,ref=1)+sex+ImagingCentreID+handedness+EstimatedTotalIntraCranialVol,data=tmp)
  }else{
    fit=glm(tmp[,i]~relevel(membership,ref=1)+sex+ImagingCentreID+handedness+EstimatedTotalIntraCranialVol,data=tmp)
  }
  
  a=summary(fit)$coefficients
  #further
  coef_p[(i-2)*3,2]=a[2,4]
  coef_p[(i-2)*3,4]=a[3,4]
  
  
  ################################late vs low
  fit=lmer(var[,i]~Age_stan*relevel(membership,ref=3)+sex+(1|SubID)+ImagingCentreID+handedness+EstimatedTotalIntraCranialVol, 
           data=var,REML=TRUE,control = lmerControl(calc.derivs = FALSE))
  
  a=summary(fit)$coefficients
  #trajectory
  coef_p[(i-2)*3-1,6]=a[15,5]
  
  tmp=var%>%filter(Age==14)
  fit=glm(tmp[,i]~relevel(membership,ref=3)+sex+ImagingCentreID+handedness+EstimatedTotalIntraCranialVol,data=tmp)
  a=summary(fit)$coefficients
  #baseline
  coef_p[(i-2)*3-2,6]=a[2,4]
  
  tmp=var%>%filter(Age==23)
  if (all(is.na(tmp[,i]))){
    tmp=var%>%filter(Age==19)
    fit=glm(tmp[,i]~relevel(membership,ref=3)+sex+ImagingCentreID+handedness+EstimatedTotalIntraCranialVol,data=tmp)
  }else{
    fit=glm(tmp[,i]~relevel(membership,ref=3)+sex+ImagingCentreID+handedness+EstimatedTotalIntraCranialVol,data=tmp)
  }

  a=summary(fit)$coefficients
  #further
  coef_p[(i-2)*3,6]=a[2,4]

}

rownames(coef_p)[seq(2,nrow(coef_p),by=3)]=paste0(rownames(coef_p)[seq(2,nrow(coef_p),by=3)],"_trajectory")
rownames(coef_p)[seq(3,nrow(coef_p),by=3)]=paste0(rownames(coef_p)[seq(3,nrow(coef_p),by=3)],"_outcome")


# padhd
i=61
tmp=var%>%filter(Age==14)
fit=glm(padhd~relevel(membership,ref=1)+sex+ImagingCentreID+handedness+EstimatedTotalIntraCranialVol,data=tmp)
a=summary(fit)$coefficients
coef_p[(i-2)*3-2,2]=a[2,4]
coef_p[(i-2)*3-2,4]=a[3,4]
fit=glm(padhd~relevel(membership,ref=3)+sex+ImagingCentreID+handedness+EstimatedTotalIntraCranialVol,data=tmp)
a=summary(fit)$coefficients
coef_p[(i-2)*3-2,6]=a[2,4]
tmp=var%>%filter(Age==16)
fit=glm(padhd~relevel(membership,ref=1)+sex+ImagingCentreID+handedness+EstimatedTotalIntraCranialVol,data=tmp)
a=summary(fit)$coefficients
coef_p[(i-2)*3,2]=a[2,4]
coef_p[(i-2)*3,4]=a[3,4]
fit=glm(padhd~relevel(membership,ref=3)+sex+ImagingCentreID+handedness+EstimatedTotalIntraCranialVol,data=tmp)
a=summary(fit)$coefficients
coef_p[(i-2)*3,6]=a[2,4]
