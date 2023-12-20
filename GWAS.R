if (!require("dplyr")){
  install.packages('dplyr')
  library(dplyr)
}


set.seed(026)

cov=read.table("./Integ_var_invar.txt",sep=" ",header=TRUE)%>%
  select(SubID,sex,ImagingCentreID,Handedness)%>%unique()
imagen=read.table("./IMAGEN_full_v4_all_Age.txt",sep=" ",header=TRUE)
group_member=read.table("./group_membership_3group.txt",sep=" ",header=TRUE)
imagen=imagen%>%left_join(group_member[,1:2])%>%mutate(membership=factor(membership))

imagen=imagen%>%filter(Time==14)%>%dplyr::select(-Time)
imagen=na.omit(imagen)

imagen_norm=imagen[,c(1:36,39:45,53)]%>%mutate(across(c(bankssts:Accumbens_area),scale))%>%
  mutate(membership=as.numeric(membership))

# Taking Group 3 vs Group 1/2 as an example
tmp=imagen_norm%>%
  mutate(membership=case_when(membership!=3~1,membership==3~membership))%>%
  mutate(membership=as.factor(membership))%>%
  left_join(cov)%>%
  mutate(across(c(sex,ImagingCentreID,Handedness),as.factor))

tmp$membership=relevel(tmp$membership,ref=1)

coef_p=matrix(NA,41,3)

for (i in 3:43){
  fit_model=glm(membership~tmp[,i]+sex+ImagingCentreID+handedness,data=tmp,family=binomial(),control=list(maxit=100))
  t=summary(fit_model)$coefficients
  coef_p[i-2,]=t[2,c(1,3:4)]
}
rownames(coef_p)=colnames(tmp)[3:43]
coef_p=data.frame(coef_p)%>%arrange(desc(abs(X1)))
coef_p=coef_p%>%select(ROI,X1,X2,X3)%>%rename(loading=X1,t=X2,p=X3)

coef_3vsall=rownames(coef_p)[1:10]

imagen=imagen%>%left_join(cov)%>%
  mutate(across(c(sex,ImagingCentreID,handedness),as.factor))

imagen_adjust=imagen

# adjust for site, sex and handedness
for (i in 3:47){
  fit_model=glm(imagen[,i]~ImagingCentreID+sex+handedness,data=imagen,na.action=na.exclude)
  imagen_adjust[,i]=resid(fit_model)+fit_model$coefficients[1]
}

tmp=imagen_adjust%>%mutate(membership=as.numeric(membership))%>%
  mutate(membership=case_when(
    membership==3~1,
    TRUE~membership
  ))%>%
  mutate(membership=factor(membership))

tmp=tmp[,c(coef_3vsall[1:10],"membership")]

tmp$membership=factor(tmp$membership,level=c("1","3"))
fit_model_latevsnonlate=glm(membership~.-membership,data=tmp,family=binomial(),control=list(maxit=100))
fit_model_latevsnonlate=step(fit_model_latevsnonlate,trace=FALSE)
tmp$pred=predict(fit_model_latevsnonlate,tmp,type="link")
# output: pred = Group3-reweighted GMV
