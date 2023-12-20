if (!require("dplyr")){
  install.packages('dplyr')
  library(dplyr)
}

if (!require("lavaan")){
  install.packages('lavaan')
  library(lavaan)
}

################################################################################
# 1. DMP analysis
################################################################################
group_member=read.table("./group_membership_3group.txt", header = T, sep=" ")
cov=read.csv("Methy_cov14_merge.csv")

load("t14_merge.rda")
rownames(data)=data[,1]
data=data[,-1]
data=data.frame(t(data))

ID=data.frame(colnames(data))%>%rename(ID=`colnames.data.`)%>%mutate(ID=as.numeric(gsub("X","",ID)))

covariate=data%>%select(ID)%>%left_join(cov)%>%
  left_join(group_member%>%rename(ID=SubID))
covariate$membership=factor(covariate$membership,levels=c("1","2"))

methy_mat=matrix(NA,nrow=372582,ncol=4)
for (i in 2:372583){
  fit_model=glm(membership~.-membership-ID+data[,i],data=covariate,family=binomial())
  methy_mat[i-1,2:3]=summary(fit_model)$coefficient[22,c(1,4)]
}
methy_mat[,1]=colnames(data)[2:372583]
methy_mat[,4]=p.adjust(methy_mat[,3],method="BH")

write.csv(methy_mat,"methy_mat.csv",row.names=F)


################################################################################
# 2. Mediation analysis
################################################################################
set.seed(026)
data=read.csv("./methylation_final.csv",header=T)%>%select(SubID,cg06064461)
peak=read.table("peak.txt",header=T,sep=" ")
child=read.csv("./child_adversity.csv",header=T)%>%
  select(SubID,starts_with("CTQ"),starts_with("FamStress"),starts_with("Childexp"))

data=data%>%left_join(group_member)%>%left_join(cov)%>%
  left_join(peak)%>%
  mutate(peak_gmv=scale(peak_gmv))%>%
  select(SubID,peak_gmv,cg06064461,sit1,sit2,sit3,sit4,sit5,sit6,
         PC1,PC2,wave1,wave2,V1,V2,PCA_1,PCA_2,PCA_3,PCA_4,sex)%>%
  left_join(child)%>%
  na.omit()

matrix_cg=matrix(NA,13,12)
colnames(matrix_cg)=c("c","c_p","a","a_p","b","b_p","ab","ab_p","total",'total_low',"total_up","total_p")
rownames(matrix_cg)=colnames(data)[21:33]
for (i in c(21:33)){ # all environment factors
  
  tmp=na.omit(data[,c(i,2,3,4:20)])
  colnames(tmp)[1:3]=c("X","Y","M") # X indicate peak GMV
  
  model = ' # direct effect
          Y~c*X+sit1+sit2+sit3+sit4+sit5+sit6+
                    PC1+PC2+wave1+wave2+V1+V2+PCA_1+PCA_2+PCA_3+PCA_4+sex
         # mediator
           M ~ a*X+sit1+sit2+sit3+sit4+sit5+sit6+
                    PC1+PC2+wave1+wave2+V1+V2+PCA_1+PCA_2+PCA_3+PCA_4+sex
           Y ~ b*M
         # indirect effect (a*b)
           ab := a*b
         # total effect
           total := c + (a*b)
        '

  fit=sem(model,data=tmp,se="bootstrap",bootstrap=1000)
  t=parameterestimates(fit,standardized = F,boot.ci.type="bca.simple",level=0.95, ci=TRUE)%>%filter(label %in% c("a","b","c","ab","total")) 
  
  matrix_cg[i-20,1:2]=unlist(t[t$label=="c",c("est","pvalue")])
  matrix_cg[i-20,3:4]=unlist(t[t$label=="a",c("est","pvalue")])
  matrix_cg[i-20,5:6]=unlist(t[t$label=="b",c("est","pvalue")])
  matrix_cg[i-20,7:8]=unlist(t[t$label=="ab",c("est","pvalue")])
  matrix_cg[i-20,9:12]=unlist(t[t$label=="total",c("est","ci.lower","ci.upper","pvalue")])
  
}

write.csv(matrix_cg,"Env~CpG~peak.csv")

