# README

Please see the google colab link [https://colab.research.google.com/drive/1TWE3ZT30fY1umL_ILoeC2iwS8oTn66aT?usp=sharing] for tutorial. 

To install the $\texttt{htlgmm}$ package, the user can follow

```R
devtools::install_github("RuzhangZhao/htlgmm")
```



One example for $p_{\mathrm{Z}}=20,p_{\mathrm{W}}=3$

```R
# Simulation Example
library(htlgmm)
pZ=20 # overlapping features
pW=3 # unmatched features 
coef<-c(rep(0,pZ+pW))
coef[1:3]<-0.5
coef[c(21,22)]<-0.5
which(coef!=0)

n=400
nE=2000
n_joint=n+nE
main_index<-1:n

set.seed(202)
Z_joint<-matrix(rnorm(n_joint*pZ),n_joint,pZ)
colnames(Z_joint)<-paste0("Z",1:pZ)
W_joint<-matrix(rnorm(n_joint*pW),n_joint,pW)
colnames(W_joint)<-paste0("W",1:pW)
Z<-Z_joint[main_index,]  # separate main and external study for Z
ZE<-Z_joint[-main_index,]

W<-W_joint[main_index,] # only need main study for W

y_joint<-cbind(Z_joint,W_joint)%*%coef+rnorm(n_joint,0,1)
y<-y_joint[main_index] # separate main study y
yE<-y_joint[-main_index] # separate external study y

Z<-scale(Z)
ZE<-scale(ZE)
W<-scale(W)

library(glmnet)
res_glmnet<-cv.glmnet(x=cbind(Z,W),y=y)
est_coef_glmnet<-coef.glmnet(res_glmnet,s="lambda.min")[-1] # without intercept

reslm<-lm(y~.,data = data.frame(y=yE,ZE))
study_external = list(
                Coeff=reslm$coefficients[-1],
                Covariance=vcov(reslm)[-1,-1],
                Sample_size=nE)
study_info<-list()
study_info[[1]] <- study_external # The study_info is also a list because we can support multiple external studies. 

y<-scale(y,center = TRUE, scale = FALSE) # only do centering for main study y. 

res_htlgmm<-cv.htlgmm(y=y,
                      Z=Z,
                      W=W,
                      study_info = study_info,
                      family = "gaussian",
                      penalty_type = "lasso",
                      use_sparseC = TRUE)

est_coef_htlgmm<-res_htlgmm$beta
names(est_coef_htlgmm)<-c(paste0('Z',1:pZ),paste0('W',1:pW))

est_coef_htlgmm<-res_htlgmm$beta
ee_htlgmm<-sum((coef-est_coef_htlgmm)^2)
ee_htlgmm<-round(ee_htlgmm,4)

ee_lasso<-sum((coef-est_coef_glmnet)^2)
ee_lasso<-round(ee_lasso,4)
print(paste0("Estimation Error: ","lasso(",ee_lasso,"); htlgmm(",ee_htlgmm,")"))

```

```R
[1] "Estimation Error: lasso(0.0414); htlgmm(0.0212)"
```
