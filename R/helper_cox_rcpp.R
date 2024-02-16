# https://medium.com/analytics-vidhya/implementing-the-cox-model-in-r-b1292d6ab6d2
# https://cran.r-project.org/web/packages/coxed/vignettes/simulating_survival_data.html
# library(survival)
# n <- 1000
# t <- rexp(100)
# c <- rbinom(100, 1, .2) ## censoring indicator (independent process)
# x <- rbinom(100, 1, exp(-t)) ## some arbitrary relationship btn x and t
# x1 <- rbinom(100, 1, exp(-t)) ## some arbitrary relationship btn x and t
# x2 <- rbinom(100, 1, exp(-t)) ## some arbitrary relationship btn x and t
# X = cbind(x,x1,x2)
# XR = cbind(x,x1)
# events=c
# times=t
# beta = c(0.5,1,2)
# hat_theta = c(0.3,0.9)
# C_half=diag(1,nrow=ncol(X)+ncol(XR))
# betamax <- coxph(Surv(t, c) ~ x)
# betamax
# beta1 <- coxph(Surv(t, c) ~ x, init = c(1), control=list('iter.max'=0))

riskset <- function(t,times,events=NULL,entry=NULL) {
    if(is.null(events)){events=rep(1,length(times))}
    if(is.null(entry)){
        return(which((times==t & events==1)|(times>t)))
    }else{
        return(which((t>=entry) &
                         ((times==t & events==1)|(times>t))))
    }
}

cox_gradient<-function(beta,
                       times,
                       X,
                       events=NULL,
                       entry=NULL){
    if(is.null(events)){events=rep(1,length(times))}
    gradient <- rowSums(sapply(which(events==1), function(i){
        ts = times[i]
        xs = X[i,]
        x_riskset = X[riskset(ts,times,events,entry),,drop=F]
        if(nrow(x_riskset)>0){
            exp_x_beta = exp(x_riskset%*%beta)
            prop=crossprod(x_riskset,exp_x_beta)/sum(exp_x_beta)
        }else{
            prop=0
        }
        -xs + prop
    }))
    return(gradient)
}
cox_hessian <- function(beta,
                        times,
                        X,
                        events=NULL,
                        entry=NULL) {

    hessian <- lapply(1:length(times), function(i){
        ts = times[i]
        xs = X[i,]
        x_riskset = X[riskset(ts,times,events,entry),,drop=F]
        if(nrow(x_riskset)>0){
            exp_x_beta = exp(x_riskset%*%beta)
            sum_exp = sum(exp_x_beta)
            mat1 = crossprod(x_riskset*c(exp_x_beta),x_riskset)
            mat2 = crossprod(x_riskset,exp_x_beta)
            mat2 = mat2%*%t(mat2)
            hes = mat1/sum_exp - mat2/sum_exp^2
        }else{
            hes = matrix(0,nrow=ncol(X),ncol=ncol(X))
        }
        hes
    })
    Reduce('+',hessian)
}

#cox_gradient(c(1,1),Ts,Xs,events = c)
#cox_hessian(c(1,1),Ts,Xs,events = c)

pseudo_Xy_cox = function(C_half,Z,W,A,times,events,
                    beta,hat_thetaZ,
                    hat_thetaA=NULL,
                    entry=NULL){
    X=cbind(A,Z,W)
    XR=cbind(A,Z)
    beta=matrix(beta,ncol=1)
    pX = ncol(X)
    pXR = ncol(XR)
    gradhessian <- lapply(which(events==1), function(i){
        ts = times[i]
        cur_riskset=riskset(ts,times,events,entry)
        x_riskset = X[cur_riskset,,drop=F]
        xr_riskset = XR[cur_riskset,,drop=F]
        hat_theta = c(hat_thetaA,hat_thetaZ)
        if(nrow(x_riskset)>0){
            exp_x_beta = exp(prod_rcpp(x_riskset,beta))
            mat2 = crossprod_rcpp(x_riskset,exp_x_beta)
            mat3 = crossprod_rcpp(xr_riskset,exp_x_beta)
            exp_x_beta = c(exp_x_beta)
            sum_exp = sum(exp_x_beta)
            mat1 = crossprod_rcpp(cbind(x_riskset,xr_riskset),exp_x_beta*x_riskset)
            mat2 = prod_rcpp(rbind(mat2,mat3),t(mat2))
            hes = mat1/sum_exp - mat2/sum_exp^2

            prop=crossprodv_rcpp(x_riskset,exp_x_beta)/sum_exp
            exp_xr_theta = exp(prodv_rcpp(xr_riskset,hat_theta))
            prop1=crossprodv_rcpp(xr_riskset,exp_xr_theta)/sum(exp_xr_theta)
            prop = c(prop,prop1)
            hes = cbind(hes,prop)
        }else{
            hes = matrix(0,nrow=pX+pXR,ncol=pX+1)
        }
        hes
    })
    gradhessian = Reduce('+',gradhessian)

    ps_X = prod_rcpp(C_half,gradhessian[,-c(pX+1)])

    # grad1 <- rowSums(sapply(1:length(times), function(i){
    #     ts = times[i]
    #     cur_riskset=riskset(ts,times,events,entry)
    #     x_riskset = X[cur_riskset,,drop=F]
    #     xr_riskset = XR[cur_riskset,,drop=F]
    #     if(nrow(x_riskset)>0){
    #         exp_x_beta = exp(x_riskset%*%beta)
    #         prop=crossprod(x_riskset,exp_x_beta)/sum(exp_x_beta)
    #         exp_xr_theta = exp(xr_riskset%*%theta)
    #         prop1=crossprod(xr_riskset,exp_xr_theta)/sum(exp_xr_beta)
    #         prop = c(prop,prop1)
    #     }else{
    #         prop=rep(0,ncol(X)+ncol(XR))
    #     }
    #     prop
    # }))
    grad0 = -crossprodv_rcpp(X,events)
    grad1 = gradhessian[1:ncol(X),(pX+1)]
    grad2 = -gradhessian[-c(1:ncol(X)),(pX+1)]
    grad3 = gradhessian[c(1:ncol(XR)),(pX+1)]
    grad = c(grad0+grad1,grad2+grad3)
    ps_y = c(prod_rcpp(ps_X,beta) - prodv_rcpp(C_half,grad))
    list("pseudo_X"=ps_X,"pseudo_y"=ps_y)
}

cv_cox_lambda_func<-function(index_fold,Z,W,A,times,events,
                             C_half,beta_initial,lambda_list,
                             hat_thetaZ,hat_thetaA=NULL,
                             type_measure = "deviance",
                             w_adaptive=NULL,final_alpha=1){
    fold_cox_lambda<-sapply(1:length(index_fold), function(cur_fold){
        index_test<-index_fold[[cur_fold]]
        Ztrain<-Z[-index_test,,drop=FALSE]
        Ztest<-Z[index_test,,drop=FALSE]
        if(!is.null(W)){
            Wtrain<-W[-index_test,,drop=FALSE]
            Wtest<-W[index_test,,drop=FALSE]
        }else{
            Wtrain<-NULL
            Wtest<-NULL}
        if(!is.null(A)){
            Atrain<-A[-index_test,,drop=FALSE]
            Atest<-A[index_test,,drop=FALSE]
        }else{
            Atrain<-NULL
            Atest<-NULL}
        timestrain<-times[-index_test]
        timestest<-times[index_test]
        eventstrain<-events[-index_test]
        eventstest<-events[index_test]

        psXy<-pseudo_Xy_cox(C_half,Ztrain,Wtrain,Atrain,timestrain,eventstrain,
                                            beta_initial,hat_thetaZ=hat_thetaZ,
                                            hat_thetaA=hat_thetaA)
        initial_sf<-nrow(Z)/sqrt(nrow(psXy$pseudo_X))
        pseudo_X_train = psXy$pseudo_X/initial_sf
        pseudo_y_train = psXy$pseudo_y/initial_sf
        pseudo_X_train<<-pseudo_X_train
        pseudo_y_train<<-pseudo_y_train
        cox_lam<-sapply(lambda_list,function(cur_lam){
            cv_fit<-glmnet(x= (pseudo_X_train),y= (pseudo_y_train),
                           standardize=F,intercept=F,
                           alpha = final_alpha,lambda = cur_lam)#penalty.factor = w_adaptive,
            cur_beta<-coef.glmnet(cv_fit)[-1]
            if(type_measure == "C"){
            tmp=Cindex(pred=prodv_rcpp(cbind(Atest,Ztest,Wtest),cur_beta),y=cbind(time=timestest,status=eventstest))
              #tmp=iAUC(timestest,eventstest,prodv_rcpp(cbind(Atest,Ztest,Wtest),cur_beta))
            }else{
              tmp=partial_likelihood(timestest,eventstest,cbind(Atest,Ztest,Wtest),cur_beta)
            }
            tmp
        })
        cox_lam
    })
    rowSums(fold_cox_lambda)
}


iAUC<-function(times,events,predict_risk){
    time_points <- seq(min(times), max(times), by = (max(times)-min(times))/length(times)/2)
    roc_results <- timeROC::timeROC(times,events,predict_risk,cause = 1,times = time_points)
    mean(roc_results$AUC,na.rm = T)
}

partial_likelihood<-function(times,events,X,beta,entry=NULL){
    #fit1 <- coxph(y ~.,data= data.frame(y=Surv(times, events),X), init = beta, iter.max=0)
    #fit1$loglik
    likelihood1=sum(X[events==1,]%*%beta)
    likelihood2=sum(sapply(which(events==1), function(i){
        ts = times[i]
        x_riskset = X[riskset(ts,times,events,entry),,drop=F]
        if(nrow(x_riskset)>0){
            exp_x_beta = exp(prodv_rcpp(x_riskset,beta))
            prop=log(sum(exp_x_beta))
        }else{
            prop=0
        }
        prop
    }))
    likelihood1-likelihood2
}

library(caret)
cv.htlgmm.cox<-function(times,events,Z,W,hat_thetaZ,
                        A=NULL,
                        hat_thetaA=NULL,
                        beta_initial=NULL,
                        type_measure = "deviance",
                        C_half=NULL,
                        nopenalty = FALSE,
                        nfolds = 10,
                        lambda_list = NULL,
                        nlambda = 100,
                        lambda.min.ratio = 0.0001){
    psXy = pseudo_Xy_cox(C_half,Z,W,A,times,events,
                         beta_initial,hat_thetaZ=hat_thetaZ,
                         hat_thetaA=NULL)
    initial_sf<-nrow(Z)/sqrt(nrow(psXy$pseudo_X))
    psX = psXy$pseudo_X/initial_sf
    psy = psXy$pseudo_y/initial_sf
    beta0=prodv_rcpp(choinv_rcpp(self_crossprod_rcpp(psX)),crossprodv_rcpp(psX,psy))
    if(nopenalty){return(list("beta_nopenalty"=beta0))}
    if(is.null(lambda_list)){
        fit_lambda<-glmnet(x= psX,y= psy,standardize=F,intercept=F)
        lambda.max<-fit_lambda$lambda[1]
        lambda_list <-exp(seq(log(lambda.max),log(lambda.max*lambda.min.ratio),
                              length.out=nlambda))
    }
    index_fold = createFolds(events,k = nfolds)
    cv_res=cv_cox_lambda_func(index_fold,Z,W,A,times,events,
                              C_half,beta_initial,lambda_list,
                              hat_thetaZ,hat_thetaA,type_measure)
    final.lambda.min=lambda_list[which.max(cv_res)]
    fit_final=glmnet(x= psX,y= psy,standardize=F,
                     intercept=F,lambda = final.lambda.min)
    beta=coef(fit_final)[-1]
    return_list<-list("beta"=beta,
                      "lambda_list"=lambda_list,
                      "lambda_min"=final.lambda.min,
                      "cv_dev"=cv_res,
                      "beta_nopenalty"=beta0
    )
}


surv.simulate<-function(foltime,
                        dist.ev, anc.ev, beta0.ev,
                        dist.cens, anc.cens, beta0.cens,
                        beta, X, z=NULL){
    nX=nrow(X)
    if (!is.null(z) && z[[1]][1] == "gamma")    az1 <- rgamma(1, as.numeric(z[[1]][2]), as.numeric(z[[1]][3]))
    if (!is.null(z) && z[[1]][1] == "exp")      az1 <- rgamma(1, 1, as.numeric(z[[1]][2]))
    if (!is.null(z) && z[[1]][1] == "weibull")  az1 <- rweibull(1, as.numeric(z[[1]][2]), as.numeric(z[[1]][3]))
    if (!is.null(z) && z[[1]][1] == "unif")     az1 <- runif(1, as.numeric(z[[1]][2]), as.numeric(z[[1]][3]))
    if (!is.null(z) && z[[1]][1] == "invgauss") az1 <- rinvgauss(1, as.numeric(z[[1]][2]), as.numeric(z[[1]][3]))
    if (is.null(z))                             az1 <- 1

    # Time to censorship
    if (dist.cens == "llogistic") {
        tc <- exp(rllogis(nX, beta0.cens, anc.cens))
    }
    else {
        if (dist.cens == "weibull") {
            a.cens <- anc.cens
            b.cens <- (exp(-anc.cens * (beta0.cens)))^(1/anc.cens)
            tc <- rweibull(nX, a.cens, b.cens)
        }else {
            if (dist.cens == "lnorm") {
                tc <- rlnorm(nX, beta0.cens, anc.cens)
            }else {
                if (dist.cens== "unif") {
                    tc <- runif(nX, beta0.cens, anc.cens)
                }
            }
        }
    }

    suma <- X%*%beta
    if (dist.ev == 'llogistic'){
        tb <- az1*exp(rllogis(nX, beta0.ev + suma, anc.ev))
    }else{
        if (dist.ev == 'weibull'){
            a.ev   <- anc.ev
            b.ev   <- (exp(-anc.ev*(beta0.ev + suma)))^(1/anc.ev)
            tb  <- az1*rweibull(nX, a.ev, b.ev)
        }else{
            if (dist.ev == 'lnorm'){
                tb  <- az1*rlnorm(nX, beta0.ev + suma, anc.ev)
            } #if
        } #if
    } #if
    it = as.numeric(tb < tc)
    start = rep(0,nX)
    stop = pmin(tb,tc)
    it[stop > foltime] = 0
    stop[stop > foltime]  = foltime
    sim.ind <- data.frame(status=it, start=start, stop=stop, z=az1)
    return(sim.ind)
}

##############
pZ=20 # overlapping features
pW=3 # unmatched features
coef<-c(rep(0,pZ+pW))
coef[1:3]<-0.2
coef[c(pZ+1,pZ+2)]<-0.2
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

eta = cbind(Z_joint,W_joint)%*%coef
base_hazard <- 0.1  # Example baseline hazard
survival_times <- -log(runif(n_joint)) / (base_hazard * exp(eta))

# Add censoring
censoring_times <- rexp(n_joint, rate = 0.05)
event_joint <- survival_times <= censoring_times
observed_times_joint <- pmin(survival_times, censoring_times)

event = event_joint[main_index]
eventE = event_joint[-main_index]
observed_times = observed_times_joint[main_index]
observed_timesE = observed_times_joint[-main_index]

library(survival)

# Create a survival object
surv_data <- Surv(time = observed_times, event = event)

# Fit Cox model
fit <- coxph(surv_data ~., data=data.frame(surv_data,Z,W))
beta_initial=fit$coefficients

#iAUC(observed_times,event,cbind(Z,W)%*%beta_final)

#betas =fit_lambda$glmnet.fit$beta
#iAUC_lambda = sapply(1:ncol(betas), function(i){
#    iAUC(observed_times,event,cbind(Z,W)%*%betas[,i])
#})
#beta_final2=betas[,which.max(iAUC_lambda)]

surv_dataE <- Surv(time = observed_timesE, event = eventE)
fitE <- coxph(surv_dataE ~., data=data.frame(surv_dataE,ZE))
hat_thetaZ=fitE$coefficients

library(glmnet)
glmnet_res =cv.glmnet(x=cbind(Z,W),y=surv_data,family="cox",alpha=1)
beta_initial_lasso=as.vector(coef(glmnet_res,s="lambda.min"))

glmnet_res =cv.glmnet(x=cbind(Z,W),y=surv_data,family="cox",alpha=0)
beta_initial_ridge=as.vector(coef(glmnet_res,s="lambda.min"))


C_half = magic::adiag(diag(1,nrow=ncol(Z)+ncol(W)),diag(1,nrow=ncol(Z)))

htlgmm_res = cv.htlgmm.cox(observed_times,event,Z,W,hat_thetaZ,beta_initial = beta_initial,C_half =C_half,type_measure = "C")
htlgmm_res2 = cv.htlgmm.cox(observed_times,event,Z,W,hat_thetaZ,beta_initial = beta_initial_ridge,C_half =C_half)

print("-------- Sum of Square between true coefficient and estimators -------------")
print(paste0("Main study coxph :",round(sum((beta_initial-coef)^2),4)))
print(paste0("Main study coxph lasso :",round(sum((beta_initial_lasso-coef)^2),4)))
print(paste0("Main study coxph ridge :",round(sum((beta_initial_ridge-coef)^2),4)))
print(paste0("HTLGMM no penalty: coxph as initial :",round(sum((htlgmm_res$beta_nopenalty-coef)^2),4)))
print(paste0("HTLGMM lasso: coxph as initial :",round(sum((htlgmm_res$beta-coef)^2),4)))
print(paste0("HTLGMM lasso: coxph ridge as initial :",round(sum((htlgmm_res2$beta-coef)^2),4)))

#################################
## from package

### Load package
library("survsim")

### 3. Simple survival data (Section 3.3, pp. 5-6)
# dist.ev <- "weibull"
# anc.ev <- 1
# beta0.ev <- 5.268
# dist.cens <- "weibull"
# anc.cens <- 1
# beta0.cens <- 5.368
# x <- list(c("bern", 0.3), c("bern", 0.4))
# beta <- list(-0.4, -0.25)
# set.seed(11092014)
# simple.dat <- simple.surv.sim(300, 365, dist.ev, anc.ev, beta0.ev, dist.cens, anc.cens,
#   beta0.cens, , beta, x)

pZ=20 # overlapping features
pW=3 # unmatched features
coef<-c(rep(0,pZ+pW))
coef[1:3]<-0.2
coef[c(pZ+1,pZ+2)]<-0.2
which(coef!=0)
n=400
nE=2000
n_joint=n+nE
main_index<-1:n
xlist <- lapply(1:(pZ+pW), function(i){c("normal",0,1)})
#xlist <- list(c("normal",0,1), c("bern", 0.4))
#coeflist <- list(-0.4, -0.25)
coeflist <- lapply(coef, function(i){i})
set.seed(11092014)
simple.dat <- simple.surv.sim(n=n_joint,foltime = 365,dist.ev =  "weibull", anc.ev = 1, beta0.ev = 5.268, dist.cens = "weibull", anc.cens = 1,
                              beta0.cens = 5.368,beta = coeflist, x = xlist)
Z_joint<-sapply(c("x",paste0("x.",1:(pZ-1))), function(i){simple.dat[[i]]})
W_joint<-sapply(paste0("x.",pZ:(pZ+pW-1)), function(i){simple.dat[[i]]})
colnames(Z_joint)<-paste0("Z",1:pZ)
colnames(W_joint)<-paste0("W",1:pW)
Z=Z_joint[main_index,]  # separate main and external study for Z
ZE=Z_joint[-main_index,]

W=W_joint[main_index,]
observed_times_joint = simple.dat$stop-simple.dat$start
event_joint =simple.dat$status



event = event_joint[main_index]
eventE = event_joint[-main_index]
observed_times = observed_times_joint[main_index]
observed_timesE = observed_times_joint[-main_index]

library(survival)

# Create a survival object
surv_data <- Surv(time = observed_times, event = event)

# Fit Cox model
fit <- coxph(surv_data ~., data=data.frame(surv_data,Z,W))
beta_initial=fit$coefficients

surv_dataE <- Surv(time = observed_timesE, event = eventE)
fitE <- coxph(surv_dataE ~., data=data.frame(surv_dataE,ZE))
hat_thetaZ=fitE$coefficients

library(glmnet)
glmnet_res =cv.glmnet(x=cbind(Z,W),y=surv_data,family="cox",alpha=1)
beta_initial_lasso=as.vector(coef(glmnet_res,s="lambda.min"))

glmnet_res =cv.glmnet(x=cbind(Z,W),y=surv_data,family="cox",alpha=0)
beta_initial_ridge=as.vector(coef(glmnet_res,s="lambda.min"))


C_half = magic::adiag(diag(1,nrow=ncol(Z)+ncol(W)),diag(1,nrow=ncol(Z)))

htlgmm_res = cv.htlgmm.cox(observed_times,event,Z,W,hat_thetaZ,beta_initial = beta_initial,C_half =C_half,type_measure = "deviance")
htlgmm_res2 = cv.htlgmm.cox(observed_times,event,Z,W,hat_thetaZ,beta_initial = beta_initial_ridge,C_half =C_half)

print("-------- Sum of Square between true coefficient and estimators -------------")
print(paste0("Main study coxph :",round(sum((-beta_initial-coef)^2),4)))
print(paste0("Main study coxph lasso :",round(sum((-beta_initial_lasso-coef)^2),4)))
print(paste0("Main study coxph ridge :",round(sum((-beta_initial_ridge-coef)^2),4)))
print(paste0("HTLGMM no penalty: coxph as initial :",round(sum((-htlgmm_res$beta_nopenalty-coef)^2),4)))
print(paste0("HTLGMM lasso: coxph as initial :",round(sum((-htlgmm_res$beta-coef)^2),4)))
print(paste0("HTLGMM lasso: coxph ridge as initial :",round(sum((-htlgmm_res2$beta-coef)^2),4)))


########   Prepare for correlation structure


### 3. Simple survival data (Section 3.3, pp. 5-6)
# dist.ev <- "weibull"
# anc.ev <- 1
# beta0.ev <- 5.268
# dist.cens <- "weibull"
# anc.cens <- 1
# beta0.cens <- 5.368
# x <- list(c("bern", 0.3), c("bern", 0.4))
# beta <- list(-0.4, -0.25)
# set.seed(11092014)
# simple.dat <- simple.surv.sim(300, 365, dist.ev, anc.ev, beta0.ev, dist.cens, anc.cens,
#   beta0.cens, , beta, x)

pZ=20 # overlapping features
pW=3 # unmatched features
coef<-c(rep(0,pZ+pW))
coef[1:3]<- -0.5
coef[c(pZ+1,pZ+2)]<- -0.5
which(coef!=0)
n=400
nE=2000
n_joint=n+nE
main_index<-1:n
xlist <- lapply(1:(pZ+pW), function(i){c("normal",0,1)})
#xlist <- list(c("normal",0,1), c("bern", 0.4))
#coeflist <- list(-0.4, -0.25)
coeflist <- lapply(coef, function(i){i})
set.seed(11092014)
#simple.dat <- simple.surv.sim(n=n_joint,foltime = 365,dist.ev =  "weibull", anc.ev = 1, beta0.ev = 5.268, dist.cens = "weibull", anc.cens = 1,
#                              beta0.cens = 5.368,beta = coeflist, x = xlist)
#Z_joint<-sapply(c("x",paste0("x.",1:(pZ-1))), function(i){simple.dat[[i]]})
#W_joint<-sapply(paste0("x.",pZ:(pZ+pW-1)), function(i){simple.dat[[i]]})
Z_joint<-matrix(rnorm(n_joint*pZ),n_joint,pZ)
colnames(Z_joint)<-paste0("Z",1:pZ)
W_joint<-matrix(rnorm(n_joint*pW),n_joint,pW)
colnames(W_joint)<-paste0("W",1:pW)
X_joint<-cbind(Z_joint,W_joint)
# simple.ev.dat = sapply(1:n_joint, function(i){
#     simple.ev.dat <- simple.ev.sim(foltime = 365,dist.ev =  "weibull", anc.ev = 1, beta0.ev = 5.268, dist.cens = "weibull", anc.cens = 1,
#                                    beta0.cens = 5.368,beta = coeflist, eff = X_joint[i,],i=i)
#     c(simple.ev.dat$status,simple.ev.dat$start,simple.ev.dat$stop)
# })
# observed_times_joint = simple.ev.dat[3,]-simple.ev.dat[2,]
# event_joint =simple.ev.dat[1,]
simple.ev.dat<-surv.simulate(foltime = 365,
              dist.ev =  "weibull", anc.ev = 1, beta0.ev = 5.268,
              dist.cens = "weibull", anc.cens = 1,beta0.cens = 5.368,
              beta = coef, X = X_joint)
observed_times_joint = simple.ev.dat$stop-simple.ev.dat$start
event_joint = simple.ev.dat$status
Z=Z_joint[main_index,]  # separate main and external study for Z
ZE=Z_joint[-main_index,]
W=W_joint[main_index,]

event = event_joint[main_index]
eventE = event_joint[-main_index]
observed_times = observed_times_joint[main_index]
observed_timesE = observed_times_joint[-main_index]

library(survival)

# Create a survival object
surv_data <- Surv(time = observed_times, event = event)

# Fit Cox model
fit <- coxph(surv_data ~., data=data.frame(surv_data,Z,W))
beta_initial=fit$coefficients

surv_dataE <- Surv(time = observed_timesE, event = eventE)
fitE <- coxph(surv_dataE ~., data=data.frame(surv_dataE,ZE))
hat_thetaZ=fitE$coefficients

library(glmnet)
glmnet_res =cv.glmnet(x=cbind(Z,W),y=surv_data,family="cox",alpha=1)
beta_initial_lasso=as.vector(coef(glmnet_res,s="lambda.min"))

glmnet_res =cv.glmnet(x=cbind(Z,W),y=surv_data,family="cox",alpha=0)
beta_initial_ridge=as.vector(coef(glmnet_res,s="lambda.min"))


C_half = magic::adiag(diag(1,nrow=ncol(Z)+ncol(W)),diag(1,nrow=ncol(Z)))

#htlgmm_res = cv.htlgmm.cox(observed_times,event,Z,W,hat_thetaZ,beta_initial = beta_initial,C_half =C_half,type_measure = "C")
#htlgmm_res2 = cv.htlgmm.cox(observed_times,event,Z,W,hat_thetaZ,beta_initial = beta_initial_ridge,C_half =C_half,type_measure = "C")

htlgmm_res_dev = cv.htlgmm.cox(observed_times,event,Z,W,hat_thetaZ,beta_initial = beta_initial,C_half =C_half,type_measure = "deviance")
htlgmm_res_dev2 = cv.htlgmm.cox(observed_times,event,Z,W,hat_thetaZ,beta_initial = beta_initial_ridge,C_half =C_half,type_measure = "deviance")


print("-------- Sum of Square between true coefficient and estimators -------------")
print(paste0("Main study coxph :",round(sum((beta_initial-coef)^2),4)))
print(paste0("Main study coxph lasso :",round(sum((beta_initial_lasso-coef)^2),4)))
print(paste0("Main study coxph ridge :",round(sum((beta_initial_ridge-coef)^2),4)))
print(paste0("HTLGMM no penalty: coxph as initial :",round(sum((htlgmm_res$beta_nopenalty-coef)^2),4)))
#print(paste0("HTLGMM lasso: coxph as initial; Cindex :",round(sum((htlgmm_res$beta-coef)^2),4)))
#print(paste0("HTLGMM lasso: coxph ridge as initial; Cindex :",round(sum((htlgmm_res2$beta-coef)^2),4)))
print(paste0("HTLGMM lasso: coxph as initial; deviance :",round(sum((htlgmm_res_dev$beta-coef)^2),4)))
print(paste0("HTLGMM lasso: coxph ridge as initial; deviance :",round(sum((htlgmm_res_dev2$beta-coef)^2),4)))



########### One parameter settings

pZ=1
pW=1
coef<-c(rep(0,pZ+pW))
coef<- c(-0.2,-0.2)
n=400
nE=2000
n_joint=n+nE
main_index<-1:n
Z_joint<-matrix(rnorm(n_joint*pZ),n_joint,pZ)
colnames(Z_joint)<-paste0("Z",1:pZ)
W_joint<-matrix(rnorm(n_joint*pW),n_joint,pW)
colnames(W_joint)<-paste0("W",1:pW)
X_joint<-cbind(Z_joint,W_joint)
simple.ev.dat<-surv.simulate(foltime = 365,
                             dist.ev =  "weibull", anc.ev = 1, beta0.ev = 5.268,
                             dist.cens = "weibull", anc.cens = 1,beta0.cens = 5.368,
                             beta = coef, X = X_joint)
observed_times_joint = simple.ev.dat$stop-simple.ev.dat$start
event_joint = simple.ev.dat$status
Z=Z_joint[main_index,,drop=F]  # separate main and external study for Z
ZE=Z_joint[-main_index,,drop=F]
W=W_joint[main_index,,drop=F]

event = event_joint[main_index]
eventE = event_joint[-main_index]
observed_times = observed_times_joint[main_index]
observed_timesE = observed_times_joint[-main_index]

library(survival)
surv_data <- Surv(time = observed_times, event = event)
fit <- coxph(surv_data ~., data=data.frame(surv_data,Z,W))
beta_initial=fit$coefficients
surv_dataE <- Surv(time = observed_timesE, event = eventE)
fitE <- coxph(surv_dataE ~., data=data.frame(surv_dataE,ZE))
hat_thetaZ=fitE$coefficients
C_half = magic::adiag(diag(1,nrow=pZ+pW),diag(1,nrow=pZ))

htlgmm_res = cv.htlgmm.cox(observed_times,event,Z,W,hat_thetaZ,beta_initial = beta_initial,nopenalty = T,C_half =C_half)
beta_htlgmm = htlgmm_res$beta_nopenalty
bias_initial = beta_initial - coef
bias_initial
bias_htlgmm = beta_htlgmm - coef
bias_htlgmm
print(paste0("Main study coxph :",round(sum((beta_initial-coef)^2),4)))
print(paste0("HTLGMM no penalty: coxph as initial :",round(sum((beta_htlgmm-coef)^2),4)))



pZ=20 # overlapping features
pW=3 # unmatched features
coef<-c(rep(0,pZ+pW))
coef[1:3]<- 0.5
coef[c(pZ+1,pZ+2)]<- 0.5
which(coef!=0)
n=400
nE=2000
n_joint=n+nE
main_index<-1:n

set.seed(11092014)
Z_joint<-matrix(rnorm(n_joint*pZ),n_joint,pZ)
colnames(Z_joint)<-paste0("Z",1:pZ)
W_joint<-matrix(rnorm(n_joint*pW),n_joint,pW)
colnames(W_joint)<-paste0("W",1:pW)
X_joint<-cbind(Z_joint,W_joint)
covs <- data.frame(id = 1:n_joint, X_joint)
names(coef)=colnames(X_joint)
s1 <- simsurv::simsurv(dist = c("weibull"),lambdas = 0.1, gammas = 1.5,
                       betas = coef,x = covs, maxt = 5)
observed_times_joint = s1$eventtime
event_joint = s1$status
Z=Z_joint[main_index,]  # separate main and external study for Z
ZE=Z_joint[-main_index,]
W=W_joint[main_index,]

event = event_joint[main_index]
eventE = event_joint[-main_index]
observed_times = observed_times_joint[main_index]
observed_timesE = observed_times_joint[-main_index]


library(survival)

# Create a survival object
surv_data <- Surv(time = observed_times, event = event)

# Fit Cox model
fit <- coxph(surv_data ~., data=data.frame(surv_data,Z,W))
beta_initial=fit$coefficients

surv_dataE <- Surv(time = observed_timesE, event = eventE)
fitE <- coxph(surv_dataE ~., data=data.frame(surv_dataE,ZE))
hat_thetaZ=fitE$coefficients

library(glmnet)
glmnet_res =cv.glmnet(x=cbind(Z,W),y=surv_data,family="cox",alpha=1)
beta_initial_lasso=as.vector(coef(glmnet_res,s="lambda.min"))

glmnet_res =cv.glmnet(x=cbind(Z,W),y=surv_data,family="cox",alpha=0)
beta_initial_ridge=as.vector(coef(glmnet_res,s="lambda.min"))


C_half = magic::adiag(diag(1,nrow=ncol(Z)+ncol(W)),diag(1,nrow=ncol(Z)))

#htlgmm_res = cv.htlgmm.cox(observed_times,event,Z,W,hat_thetaZ,beta_initial = beta_initial,C_half =C_half,type_measure = "C")
#htlgmm_res2 = cv.htlgmm.cox(observed_times,event,Z,W,hat_thetaZ,beta_initial = beta_initial_ridge,C_half =C_half,type_measure = "C")

htlgmm_res_dev = cv.htlgmm.cox(observed_times,event,Z,W,hat_thetaZ,beta_initial = beta_initial,C_half =C_half,type_measure = "deviance")
htlgmm_res_dev2 = cv.htlgmm.cox(observed_times,event,Z,W,hat_thetaZ,beta_initial = beta_initial_ridge,C_half =C_half,type_measure = "deviance")


print("-------- Sum of Square between true coefficient and estimators -------------")
print(paste0("Main study coxph :",round(sum((beta_initial-coef)^2),4)))
print(paste0("Main study coxph lasso :",round(sum((beta_initial_lasso-coef)^2),4)))
print(paste0("Main study coxph ridge :",round(sum((beta_initial_ridge-coef)^2),4)))
print(paste0("HTLGMM no penalty: coxph as initial :",round(sum((htlgmm_res_dev$beta_nopenalty-coef)^2),4)))
print(paste0("HTLGMM lasso: coxph as initial; deviance :",round(sum((htlgmm_res_dev$beta-coef)^2),4)))
print(paste0("HTLGMM lasso: coxph ridge as initial; deviance :",round(sum((htlgmm_res_dev2$beta-coef)^2),4)))






if(0){

    bias_vec=c()
    for( iii in 1:100){
        message(iii)
        set.seed(iii)
        pZ=1
        pW=1
        coef<-c(rep(0,pZ+pW))
        coef<- c(0.5,0.5)
        n=400
        nE=2000
        n_joint=n+nE
        main_index<-1:n
        Z_joint<-matrix(rnorm(n_joint*pZ),n_joint,pZ)
        colnames(Z_joint)<-paste0("Z",1:pZ)
        W_joint<-matrix(rnorm(n_joint*pW),n_joint,pW)
        W_joint = W_joint+0.3*Z_joint
        W_joint = scale(W_joint)
        colnames(W_joint)<-paste0("W",1:pW)
        X_joint<-cbind(Z_joint,W_joint)
        covs <- data.frame(id = 1:n_joint, X_joint)
        names(coef)=colnames(X_joint)
        s1 <- simsurv::simsurv(dist = c("weibull"),lambdas = 0.1, gammas = 1.5,
                               betas = coef,x = covs, maxt = 5)
        observed_times_joint = s1$eventtime
        event_joint = s1$status
        Z=Z_joint[main_index,,drop=F]  # separate main and external study for Z
        ZE=Z_joint[-main_index,,drop=F]
        W=W_joint[main_index,,drop=F]

        event = event_joint[main_index]
        eventE = event_joint[-main_index]
        observed_times = observed_times_joint[main_index]
        observed_timesE = observed_times_joint[-main_index]

        library(survival)
        surv_data <- Surv(time = observed_times, event = event)
        fit <- coxph(surv_data ~., data=data.frame(surv_data,Z,W))
        beta_initial=fit$coefficients
        surv_dataE <- Surv(time = observed_timesE, event = eventE)
        fitE <- coxph(surv_dataE ~., data=data.frame(surv_dataE,ZE))
        hat_thetaZ=fitE$coefficients
        C_half = magic::adiag(diag(1,nrow=pZ+pW),diag(1,nrow=pZ))

        htlgmm_res = cv.htlgmm.cox(observed_times,event,Z,W,hat_thetaZ,beta_initial = beta_initial,nopenalty = T,C_half =C_half)
        beta_htlgmm = htlgmm_res$beta_nopenalty
        bias_initial = beta_initial - coef
        bias_htlgmm = beta_htlgmm - coef
        bias_vec<-rbind(bias_vec,c(bias_initial,bias_htlgmm))
    }

}
