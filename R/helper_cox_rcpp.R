if(0){
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

## functions to compute riskset
right_id<-function(times){
    right_equal_id = c(NA,length(times))
    record_i = length(times)
    record_time = times[record_i]
    for(i in length(times):1){
        if (times[i]<record_time){
            record_time=times[i]
            record_i=i
        }
        right_equal_id[i]=record_i
    }
    right_equal_id
}

left_id<-function(times){
    left_equal_id = c(NA,length(times))
    record_i = 1
    record_time = times[1]
    for(i in 1:length(times)){
        if (times[i]>record_time){
            record_time=times[i]
            record_i=i
        }
        left_equal_id[i]=record_i
    }
    left_equal_id
}

reorder_U1U2U3 = function(Z,W,A,times,events,
                          beta,tilde_thetaZ,
                          hat_thetaA,
                          left_equal_id,
                          right_equal_id,X=NULL,XR=NULL){
    if(is.null(X)){X=cbind(A,Z,W)}
    if(is.null(XR)){XR=cbind(A,Z)}
    pX = ncol(X)
    if(!is.null(A)){pA=ncol(A)}else{pA=0}
    if(!is.null(W)){pW=ncol(W)}else{pW=0}
    pZ=ncol(Z)
    pXR = ncol(XR) # XR is the first several cols of X
    nX = length(times)
    tilde_theta = c(hat_thetaA,tilde_thetaZ)
    # U1 section
    s01 = sapply(which(events==1), function(i){
        x_riskset = X[left_equal_id[i]:nX,,drop=F]
        exp_x_beta = exp(prodv_rcpp(x_riskset,beta))
        inv_sum_exp = 1/sum(exp_x_beta)
        s1=crossprodv_rcpp(x_riskset,exp_x_beta)
        c(inv_sum_exp,s1)
    })
    H_T_x_beta=sapply(1:nX, function(j){
        ind=sum(events[1:right_equal_id[j]])
        if(ind==0){inv_s0list=0}else{
            inv_s0list=s01[1,1:ind]}
        s1list=s01[-1,1:ind,drop=F]
        X[j,]*sum(inv_s0list)-colSums(t(s1list)*inv_s0list^2)
    })
    if(is.null(dim(H_T_x_beta)[1])){
        H_T_x_beta = matrix(H_T_x_beta,nrow=1)
    }
    exp_H = t(H_T_x_beta)*exp(prodv_rcpp(X,beta))
    delta_x_s = matrix(0,nrow=nX,ncol = ncol(X))
    delta_x_s[events==1,] = t(s01[-1,])*s01[1,]
    U1mat=-X*events+delta_x_s+exp_H

    ## General format for U2 and U3
    s01R <- sapply(which(events==1), function(i){
        xr_riskset = XR[left_equal_id[i]:nX,,drop=F]
        exp_xr_theta = exp(prodv_rcpp(xr_riskset,tilde_theta))
        inv_sum_exp_xr = 1/sum(exp_xr_theta)
        s1_r=crossprodv_rcpp(xr_riskset,exp_xr_theta)
        c(inv_sum_exp_xr,s1_r)
    })

    HR_T_x_theta<-sapply(1:nX, function(j){
        ind=sum(events[1:right_equal_id[j]])
        if(ind==0){inv_s0Rlist=0}else{
            inv_s0Rlist=s01R[1,1:ind]}
        s1Rlist=s01R[-1,1:ind,drop=F]
        XR[j,]*sum(inv_s0Rlist)-colSums(t(s1Rlist)*inv_s0Rlist^2)
    })
    if(is.null(dim(HR_T_x_theta)[1])){
        HR_T_x_theta = matrix(HR_T_x_theta,nrow=1)
    }
    exp_HR = t(HR_T_x_theta)*exp(prodv_rcpp(XR,tilde_theta))
    delta_xr_s = matrix(0,nrow=nX,ncol = ncol(XR))
    delta_xr_s[events==1,] = t(s01R[-1,])*s01R[1,]

    ## First section of U2 : truncate from above
    exp_H_XR_Z = exp_HR[,(pA+1):(pA+pZ),drop=F]
    delta_xr_z_s = delta_xr_s[,(pA+1):(pA+pZ),drop=F]

    ## Second section of U2 : truncate from U1
    exp_H_X_Z = exp_H[,(pA+1):(pA+pZ),drop=F]
    delta_x_z_s = delta_x_s[,(pA+1):(pA+pZ),drop=F]
    U2mat = -(delta_xr_z_s+exp_H_XR_Z)+delta_x_z_s+exp_H_X_Z
    Umat = cbind(U1mat,U2mat)
    return_list=list("Umat"=Umat)
    if(pA>0){
        exp_H_XR_A = exp_HR[,1:pA,drop=F]
        delta_xr_a_s = delta_xr_s[,1:pA,drop=F]
        U3mat = -A*events+delta_xr_a_s+exp_H_XR_A
        return_list=c(return_list,list("U3mat"=U3mat))
    }
    return_list
}

Delta_opt_cox_rcpp<-function(Z,W,A,times,events,
                             beta,tilde_thetaZ,
                             V_thetaZ,
                             left_equal_id,
                             right_equal_id,
                             hat_thetaA=NULL,
                             V_thetaA=NULL,
                             X=NULL,XR=NULL){
    nX=length(times)
    if(is.null(X)){X=cbind(A,Z,W)}
    if(is.null(XR)){XR=cbind(A,Z)}
    pX=ncol(X)
    pXR=ncol(XR)
    if(!is.null(A)){pA=ncol(A)}else{pA=0}
    if(!is.null(W)){pW=ncol(W)}else{pW=0}
    pZ=ncol(Z)
    tilde_theta=c(hat_thetaA,tilde_thetaZ)
    Umatlist = reorder_U1U2U3(Z,W,A,times,events,
                              beta,tilde_thetaZ,
                              hat_thetaA,
                              left_equal_id,
                              right_equal_id,X,XR)
    Umat=Umatlist$Umat
    Delta_U=(1/nX)*self_crossprod_rcpp(Umat)
    tilde_theta = as.matrix(tilde_theta,ncol=1)
    dU2U3dtheta <- lapply(which(events==1), function(i){
        xr_riskset = XR[left_equal_id[i]:nX,,drop=F]
        exp_xr_theta = exp(prod_rcpp(xr_riskset,tilde_theta))
        mat2r = crossprod_rcpp(xr_riskset,exp_xr_theta)
        exp_xr_theta = c(exp_xr_theta)
        sum_exp_xr = sum(exp_xr_theta)
        mat1r = crossprod_rcpp(xr_riskset,exp_xr_theta*xr_riskset)
        mat2r = prod_rcpp(mat2r,t(mat2r))
        hes = mat1r/sum_exp_xr - mat2r/sum_exp_xr^2
    })
    dU2U3dtheta = Reduce('+',dU2U3dtheta)/nX
    GammaZZ=dU2U3dtheta[(pA+1):(pA+pZ),(pA+1):(pA+pZ),drop=F]

    if(is.null(dim(V_thetaZ)[1])){V_thetaZ=as.matrix(V_thetaZ,nrow=pZ,ncol=pZ)}
    Delta22_theta=prod_rcpp(prod_rcpp(GammaZZ,(nX*V_thetaZ)),t(GammaZZ))
    Delta_theta = rbind(matrix(0,nrow=pX,ncol=pX+pZ),
                        cbind(matrix(0,nrow=pZ,ncol=pX),Delta22_theta))
    if(pA>0){
        U3mat=Umatlist$U3mat
        GammaZA=dU2U3dtheta[(pA+1):(pA+pZ),1:pA,drop=F]
        inv_GammaAA=choinv_rcpp(dU2U3dtheta[1:pA,1:pA,drop=F]+diag(1e-15,pA))
        DDA=prod_rcpp(U3mat,prod_rcpp(inv_GammaAA,t(GammaZA)))
        Cov_Utheta=(1/nX)*crossprod_rcpp(Umat,DDA)
        Cov_U1theta=Cov_Utheta[1:pX,,drop=F]
        Cov_U2theta=Cov_Utheta[-c(1:pX),,drop=F]
        Delta22_thetaA = prod_rcpp(prod_rcpp(GammaZA,(nX*V_thetaA)),t(GammaZA))
        +Cov_U2theta+t(Cov_U2theta)
        Delta12_thetaA = Cov_U1theta
        Delta_theta = Delta_theta+rbind(cbind(matrix(0,nrow=pX,ncol=pX),Delta12_thetaA),
                                        cbind(t(Delta12_thetaA),Delta22_thetaA))
    }
    Delta_U+Delta_theta
}


cv_cox_lambda_func<-function(index_fold,Z,W,A,times,events,
                             C_half,beta_initial,lambda_list,
                             tilde_thetaZ,hat_thetaA=NULL,
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

        left_equal_id_train = left_id(timestrain)

        psXy<-pseudo_Xy_cox(C_half,Ztrain,Wtrain,Atrain,timestrain,
                            eventstrain,beta_initial,
                            tilde_thetaZ=tilde_thetaZ,
                            hat_thetaA=hat_thetaA,
                            left_equal_id=left_equal_id_train)
        initial_sf<-nrow(Z)/sqrt(nrow(psXy$pseudo_X))
        pseudo_X_train = psXy$pseudo_X/initial_sf
        pseudo_y_train = psXy$pseudo_y/initial_sf
        cox_lam<-sapply(lambda_list,function(cur_lam){
            cv_fit<-glmnet(x= (pseudo_X_train),y= (pseudo_y_train),
                           standardize=F,intercept=F,
                           alpha = final_alpha,lambda = cur_lam)#penalty.factor = w_adaptive,
            cur_beta<-coef.glmnet(cv_fit)[-1]
            if(type_measure == "C"){
              tmp=Cindex(pred=prodv_rcpp(cbind(Atest,Ztest,Wtest),cur_beta),y=cbind(time=timestest,status=eventstest))
            }else{
              tmp=partial_likelihood(timestest,eventstest,cbind(Atest,Ztest,Wtest),cur_beta,left_equal_id_train)
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

partial_likelihood<-function(times,events,X,beta,left_equal_id){
    likelihood1=sum(X[events==1,]%*%beta)
    likelihood2=sum(sapply(which(events==1), function(i){
        x_riskset = X[left_equal_id[i]:nrow(X),,drop=F]
        exp_x_beta = exp(prodv_rcpp(x_riskset,beta))
        prop=log(sum(exp_x_beta))
        prop
    }))
    likelihood1-likelihood2
}


pseudo_Xy_cox = function(C_half,Z,W,A,times,events,
                         beta,tilde_thetaZ,
                         hat_thetaA=NULL,
                         X=NULL,XR=NULL,
                         left_equal_id=NULL){
    if(is.null(X)) X=cbind(A,Z,W)
    if(is.null(XR)) XR=cbind(A,Z)
    beta=matrix(beta,ncol=1)
    tilde_theta = c(hat_thetaA,tilde_thetaZ)
    pX = ncol(X)
    if(!is.null(A)){pA=ncol(A)}else{pA=0}
    if(!is.null(W)){pW=ncol(W)}else{pW=0}
    pZ=ncol(Z)
    pXR = ncol(XR)
    nX = nrow(X)
    gradhessian <- lapply(which(events==1), function(i){
        cur_riskset=left_equal_id[i]:nX
        x_riskset = X[cur_riskset,,drop=F]
        xr_riskset = XR[cur_riskset,,drop=F]
        z_riskset = Z[cur_riskset,,drop=F]
        exp_x_beta = exp(prod_rcpp(x_riskset,beta))
        mat2 = crossprod_rcpp(x_riskset,exp_x_beta)
        mat3 = crossprod_rcpp(z_riskset,exp_x_beta)
        exp_x_beta = c(exp_x_beta)
        sum_exp = sum(exp_x_beta)
        mat1 = crossprod_rcpp(cbind(x_riskset,z_riskset),exp_x_beta*x_riskset)
        mat2 = prod_rcpp(rbind(mat2,mat3),t(mat2))
        hes = mat1/sum_exp - mat2/sum_exp^2

        prop=crossprodv_rcpp(x_riskset,exp_x_beta)/sum_exp
        exp_xr_theta = exp(prodv_rcpp(xr_riskset,tilde_theta))
        prop1=crossprodv_rcpp(z_riskset,exp_xr_theta)/sum(exp_xr_theta)
        prop = c(prop,prop1)
        cbind(hes,prop)
    })
    gradhessian = Reduce('+',gradhessian)

    ps_X = prod_rcpp(C_half,gradhessian[,-c(pX+1)])

    grad0 = -crossprodv_rcpp(X,events) ## sum -delta_i*x_i
    grad1 = gradhessian[1:pX,(pX+1)] ## sum delta_i*sum_j exp(x_j beta)x_j / sum_j exp(x_j beta)
    grad2 = -gradhessian[-c(1:pX),(pX+1)] ## sum delta_i*sum_j exp(xr_j theta)z_j / sum_j exp(xr_j theta)
    grad3 = gradhessian[(1+pA):(pA+pZ),(pX+1)] ## sum delta_i*sum_j exp(x_j beta)z_j / sum_j exp(x_j beta)
    grad = c(grad0+grad1,grad2+grad3)
    ps_y = c(prod_rcpp(ps_X,beta) - prodv_rcpp(C_half,grad))
    list("pseudo_X"=ps_X,"pseudo_y"=ps_y)
}

library(caret)
library(survival)
library(glmnet)
cv.htlgmm.cox<-function(times,events,Z,W,tilde_thetaZ,
                        A=NULL,
                        hat_thetaA=NULL,
                        beta_initial=NULL,
                        V_thetaZ=NULL,
                        V_thetaA=NULL,
                        type_measure = "deviance",
                        C_half=NULL,
                        nopenalty = FALSE,
                        nfolds = 10,
                        lambda_list = NULL,
                        nlambda = 100,
                        lambda.min.ratio = 0.0001,
                        seed.use = 97){
    set.seed(seed.use)
    timesorder=order(times)
    times=times[timesorder]
    events=events[timesorder]
    Z=Z[timesorder,,drop=F]
    if(!is.null(W)){W=W[timesorder,,drop=F]}
    if(!is.null(A)){A=A[timesorder,,drop=F]}
    X=cbind(A,Z,W)
    XR=cbind(A,Z)

    nX = length(times)

    right_equal_id = right_id(times)
    left_equal_id = left_id(times)

    if(is.null(C_half)){
        inv_C = Delta_opt_cox_rcpp(Z,W,A,times,events,
                                   beta=beta_initial,
                                   tilde_thetaZ=tilde_thetaZ,
                                   V_thetaZ=V_thetaZ,
                                   left_equal_id,
                                   right_equal_id,
                                   hat_thetaA=hat_thetaA,
                                   V_thetaA=V_thetaA,
                                   X=X,XR=XR)
        C_half<-sqrtchoinv_rcpp(inv_C+diag(1e-15,nrow(inv_C)))
    }

    psXy = pseudo_Xy_cox(C_half,Z,W,A,times,events,
                         beta_initial,tilde_thetaZ=tilde_thetaZ,
                         hat_thetaA=hat_thetaA,X=X,XR=XR,
                         left_equal_id=left_equal_id)
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
                              tilde_thetaZ,hat_thetaA,type_measure)
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
tilde_thetaZ=fitE$coefficients

library(glmnet)
glmnet_res =cv.glmnet(x=cbind(Z,W),y=surv_data,family="cox",alpha=1)
beta_initial_lasso=as.vector(coef(glmnet_res,s="lambda.min"))

glmnet_res =cv.glmnet(x=cbind(Z,W),y=surv_data,family="cox",alpha=0)
beta_initial_ridge=as.vector(coef(glmnet_res,s="lambda.min"))


C_half = magic::adiag(diag(1,nrow=ncol(Z)+ncol(W)),diag(1,nrow=ncol(Z)))

htlgmm_res = cv.htlgmm.cox(observed_times,event,Z,W,tilde_thetaZ,beta_initial = beta_initial,C_half =C_half,type_measure = "C")
htlgmm_res2 = cv.htlgmm.cox(observed_times,event,Z,W,tilde_thetaZ,beta_initial = beta_initial_ridge,C_half =C_half)

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
tilde_thetaZ=fitE$coefficients

library(glmnet)
glmnet_res =cv.glmnet(x=cbind(Z,W),y=surv_data,family="cox",alpha=1)
beta_initial_lasso=as.vector(coef(glmnet_res,s="lambda.min"))

glmnet_res =cv.glmnet(x=cbind(Z,W),y=surv_data,family="cox",alpha=0)
beta_initial_ridge=as.vector(coef(glmnet_res,s="lambda.min"))


C_half = magic::adiag(diag(1,nrow=ncol(Z)+ncol(W)),diag(1,nrow=ncol(Z)))

htlgmm_res = cv.htlgmm.cox(observed_times,event,Z,W,tilde_thetaZ,beta_initial = beta_initial,C_half =C_half,type_measure = "deviance")
htlgmm_res2 = cv.htlgmm.cox(observed_times,event,Z,W,tilde_thetaZ,beta_initial = beta_initial_ridge,C_half =C_half)

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
tilde_thetaZ=fitE$coefficients

library(glmnet)
glmnet_res =cv.glmnet(x=cbind(Z,W),y=surv_data,family="cox",alpha=1)
beta_initial_lasso=as.vector(coef(glmnet_res,s="lambda.min"))

glmnet_res =cv.glmnet(x=cbind(Z,W),y=surv_data,family="cox",alpha=0)
beta_initial_ridge=as.vector(coef(glmnet_res,s="lambda.min"))


C_half = magic::adiag(diag(1,nrow=ncol(Z)+ncol(W)),diag(1,nrow=ncol(Z)))

#htlgmm_res = cv.htlgmm.cox(observed_times,event,Z,W,tilde_thetaZ,beta_initial = beta_initial,C_half =C_half,type_measure = "C")
#htlgmm_res2 = cv.htlgmm.cox(observed_times,event,Z,W,tilde_thetaZ,beta_initial = beta_initial_ridge,C_half =C_half,type_measure = "C")

htlgmm_res_dev = cv.htlgmm.cox(observed_times,event,Z,W,tilde_thetaZ,beta_initial = beta_initial,C_half =C_half,type_measure = "deviance")
htlgmm_res_dev2 = cv.htlgmm.cox(observed_times,event,Z,W,tilde_thetaZ,beta_initial = beta_initial_ridge,C_half =C_half,type_measure = "deviance")


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
tilde_thetaZ=fitE$coefficients
C_half = magic::adiag(diag(1,nrow=pZ+pW),diag(1,nrow=pZ))

htlgmm_res = cv.htlgmm.cox(observed_times,event,Z,W,tilde_thetaZ,beta_initial = beta_initial,nopenalty = T,C_half =C_half)
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

set.seed(2)
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
tilde_thetaZ=fitE$coefficients
V_thetaZ = fitE$var
library(glmnet)
glmnet_res =cv.glmnet(x=cbind(Z,W),y=surv_data,family="cox",alpha=1)
beta_initial_lasso=as.vector(coef(glmnet_res,s="lambda.min"))

#glmnet_res =cv.glmnet(x=cbind(Z,W),y=surv_data,family="cox",alpha=0)
#beta_initial_ridge=as.vector(coef(glmnet_res,s="lambda.min"))

C_half = magic::adiag(diag(1,nrow=ncol(Z)+ncol(W)),diag(1,nrow=ncol(Z)))

#htlgmm_res = cv.htlgmm.cox(observed_times,event,Z,W,tilde_thetaZ,beta_initial = beta_initial,C_half =C_half,type_measure = "C")
#htlgmm_res2 = cv.htlgmm.cox(observed_times,event,Z,W,tilde_thetaZ,beta_initial = beta_initial_ridge,C_half =C_half,type_measure = "C")

htlgmm_res_dev = cv.htlgmm.cox(observed_times,event,Z,W,tilde_thetaZ,beta_initial = beta_initial,C_half =C_half,type_measure = "deviance")
#htlgmm_res_dev2 = cv.htlgmm.cox(observed_times,event,Z,W,tilde_thetaZ,beta_initial = beta_initial_ridge,C_half =C_half,type_measure = "deviance")

print("-------- Sum of Square between true coefficient and estimators -------------")
print(paste0("Main study coxph :",round(sum((beta_initial-coef)^2),4)))
print(paste0("Main study coxph lasso :",round(sum((beta_initial_lasso-coef)^2),4)))
#print(paste0("Main study coxph ridge :",round(sum((beta_initial_ridge-coef)^2),4)))
print(paste0("HTLGMM no penalty:",round(sum((htlgmm_res_dev$beta_nopenalty-coef)^2),4)))
print(paste0("HTLGMM lasso:",round(sum((htlgmm_res_dev$beta-coef)^2),4)))
#print(paste0("HTLGMM lasso: coxph as initial:",round(sum((htlgmm_res_dev$beta-coef)^2),4)))
#print(paste0("HTLGMM lasso: coxph ridge as initial; deviance :",round(sum((htlgmm_res_dev2$beta-coef)^2),4)))

htlgmm_res_dev_opt = cv.htlgmm.cox(observed_times,event,Z,W,tilde_thetaZ,beta_initial = beta_initial,
                                       V_thetaZ=V_thetaZ,type_measure = "deviance")

print(paste0("HTLGMM opt no penalty:",round(sum((htlgmm_res_dev_opt$beta_nopenalty-coef)^2),4)))
print(paste0("HTLGMM opt lasso:",round(sum((htlgmm_res_dev_opt$beta-coef)^2),4)))


if(0){

    bias_vec=c()
    for( iii in 1:100){
        message(iii)
        set.seed(iii)
        pZ=1
        pW=1
        coef<-c(rep(0,pZ+pW))
        coef<- c(0.2,0.2)
        n=400
        nE=2000
        n_joint=n+nE
        main_index<-1:n
        Z_joint<-matrix(rnorm(n_joint*pZ),n_joint,pZ)
        colnames(Z_joint)<-paste0("Z",1:pZ)
        W_joint<-matrix(rnorm(n_joint*pW),n_joint,pW)
        W_joint = W_joint+0*Z_joint
        W_joint = scale(W_joint)
        colnames(W_joint)<-paste0("W",1:pW)
        X_joint<-cbind(Z_joint,W_joint)
        covs <- data.frame(id = 1:n_joint, X_joint)
        names(coef)=colnames(X_joint)
        s1 <- simsurv::simsurv(dist = c("weibull"),lambdas = 0.1, gammas = 1.5,
                              betas = coef,x = covs, maxt = 5)
        observed_times_joint = s1$eventtime
        event_joint = s1$status
        # simple.ev.dat<-surv.simulate(foltime = 365,
        #                              dist.ev =  "weibull", anc.ev = 1, beta0.ev = 4.268,
        #                              dist.cens = "weibull", anc.cens = 1,beta0.cens = 5.368,
        #                              beta = coef, X = X_joint)
        # observed_times_joint = simple.ev.dat$stop-simple.ev.dat$start
        # event_joint = simple.ev.dat$status
        Z=Z_joint[main_index,,drop=F]  # separate main and external study for Z
        ZE=Z_joint[-main_index,,drop=F]
        W=W_joint[main_index,,drop=F]

        event = event_joint[main_index]
        eventE = event_joint[-main_index]
        observed_times = observed_times_joint[main_index]
        observed_timesE = observed_times_joint[-main_index]

        library(survival)
        surv_data <- Surv(time = observed_times, event = event)
        fit <- coxph(surv_data ~., data=data.frame(surv_data,Z,W),robust=T)
        beta_initial=fit$coefficients
        fitvar=fit$var
        surv_dataE <- Surv(time = observed_timesE, event = eventE)
        fitE <- coxph(surv_dataE ~., data=data.frame(surv_dataE,ZE))
        tilde_thetaZ=fitE$coefficients
        C_half = magic::adiag(diag(1,nrow=pZ+pW),diag(1,nrow=pZ))

        htlgmm_res = cv.htlgmm.cox(observed_times,event,Z,W,tilde_thetaZ,beta_initial = beta_initial,nopenalty = T,C_half =C_half)
        htlgmm_res2 = cv.htlgmm.cox(observed_times,event,Z,W,tilde_thetaZ,V_thetaZ=fitvar,beta_initial = beta_initial,nopenalty = T)
        beta_htlgmm = htlgmm_res$beta#_nopenalty
        beta_htlgmm2 = htlgmm_res2$beta#_nopenalty
        bias_initial = beta_initial - coef
        bias_htlgmm = beta_htlgmm - coef
        bias_htlgmm2 = beta_htlgmm2 - coef
        bias_vec<-rbind(bias_vec,c(bias_initial,bias_htlgmm,bias_htlgmm2))
    }

}
}
