
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
    GammaZZ=dU2U3dtheta[(pA+1):(pA+pZ),(pA+1):(pA+pZ),drop=F]*(-1)

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
        left_equal_id_test = left_id(timestest)

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
                           penalty.factor = w_adaptive,
                           alpha = final_alpha,lambda = cur_lam)
            cur_beta<-coef.glmnet(cv_fit)[-1]
            if(type_measure == "C"){
              tmp=Cindex(pred=prodv_rcpp(cbind(Atest,Ztest,Wtest),cur_beta),y=cbind(time=timestest,status=eventstest))
            }else{
              tmp=partial_likelihood(timestest,eventstest,cbind(Atest,Ztest,Wtest),cur_beta,left_equal_id_test)
            }
            tmp
        })
        cox_lam
    })
    rowSums(fold_cox_lambda)
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

beta_initial_cox_func<-function(surv_data,X,
                                initial_with_type,
                                fix_penalty){
    if(initial_with_type %in% c("coxph","ridge","lasso")){
        if(initial_with_type == "ridge"){initial_alpha=0}else{initial_alpha=1}
        if(initial_with_type == "coxph"){
            fit_initial <- coxph(surv_data ~., data=data.frame(surv_data,X))
            beta_initial=fit_initial$coefficients
        }else{
            fit_initial=cv.glmnet(x=X,y=surv_data,
                                  family="cox",alpha=initial_alpha,
                                  penalty.factor = fix_penalty)
            beta_initial=as.vector(coef(fit_initial,s="lambda.min"))
        }
    }else{stop("Select Initial Type from c('coxph','ridge','lasso')")}
    list("beta_initial"=beta_initial)
}

thetaA_cox_func<-function(pA,Z,A,surv_data,tilde_thetaZ,robust,hat_thetaA=NULL,V_thetaA=NULL){
    if(pA!=0){
        if(is.null(hat_thetaA)){
            if(!is.null(V_thetaA)){
                stop("With customized hat_thetaA input, V_thetaA is also needed.")
            }
            offset_term = prodv_rcpp(Z,tilde_thetaZ)

            hat_thetaA_coxph = coxph(surv_data~.+offset(offset_term),
                                     data=data.frame(surv_data,A),robust=robust)
            hat_thetaA=hat_thetaA_coxph$coefficients
            V_thetaA=hat_thetaA_coxph$var
        }
    }
    return(list("hat_thetaA"=hat_thetaA,
                "V_thetaA"=V_thetaA))
}


cv_C_lambda_Cweight_func<-function(tune_weight_method,
                                   index_fold,Z,W,A,times,events,
                                   C_half,inv_C,beta_initial,
                                   initial_with_type,
                                   hat_thetaA,
                                   tilde_thetaZ,
                                   V_thetaZ,
                                   weight_list,pZ,pW,pA,
                                   w_adaptive,final_alpha,
                                   pseudo_Xy_cox,lambda.min.ratio,
                                   nlambda,V_thetaA,use_sparseC,
                                   robust,
                                   fold_self_beta,
                                   X=NULL,XR=NULL,
                                   fix_lambda_list=NULL,sC_half=NULL){
    sample_p<-0.3
    index_test<-c(sample(which(events==1),round(sum(events)*sample_p)),
                  sample(which(events==0),round(sum(1-events)*sample_p)))

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
    left_equal_id_test = left_id(timestest)
    right_equal_id_train = right_id(timestrain)

    surv_data_train = Surv(time=timestrain,event=eventstrain)

    # renew the initial value
    if(fold_self_beta){
        initial_res1<-beta_initial_cox_func(surv_data_train,cbind(Atrain,Ztrain,Wtrain),
                                            initial_with_type,w_adaptive)
        beta_initial1<-initial_res1$beta_initial
    }else{beta_initial1<-beta_initial}


    thetaA_list<-thetaA_cox_func(pA,Ztrain,Atrain,surv_data_train,
                                 tilde_thetaZ,robust)
    hat_thetaA1<-thetaA_list$hat_thetaA
    V_thetaA1<-thetaA_list$V_thetaA

    inv_C_train = Delta_opt_cox_rcpp(Z=Ztrain,W=Wtrain,A=Atrain,
                                     times=timestrain,events=eventstrain,
                                     beta=beta_initial1,
                                     tilde_thetaZ=tilde_thetaZ,
                                     V_thetaZ=V_thetaZ,
                                     left_equal_id = left_equal_id_train,
                                     right_equal_id = right_equal_id_train,
                                     hat_thetaA=hat_thetaA1,
                                     V_thetaA=V_thetaA1)

    if(use_sparseC){inv_C_train<-diag(diag(inv_C_train))}
    C_half_train<-sqrtchoinv_rcpp2(inv_C_train)

    dev_lam_weight_fold<-lapply(1:length(weight_list),function(weight_id){
        cur_weight=weight_list[weight_id]

        C_half_weight<-weighted_C_half_func(inv_C_train,cur_weight,pA+pZ+pW,pZ,tune_weight_method,C_half_train)


        pseudo_Xy_list_train<-pseudo_Xy_cox(C_half_weight,Ztrain,Wtrain,Atrain,timestrain,
                                            eventstrain,beta_initial1,
                                            tilde_thetaZ=tilde_thetaZ,
                                            hat_thetaA=hat_thetaA1,
                                            left_equal_id=left_equal_id_train)
        initial_sf_train<-nrow(Ztrain)/sqrt(nrow(pseudo_Xy_list_train$pseudo_X))
        pseudo_X_train<-pseudo_Xy_list_train$pseudo_X/initial_sf_train
        pseudo_y_train<-pseudo_Xy_list_train$pseudo_y/initial_sf_train

        if(!is.null(fix_lambda_list)){
            lambda_list_weight<-fix_lambda_list
        }else{
            innerprod<-crossprodv_rcpp(pseudo_X_train,pseudo_y_train)[which(w_adaptive!=0)]
            lambda.max<-max(abs(innerprod))/nrow(pseudo_X_train)
            lambda_list_weight <-exp(seq(log(lambda.max),log(lambda.max*lambda.min.ratio),
                                         length.out=nlambda))
        }

        sapply(lambda_list_weight,function(cur_lam){
            cv_fit<-glmnet(x=pseudo_X_train,y=pseudo_y_train,
                           standardize=F,intercept=F,alpha = final_alpha,
                           penalty.factor = w_adaptive,
                           lambda = cur_lam)
            cur_beta<-coef.glmnet(cv_fit)[-1]
            cur_C<-Cindex(pred=prodv_rcpp(cbind(Atest,Ztest,Wtest),cur_beta),y=cbind(time=timestest,status=eventstest))
            cur_dev<-partial_likelihood(timestest,eventstest,cbind(Atest,Ztest,Wtest),cur_beta,left_equal_id_test)
            c(cur_dev,cur_C)
        })
    })
    # row is weight_list & col is lambda_list
    cv_dev1<-do.call(rbind, lapply(dev_lam_weight_fold, function(m) m[1,]))
    cv_auc1<-do.call(rbind, lapply(dev_lam_weight_fold, function(m) m[2,]))

    ids_dev<-which(cv_dev1==max(cv_dev1),arr.ind = TRUE)
    ids_auc<-which(cv_auc1==max(cv_auc1),arr.ind = TRUE)

    final_weight<-weight_list[ids_auc[1,1]] #nrow(ids_auc)
    final_weight_dev<-weight_list[ids_dev[1,1]] #nrow(ids_dev)

    return(list("final_weight"=final_weight,
                "final_weight_dev"=final_weight_dev,
                "deviance"=cv_dev1,
                "auc"=cv_auc1))
}


htlgmm.cox.default<-function(y,Z,W=NULL,
                             ext_study_info=NULL,
                             A=NULL,
                             penalty_type = "lasso",
                             initial_with_type = "lasso",
                             beta_initial = NULL,
                             weight_adaptivelasso = NULL,
                             alpha = NULL,
                             hat_thetaA = NULL,
                             V_thetaA = NULL,
                             robust = TRUE,
                             remove_penalty_Z = FALSE,
                             remove_penalty_W = FALSE,
                             inference = "default",
                             fix_C = NULL,
                             fix_inv_C = NULL,
                             refine_C = FALSE,
                             sqrt_matrix ="cholesky",
                             use_cv = TRUE,
                             type_measure = "C",
                             nfolds = 10,
                             foldid = NULL,
                             fix_lambda = NULL,
                             lambda_list = NULL,
                             nlambda = 100,
                             lambda.min.ratio = 0.0001,
                             tune_weight = FALSE,
                             fix_weight = NULL,
                             weight_list = NULL,
                             tune_weight_method = 1,
                             gamma_adaptivelasso = 1/2,
                             use_sparseC = FALSE,
                             seed.use = 97){

    set.seed(seed.use)
    if(!is.list(y)){
        stop("Please input y in the form of list with 'time' and 'event'.")
    }
    times=y[[1]]
    events=y[[2]]
    if(sum(is.na(times))>0){stop("time includes NA")}
    if(sum(is.na(events))>0){stop("event includes NA")}
    if(length(times)!=length(events)){stop("time and event should be the same length")}
    if(sum(is.na(Z))>0){stop("Z includes NA")}
    if(!is.null(W)&sum(is.na(W))>0){stop("W includes NA")}
    if(!is.null(A)&sum(is.na(A))>0){stop("A includes NA")}
    if (is.null(ext_study_info)){stop("Please input ext_study_info as trained model")}
    if(!penalty_type %in% c("none","adaptivelasso","lasso","ridge")){
        stop("Select penalty type from c('none','adaptivelasso','lasso','ridge').")}
    if(!type_measure%in% c("default", "deviance", "C")){
        stop("Select type_measure from c('default','deviance','C')")}
    if(is.null(dim(Z)[1])){
        warning("Z is input as a vector, convert Z into matrix with size nZ*1")
        Z=matrix(Z,ncol=1)}
    final_alpha = 1
    if(penalty_type == "ridge"){final_alpha = 0}

    if(length(ext_study_info[[1]]$Coeff) != ncol(Z)){
        stop("Please match ext_study_info with Z. (Current input length not match!)")
    }

    if(!penalty_type %in% c("none","adaptivelasso","lasso","ridge","elasticnet")){
        stop("Select penalty type from c('none','adaptivelasso','lasso','ridge','elasticnet').")
    }
    if(!sqrt_matrix %in% c("cholesky","svd")){
        stop("Select penalty type from c('cholesky','svd').")
    }

    if(inference == 'default'){
        if(penalty_type%in%c("none","adaptivelasso")){
            inference = TRUE
        }else{inference = FALSE}
    }

    ###########--------------###########
    # determine which kind of penalty for the final model
    if(!is.null(alpha) & penalty_type != "elasticnet"){
        warning("When using alpha between 0 and 1, please set penalty_type to be 'elasticnet'.")
    }
    if(penalty_type%in%c("adaptivelasso","lasso")){final_alpha = 1}
    if(penalty_type == "ridge"){final_alpha = 0}
    if(penalty_type == "elasticnet"){
        if(is.null(alpha)|!is.numeric(alpha)){
            stop("When using penalty_type is 'elasticnet', need the input of alpha between 0 and 1.")
        }else if(alpha<0|alpha>1){
            stop("When using penalty_type is 'elasticnet', need the input of alpha between 0 and 1.")
        }else{final_alpha=alpha}
    }


    tilde_thetaZ = ext_study_info[[1]]$Coeff
    V_thetaZ = ext_study_info[[1]]$Covariance
    surv_data = Surv(time=times,event=events)

    pZ=ncol(Z)
    if(is.null(W)){pW=0}else{pW=ncol(W)}
    if(is.null(A)){pA=0}else{pA=ncol(A)}


    thetaA_list<-thetaA_cox_func(pA,Z,A,surv_data,tilde_thetaZ,
                                 robust,hat_thetaA,V_thetaA)
    hat_thetaA<-thetaA_list$hat_thetaA
    V_thetaA<-thetaA_list$V_thetaA

    # sort the time and event
    set.seed(seed.use)
    timesorder=order(times)
    times=times[timesorder]
    events=events[timesorder]
    Z=Z[timesorder,,drop=F]
    if(!is.null(W)){W=W[timesorder,,drop=F]}
    if(!is.null(A)){A=A[timesorder,,drop=F]}
    X=cbind(A,Z,W)
    XR=cbind(A,Z)

    nX=nZ=length(times)

    right_equal_id = right_id(times)
    left_equal_id = left_id(times)


    Zid<-(pA+1):(pA+pZ)
    if(pW>0){Wid<-(pA+pZ+1):(pA+pZ+pW)}else{Wid=NULL}
    Acolnames=NULL
    if(pA>0){
        Acolnames=colnames(A)
        if(is.null(Acolnames[1])){
            Acolnames=paste0('A',1:pA)
        }
        if(length(unique(A[,1])) == 1){
            if(unique(A[,1]) == 1){
                Acolnames[1]='intercept'
            }
        }
        colnames(A)=Acolnames
    }
    Zcolnames=colnames(Z)
    if(pW>0){Wcolnames=colnames(W)}else{Wcolnames=NULL}
    if(is.null(Zcolnames[1])){
        Zcolnames=paste0('Z',1:pZ)
        colnames(Z)=Zcolnames
    }
    if(!is.null(W) & is.null(Wcolnames[1])){
        Wcolnames=paste0('W',1:pW)
        colnames(W)=Wcolnames
    }
    Xcolnames<-c(Acolnames,Zcolnames,Wcolnames)

    fix_penalty<-c(rep(0,pA),rep(1,pZ+pW))
    if(remove_penalty_Z){fix_penalty[Zid]<-0}
    if(remove_penalty_W){fix_penalty[Wid]<-0}
    if(remove_penalty_Z & remove_penalty_W){
        penalty_type = "none"
        warning("All penalties are removed, turn to no penalties!")
    }
    if(penalty_type == "none"){
        initial_with_type = "coxph"
        use_cv = FALSE
    }
    if(!is.null(beta_initial) & length(beta_initial)!=pA+pZ+pW){
        warning("beta_initial should be from A,Z,W.\n Length not match, compute default initial instead.")
        beta_initial=NULL
    }
    if(is.null(beta_initial)){
        beta_initial<-beta_initial_cox_func(surv_data,X,initial_with_type,
                                        fix_penalty)
    }

    if (penalty_type == "adaptivelasso"){
        if(is.null(weight_adaptivelasso)){
            w_adaptive<-1/abs(beta_initial)^gamma_adaptivelasso
            w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*100
            w_adaptive<-w_adaptive*fix_penalty
        }else{
            if(length(weight_adaptivelasso) != length(fix_penalty) & length(weight_adaptivelasso) != length(fix_penalty)-1 ){
                stop("The length of 'weight_adaptivelasso' is invalid. Either the same length as beta_initial or without intercept is accepted.")
            }else{
                if(length(weight_adaptivelasso) == length(fix_penalty)){
                    w_adaptive<-weight_adaptivelasso*fix_penalty
                }else if(Acolnames[1] == "intercept"){
                    w_adaptive<-c(0,weight_adaptivelasso)*fix_penalty
                }else{
                    stop("The length of 'weight_adaptivelasso' is invalid. Either the same length as beta_initial or without intercept is accepted.")
                }
            }

        }
    }else{w_adaptive<-fix_penalty}

    if(is.null(fix_C) & is.null(fix_inv_C)){

        inv_C = Delta_opt_cox_rcpp(Z,W,A,times,events,
                                   beta=beta_initial,
                                   tilde_thetaZ=tilde_thetaZ,
                                   V_thetaZ=V_thetaZ,
                                   left_equal_id,
                                   right_equal_id,
                                   hat_thetaA=hat_thetaA,
                                   V_thetaA=V_thetaA,
                                   X=X,XR=XR)
        if(refine_C){
            return(list("Delta_opt"=inv_C))
        }

        sC_half<-diag(1/sqrt(diag(inv_C)))
        if(use_sparseC){
            C_half<-sC_half
        }else{
            if(sqrt_matrix =="svd"){
                inv_C_svd=fast.svd(inv_C+diag(1e-15,nrow(inv_C)))
                C_half=prod_rcpp(inv_C_svd$v,(t(inv_C_svd$u)*1/sqrt(inv_C_svd$d)))
            }else if(sqrt_matrix =="cholesky"){
                C_half<-sqrtchoinv_rcpp2(inv_C)
            }
        }
    }else if(is.null(fix_inv_C)){
        if(nrow(fix_C)!=pA+pZ+pW+pZ){
            stop("Input fix_C dimension is wrong!")}
        C_half<-sqrtcho_rcpp2(fix_C)
    }else{
        if(nrow(fix_inv_C)!=pA+pZ+pW+pZ){
            stop("Input fix_inv_C dimension is wrong!")}
        C_half<-sqrtchoinv_rcpp2(fix_inv_C)
    }

    psXy = pseudo_Xy_cox(C_half,Z,W,A,times,events,
                         beta_initial,tilde_thetaZ=tilde_thetaZ,
                         hat_thetaA=hat_thetaA,X=X,XR=XR,
                         left_equal_id=left_equal_id)
    initial_sf<-nrow(Z)/sqrt(nrow(psXy$pseudo_X))
    psX = psXy$pseudo_X/initial_sf
    psy = psXy$pseudo_y/initial_sf

    if(!is.null(lambda_list)){fix_lambda_list<-lambda_list}else{
        fix_lambda_list<-NULL}
    if(penalty_type != "none" & is.null(fix_lambda)&is.null(lambda_list)){
        # if(penalty_type == "adaptivelasso"){
        #     fit_final<-glmnet(x= pseudo_X,y= pseudo_y,standardize=F,
        #                       intercept=F,alpha = final_alpha,penalty.factor = w_adaptive)
        #     lambda_list<-fit_final$lambda
        #     lambda_list<-lambda_list[!is.na(lambda_list)]
        # }else{
        innerprod<-crossprodv_rcpp(psX,psy)[which(fix_penalty!=0)]
        lambda.max<-max(abs(innerprod))/nrow(psX)
        lambda_list <-exp(seq(log(lambda.max),log(lambda.max*lambda.min.ratio),
                              length.out=nlambda))
        # }
    }
    if(!is.null(fix_lambda)){
        use_cv = FALSE
        if(!is.numeric(fix_lambda)|fix_weight<0){
            stop("fix_lambda should be >= 0.")
        }
    }

    if(tune_weight){
        if(is.null(weight_list)){
            if(tune_weight_method %in% c(1,2,3,6)){
                weight_list<-c(1,2,4,8,16)
            }else{
                weight_list<-c(1,2,4,8,16)-1
            }
        }
    }

    if(!use_cv){
        if(penalty_type == "none"){
            beta=prodv_rcpp(choinv_rcpp(self_crossprod_rcpp(psX)),crossprodv_rcpp(psX,psy))
            return_list<-list("beta"=beta)
        }else{
            if(!is.null(fix_lambda)){
                fit_final_fixed_lambda=glmnet(x= psX,y= psy,standardize=F,
                                              intercept=F,alpha = final_alpha,
                                              penalty.factor = w_adaptive,
                                              lambda = fix_lambda)
                beta<-coef.glmnet(fit_final_fixed_lambda)[-1]
                return_list<-list("beta"=beta,
                                  "fix_lambda"=fix_lambda)
            }else{
                fit_final_lambda_list=glmnet(x= psX,y= psy,standardize=F,
                                             intercept=F,alpha = final_alpha,
                                             penalty.factor = w_adaptive,
                                             lambda = lambda_list)
                return_list<-list("beta"=fit_final_lambda_list$beta,
                                  "lambda_list"=fit_final_lambda_list$lambda)
                if(inference){warning("When use_cv=F,fix_lambda is NULL, no inference will be done")}
                inference=FALSE
            }
        }
    }

    if(use_cv){
        if(is.null(foldid)){
            index_fold = createFolds(events,k = nfolds)
        }else{
            uni_fold=sort(unique(foldid))
            nfolds=length(uni_fold)
            index_fold<-lapply(1:length(uni_fold),function(i){which(foldid==uni_fold[i])})
        }

        if(tune_weight){
            fold_self_beta = TRUE
            cv_res<-cv_C_lambda_Cweight_func(tune_weight_method,
                                             index_fold,Z,W,A,times,events,
                                             C_half,inv_C,beta_initial,
                                             initial_with_type,
                                             hat_thetaA,
                                             tilde_thetaZ,
                                             V_thetaZ,
                                             weight_list,pZ,pW,pA,
                                             w_adaptive,final_alpha,
                                             pseudo_Xy_cox,lambda.min.ratio,
                                             nlambda,V_thetaA,use_sparseC,
                                             robust,
                                             fold_self_beta,
                                             X,XR,fix_lambda_list,sC_half)

            cv_auc<-cv_res$auc
            ids_auc<-which(cv_auc==max(cv_auc),arr.ind = TRUE)
            ids_auc<-ids_auc[1,] #nrow(ids_auc)
            # cv_dev<-cv_res$deviance
            # ids_dev<-which(cv_dev==min(cv_dev),arr.ind = TRUE)
            # ids_dev<-ids_dev[1,] #nrow(ids_dev)
            print(paste0("auc_weightid:",ids_auc[1]))

            final.ratio.min<-1
            weight<-cv_res$final_weight
            final.weight.min<-weight
            C_half<-weighted_C_half_func(inv_C,weight,pA+pZ+pW,pZ,tune_weight_method,C_half)
            if( !((tune_weight_method%in%c(1,4,5,6) & weight == 1)|
                  (tune_weight_method%in%c(2,3) & weight == 0)) ){

                pseudo_Xy_list<-pseudo_Xy_cox(C_half,Z,W,A,times,events,
                                              beta_initial,tilde_thetaZ=tilde_thetaZ,
                                              hat_thetaA=hat_thetaA,X=X,XR=XR,
                                              left_equal_id=left_equal_id)

                initial_sf<-nrow(Z)/sqrt(nrow(pseudo_Xy_list$pseudo_X))
                psX<-pseudo_Xy_list$pseudo_X/initial_sf
                psy<-pseudo_Xy_list$pseudo_y/initial_sf
                if(!is.null(fix_lambda_list)){
                    lambda_list<-fix_lambda_list
                }else{
                    innerprod<-crossprodv_rcpp(psX,psy)[which(w_adaptive!=0)]
                    lambda.max<-max(abs(innerprod))/nrow(psX)
                    lambda_list <-exp(seq(log(lambda.max),log(lambda.max*lambda.min.ratio),
                                          length.out=nlambda))
                }
            }

            res_weight<-cv_cox_lambda_func(index_fold,Z,W,A,times,events,
                               C_half,beta_initial,lambda_list,
                               tilde_thetaZ,hat_thetaA,type_measure,
                               w_adaptive,final_alpha)

            cv_auc1<-res_weight
            max_id<-which.max(cv_auc1)
            final.lambda.min<-lambda_list[max_id]
            return_list<-list("cv_Cindex"=cv_auc1)

            return_list<-c(return_list,
                           list("lambda_list"=lambda_list,
                                "weight_list"=weight_list))
        }


        if(!tune_weight){
            cv_res=cv_cox_lambda_func(index_fold,Z,W,A,times,events,
                                      C_half,beta_initial,lambda_list,
                                      tilde_thetaZ,hat_thetaA,type_measure,
                                      w_adaptive,final_alpha)
            final.lambda.min=lambda_list[which.max(cv_res)]
            final.weight.min = 1
            return_list<-list("lambda_list"=lambda_list)
            if(type_measure == "C"){
                return_list = c(return_list,list("cv_Cindex"=cv_res))
            }else{return_list = c(return_list,list("cv_dev"=cv_res))}
        }

        ###########--------------###########
        # Model final: use best lambda, best weight to build model

        fit_final_lam=glmnet(x= psX,y= psy,standardize=F,
                         intercept=F,alpha = final_alpha,
                         penalty.factor = w_adaptive,
                         lambda = final.lambda.min)

        beta<-coef.glmnet(fit_final_lam)[-1]

        return_list<-c(list("beta"=beta,
                            "lambda_min"=final.lambda.min,
                            "weight_min"=final.weight.min),
                       return_list)

    }


    ###########--------------###########
    # perform inference

    if(inference){
        index_nonzero<-which(beta!=0)
        if(length(index_nonzero) > 1){
            if(penalty_type == "lasso"){
                warning("Current penalty is lasso, please turn to adaptivelasso for inference")
            }

            ###########--------------###########
            # refine C will cover the previously used C
            inv_C = Delta_opt_cox_rcpp(Z,W,A,times,events,
                                       beta=beta,
                                       tilde_thetaZ=tilde_thetaZ,
                                       V_thetaZ=V_thetaZ,
                                       left_equal_id,
                                       right_equal_id,
                                       hat_thetaA=hat_thetaA,
                                       V_thetaA=V_thetaA,
                                       X=X,XR=XR)
            if(refine_C){
                C_half<-sqrtchoinv_rcpp2(inv_C)
                pseudo_Xy_list<-pseudo_Xy_cox(C_half,Z,W,A,times,events,
                                     beta=beta,tilde_thetaZ=tilde_thetaZ,
                                     hat_thetaA=hat_thetaA,X=X,XR=XR,
                                     left_equal_id=left_equal_id)
                psX<-pseudo_Xy_list$pseudo_X/nZ
                psXtX<-self_crossprod_rcpp(psX)
                psXtX_non0<-psXtX[index_nonzero,index_nonzero,drop=F]
                inv_psXtX_non0<-choinv_rcpp2(psXtX_non0)
                inv_psXtX_final<-inv_psXtX_non0

            }else{

                sC_half<-diag(1/sqrt(diag(inv_C)))
                if(use_sparseC){
                    C_half<-sC_half
                }else{
                    if(sqrt_matrix =="svd"){
                        inv_C_svd=fast.svd(inv_C+diag(1e-15,nrow(inv_C)))
                        C_half=prod_rcpp(inv_C_svd$v,(t(inv_C_svd$u)*1/sqrt(inv_C_svd$d)))
                    }else if(sqrt_matrix =="cholesky"){
                        C_half<-sqrtchoinv_rcpp2(inv_C)
                    }
                }
                runsandwich<-F
                if(tune_weight){
                    weight<-final.weight.min
                    C_half<-weighted_C_half_func(inv_C,weight,pA+pZ+pW,pZ,tune_weight_method,C_half)
                    if(!((tune_weight_method%in%c(1,4,5,6) & weight == 1)|
                         (tune_weight_method%in%c(2,3) & weight == 0))){
                        runsandwich<-T
                    }
                }
                if(!is.null(fix_C)|!is.null(fix_inv_C)|use_sparseC){runsandwich<-T}
                ###########--------------###########
                # Compute new pseudo_X

                pseudo_Xy_list<-pseudo_Xy_cox(C_half,Z,W,A,times,events,
                                              beta=beta,tilde_thetaZ=tilde_thetaZ,
                                              hat_thetaA=hat_thetaA,X=X,XR=XR,
                                              left_equal_id=left_equal_id)
                psX<-pseudo_Xy_list$pseudo_X/nZ

                psXtX<-self_crossprod_rcpp(psX)
                psXtX_non0<-psXtX[index_nonzero,index_nonzero,drop=F]
                inv_psXtX_non0<-choinv_rcpp2(psXtX_non0)
                inv_psXtX_final<-inv_psXtX_non0
                ###########--------------###########
                # When the C using is not optimal C,

                if(runsandwich){
                    inv_C_half<-sqrtcho_rcpp2(inv_C+diag(1e-15,nrow(inv_C)))
                    psX_mid<-prod_rcpp(inv_C_half,crossprod_rcpp(C_half,psX))
                    psXtX_mid<-self_crossprod_rcpp(psX_mid)
                    psXtX_mid_non0<-psXtX_mid[index_nonzero,index_nonzero,drop=F]
                    inv_psXtX_final<-prod_rcpp(prod_rcpp(inv_psXtX_non0,psXtX_mid_non0),inv_psXtX_non0)

                    if(sum(diag(inv_psXtX_non0)<=0)>0){
                        psXtX_mid_non0_half<-sqrtcho_rcpp2(psXtX_mid_non0+diag(1e-15,nrow(psXtX_mid_non0)))
                        inv_psXtX_final<-self_crossprod_rcpp(prod_rcpp(psXtX_mid_non0_half,inv_psXtX_non0))
                    }
                }

            }
            final_vcov<-inv_psXtX_final/nZ
            final_v<-diag(final_vcov)

            pval_final<-pchisq(beta[index_nonzero]^2/final_v,1,lower.tail = F)
            pval_final1<-p.adjust(pval_final,method = "BH")

            selected_pos<-index_nonzero[which(pval_final1<0.05)]
            return_list<-c(return_list,
                           list("selected_vars"=
                                    list("position"=index_nonzero,
                                         "name"=Xcolnames[index_nonzero],
                                         "coef"=beta[index_nonzero],
                                         "variance"=final_v,
                                         "variance_covariance"=final_vcov,
                                         "pval"=pval_final,
                                         "FDR_adjust_position"=selected_pos,
                                         "FDR_adjust_name"=Xcolnames[selected_pos])
                           ))

        }}

    if(is.null(fix_C) & is.null(fix_inv_C)){
        return_list<-c(return_list,list("Delta_opt"=inv_C))
    }
    return(return_list)
}


if(0){
    iAUC<-function(times,events,predict_risk){
        time_points <- seq(min(times), max(times), by = (max(times)-min(times))/length(times)/2)
        roc_results <- timeROC::timeROC(times,events,predict_risk,cause = 1,times = time_points)
        mean(roc_results$AUC,na.rm = T)
    }

    pacman::p_load(glmnet,Rcpp,caret,survival)
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
        s1 <- simsurv::simsurv(dist = c("weibull"),lambdas = .2, gammas = .1,
                              betas = coef,x = covs, maxt = 365)
        observed_times_joint = s1$eventtime
        event_joint = s1$status
        mean(event_joint)
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
