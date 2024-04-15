
Delta_opt_rcpp<-function(y,Z,W,family,
                    study_info,A=NULL,pA=NULL,pZ=NULL,
                    beta=NULL,hat_thetaA=NULL,
                    V_thetaA=NULL,use_offset=TRUE,X=NULL,XR=NULL){
    n_main=length(y)
    tilde_thetaZ=study_info[[1]]$Coeff
    tilde_theta=c(hat_thetaA,tilde_thetaZ)
    if(is.null(X)){X=cbind(A,Z,W)}
    if(is.null(XR)){XR=cbind(A,Z)}
    mu_X_beta=prodv_rcpp(X,beta)
    mu_XR_theta=prodv_rcpp(XR,tilde_theta)
    if(family == "binomial"){
        mu_X_beta=expit_rcpp(mu_X_beta)
        mu_XR_theta=expit_rcpp(mu_XR_theta)
        mu_prime_XR_theta=mu_XR_theta*(1-mu_XR_theta)
    }else{mu_prime_XR_theta=1}
    DX=X*c(mu_X_beta-y)
    DZ=Z*c(mu_X_beta-mu_XR_theta)
    DZ2=Z*c(mu_prime_XR_theta)
    V_U1=(1/n_main)*self_crossprod_rcpp(DX)
    V_U2=(1/n_main)*self_crossprod_rcpp(DZ)
    Cov_U1U2=(1/n_main)*crossprod_rcpp(DX,DZ)
    GammaZZ=(1/n_main)*crossprod_rcpp(DZ2,Z)
    V_thetaZ=study_info[[1]]$Covariance
    if(is.null(dim(V_thetaZ)[1])){V_thetaZ=as.matrix(V_thetaZ,nrow=pZ,ncol=pZ)}
    Delta22=V_U2+prod_rcpp(prod_rcpp(GammaZZ,(n_main*V_thetaZ)),t(GammaZZ))
    Delta12=Cov_U1U2
    if(pA!=0){
        GammaZA=(1/n_main)*crossprod_rcpp(DZ2,A)
        if(use_offset){
            inv_GammaAA=choinv_rcpp((1/n_main)*crossprod_rcpp(A*c(mu_prime_XR_theta),A)+diag(1e-15,pA))
            DDA=prod_rcpp(A*c(mu_XR_theta-y),prod_rcpp(inv_GammaAA,t(GammaZA)))
            Cov_U1theta=(1/n_main)*crossprod_rcpp(DX,DDA)
            Cov_U2theta=(1/n_main)*crossprod_rcpp(DZ,DDA)
        }else{
            inv_GammaXRXR=choinv_rcpp((1/n_main)*crossprod_rcpp(XR*c(mu_prime_XR_theta),XR)+diag(1e-15,pA+pZ))
            #V_thetaA = (1/n_main)*inv_GammaXRXR[1:pA,]%*%((1/n_main)* crossprod(XR*c(mu_XR_theta-y)) )%*%inv_GammaXRXR[,1:pA]
            DDXR=prod_rcpp(XR*c(mu_XR_theta-y),prod_rcpp(inv_GammaXRXR[,1:pA,drop=F],t(GammaZA)))
            Cov_U1theta=(1/n_main)*crossprod_rcpp(DX,DDXR)
            Cov_U2theta=(1/n_main)*crossprod_rcpp(DZ,DDXR)
        }
        Delta22 = Delta22 + prod_rcpp(prod_rcpp(GammaZA,(n_main*V_thetaA)),t(GammaZA))+Cov_U2theta+t(Cov_U2theta)
        Delta12 = Delta12 + Cov_U1theta
    }
    Delta = rbind(cbind(V_U1,Delta12),cbind(t(Delta12),Delta22))
    Delta
}

vcov_sandwich_rcpp<-function(y,A,Z,family,study_info,pA,
                        hat_thetaA,use_offset,XR=NULL){
    n_main = length(y)
    if(is.null(dim(XR)[1])){XR=cbind(A,Z)}
    tilde_thetaZ=study_info[[1]]$Coeff
    tilde_theta=c(hat_thetaA,tilde_thetaZ)
    mu_XR_theta=prodv_rcpp(cbind(A,Z),tilde_theta)
    if(family == "binomial"){
        mu_XR_theta = expit_rcpp(mu_XR_theta)
        mu_prime_XR_theta=mu_XR_theta*(1-mu_XR_theta)
    }else{mu_prime_XR_theta=1}
    if(use_offset){
        inv_GammaAA=choinv_rcpp((1/n_main)*crossprod_rcpp(A*c(mu_prime_XR_theta),A)+diag(1e-15,pA))
        V_thetaA = (1/n_main)*prod_rcpp(prod_rcpp(inv_GammaAA,((1/n_main)* self_crossprod_rcpp(A*c(mu_XR_theta-y)))),inv_GammaAA)
    }else{
        inv_GammaXRXR=choinv_rcpp((1/n_main)*crossprod_rcpp(XR*c(mu_prime_XR_theta),XR)+diag(1e-15,ncol(XR)))
        V_thetaA = (1/n_main)*prod_rcpp(prod_rcpp(inv_GammaXRXR[1:pA,,drop=F],((1/n_main)* self_crossprod_rcpp(XR*c(mu_XR_theta-y)))),inv_GammaXRXR[,1:pA,drop=F])
    }
}
pseudo_Xy_gaussian_rcpp<-function(
        C_half,Z,W,A,y,beta=NULL,hat_thetaA=NULL,study_info=NULL,X=NULL,XR=NULL){
    tilde_thetaZ=study_info[[1]]$Coeff
    tilde_theta=c(hat_thetaA,tilde_thetaZ)
    if(is.null(X)){X=cbind(A,Z,W)}
    if(is.null(XR)){XR=cbind(A,Z)}
    pseudo_X=prod_rcpp(C_half,crossprod_rcpp(cbind(X,Z),X))
    pseudo_y1=crossprodv_rcpp(X,y)
    pseudo_y2=prodv_rcpp(crossprod_rcpp(Z,XR),tilde_theta)
    pseudo_y=prodv_rcpp(C_half,c(pseudo_y1,pseudo_y2))
    list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
}

pseudo_Xy_binomial_rcpp<-function(
        C_half,Z,W,A,y,beta=NULL,hat_thetaA=NULL,study_info=NULL,X=NULL,XR=NULL){
    tilde_thetaZ=study_info[[1]]$Coeff
    tilde_theta=c(hat_thetaA,tilde_thetaZ)
    if(is.null(dim(X)[1])){X=cbind(A,Z,W)}
    if(is.null(dim(XR)[1])){XR=cbind(A,Z)}
    expit_beta=expit_rcpp(prodv_rcpp(X,beta))
    dexpit_beta=expit_beta*(1-expit_beta)
    pseudo_X=prod_rcpp(C_half,crossprod_rcpp(cbind(X,Z),X*dexpit_beta))
    u1=crossprodv_rcpp(X,(expit_beta-y))
    u2=crossprodv_rcpp(Z,c(expit_beta-expit_rcpp(prodv_rcpp(XR,tilde_theta))))
    pseudo_y= -prodv_rcpp(C_half,c(u1,u2)) + prodv_rcpp(pseudo_X,beta)
    list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
}

beta_initial_func<-function(y,X,A,pA,
                            family,
                            initial_with_type,
                            fix_penalty){
    if(initial_with_type %in% c("glm","ridge","lasso")){
        if(initial_with_type == "ridge"){initial_alpha=0}else{initial_alpha=1}
        if(initial_with_type == "glm"){
            df=data.frame(y,X)
            if(family=="binomial"){
                fit_initial=speedglm(y~0+.,data = df,family = binomial())
            }else if(family=="gaussian"){
                fit_initial=speedlm(y~0+.,data = df)
            }
            beta_initial=fit_initial$coefficients
            beta_initial1se=beta_initial
        }else if(pA == 0){
            fit_initial=cv.glmnet(x=X,y= y,
                                  alpha = initial_alpha,
                                  penalty.factor = fix_penalty,
                                  family=family)
            beta_initial=c(coef.glmnet(fit_initial,s="lambda.min")[-1])
            beta_initial1se=c(coef.glmnet(fit_initial,s="lambda.1se")[-1])
        }else if(length(unique(A[,1]))==1){
            if(unique(A[,1])==1){
                fit_initial=cv.glmnet(x=X[,-1,drop=F],y=y,
                                      alpha = initial_alpha,
                                      penalty.factor = fix_penalty[-1],
                                      family=family)
                beta_initial=as.vector(coef.glmnet(fit_initial,s="lambda.min"))
                beta_initial1se=as.vector(coef.glmnet(fit_initial,s="lambda.1se"))
            }else{
                stop("The first column of A is constant, then it should be 1 for intercept.")
            }
        }else{
            fit_initial=cv.glmnet(x=X,y= y,
                                  alpha = initial_alpha,
                                  penalty.factor = fix_penalty,
                                  family=family)
            beta_initial=c(coef.glmnet(fit_initial,s="lambda.min")[-1])
            beta_initial1se=c(coef.glmnet(fit_initial,s="lambda.1se")[-1])
        }
    }else{stop("Select Initial Type from c('glm','ridge','lasso')")}

    list("beta_initial"=beta_initial,"beta_initial1se"=beta_initial1se)

}

thetaA_func<-function(pA,Z,A,y,study_info,family,use_offset,V_thetaA_sandwich,hat_thetaA=NULL,V_thetaA=NULL){
    if(pA!=0){
        if(is.null(hat_thetaA)){
            if(!is.null(V_thetaA)){
                stop("With customized hat_thetaA input, V_thetaA is also needed")
            }
            if(use_offset){
                offset_term = prodv_rcpp(Z,study_info[[1]]$Coeff)
                df=data.frame(y,A)
                if(family=="binomial"){
                    hat_thetaA_glm=speedglm(y~0+.,data = df,offset = offset_term,family = binomial())
                }else if(family=="gaussian"){
                    hat_thetaA_glm=speedlm(y~0+.,data = df,offset = offset_term)
                    #hat_thetaA_glm=lm(y~0+.,data = df,offset = offset_term)
                }
                hat_thetaA=hat_thetaA_glm$coefficients
                if(V_thetaA_sandwich){
                    V_thetaA=vcov_sandwich_rcpp(y=y,A=A,Z=Z,family=family,
                                                study_info=study_info,pA=pA,
                                                hat_thetaA=hat_thetaA,
                                                use_offset=use_offset)
                }else{V_thetaA=vcov(hat_thetaA_glm)}
            }else{
                df=data.frame(y,A,Z)
                if(family=="binomial"){
                    hat_thetaA_glm=speedglm(y~0+.,data = df,family = binomial())
                }else if(family=="gaussian"){
                    hat_thetaA_glm=speedlm(y~0+.,data = df)
                }
                hat_thetaA=hat_thetaA_glm$coefficients[1:pA]
                if(V_thetaA_sandwich){
                    V_thetaA=vcov_sandwich_rcpp(y,A,Z,family,study_info,pA,
                                                hat_thetaA,use_offset)
                }else{V_thetaA=vcov(hat_thetaA_glm)[1:pA,1:pA,drop=F]}
            }
        }

        if(is.null(dim(V_thetaA)[1])){
            V_thetaA = as.matrix(V_thetaA,nrow=pA,ncol=pA)
        }
    }
    return(list("hat_thetaA"=hat_thetaA,
                "V_thetaA"=V_thetaA))
}

## cross validation function for continuous y with lambda only
cv_mse_lambda_func<-function(index_fold,Z,W,A,y,
                             C_half,beta_initial,hat_thetaA,
                             study_info,lambda_list,
                             w_adaptive,final_alpha,
                             pseudo_Xy){
    fold_mse_lambda<-sapply(1:length(index_fold), function(cur_fold){
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
        ytrain<-y[-index_test]
        ytest<-y[index_test]
        pseudo_Xy_list_train<-pseudo_Xy(C_half,Ztrain,Wtrain,Atrain,
                                        ytrain,beta = beta_initial,hat_thetaA = hat_thetaA,
                                        study_info=study_info)
        initial_sf_train<-nrow(Ztrain)/sqrt(nrow(pseudo_Xy_list_train$pseudo_X))
        pseudo_X_train<-pseudo_Xy_list_train$pseudo_X/initial_sf_train
        pseudo_y_train<-pseudo_Xy_list_train$pseudo_y/initial_sf_train
        mse_lam<-sapply(lambda_list,function(cur_lam){
            cv_fit<-glmnet(x= (pseudo_X_train),y= (pseudo_y_train),
                           standardize=F,intercept=F,
                           alpha = final_alpha,penalty.factor = w_adaptive,lambda = cur_lam)
            cur_beta<-coef.glmnet(cv_fit)[-1]
            mean(( prodv_rcpp(cbind(Atest,Ztest,Wtest),cur_beta) - ytest)^2)
        })
        mse_lam
    })
    rowMeans(fold_mse_lambda)
}

## cross validation function for continuous y with lambda and ratio
cv_mse_lambda_ratio_func<-function(index_fold,Z,W,A,y,
                                   C_half,beta_initial,hat_thetaA,
                                   study_info,lambda_list,
                                   ratio_range,pZ,pW,pA,
                                   w_adaptive,final_alpha,
                                   pseudo_Xy){
    fold_mse_lambda_ratio<-lapply(1:length(index_fold), function(cur_fold){
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
        ytrain<-y[-index_test]
        ytest<-y[index_test]
        pseudo_Xy_list_train<-pseudo_Xy(C_half,Ztrain,Wtrain,Atrain,
                                        ytrain,beta = beta_initial,hat_thetaA = hat_thetaA,
                                        study_info=study_info)
        initial_sf_train<-nrow(Ztrain)/sqrt(nrow(pseudo_Xy_list_train$pseudo_X))
        pseudo_X_train<-pseudo_Xy_list_train$pseudo_X/initial_sf_train
        pseudo_y_train<-pseudo_Xy_list_train$pseudo_y/initial_sf_train
        mse_lam_ratio_fold<-sapply(lambda_list,function(cur_lam){
            sapply(ratio_range,function(cur_ratio){
                ratio_vec<-c(rep(1,pA),rep(cur_ratio,pZ),rep(1,pW))
                w_adaptive_ratio<-w_adaptive*ratio_vec
                cv_fit<-glmnet(x=pseudo_X_train,y=pseudo_y_train,
                               standardize=F,intercept=F,alpha = final_alpha,
                               penalty.factor = w_adaptive_ratio,
                               lambda = cur_lam)
                cur_beta<-coef.glmnet(cv_fit)[-1]
                mean(( prodv_rcpp(cbind(Atest,Ztest,Wtest),cur_beta) - ytest)^2)
            })
        }) # row is ratio_range & col is lambda_list
        mse_lam_ratio_fold
    })
    cv_mse_lambda_ratio<-Reduce(`+`, fold_mse_lambda_ratio)/length(index_fold)
    cv_mse_lambda_ratio
}


## cross validation function for binary y with lambda only
cv_dev_lambda_func<-function(index_fold,Z,W,A,y,
                             C_half,beta_initial,hat_thetaA,
                             study_info,lambda_list,
                             w_adaptive,final_alpha,
                             pseudo_Xy){
    dev_fold<-lapply(1:length(index_fold), function(cur_fold){
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
        ytrain<-y[-index_test]
        ytest<-y[index_test]
        pseudo_Xy_list_train<-pseudo_Xy(C_half,Ztrain,Wtrain,Atrain,
                                        ytrain,beta = beta_initial,hat_thetaA = hat_thetaA,
                                        study_info=study_info)
        initial_sf_train<-nrow(Ztrain)/sqrt(nrow(pseudo_Xy_list_train$pseudo_X))
        pseudo_X_train<-pseudo_Xy_list_train$pseudo_X/initial_sf_train
        pseudo_y_train<-pseudo_Xy_list_train$pseudo_y/initial_sf_train
        dev_lam<-sapply(lambda_list,function(cur_lam){
            cv_fit<-glmnet(x= (pseudo_X_train),y= (pseudo_y_train),
                           standardize=F,intercept=F,
                           alpha = final_alpha,penalty.factor = w_adaptive,lambda = cur_lam)
            cur_beta<-coef.glmnet(cv_fit)[-1]
            probtest <- expit_rcpp(prodv_rcpp(cbind(Atest,Ztest,Wtest),cur_beta))
            cur_dev <- -2*mean( ytest * log(probtest) + (1 - ytest) * log(1 - probtest) )
            suppressMessages(cur_auc<-c(auc(ytest,probtest,direction = "<")))

            c(cur_dev,cur_auc)
        })
        dev_lam
    })
    sum_dev_lam<-Reduce(`+`, dev_fold)
    sum_dev_lam=sum_dev_lam/length(index_fold)
    for(i in 1:length(dev_fold)){dev_fold[[i]]=dev_fold[[i]]^2}
    sum_dev_lam_sq<-Reduce(`+`, dev_fold)
    sum_dev_lam_sq=sum_dev_lam_sq/length(index_fold)
    sum_dev_lam_sq=sum_dev_lam_sq-sum_dev_lam^2
    sum_dev_lam_sd=sqrt(sum_dev_lam_sq)
    list("deviance"=sum_dev_lam[1,],"auc"=sum_dev_lam[2,],
         "deviance_sd"=sum_dev_lam_sd[1,],"auc_sd"=sum_dev_lam_sd[2,])
}

## cross validation function for continuous y with lambda and ratio
cv_dev_lambda_ratio_func<-function(index_fold,Z,W,A,y,
                                   C_half,beta_initial,hat_thetaA,
                                   study_info,lambda_list,
                                   ratio_range,pZ,pW,pA,
                                   w_adaptive,final_alpha,
                                   pseudo_Xy){
    dev_lam_ratio<-lapply(1:length(index_fold), function(cur_fold){
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
        ytrain<-y[-index_test]
        ytest<-y[index_test]
        pseudo_Xy_list_train<-pseudo_Xy(C_half,Ztrain,Wtrain,Atrain,
                                        ytrain,beta = beta_initial,hat_thetaA = hat_thetaA,
                                        study_info=study_info)
        initial_sf_train<-nrow(Ztrain)/sqrt(nrow(pseudo_Xy_list_train$pseudo_X))
        pseudo_X_train<-pseudo_Xy_list_train$pseudo_X/initial_sf_train
        pseudo_y_train<-pseudo_Xy_list_train$pseudo_y/initial_sf_train
        dev_lam_ratio_fold<-lapply(ratio_range,function(cur_ratio){
            sapply(lambda_list,function(cur_lam){
                ratio_vec<-c(rep(1,pA),rep(cur_ratio,pZ),rep(1,pW))
                w_adaptive_ratio<-w_adaptive*ratio_vec
                cv_fit<-glmnet(x=pseudo_X_train,y=pseudo_y_train,
                               standardize=F,intercept=F,alpha = final_alpha,
                               penalty.factor = w_adaptive_ratio,
                               lambda = cur_lam)
                cur_beta<-coef.glmnet(cv_fit)[-1]
                probtest <- expit_rcpp(prodv_rcpp(cbind(Atest,Ztest,Wtest),cur_beta))
                cur_dev <- -2*sum( ytest * log(probtest) + (1 - ytest) * log(1 - probtest) )
                suppressMessages(cur_auc<-c(auc(ytest,probtest,direction = "<")))
                c(cur_dev,cur_auc)
            })
        }) # row is ratio_range & col is lambda_list

        list("deviance"=do.call(rbind, lapply(dev_lam_ratio_fold, function(m) m[1,])),
             "auc"=do.call(rbind, lapply(dev_lam_ratio_fold, function(m) m[2,])))
    })
    dev_lam_ratio1<-lapply(1:length(index_fold), function(cur_fold){
        dev_lam_ratio[[cur_fold]]$deviance
    })
    dev_lam_ratio2<-lapply(1:length(index_fold), function(cur_fold){
        dev_lam_ratio[[cur_fold]]$auc
    })
    list("deviance"=Reduce(`+`, dev_lam_ratio1)/length(index_fold),
         "auc"=Reduce(`+`, dev_lam_ratio2)/length(index_fold))
}


cv_mse_lambda_Cweight_func<-function(index_fold,Z,W,A,y,family,
                                     C_half,beta_initial,hat_thetaA,
                                     study_info,
                                     weight_list,pZ,pW,pA,
                                     w_adaptive,final_alpha,
                                     pseudo_Xy,lambda.min.ratio,
                                     nlambda,X=NULL,XR=NULL,
                                     fix_lambda_list=NULL,sC_half=NULL){
    if(is.null(X)){X=cbind(A,Z,W)}
    if(is.null(XR)){XR=cbind(A,Z)}
    item_weight_list<-lapply(weight_list,function(weight){
        if(weight<0){
            C_half_weight<-sC_half
        }else{
            C_half_weight<-C_half
        }
        C_half_weight[(pZ+pA+pW+1):nrow(C_half),(pZ+pA+pW+1):nrow(C_half)]=
            C_half[(pZ+pA+pW+1):nrow(C_half),(pZ+pA+pW+1):nrow(C_half)]*sqrt(abs(weight))
        pseudo_Xy_list<-pseudo_Xy(C_half=C_half_weight,Z=Z,W=W,A=A,y=y,
                                  beta=beta_initial,hat_thetaA=hat_thetaA,
                                  study_info=study_info,X=X,XR=XR)

        initial_sf<-nrow(Z)/sqrt(nrow(pseudo_Xy_list$pseudo_X))
        pseudo_X_weight<-pseudo_Xy_list$pseudo_X/initial_sf
        pseudo_y_weight<-pseudo_Xy_list$pseudo_y/initial_sf
        if(!is.null(fix_lambda_list)){
            lambda_list_weight<-fix_lambda_list
        }else{
            innerprod<-crossprodv_rcpp(pseudo_X_weight,pseudo_y_weight)[which(w_adaptive!=0)]
            lambda.max<-max(abs(innerprod))/nrow(pseudo_X_weight)
            lambda_list_weight <-exp(seq(log(lambda.max),log(lambda.max*lambda.min.ratio),
                                         length.out=nlambda))
        }
        list("C_half"=C_half_weight,
             "pseudo_X"=pseudo_X_weight,
             "pseudo_y"=pseudo_y_weight,
             "lambda_list"=lambda_list_weight)
    })


    mse_lam_weight<-lapply(1:length(index_fold), function(cur_fold){
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
        ytrain<-y[-index_test]
        ytest<-y[index_test]
        mse_lam_weight_fold<-sapply(1:length(weight_list),function(weight_id){
            cur_weight=weight_list[weight_id]
            C_half_weight = item_weight_list[[weight_id]]$C_half
            lambda_list_weight = item_weight_list[[weight_id]]$lambda_list
            pseudo_Xy_list_train<-pseudo_Xy(C_half_weight,Ztrain,Wtrain,Atrain,
                                            ytrain,beta = beta_initial,hat_thetaA = hat_thetaA,
                                            study_info=study_info)
            initial_sf_train<-nrow(Ztrain)/sqrt(nrow(pseudo_Xy_list_train$pseudo_X))
            pseudo_X_train<-pseudo_Xy_list_train$pseudo_X/initial_sf_train
            pseudo_y_train<-pseudo_Xy_list_train$pseudo_y/initial_sf_train

            sapply(lambda_list_weight,function(cur_lam){
                cv_fit<-glmnet(x=pseudo_X_train,y=pseudo_y_train,
                               standardize=F,intercept=F,alpha = final_alpha,
                               penalty.factor = w_adaptive,
                               lambda = cur_lam)
                cur_beta<-coef.glmnet(cv_fit)[-1]
                mean(( prodv_rcpp(cbind(Atest,Ztest,Wtest),cur_beta) - ytest)^2)
            })
        }) # row is weight_list & col is lambda_list
        mse_lam_weight_fold
    })
    list("mse"=t(Reduce(`+`, mse_lam_weight)/length(index_fold)),
         "items"=item_weight_list)
}

cv_dev_lambda_Cweight_func<-function(index_fold,Z,W,A,y,family,
                                     C_half,beta_initial,
                                     initial_with_type,
                                     hat_thetaA,
                                     study_info,
                                     weight_list,pZ,pW,pA,
                                     w_adaptive,final_alpha,
                                     pseudo_Xy,lambda.min.ratio,
                                     nlambda,V_thetaA,use_sparseC,
                                     use_offset,V_thetaA_sandwich,
                                     fold_self_beta,
                                     X=NULL,XR=NULL,
                                     fix_lambda_list=NULL,sC_half=NULL){
    if(is.null(X)){X=cbind(A,Z,W)}
    if(is.null(XR)){XR=cbind(A,Z)}
    item_weight_list<-lapply(weight_list,function(weight){
        if(weight<0){
            C_half_weight<-sC_half
        }else{
            C_half_weight<-C_half
        }
        C_half_weight[(pZ+pA+pW+1):nrow(C_half),(pZ+pA+pW+1):nrow(C_half)]<-
            C_half[(pZ+pA+pW+1):nrow(C_half),(pZ+pA+pW+1):nrow(C_half)]*sqrt(abs(weight))
        pseudo_Xy_list<-pseudo_Xy(C_half=C_half_weight,Z=Z,W=W,A=A,y=y,
                                  beta=beta_initial,hat_thetaA=hat_thetaA,
                                  study_info=study_info,X=X,XR=XR)

        initial_sf<-nrow(Z)/sqrt(nrow(pseudo_Xy_list$pseudo_X))
        pseudo_X_weight<-pseudo_Xy_list$pseudo_X/initial_sf
        pseudo_y_weight<-pseudo_Xy_list$pseudo_y/initial_sf

        if(!is.null(fix_lambda_list)){
            lambda_list_weight<-fix_lambda_list
        }else{
            innerprod<-crossprodv_rcpp(pseudo_X_weight,pseudo_y_weight)[which(w_adaptive!=0)]
            lambda.max<-max(abs(innerprod))/nrow(pseudo_X_weight)
            lambda_list_weight <-exp(seq(log(lambda.max),log(lambda.max*lambda.min.ratio),
                                         length.out=nlambda))
        }
        list("C_half"=C_half_weight,
             "pseudo_X"=pseudo_X_weight,
             "pseudo_y"=pseudo_y_weight,
             "lambda_list"=lambda_list_weight)
    })

    dev_lam_weight<-lapply(1:length(index_fold), function(cur_fold){
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
        ytrain<-y[-index_test]
        ytest<-y[index_test]
        dev_lam_weight_fold<-lapply(1:length(weight_list),function(weight_id){
            cur_weight=weight_list[weight_id]
            C_half_weight = item_weight_list[[weight_id]]$C_half
            lambda_list_weight = item_weight_list[[weight_id]]$lambda_list
            pseudo_Xy_list_train<-pseudo_Xy(C_half_weight,Ztrain,Wtrain,Atrain,
                                            ytrain,beta = beta_initial,hat_thetaA = hat_thetaA,
                                            study_info=study_info)
            initial_sf_train<-nrow(Ztrain)/sqrt(nrow(pseudo_Xy_list_train$pseudo_X))
            pseudo_X_train<-pseudo_Xy_list_train$pseudo_X/initial_sf_train
            pseudo_y_train<-pseudo_Xy_list_train$pseudo_y/initial_sf_train

            sapply(lambda_list_weight,function(cur_lam){
                cv_fit<-glmnet(x=pseudo_X_train,y=pseudo_y_train,
                               standardize=F,intercept=F,alpha = final_alpha,
                               penalty.factor = w_adaptive,
                               lambda = cur_lam)
                cur_beta<-coef.glmnet(cv_fit)[-1]
                probtest <- expit_rcpp(prodv_rcpp(cbind(Atest,Ztest,Wtest),cur_beta))
                cur_dev <- -2*sum( ytest * log(probtest) + (1 - ytest) * log(1 - probtest) )
                suppressMessages(cur_auc<-c(auc(ytest,probtest,direction = "<")))
                c(cur_dev,cur_auc)
            })
        }) # row is weight_list & col is lambda_list
        list("deviance"=do.call(rbind, lapply(dev_lam_weight_fold, function(m) m[1,])),
             "auc"=do.call(rbind, lapply(dev_lam_weight_fold, function(m) m[2,])))
    })
    dev_lam_weight1<-lapply(1:length(index_fold), function(cur_fold){
        dev_lam_weight[[cur_fold]]$deviance
    })
    dev_lam_weight2<-lapply(1:length(index_fold), function(cur_fold){
        dev_lam_weight[[cur_fold]]$auc
    })
    list("deviance"=Reduce(`+`, dev_lam_weight1)/length(index_fold),
         "auc"=Reduce(`+`, dev_lam_weight2)/length(index_fold),
         "items"=item_weight_list)
}



cv_dev_lambda_Cweight_func2<-function(index_fold,Z,W,A,y,family,
                                      C_half,beta_initial,
                                      initial_with_type,
                                      hat_thetaA,
                                      study_info,
                                      weight_list,pZ,pW,pA,
                                      w_adaptive,final_alpha,
                                      pseudo_Xy,lambda.min.ratio,
                                      nlambda,V_thetaA,use_sparseC,
                                      use_offset,V_thetaA_sandwich,
                                      fold_self_beta,
                                      X=NULL,XR=NULL,
                                      fix_lambda_list=NULL,sC_half=NULL){
    if(is.null(X)){X=cbind(A,Z,W)}
    if(is.null(XR)){XR=cbind(A,Z)}
    item_weight_list<-lapply(weight_list,function(weight){
        if(weight<0){
            C_half_weight<-sC_half
        }else{
            C_half_weight=C_half
            C_half_weight[(pZ+pA+pW+1):nrow(C_half),(pZ+pA+pW+1):nrow(C_half)]<-
                C_half[(pZ+pA+pW+1):nrow(C_half),(pZ+pA+pW+1):nrow(C_half)]*sqrt(weight)
        }
        pseudo_Xy_list<-pseudo_Xy(C_half=C_half_weight,Z=Z,W=W,A=A,y=y,
                                  beta=beta_initial,hat_thetaA=hat_thetaA,
                                  study_info=study_info,X=X,XR=XR)

        initial_sf<-nrow(Z)/sqrt(nrow(pseudo_Xy_list$pseudo_X))
        pseudo_X_weight<-pseudo_Xy_list$pseudo_X/initial_sf
        pseudo_y_weight<-pseudo_Xy_list$pseudo_y/initial_sf

        if(!is.null(fix_lambda_list)){
            lambda_list_weight<-fix_lambda_list
        }else{
            innerprod<-crossprodv_rcpp(pseudo_X_weight,pseudo_y_weight)[which(w_adaptive!=0)]
            lambda.max<-max(abs(innerprod))/nrow(pseudo_X_weight)
            lambda_list_weight <-exp(seq(log(lambda.max),log(lambda.max*lambda.min.ratio),
                                         length.out=nlambda))
        }
        list("C_half"=C_half_weight,
             "pseudo_X"=pseudo_X_weight,
             "pseudo_y"=pseudo_y_weight,
             "lambda_list"=lambda_list_weight)
    })

    dev_lam_weight<-lapply(1:length(index_fold), function(cur_fold){
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
        ytrain<-y[-index_test]
        ytest<-y[index_test]
        if(fold_self_beta){
            initial_res1<-beta_initial_func(ytrain,cbind(Atrain,Ztrain,Wtrain),
                                            Atrain,pA,family,initial_with_type,w_adaptive)
            beta_initial1<-initial_res1$beta_initial
        }else{
            beta_initial1<-beta_initial
        }

        thetaA_list<-thetaA_func(pA,Ztrain,Atrain,ytrain,study_info,
                                 family,use_offset,
                                 V_thetaA_sandwich)
        hat_thetaA1<-thetaA_list$hat_thetaA
        V_thetaA1<-thetaA_list$V_thetaA

        inv_C_train = Delta_opt_rcpp(y=ytrain,Z=Ztrain,W=Wtrain,
                                     family=family,
                                     study_info=study_info,
                                     A=Atrain,pA=pA,pZ=pZ,beta=beta_initial1,
                                     hat_thetaA=hat_thetaA1,
                                     V_thetaA=V_thetaA1,
                                     use_offset=use_offset)
        sC_half_train<-diag(1/sqrt(diag(inv_C_train)))
        if(use_sparseC){C_half_train<-sC_half_train
        }else{C_half_train<-sqrtchoinv_rcpp(inv_C_train+diag(1e-15,nrow(inv_C_train)))}

        dev_lam_weight_fold<-lapply(1:length(weight_list),function(weight_id){
            cur_weight=weight_list[weight_id]
            if(cur_weight<0){
                C_half_weight<-sC_half_train
            }else{
                C_half_weight<-C_half_train
            }
            C_half_weight[(pZ+pA+pW+1):nrow(C_half_weight),(pZ+pA+pW+1):nrow(C_half_weight)]<-
                C_half_weight[(pZ+pA+pW+1):nrow(C_half_weight),(pZ+pA+pW+1):nrow(C_half_weight)]*sqrt(abs(cur_weight))


            #C_half_weight = item_weight_list[[weight_id]]$C_half
            lambda_list_weight = item_weight_list[[weight_id]]$lambda_list
            pseudo_Xy_list_train<-pseudo_Xy(C_half_weight,Ztrain,Wtrain,Atrain,
                                            ytrain,beta = beta_initial1,hat_thetaA = hat_thetaA,
                                            study_info=study_info)
            initial_sf_train<-nrow(Ztrain)/sqrt(nrow(pseudo_Xy_list_train$pseudo_X))
            pseudo_X_train<-pseudo_Xy_list_train$pseudo_X/initial_sf_train
            pseudo_y_train<-pseudo_Xy_list_train$pseudo_y/initial_sf_train

            sapply(lambda_list_weight,function(cur_lam){
                cv_fit<-glmnet(x=pseudo_X_train,y=pseudo_y_train,
                               standardize=F,intercept=F,alpha = final_alpha,
                               penalty.factor = w_adaptive,
                               lambda = cur_lam)
                cur_beta<-coef.glmnet(cv_fit)[-1]
                probtest <- expit_rcpp(prodv_rcpp(cbind(Atest,Ztest,Wtest),cur_beta))
                cur_dev <- -2*sum( ytest * log(probtest) + (1 - ytest) * log(1 - probtest) )
                suppressMessages(cur_auc<-c(auc(ytest,probtest,direction = "<")))
                c(cur_dev,cur_auc)
            })
        }) # row is weight_list & col is lambda_list
        list("deviance"=do.call(rbind, lapply(dev_lam_weight_fold, function(m) m[1,])),
             "auc"=do.call(rbind, lapply(dev_lam_weight_fold, function(m) m[2,])))
    })
    dev_lam_weight1<-lapply(1:length(index_fold), function(cur_fold){
        dev_lam_weight[[cur_fold]]$deviance
    })
    dev_lam_weight2<-lapply(1:length(index_fold), function(cur_fold){
        dev_lam_weight[[cur_fold]]$auc
    })
    list("deviance"=Reduce(`+`, dev_lam_weight1)/length(index_fold),
         "auc"=Reduce(`+`, dev_lam_weight2)/length(index_fold),
         "items"=item_weight_list)
}

cv_dev_lambda_Cweight_func4<-function(index_fold,Z,W,A,y,family,
                                      C_half,beta_initial,
                                      initial_with_type,
                                      hat_thetaA,
                                      study_info,
                                      weight_list,pZ,pW,pA,
                                      w_adaptive,final_alpha,
                                      pseudo_Xy,lambda.min.ratio,
                                      nlambda,V_thetaA,use_sparseC,
                                      use_offset,V_thetaA_sandwich,
                                      fold_self_beta,
                                      X=NULL,XR=NULL,
                                      fix_lambda_list=NULL,sC_half=NULL){
    sample_p<-0.3
    index_test<-c(sample(which(y==1),round(sum(y)*sample_p)),
                  sample(which(y==0),round(sum(1-y)*sample_p)))

    #index_test<-Reduce(c,lapply(1:round(length(index_fold)*0.3), function(i){index_fold[[i]]}))

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
    ytrain<-y[-index_test]
    ytest<-y[index_test]

    # renew the initial value
    if(fold_self_beta){
        initial_res1<-beta_initial_func(ytrain,cbind(Atrain,Ztrain,Wtrain),
                                        Atrain,pA,family,initial_with_type,w_adaptive)
        beta_initial1<-initial_res1$beta_initial
    }else{
        beta_initial1<-beta_initial
    }

    thetaA_list<-thetaA_func(pA,Ztrain,Atrain,ytrain,study_info,
                             family,use_offset,
                             V_thetaA_sandwich)
    hat_thetaA1<-thetaA_list$hat_thetaA
    V_thetaA1<-thetaA_list$V_thetaA

    inv_C_train = Delta_opt_rcpp(y=ytrain,Z=Ztrain,W=Wtrain,
                                 family=family,
                                 study_info=study_info,
                                 A=Atrain,pA=pA,pZ=pZ,beta=beta_initial1,
                                 hat_thetaA=hat_thetaA1,
                                 V_thetaA=V_thetaA1,
                                 use_offset=use_offset)
    sC_half_train<-diag(1/sqrt(diag(inv_C_train)))
    if(use_sparseC){C_half_train<-sC_half_train
    }else{C_half_train<-sqrtchoinv_rcpp(inv_C_train+diag(1e-15,nrow(inv_C_train)))}

    dev_lam_weight_fold<-lapply(1:length(weight_list),function(weight_id){
        cur_weight=weight_list[weight_id]
        if(cur_weight<0){
            C_half_weight<-sC_half_train
        }else{
            C_half_weight<-C_half_train
        }
        C_half_weight[(pZ+pA+pW+1):nrow(C_half_weight),(pZ+pA+pW+1):nrow(C_half_weight)]<-
            C_half_weight[(pZ+pA+pW+1):nrow(C_half_weight),(pZ+pA+pW+1):nrow(C_half_weight)]*sqrt(abs(cur_weight))

        pseudo_Xy_list_train<-pseudo_Xy(C_half_weight,Ztrain,Wtrain,Atrain,
                                        ytrain,beta = beta_initial1,hat_thetaA = hat_thetaA1,
                                        study_info=study_info)
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
            probtest <- expit_rcpp(prodv_rcpp(cbind(Atest,Ztest,Wtest),cur_beta))
            cur_dev <- -2*sum( ytest * log(probtest) + (1 - ytest) * log(1 - probtest) )
            suppressMessages(cur_auc<-c(auc(ytest,probtest,direction = "<")))
            c(cur_dev,cur_auc)
        })
    })
    # row is weight_list & col is lambda_list
    cv_dev1<-do.call(rbind, lapply(dev_lam_weight_fold, function(m) m[1,]))
    cv_auc1<-do.call(rbind, lapply(dev_lam_weight_fold, function(m) m[2,]))


    ids_dev<-which(cv_dev1==min(cv_dev1),arr.ind = TRUE)
    ids_auc<-which(cv_auc1==max(cv_auc1),arr.ind = TRUE)

    final_weight<-weight_list[ids_auc[1,1]] #nrow(ids_auc)
    final_weight_dev<-weight_list[ids_dev[1,1]] #nrow(ids_dev)

    return(list("final_weight"=final_weight,
                "final_weight_dev"=final_weight_dev,
                "deviance"=cv_dev1,
                "auc"=cv_auc1))
}






## cross validation function for continuous y with lambda and ratio
cv_auc_ext<-function(index_fold,Z,W,A,y,study_info,hat_thetaA){
    sapply(1:length(index_fold), function(cur_fold){
        index_test<-index_fold[[cur_fold]]
        Ztrain<-Z[-index_test,,drop=FALSE]
        Ztest<-Z[index_test,,drop=FALSE]
        if(!is.null(A)){
            Atrain<-A[-index_test,,drop=FALSE]
            Atest<-A[index_test,,drop=FALSE]
        }else{
            Atrain<-NULL
            Atest<-NULL}
        ytrain<-y[-index_test]
        ytest<-y[index_test]
        hat_theta<-c(hat_thetaA,study_info[[1]]$Coeff)
        probtest <- expit_rcpp(prodv_rcpp(cbind(Atest,Ztest),hat_theta))
        suppressMessages(cur_auc<-c(auc(ytest,probtest,direction = "<")))
        cur_auc
    })
}





htlgmm.default<-function(
        y,Z,W=NULL,
        study_info=NULL,
        A=1,
        penalty_type = "lasso",
        family = "gaussian",
        initial_with_type = "ridge",
        beta_initial = NULL,
        alpha = NULL,
        hat_thetaA = NULL,
        V_thetaA = NULL,
        use_offset = TRUE,
        V_thetaA_sandwich = FALSE,
        remove_penalty_Z = FALSE,
        remove_penalty_W = FALSE,
        inference = TRUE,
        fix_C = NULL,
        refine_C = FALSE,
        sqrt_matrix ="cholesky",
        use_cv = TRUE,
        type_measure = "default",
        nfolds = 10,
        foldid = NULL,
        fix_lambda = NULL,
        lambda_list = NULL,
        nlambda = 100,
        lambda.min.ratio = 0.0001,
        tune_ratio = FALSE,
        fix_ratio = NULL,
        ratio_list = NULL,
        tune_weight = FALSE,
        fix_weight = NULL,
        weight_list = NULL,
        tune_weight_method = "holdout",
        gamma_adaptivelasso = 1/2,
        use_sparseC = TRUE,
        seed.use = 97,
        output_all_betas=FALSE
){

    ## Initial Preparation
    set.seed(seed.use)
    if(sum(is.na(y))>0){stop("y includes NA")}
    if(sum(is.na(Z))>0){stop("Z includes NA")}
    if(!is.null(W)&sum(is.na(W))>0){stop("W includes NA")}
    if(!is.null(A)&sum(is.na(A))>0){stop("A includes NA")}
    if (is.null(study_info)){stop("Please input study_info as trained model")}
    if(!penalty_type %in% c("none","adaptivelasso","lasso","ridge","elasticnet")){
        stop("Select penalty type from c('none','adaptivelasso','lasso','ridge','elasticnet').")
    }
    if(!sqrt_matrix %in% c("cholesky","svd")){
        stop("Select penalty type from c('cholesky','svd').")
    }
    if(!type_measure%in% c("default", "mse", "deviance", "auc")){
        stop("Select type_measure from c('default','mse','deviance','auc'). When family == 'gaussian', type_measure is 'mse' no matter which input is given. When family == 'binomial', default is set to be 'auc'.")
    }
    if(tune_weight){
        if(!tune_weight_method%in%c("raw","holdout","exact","1se")){
            stop("Select weight tuning method from the following ones: 'raw','holdout','exact','1se'. ")
        }
    }
    if(is.null(dim(Z)[1])){
        warning("Z is input as a vector, convert Z into matrix with size nZ*1")
        Z=matrix(Z,ncol=1)
    }

    ###########--------------###########
    # determine which kind of penalty for the final model
    if(!is.null(alpha) & penalty_type != "elasticnet"){
        stop("When using alpha between 0 and 1, please set penalty_type to be 'elasticnet'.")
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

    ###########--------------###########
    # deal with in consistence for ratio and weight
    if(tune_weight & tune_ratio){
        stop("Not support tune weight and ratio together. One can set fix_weight or fix_ratio.")
    }
    if(!is.null(fix_ratio)){
        if(use_cv&tune_ratio){
            tune_ratio<-FALSE
            warning("Ratio is fixed, set tune_ratio as FALSE")
        }else if(remove_penalty_Z | remove_penalty_W){
            remove_penalty_Z<-FALSE
            remove_penalty_W<-FALSE
            warning("Ratio is fixed, set remove_penalty's as FALSE")
        }
        if(!is.numeric(fix_ratio)|fix_ratio<0){
            stop("fix_ratio should be >=0.")
        }
    }

    if(!is.null(fix_weight)){
        if(use_cv&tune_weight){
            tune_ratio<-FALSE
            warning("Weight is fixed, set tune_weight as FALSE")
        }
    }

    ###########--------------###########
    # define dimensions for A,Z,W
    nZ=nrow(Z)
    if(length(study_info[[1]])==3){
        nZext=study_info[[1]]$Sample_size
    }else{nZext=NULL}
    pZ=ncol(Z)
    if(is.null(W)){pW=0}else{pW=ncol(W)}
    if(is.null(A)){pA=0}else{
        if(is.null(dim(A)[1])){
            if(length(A)==1){
                if(A==1){A=matrix(1,nrow=nZ,ncol=1)}
            }
        }
        pA=ncol(A)
    }
    if(nZ<2*pZ+pW+pA){use_sparseC=TRUE}

    if(family == "gaussian"){pseudo_Xy=pseudo_Xy_gaussian_rcpp
    }else if(family == "binomial"){pseudo_Xy=pseudo_Xy_binomial_rcpp}

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
    X=cbind(A,Z,W)
    XR=cbind(A,Z)

    ###########--------------###########
    # compute hat_thetaA and V_thetaA

    # if(pA!=0){
    #     if(is.null(hat_thetaA)){
    #         if(!is.null(V_thetaA)){
    #             stop("With customized hat_thetaA input, V_thetaA is also needed")
    #         }
    #         if(use_offset){
    #             offset_term = prodv_rcpp(Z,study_info[[1]]$Coeff)
    #             df=data.frame(y,A)
    #             if(family=="binomial"){
    #                 hat_thetaA_glm=speedglm(y~0+.,data = df,offset = offset_term,family = binomial())
    #             }else if(family=="gaussian"){
    #                 hat_thetaA_glm=speedlm(y~0+.,data = df,offset = offset_term)
    #                 #hat_thetaA_glm=lm(y~0+.,data = df,offset = offset_term)
    #             }
    #             hat_thetaA=hat_thetaA_glm$coefficients
    #             if(V_thetaA_sandwich){
    #                 V_thetaA=vcov_sandwich_rcpp(y=y,A=A,Z=Z,family=family,
    #                                             study_info=study_info,pA=pA,
    #                                             hat_thetaA=hat_thetaA,
    #                                             use_offset=use_offset,
    #                                             XR=XR)
    #             }else{V_thetaA=vcov(hat_thetaA_glm)}
    #         }else{
    #             df=data.frame(y,A,Z)
    #             if(family=="binomial"){
    #                 hat_thetaA_glm=speedglm(y~0+.,data = df,family = binomial())
    #             }else if(family=="gaussian"){
    #                 hat_thetaA_glm=speedlm(y~0+.,data = df)
    #             }
    #             hat_thetaA=hat_thetaA_glm$coefficients[1:pA]
    #             if(V_thetaA_sandwich){
    #                 V_thetaA=vcov_sandwich_rcpp(y,A,Z,family,study_info,pA,
    #                                        hat_thetaA,use_offset)
    #             }else{V_thetaA=vcov(hat_thetaA_glm)[1:pA,1:pA,drop=F]}
    #         }
    #     }
    #
    #     if(is.null(dim(V_thetaA)[1])){
    #         V_thetaA = as.matrix(V_thetaA,nrow=pA,ncol=pA)
    #     }
    # }

    thetaA_list<-thetaA_func(pA,Z,A,y,study_info,
                             family,use_offset,
                             V_thetaA_sandwich,
                             hat_thetaA,V_thetaA)
    hat_thetaA<-thetaA_list$hat_thetaA
    V_thetaA<-thetaA_list$V_thetaA

    ###########--------------###########
    # define penalty assignment

    fix_penalty<-c(rep(0,pA),rep(1,pZ+pW))
    if(remove_penalty_Z){fix_penalty[Zid]<-0}else{
        if(!is.null(fix_ratio)){fix_penalty[Zid]<-fix_ratio}
    }
    if(remove_penalty_W){fix_penalty[Wid]<-0}
    if((remove_penalty_Z & remove_penalty_W)|(length(unique(fix_penalty))==1 & unique(fix_penalty)[1] == 0) ){
        penalty_type = "none"
        warning("All penalties are removed, turn to no penalties!")
    }

    ###########--------------###########
    # compute initial beta if not given
    beta_initial1se<-NULL
    if(penalty_type == "none"){
        initial_with_type = "glm"
        use_cv = FALSE
        tune_ratio = FALSE
    }
    if(!is.null(beta_initial) & length(beta_initial)!=pA+pZ+pW){
        warning("beta_initial should be from A,Z,W.\n Length not match, compute default initial instead.")
        beta_initial=NULL
    }
    if(is.null(beta_initial)){
        initial_res<-beta_initial_func(y,X,A,pA,
                                    family,
                                    initial_with_type,
                                    fix_penalty)
        beta_initial<-initial_res$beta_initial
        beta_initial1se<-initial_res$beta_initial1se
    }
    if (penalty_type == "adaptivelasso"){
        w_adaptive<-1/abs(beta_initial)^gamma_adaptivelasso
        w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*100
        w_adaptive<-w_adaptive*fix_penalty
    }else{w_adaptive<-fix_penalty}

    ###########--------------###########
    # estimation of C

    if(is.null(fix_C)){

        inv_C = Delta_opt_rcpp(y=y,Z=Z,W=W,
                               family=family,
                               study_info=study_info,
                               A=A,pA=pA,pZ=pZ,beta=beta_initial,
                               hat_thetaA=hat_thetaA,
                               V_thetaA=V_thetaA,
                               use_offset = use_offset,
                               X=X,XR=XR)

        sC_half<-diag(1/sqrt(diag(inv_C)))
        if(use_sparseC){
            C_half<-sC_half
        }else{
            if(sqrt_matrix =="svd"){
                inv_C_svd=fast.svd(inv_C+diag(1e-15,nrow(inv_C)))
                C_half=prod_rcpp(inv_C_svd$v,(t(inv_C_svd$u)*1/sqrt(inv_C_svd$d)))
            }else if(sqrt_matrix =="cholesky"){
                C_half<-sqrtchoinv_rcpp(inv_C+diag(1e-15,nrow(inv_C)))
            }
        }
    }else{
        if(nrow(fix_C)!=pA+pZ+pW+pZ){
            stop("Input fix_C dimension is wrong!")}
        C_half<-sqrtcho_rcpp(fix_C+diag(1e-15,nrow(fix_C)))
    }

    if(!is.null(fix_weight)){
        if(!is.numeric(fix_weight)|fix_weight<0){
            stop("fix_weight should be >=0.")
        }
        C_half[(pZ+pA+pW+1):nrow(C_half),(pZ+pA+pW+1):nrow(C_half)]=
            C_half[(pZ+pA+pW+1):nrow(C_half),(pZ+pA+pW+1):nrow(C_half)]*sqrt(fix_weight)
    }
    ###########--------------###########
    # Prepare for final model

    pseudo_Xy_list<-pseudo_Xy(C_half=C_half,Z=Z,W=W,A=A,y=y,
                              beta=beta_initial,hat_thetaA=hat_thetaA,
                              study_info=study_info,X=X,XR=XR)

    initial_sf<-nZ/sqrt(nrow(pseudo_Xy_list$pseudo_X))
    pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
    pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf

    ###########--------------###########
    # generate lambda list from glmnet

    if(!is.null(lambda_list)){fix_lambda_list<-lambda_list}else{
        fix_lambda_list<-NULL}
    if(penalty_type != "none" & is.null(fix_lambda)&is.null(lambda_list)){
        # if(penalty_type == "adaptivelasso"){
        #     fit_final<-glmnet(x= pseudo_X,y= pseudo_y,standardize=F,
        #                       intercept=F,alpha = final_alpha,penalty.factor = w_adaptive)
        #     lambda_list<-fit_final$lambda
        #     lambda_list<-lambda_list[!is.na(lambda_list)]
        # }else{
            innerprod<-crossprodv_rcpp(pseudo_X,pseudo_y)[which(fix_penalty!=0)]
            lambda.max<-max(abs(innerprod))/nrow(pseudo_X)
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

    ###########--------------###########
    # generate ratio list from glmnet

    if(tune_ratio & !remove_penalty_Z & !remove_penalty_W){
        if(is.null(ratio_list)){
            if(is.null(nZext)){
                warning("No sample size is in study_info, ratio estimation will be bad.")
                nZext=nZ
            }
            # ratio_lower<-sqrt(nZ/(nZ+nZext))/2
            # ratio_upper<-(nZ)^(1/3)/2
            # ratio_count<-10
            # ratio_list<-(seq(sqrt(ratio_lower),sqrt(ratio_upper),(sqrt(ratio_upper)-sqrt(ratio_lower))/ratio_count)^2)
            # ratio_list<-c(1,ratio_list)
            ratio_list<-c(1,1/2,2,4)
        }
    }else{tune_ratio<-FALSE}

    ###########--------------###########
    # generate weight list from glmnet
    if(tune_weight){
        if(is.null(weight_list)){
            weight_list<-c(1,2,4,8,16)
            if(!use_sparseC){weight_list<-c(weight_list,-weight_list)}
        }
    }

    ###########--------------###########
    # No Cross Validation

    if(!use_cv){
        if(penalty_type == "none"){
            fit_final_ols=lm(y~0+.,data = data.frame(y= pseudo_y,pseudo_X))
            beta=fit_final_ols$coefficients
            return_list<-list("beta"=beta)
        }else{
            if(!is.null(fix_lambda)){
                fit_final_fixed_lambda<-glmnet(x= pseudo_X,y= pseudo_y,standardize=F,
                                               intercept=F,alpha = final_alpha,
                                               penalty.factor = w_adaptive,
                                               lambda = fix_lambda)
                beta<-coef.glmnet(fit_final_fixed_lambda)[-1]
                return_list<-list("beta"=beta,
                                  "fix_lambda"=fix_lambda)
            }else{
                fit_final_lambda_list<-glmnet(x= pseudo_X,y= pseudo_y,standardize=F,
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

    ###########--------------###########
    # Cross Validation: Just for lambda; lambda and weight; lambda and ratio
    if(use_cv){
        if(is.null(foldid)){
            if(length(unique(y)) <= 2){
                index_fold<-createFolds(as.numeric(y>0),k = nfolds)
            }else{index_fold<-createFolds(y,k = nfolds)}
        }else{
            uni_fold=sort(unique(foldid))
            nfolds=length(uni_fold)
            index_fold<-lapply(1:length(uni_fold),function(i){which(foldid==uni_fold[i])})
        }

        ###########--------------###########
        # Mode 1: lambda and weight

        if(tune_weight & !tune_ratio){
            if(family == "gaussian"){
                cv_res<-cv_mse_lambda_Cweight_func(index_fold,Z,W,A,y,family,
                                                   C_half,beta_initial,hat_thetaA,
                                                   study_info,
                                                   weight_list,pZ,pW,pA,
                                                   w_adaptive,final_alpha,
                                                   pseudo_Xy,lambda.min.ratio,
                                                   nlambda,X,XR,
                                                   fix_lambda_list,sC_half)
                cv_mse<-cv_res$mse
                ids<-which(cv_mse==min(cv_mse),arr.ind = TRUE)[1,]
                final.weight.min<-weight_list[ids[1]]
                final.lambda.min<-cv_res$items[[ids[1]]]$lambda_list[ids[2]]
                return_list<-list("cv_mse"=cv_mse)
            }else if(family == "binomial"){
                if(tune_weight_method == "holdout"){
                    cv_dev_lambda_Cweight_func = cv_dev_lambda_Cweight_func4
                }else if(tune_weight_method == "exact"){
                    cv_dev_lambda_Cweight_func = cv_dev_lambda_Cweight_func2
                }else if(tune_weight_method == "1se"){
                    if(!is.null(beta_initial1se)){beta_initial<-beta_initial1se}
                }
                fold_self_beta = TRUE
                cv_res<-cv_dev_lambda_Cweight_func(index_fold,Z,W,A,y,family,
                                                   C_half,beta_initial,
                                                   initial_with_type,
                                                   hat_thetaA,
                                                   study_info,
                                                   weight_list,pZ,pW,pA,
                                                   w_adaptive,final_alpha,
                                                   pseudo_Xy,lambda.min.ratio,
                                                   nlambda,V_thetaA,use_sparseC,
                                                   use_offset,V_thetaA_sandwich,
                                                   fold_self_beta,
                                                   X,XR,fix_lambda_list,sC_half)
                cv_auc<-cv_res$auc
                ids_auc<-which(cv_auc==max(cv_auc),arr.ind = TRUE)
                ids_auc<-ids_auc[1,] #nrow(ids_auc)
                cv_dev<-cv_res$deviance
                ids_dev<-which(cv_dev==min(cv_dev),arr.ind = TRUE)
                ids_dev<-ids_dev[1,] #nrow(ids_dev)
                print(paste0("auc_weightid:",ids_auc[1]))

                final.ratio.min<-1

                if(tune_weight_method == "holdout"){
                    weight<-cv_res$final_weight

                    if(weight<0){
                        C_half<-sC_half
                    }

                    C_half[(pZ+pA+pW+1):nrow(C_half),(pZ+pA+pW+1):nrow(C_half)]<-
                        C_half[(pZ+pA+pW+1):nrow(C_half),(pZ+pA+pW+1):nrow(C_half)]*sqrt(abs(weight))

                    if(weight!=1){
                        pseudo_Xy_list<-pseudo_Xy(C_half=C_half,Z=Z,W=W,A=A,y=y,
                                                  beta=beta_initial,hat_thetaA=hat_thetaA,
                                                  study_info=study_info,X=X,XR=XR)

                        initial_sf<-nrow(Z)/sqrt(nrow(pseudo_Xy_list$pseudo_X))
                        pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
                        pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
                        if(!is.null(fix_lambda_list)){
                            lambda_list<-fix_lambda_list
                        }else{
                            innerprod<-crossprodv_rcpp(pseudo_X,pseudo_y)[which(w_adaptive!=0)]
                            lambda.max<-max(abs(innerprod))/nrow(pseudo_X)
                            lambda_list <-exp(seq(log(lambda.max),log(lambda.max*lambda.min.ratio),
                                                         length.out=nlambda))
                        }
                    }
                    pseudo_X_dev<-pseudo_X
                    pseudo_y_dev<-pseudo_y
                    final.weight.min<-weight

                    res_weight<-cv_dev_lambda_func(index_fold,Z,W,A,y,
                                       C_half,beta_initial,hat_thetaA,
                                       study_info,lambda_list,
                                       w_adaptive,final_alpha,pseudo_Xy)

                    cv_auc1<-res_weight$auc
                    max_id<-which.max(cv_auc1)
                    final.lambda.min<-lambda_list[max_id]

                    cv_dev1<-res_weight$deviance
                    min_id<-which.min(cv_dev1)
                    final.lambda.dev.min<-lambda_list[min_id]
                    final.weight.dev.min<-weight

                    return_list<-list("cv_auc"=cv_auc,"cv_dev"=cv_dev)
                    final.ratio.dev.min<-1
                    return_list<-c(return_list,
                                   list("lambda_list"=lambda_list,
                                        "weight_list"=weight_list))


                }else{
                    pseudo_X<-cv_res$items[[ids_auc[1]]]$pseudo_X
                    pseudo_y<-cv_res$items[[ids_auc[1]]]$pseudo_y

                    pseudo_X_dev<-cv_res$items[[ids_dev[1]]]$pseudo_X
                    pseudo_y_dev<-cv_res$items[[ids_dev[1]]]$pseudo_y

                    final.weight.min<-weight_list[ids_auc[1]]
                    final.lambda.min<-cv_res$items[[ids_auc[1]]]$lambda_list[ids_auc[2]]
                    final.weight.dev.min<-weight_list[ids_dev[1]]
                    final.lambda.dev.min<-cv_res$items[[ids_dev[1]]]$lambda_list[ids_dev[2]]

                    return_list<-list("cv_auc"=cv_auc,"cv_dev"=cv_dev)

                    final.ratio.dev.min<-1
                    return_list<-c(return_list,
                                   list("lambda_list"=
                                            sapply(1:length(weight_list), function(i){
                                                cv_res$items[[i]]$lambda_list}),
                                        "weight_list"=weight_list))

                    if(output_all_betas){
                        all_betas<-lapply(1:length(weight_list), function(i){
                            if(family == "gaussian"){cv_here<-cv_mse[i,]
                            }else{cv_here<- -cv_auc[i,]}
                            ids1<-which.min(cv_here)[1]
                            final.lambda.min1<-cv_res$items[[i]]$lambda_list[ids1]

                            fit_final_lam_weight1<-glmnet(x= cv_res$items[[i]]$pseudo_X,
                                                          y= cv_res$items[[i]]$pseudo_y,
                                                          standardize=F,intercept=F,
                                                          alpha = final_alpha,
                                                          penalty.factor = w_adaptive,
                                                          lambda = final.lambda.min1)
                            beta1<-coef.glmnet(fit_final_lam_weight1)[-1]
                        })
                        return_list<-c(return_list,list("all_betas"=all_betas))
                    }
                }
            }

        }

        ###########--------------###########
        # Mode 2: lambda and ratio


        if(!tune_weight & tune_ratio){
            if(family == "gaussian"){
                cv_mse<-cv_mse_lambda_ratio_func(index_fold,Z,W,A,y,
                                                 C_half,beta_initial,hat_thetaA,
                                                 study_info,lambda_list,
                                                 ratio_list,pZ,pW,pA,
                                                 w_adaptive,final_alpha,pseudo_Xy)
                ids<-which(cv_mse==min(cv_mse),arr.ind = TRUE)
                final.ratio.min<-ratio_list[ids[1]]
                final.lambda.min<-lambda_list[ids[2]]
                return_list<-list("cv_mse"=cv_mse)
            }else if(family == "binomial"){
                cv_res<-cv_dev_lambda_ratio_func(index_fold,Z,W,A,y,
                                                 C_half,beta_initial,hat_thetaA,
                                                 study_info,lambda_list,
                                                 ratio_list,pZ,pW,pA,
                                                 w_adaptive,final_alpha,pseudo_Xy)
                cv_auc<-cv_res$auc
                ids_auc<-which(cv_auc==max(cv_auc),arr.ind = TRUE)[1,]
                final.ratio.min<-ratio_list[ids_auc[1]]
                final.lambda.min<-lambda_list[ids_auc[2]]
                cv_dev<-cv_res$deviance
                ids_dev<-which(cv_dev==min(cv_dev),arr.ind = TRUE)[1,]
                final.ratio.dev.min<-ratio_list[ids_dev[1]]
                final.lambda.dev.min<-lambda_list[ids_dev[2]]
                return_list<-list("cv_auc"=cv_auc,"cv_dev"=cv_dev)
                final.weight.dev.min<-1
            }
            final.weight.min<-1
            return_list<-c(return_list,
                           list("lambda_list"=lambda_list,
                                "ratio_list"=ratio_list))

            if(output_all_betas){
                all_betas<-lapply(1:length(ratio_list), function(i){
                    if(family == "gaussian"){cv_here<-cv_mse[i,]
                    }else{cv_here<- -cv_auc[i,]}
                    ids1<-which.min(cv_here)[1]
                    final.lambda.min1<-lambda_list[ids1]
                    ratio_vec1<-c(rep(1,pA),rep(ratio_list[i],pZ),rep(1,pW))
                    w_adaptive_ratio1<-w_adaptive*ratio_vec1

                    fit_final_lam_ratio1<-glmnet(x= pseudo_X,
                                                  y= pseudo_y,
                                                  standardize=F,intercept=F,
                                                  alpha = final_alpha,
                                                  penalty.factor = w_adaptive_ratio1,
                                                  lambda = final.lambda.min1)
                    beta1<-coef.glmnet(fit_final_lam_ratio1)[-1]
                })
                return_list<-c(return_list,list("all_betas"=all_betas))
            }

        }


        ###########--------------###########
        # Mode 3: lambda only

        if(!tune_ratio & !tune_weight){
            if(family == "gaussian"){
                cv_mse<-cv_mse_lambda_func(index_fold,Z,W,A,y,
                                           C_half,beta_initial,hat_thetaA,
                                           study_info,lambda_list,
                                           w_adaptive,final_alpha,pseudo_Xy)
                final.lambda.min<-lambda_list[which.min(cv_mse)]
                return_list<-list("cv_mse"=cv_mse)
            }else if(family == "binomial"){
                cv_res<-cv_dev_lambda_func(index_fold,Z,W,A,y,
                                           C_half,beta_initial,hat_thetaA,
                                           study_info,lambda_list,
                                           w_adaptive,final_alpha,pseudo_Xy)
                cv_auc<-cv_res$auc
                #cv_auc_sd<-cv_dev$auc_sd
                max_id<-which.max(cv_auc)
                final.lambda.min<-lambda_list[max_id]

                #cv_dev_sd<-cv_dev$deviance_sd
                cv_dev<-cv_res$deviance
                min_id<-which.min(cv_dev)
                final.lambda.dev.min<-lambda_list[min_id]
                return_list<-list("cv_auc"=cv_auc,"cv_dev"=cv_dev)

                final.ratio.dev.min<-1
                final.weight.dev.min<-1
            }
            final.ratio.min<-1
            final.weight.min<-1

            return_list<-c(return_list,
                           list("lambda_list"=lambda_list))

        }


        ratio_vec<-c(rep(1,pA),rep(final.ratio.min,pZ),rep(1,pW))
        w_adaptive_ratio<-w_adaptive*ratio_vec

        ###########--------------###########
        # Model final: use best lambda, best ratio, best weight to build model

        fit_final_lam<-glmnet(x= pseudo_X,y= pseudo_y,standardize=F,
                                    intercept=F,alpha = final_alpha,
                                    penalty.factor = w_adaptive_ratio,
                                    lambda = final.lambda.min)
        beta<-coef.glmnet(fit_final_lam)[-1]

        return_list<-c(list("beta"=beta,
                            "lambda_min"=final.lambda.min,
                            "ratio_min"=final.ratio.min,
                            "weight_min"=final.weight.min),
                       return_list)

        if(family == "binomial"){
            if(tune_weight){
                fit_final_lam_dev<-glmnet(x= pseudo_X_dev,y= pseudo_y_dev,standardize=F,
                                          intercept=F,alpha = final_alpha,
                                          penalty.factor = w_adaptive_ratio,
                                          lambda = final.lambda.dev.min)
            }else{
                fit_final_lam_dev<-glmnet(x= pseudo_X,y= pseudo_y,standardize=F,
                                          intercept=F,alpha = final_alpha,
                                          penalty.factor = w_adaptive_ratio,
                                          lambda = final.lambda.dev.min)
            }
            beta_dev=coef.glmnet(fit_final_lam_dev)[-1]
            return_list<-c(return_list,
                           list("beta_dev"=beta_dev,
                                "lambda_dev_min"=final.lambda.dev.min,
                                "ratio_dev_min"=final.ratio.dev.min,
                                "weight_dev_min"=final.weight.dev.min))
        }

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
            inv_C = Delta_opt_rcpp(y=y,Z=Z,W=W,
                                   family=family,
                                   study_info=study_info,
                                   A=A,pA=pA,pZ=pZ,beta=beta,
                                   hat_thetaA=hat_thetaA,
                                   V_thetaA = V_thetaA,
                                   use_offset = use_offset,
                                   X=X,XR=XR)
            if(output_all_betas){
                C_half<-sqrtchoinv_rcpp(inv_C+diag(1e-15,nrow(inv_C)))
                pseudo_Xy_list<-pseudo_Xy(C_half=C_half,Z=Z,W=W,A=A,y=y,
                                          beta=beta,hat_thetaA=hat_thetaA,
                                          study_info=study_info,X=X,XR=XR)
                psX<-pseudo_Xy_list$pseudo_X/nZ

                psXtX<-self_crossprod_rcpp(psX)
                psXtX_non0<-psXtX[index_nonzero,index_nonzero,drop=F]
                inv_psXtX_non0<-choinv_rcpp(psXtX_non0)
                inv_psXtX_final<-inv_psXtX_non0

            }else{

                if(refine_C){
                    sC_half<-diag(1/sqrt(diag(inv_C)))
                    if(use_sparseC){
                        C_half<-sC_half
                    }else{
                        if(sqrt_matrix =="svd"){
                            inv_C_svd=fast.svd(inv_C+diag(1e-15,nrow(inv_C)))
                            C_half=prod_rcpp(inv_C_svd$v,(t(inv_C_svd$u)*1/sqrt(inv_C_svd$d)))
                        }else if(sqrt_matrix =="cholesky"){
                            C_half<-sqrtchoinv_rcpp(inv_C+diag(1e-15,nrow(inv_C)))
                        }
                    }
                    if(final.weight.min<0){
                        C_half<-sC_half
                    }
                    C_half[(pZ+pA+pW+1):nrow(C_half),(pZ+pA+pW+1):nrow(C_half)]=
                        C_half[(pZ+pA+pW+1):nrow(C_half),(pZ+pA+pW+1):nrow(C_half)]*sqrt(abs(final.weight.min))
                    use_sparseC<-T
                }

            ###########--------------###########
            # Compute new pseudo_X

            pseudo_Xy_list<-pseudo_Xy(C_half=C_half,Z=Z,W=W,A=A,y=y,
                                      beta=beta,hat_thetaA=hat_thetaA,
                                      study_info=study_info,X=X,XR=XR)
            psX<-pseudo_Xy_list$pseudo_X/nZ

            psXtX<-self_crossprod_rcpp(psX)
            psXtX_non0<-psXtX[index_nonzero,index_nonzero,drop=F]
            inv_psXtX_non0<-choinv_rcpp(psXtX_non0)
            inv_psXtX_final<-inv_psXtX_non0
            ###########--------------###########
            # When the C using is not optimal C,

            if(!is.null(fix_C)|final.weight.min!=1 |use_sparseC){
                inv_C_half<-sqrtcho_rcpp(inv_C+diag(1e-15,nrow(inv_C)))
                psX_mid<-prod_rcpp(inv_C_half,crossprod_rcpp(C_half,psX))
                psXtX_mid<-self_crossprod_rcpp(psX_mid)
                psXtX_mid_non0<-psXtX_mid[index_nonzero,index_nonzero,drop=F]
                inv_psXtX_final<-prod_rcpp(prod_rcpp(inv_psXtX_non0,psXtX_mid_non0),inv_psXtX_non0)

                if(sum(diag(inv_psXtX_non0)<=0)>0){
                    psXtX_mid_non0_half<-sqrtcho_rcpp(psXtX_mid_non0+diag(1e-15,nrow(psXtX_mid_non0)))
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
    return(return_list)
}
