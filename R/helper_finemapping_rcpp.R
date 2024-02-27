Delta_opt_finemapping<-function(y,Z,W,family,
                                study_info,A=NULL,pA=NULL,pZ=NULL,
                                beta=NULL,hat_thetaA=NULL,
                                V_thetaA=NULL,X=NULL,corZ=NULL){
    n_main=length(y)
    if(is.null(dim(X)[1])){X=cbind(A,Z,W)}
    X_beta = prodv_rcpp(X,beta)
    if(pA==0){ A_thetaA=0
    }else{A_thetaA = prodv_rcpp(A,hat_thetaA)}
    if(family == "binomial"){
        X_beta=expit_rcpp(X_beta)
        DZ<-sapply(1:ncol(Z), function(id){
            XR_thetaid<-expit_rcpp(c(A_thetaA+Z[,id]*study_info[[id]]$Coeff))
            X_beta-XR_thetaid})*Z #the same shape with Z
        DDZ<-sapply(1:ncol(Z), function(id){
            dexpit_rcpp(A_thetaA+Z[,id]*study_info[[id]]$Coeff)})*Z
        GammaZZ_vec =colMeans(timesm_rcpp(DDZ,Z))
        A_thetaA = expit_rcpp(A_thetaA)
        mu_prime_A_thetaA = A_thetaA*(1-A_thetaA)
    }else{
        DZ<-sapply(1:ncol(Z), function(id){
            XR_thetaid<-c(A_thetaA+Z[,id]*study_info[[id]]$Coeff)
            X_beta-XR_thetaid
        }) #the same shape with Z
        DDZ = Z
        GammaZZ_vec = colMeans(square_rcpp(Z))
        mu_prime_A_thetaA = 1
    }
    DX = X*c(X_beta-y)
    V_U1=(1/n_main)*self_crossprod_rcpp(DX)
    V_U2=(1/n_main)*self_crossprod_rcpp(DZ)
    Cov_U1U2=(1/n_main)*crossprod_rcpp(DX,DZ)

    V_thetaZ_vec<-sapply(1:ncol(Z), function(id){
        sqrt(study_info[[id]]$Covariance)})

    Delta22=V_U2+n_main*t((GammaZZ_vec*V_thetaZ_vec)*corZ)*(GammaZZ_vec*V_thetaZ_vec)
    #Delta22=V_U2+GammaZZ%*%(n_main*V_thetaZ)%*%t(GammaZZ)
    Delta12=Cov_U1U2
    if(pA!=0){
        GammaZA=(1/n_main)*crossprod_rcpp(DDZ,A)
        inv_GammaAA=ginv((1/n_main)*crossprod_rcpp(A*c(mu_prime_A_thetaA),A))
        Cov_U1theta=(1/n_main)*crossprod_rcpp(DX,A*c(A_thetaA-y))%*%(inv_GammaAA%*%t(GammaZA))
        Cov_U2theta=(1/n_main)*crossprod_rcpp(DZ,A*c(A_thetaA-y))%*%(inv_GammaAA%*%t(GammaZA))
        Delta22 = Delta22 + prod_rcpp(prod_rcpp(GammaZA,(n_main*V_thetaA)),t(GammaZA))
        + Cov_U2theta+t(Cov_U2theta)
        Delta12 = Delta12 + Cov_U1theta
    }
    Delta = rbind(cbind(V_U1,Delta12),cbind(t(Delta12),Delta22))
    Delta
}

pseudo_Xy_gaussian_finemapping<-function(
        C_half,Z,W,A,y,beta=NULL,hat_thetaA=NULL,study_info=NULL,X=NULL){
    if(is.null(dim(X)[1])){X=cbind(A,Z,W)}
    if(is.null(A)){A_thetaA = 0
    }else{A_thetaA = prodv_rcpp(A,hat_thetaA)}
    pseudo_X=prod_rcpp(C_half,crossprod_rcpp(cbind(X,Z),X))
    pseudo_y1=crossprodv_rcpp(X,y)
    pseudo_y2<-sapply(1:ncol(Z), function(id){
        c(c(A_thetaA+Z[,id]*study_info[[id]]$Coeff)%*%Z[,id])
    })
    pseudo_y=prodv_rcpp(C_half,c(pseudo_y1,pseudo_y2))
    list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
}
pseudo_Xy_binomial_finemapping<-function(
        C_half,Z,W,A,y,beta=NULL,hat_thetaA=NULL,study_info=NULL,X=NULL){
    if(is.null(dim(X)[1])){X=cbind(A,Z,W)}
    X_beta=prodv_rcpp(X,beta)
    if(is.null(A)){A_thetaA = 0
    }else{A_thetaA = prodv_rcpp(A,hat_thetaA)}
    expit_beta=expit_rcpp(X_beta)
    dexpit_beta=expit_beta*(1-expit_beta)
    #pseudo_X=C_half%*%Rfast::Crossprod(cbind(X,Z),X*dexpit_beta)
    pseudo_X=prod_rcpp(C_half,crossprod_rcpp(cbind(X,Z),X*dexpit_beta))

    ps_y0=crossprodv_rcpp(cbind(X,Z),c(dexpit_beta*X_beta-expit_beta))
    ps_y1=crossprodv_rcpp(X,y)
    ps_y2<-sapply(1:ncol(Z), function(id){
        c(expit_rcpp(c(A_thetaA+Z[,id]*study_info[[id]]$Coeff))%*%Z[,id])
    })
    pseudo_y=prodv_rcpp(C_half,c(ps_y0+c(ps_y1,ps_y2)))
    list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
}




fm.htlgmm.default<-function(
        y,Z,W=NULL,
        study_info=NULL,
        A=1,
        penalty_type = "adaptivelasso",
        family = "gaussian",
        initial_with_type = "ridge",
        beta_initial = NULL,
        hat_thetaA = NULL,
        V_thetaA = NULL,
        remove_penalty_Z = FALSE,
        remove_penalty_W = FALSE,
        inference = TRUE,
        refine_C = TRUE,
        sqrt_matrix ="cholesky",
        use_cv = TRUE,
        type_measure = "default",
        nfolds = 10,
        fix_lambda = NULL,
        lambda_list = NULL,
        nlambda = 100,
        lambda.min.ratio = 0.0001,
        tune_ratio = FALSE,
        fix_ratio = NULL,
        ratio_list = NULL,
        gamma_adaptivelasso = 1/2,
        use_sparseC = TRUE,
        seed.use = 97,
        corZ = NULL
){
    set.seed(seed.use)
    if (is.null(study_info)){stop("Please input study_info as trained model")}
    if(!penalty_type %in% c("none","adaptivelasso","lasso","ridge")){
        stop("Select penalty type from c('none','adaptivelasso','lasso','ridge').")
    }
    if(!type_measure%in% c("default", "mse", "deviance", "auc")){
        stop("Select type_measure from c('default','deviance','auc')")
    }
    if(!sqrt_matrix %in% c("cholesky","svd")){
        stop("Select penalty type from c('cholesky','svd').")
    }
    if(is.null(dim(Z)[1])){
        warning("Z is input as a vector, convert Z into matrix with size nZ*1")
        Z=matrix(Z,ncol=1)
    }
    final_alpha = 1
    if(penalty_type == "ridge"){final_alpha = 0}

    if(!is.null(fix_ratio)){
        if(tune_ratio){
            stop("If ratio is fixed, please set tune_ratio as FALSE")
        }else if(remove_penalty_Z | remove_penalty_W){
            stop("If ratio is fixed, please set remove_penalty's as FALSE")
        }
    }

    nZ=nrow(Z)
    nZext=study_info[[1]]$Sample_size
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

    if(family == "gaussian"){pseudo_Xy=pseudo_Xy_gaussian_finemapping
    }else if(family == "binomial"){pseudo_Xy=pseudo_Xy_binomial_finemapping}


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
    X<-cbind(A,Z,W)

    if(pA!=0){
        if(is.null(hat_thetaA)){
            if(!is.null(V_thetaA)){
                stop("With customized hat_thetaA input, V_thetaA is also needed")
            }
            df=data.frame(y,A)
            if(family=="binomial"){
                hat_thetaA_glm=speedglm(y~0+.,data = df,family = binomial())
            }else if(family=="gaussian"){
                hat_thetaA_glm=speedlm(y~0+.,data = df)
            }
            hat_thetaA=hat_thetaA_glm$coefficients
            V_thetaA=vcov(hat_thetaA_glm)
        }
    }


    fix_penalty<-c(rep(0,pA),rep(1,pZ+pW))
    if(remove_penalty_Z){fix_penalty[Zid]<-0}else{
        if(!is.null(fix_ratio)){fix_penalty[Zid]<-fix_ratio}
    }
    if(remove_penalty_W){fix_penalty[Wid]<-0}
    if((remove_penalty_Z & remove_penalty_W)|(length(unique(fix_penalty))==1 & unique(fix_penalty)[1] == 0) ){
        penalty_type = "none"
        warning("All penalties are removed, turn to no penalties!")
    }
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
            }else if(pA == 0){
                #res <- susie(X,y,L=10)
                #beta_initial = coef(res)[-1]
                #w_adaptive<-1/abs(beta_initial)^(1/2)
                #w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*100
                #fit_initial=cv.glmnet(x= X,y= y,alpha = initial_alpha,
                #            penalty.factor = w_adaptive,family=family)

                fit_initial=cv.glmnet(x=X,y= y,
                                      alpha = initial_alpha,
                                      penalty.factor = fix_penalty,
                                      family=family)
                beta_initial=c(coef.glmnet(fit_initial,s="lambda.min")[-1])
            }else if(length(unique(A[,1]))==1){
                if(unique(A[,1])==1){
                    fit_initial=cv.glmnet(x=X[,-1],y=y,
                                          alpha = initial_alpha,
                                          penalty.factor = fix_penalty[-1],
                                          family=family)
                    beta_initial=as.vector(coef.glmnet(fit_initial,s="lambda.min"))
                }else{
                    stop("The first column of A is constant, then it should be 1 for intercept.")
                }
            }else{
                fit_initial=cv.glmnet(x=X,y=y,
                                      alpha = initial_alpha,
                                      penalty.factor = fix_penalty,
                                      family=family)
                beta_initial=c(coef.glmnet(fit_initial,s="lambda.min")[-1])
            }
        }else{stop("Select Initial Type from c('glm','ridge','lasso')")}
    }
    if (penalty_type == "adaptivelasso"){
        w_adaptive<-1/abs(beta_initial)^gamma_adaptivelasso
        w_adaptive[is.infinite(w_adaptive)]<-max(w_adaptive[!is.infinite(w_adaptive)])*100
        w_adaptive<-w_adaptive*fix_penalty
    }else{w_adaptive<-fix_penalty}
    # Estimation of C
    if(is.null(corZ)){corZ = cor(Z)}
    inv_C = Delta_opt_finemapping(y=y,Z=Z,W=W,
                                  family=family,
                                  study_info=study_info,
                                  A=A,pA=pA,pZ=pZ,beta=beta_initial,
                                  hat_thetaA=hat_thetaA,
                                  V_thetaA = V_thetaA,
                                  X=X,
                                  corZ=corZ)
    if(use_sparseC){
        C_half<-diag(1/sqrt(diag(inv_C)))
    }else{
        if(sqrt_matrix =="svd"){
            inv_C_svd=fast.svd(inv_C+diag(1e-15,nrow(inv_C)))
            C_half=prod_rcpp(inv_C_svd$v,(t(inv_C_svd$u)*1/sqrt(inv_C_svd$d)))
            #C_half<-inv_C_svd$v%*%diag(1/sqrt(inv_C_svd$d))%*%t(inv_C_svd$u)
        }else if(sqrt_matrix =="cholesky"){
            C_half<-sqrtchoinv_rcpp(inv_C+diag(1e-15,nrow(inv_C)))
        }
    }

    # Prepare for final model
    pseudo_Xy_list<-pseudo_Xy(C_half=C_half,Z=Z,W=W,A=A,y=y,
                              beta=beta_initial,hat_thetaA=hat_thetaA,
                              study_info=study_info,X=X)
    initial_sf<-nZ/sqrt(nrow(pseudo_Xy_list$pseudo_X))
    pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
    pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf

    # generate lambda list from glmnet
    if(penalty_type != "none"){
        if(is.null(fix_lambda)&is.null(lambda_list)){
            fit_final<-glmnet(x= pseudo_X,y= pseudo_y,standardize=F,
                              intercept=F,alpha = final_alpha,penalty.factor = fix_penalty)
            #lambda_list<-fit_final$lambda
            lambda.max<-fit_final$lambda[1]
            lambda_list <-exp(seq(log(lambda.max),log(lambda.max*lambda.min.ratio),length.out=nlambda))
        }
    }
    if(!is.null(fix_lambda)){
        use_cv = FALSE
        if(fix_lambda<0){stop("The fixed lambda should be nonnegative.")}
    }

    if(tune_ratio & !remove_penalty_Z & !remove_penalty_W){
        if(is.null(ratio_list)){
            ratio_lower<-sqrt(nZ/(nZ+nZext))/2
            ratio_upper<-(nZ)^(1/3)/2
            ratio_count<-10
            ratio_list<-(seq(sqrt(ratio_lower),sqrt(ratio_upper),(sqrt(ratio_upper)-sqrt(ratio_lower))/ratio_count)^2)
            ratio_list<-c(1,ratio_list)
        }
    }else{tune_ratio<-FALSE}
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
            if(!is.null(fix_ratio)){
                return_list<-c(return_list,
                               list("fix_ratio"=fix_ratio))
            }
        }
    }else{
        if(length(unique(y)) <= 2){
            set.seed(seed.use)
            index_fold<-createFolds(as.numeric(y>0),k = nfolds)
        }else{
            set.seed(seed.use)
            index_fold<-createFolds(y,k = nfolds)
        }

        if(tune_ratio){
            if(family == "gaussian"){
                cv_mse<-cv_mse_lambda_ratio_func(index_fold,Z,W,A,y,
                                                 C_half,beta_initial,hat_thetaA,
                                                 study_info,lambda_list,
                                                 ratio_list,pZ,pW,pA,
                                                 w_adaptive,final_alpha,pseudo_Xy)
                ids<-which(cv_mse==min(cv_mse),arr.ind = TRUE)
            }else if(family == "binomial"){
                cv_dev<-cv_dev_lambda_ratio_func(index_fold,Z,W,A,y,
                                                 C_half,beta_initial,hat_thetaA,
                                                 study_info,lambda_list,
                                                 ratio_list,pZ,pW,pA,
                                                 w_adaptive,final_alpha,pseudo_Xy)
                if(type_measure == "auc"){
                    cv_dev1<-cv_dev$auc
                    ids<-which(cv_dev1==max(cv_dev1),arr.ind = TRUE)
                }else{
                    cv_dev1<-cv_dev$deviance
                    ids<-which(cv_dev1==min(cv_dev1),arr.ind = TRUE)
                }
            }
            final.ratio.min<-ratio_list[ids[1]]
            final.lambda.min<-lambda_list[ids[2]]
        }else{
            if(family == "gaussian"){
                cv_mse<-cv_mse_lambda_func(index_fold,Z,W,A,y,
                                           C_half,beta_initial,hat_thetaA,
                                           study_info,lambda_list,
                                           w_adaptive,final_alpha,pseudo_Xy)
                final.lambda.min<-lambda_list[which.min(cv_mse)]
            }else if(family == "binomial"){
                cv_dev<-cv_dev_lambda_func(index_fold,Z,W,A,y,
                                           C_half,beta_initial,hat_thetaA,
                                           study_info,lambda_list,
                                           w_adaptive,final_alpha,pseudo_Xy)
                if(type_measure == "auc"){
                    cv_dev1<-cv_dev$auc
                    final.lambda.min<-lambda_list[which.max(cv_dev1)]
                }else{
                    cv_dev1<-cv_dev$deviance
                    final.lambda.min<-lambda_list[which.min(cv_dev1)]
                }
            }
            final.ratio.min<-1
        }
        ratio_vec<-c(rep(final.ratio.min,pZ),rep(1,pW+pA))
        w_adaptive_ratio<-w_adaptive*ratio_vec
        fit_final_lam_ratio<-glmnet(x= pseudo_X,y= pseudo_y,standardize=F,
                                    intercept=F,alpha = final_alpha,
                                    penalty.factor = w_adaptive_ratio,
                                    lambda = final.lambda.min)

        beta<-coef.glmnet(fit_final_lam_ratio)[-1]
        return_list<-list("beta"=beta,
                          "lambda_list"=lambda_list,
                          "ratio_list"=ratio_list,
                          "lambda_min"=final.lambda.min,
                          "ratio_min"=final.ratio.min)
        if(family == "gaussian"){
            return_list<-c(return_list,
                           list("cv_mse"=cv_mse/nfolds))
        }else if(family == "binomial"){
            return_list<-c(return_list,
                           list("cv_dev"=cv_dev$deviance/nfolds,
                                "cv_auc"=cv_dev$auc/nfolds))
        }
    }

    index_nonzero<-which(beta!=0)
    # remove intercept term related
    #if(inference){
    #    if(pA>0 & penalty_type != "none"){
    #        if(Acolnames[1]=='intercept' & index_nonzero[1] == 1){
    #            index_nonzero = index_nonzero[-1]
    #        }
    #    }
    #}
    if(length(index_nonzero) > 0){
        if(penalty_type == "lasso"){
            warning("Current penalty is lasso, please turn to adaptivelasso for inference")
        }
        # refine C
        if(refine_C){
            inv_C = Delta_opt_finemapping(y=y,Z=Z,W=W,
                                          family=family,
                                          study_info=study_info,
                                          A=A,pA=pA,pZ=pZ,beta=beta,
                                          hat_thetaA=hat_thetaA,
                                          V_thetaA = V_thetaA,
                                          X=X,
                                          corZ=corZ)
            if(sqrt_matrix =="svd"){
                inv_C_svd=fast.svd(inv_C+diag(1e-15,nrow(inv_C)))
                C_half=prod_rcpp(inv_C_svd$v,(t(inv_C_svd$u)*1/sqrt(inv_C_svd$d)))
                #C_half<-inv_C_svd$v%*%diag(1/sqrt(inv_C_svd$d))%*%t(inv_C_svd$u)
            }else if(sqrt_matrix =="cholesky"){
                C_half<-sqrtchoinv_rcpp(inv_C+diag(1e-15,nrow(inv_C)))
            }
        }

        pseudo_Xy_list<-pseudo_Xy(C_half=C_half,Z=Z,W=W,A=A,y=y,
                                  beta=beta,hat_thetaA=hat_thetaA,
                                  study_info=study_info,X=X)
        Sigsum_half<-pseudo_Xy_list$pseudo_X/nZ

        Sigsum_scaled<-self_crossprod_rcpp(Sigsum_half)
        Sigsum_scaled_nonzero<-Sigsum_scaled[index_nonzero,index_nonzero,drop=F]
        inv_Sigsum_scaled_nonzero<-choinv_rcpp(Sigsum_scaled_nonzero)
        final_v<-diag(inv_Sigsum_scaled_nonzero)/nZ

        pval_final<-pchisq(beta[index_nonzero]^2/final_v,1,lower.tail = F)
        pval_final1<-p.adjust(pval_final,method = "BH")
        selected_pos<-index_nonzero[which(pval_final1<0.05)]
        return_list<-c(return_list,
                       list("selected_vars"=
                                list("position"=index_nonzero,
                                     "name"=Xcolnames[index_nonzero],
                                     "coef"=beta[index_nonzero],
                                     "variance"=final_v,
                                     "pval"=pval_final,
                                     "FDR_adjust_position"=selected_pos,
                                     "FDR_adjust_name"=Xcolnames[selected_pos])
                       ))
        if(pA>0 & penalty_type != "none"){
            if(Acolnames[1]=='intercept' & index_nonzero[1] == 1){
                index_nonzero = index_nonzero[-1]
                Sigsum_scaled_nonzero<-Sigsum_scaled[index_nonzero,index_nonzero,drop=F]
                inv_Sigsum_scaled_nonzero<-choinv_rcpp(Sigsum_scaled_nonzero)
                final_v<-diag(inv_Sigsum_scaled_nonzero)/nZ
                pval_final<-pchisq(beta[index_nonzero]^2/final_v,1,lower.tail = F)
                pval_final1<-p.adjust(pval_final,method = "BH")
                selected_pos<-index_nonzero[which(pval_final1<0.05)]
                return_list<-c(return_list,
                               list("selected_vars_nointercept"=
                                        list("position"=index_nonzero,
                                             "name"=Xcolnames[index_nonzero],
                                             "coef"=beta[index_nonzero],
                                             "variance"=final_v,
                                             "pval"=pval_final,
                                             "FDR_adjust_position"=selected_pos,
                                             "FDR_adjust_name"=Xcolnames[selected_pos])
                               ))
            }}
    }
    return(return_list)
}


group.fm.htlgmm.default<-function(
        y,Z,W=NULL,
        study_info=NULL,
        A=1,
        penalty_type = "adaptivelasso",
        family = "gaussian",
        decor_method = "pca",
        max_cor = 0.9,
        min_cor = 0.5,
        ncor = 5,
        initial_with_type = "ridge",
        hat_thetaA = NULL,
        V_thetaA = NULL,
        remove_penalty_Z = FALSE,
        remove_penalty_W = FALSE,
        inference = TRUE,
        refine_C = TRUE,
        sqrt_matrix ="cholesky",
        use_cv = TRUE,
        type_measure = "default",
        nfolds = 10,
        fix_lambda = NULL,
        lambda_list = NULL,
        nlambda = 20,
        lambda.min.ratio = 0.0001,
        tune_ratio = FALSE,
        fix_ratio = NULL,
        ratio_list = NULL,
        gamma_adaptivelasso = 1/2,
        use_sparseC = TRUE,
        seed.use = 97){
    if(!decor_method%in%c("pca","top","pip","pippca")){
        stop("Select decor_method from c('pca','top','pip','pippca')")
    }
    set.seed(seed.use)
    if (is.null(study_info)){stop("Please input study_info as trained model")}
    if(!penalty_type %in% c("none","adaptivelasso","lasso","ridge")){
        stop("Select penalty type from c('none','adaptivelasso','lasso','ridge').")
    }
    if(!type_measure%in% c("default", "mse", "deviance", "auc")){
        stop("Select type_measure from c('default','deviance','auc')")
    }
    if(!sqrt_matrix %in% c("cholesky","svd")){
        stop("Select penalty type from c('cholesky','svd').")
    }
    if(is.null(dim(Z)[1])){
        warning("Z is input as a vector, convert Z into matrix with size nZ*1")
        Z=matrix(Z,ncol=1)
    }
    final_alpha = 1
    if(penalty_type == "ridge"){final_alpha = 0}

    if(!is.null(fix_ratio)){
        if(tune_ratio){
            stop("If ratio is fixed, please set tune_ratio as FALSE")
        }else if(remove_penalty_Z | remove_penalty_W){
            stop("If ratio is fixed, please set remove_penalty's as FALSE")
        }
    }

    corZ = cor(Z)
    abscorZ = abs(corZ)
    if(max_cor<=min_cor){
        stop("We need max_cor > min_cor")
    }
    max_cor=min(max_cor,max(abscorZ[row(corZ)!=col(corZ)]))
    min_cor=max(min_cor,min(abscorZ[row(corZ)!=col(corZ)]))
    if(max_cor<=min_cor){
        stop("The max cor in Z is smaller than preset min_cor or the min cor in Z is larger than preset max_cor")
    }
    cor_seq = c(seq(min_cor,max_cor,length.out=ncor))
    if(decor_method %in%c("pip","pippca") ){
        pZ = nrow(corZ)
        qval_cutoff = qchisq(0.05/sqrt(pZ),1,lower.tail = F)
        #diag(abscorZ) = 0
        s<-sapply(1:length(study_info), function(i){
            c(study_info[[i]]$Coeff^2/study_info[[i]]$Covariance,i,T)
        })
        s = s[,order(s[1,],decreasing = T)]
        s = rbind(s,1:ncol(s))
if(0){
        res_clu_bycor=lapply(cor_seq, function(cur_cor){
            #message(cur_cor)
            clusters_list=list()
            cur_id = 1
            cur_clu = 1
            remain_ct = pZ
            while(remain_ct>0 & cur_id < pZ){
                if(qval_cutoff>s[1,cur_id]){break}
                cur_cutoff = sqrt(1-qval_cutoff/s[1,cur_id])
                if(cur_cutoff<cur_cor){break}
                stmp = s[,s[3,]==1]
                cur_ids=which(abscorZ[s[2,cur_id],stmp[2,]]>cur_cutoff)
                clusters_list[[cur_clu]] = stmp[2,cur_ids]
                s[3,stmp[4,cur_ids]]=0
                cur_clu = cur_clu+1
                remain_ct=remain_ct-length(cur_ids)
                while (s[3,cur_id] == 0 & cur_id < pZ) {
                    cur_id = cur_id+1
                }
            }
            if(remain_ct>0 & cur_id < pZ){
                cur_clu = cur_clu-1
                s=s[,s[3,]==1]
                for(i in 1:ncol(s)){
                    clusters_list[[cur_clu+i]]=s[2,i]
                }
            }
            Z_clu<-sapply(1:length(clusters_list), function(j){
                cur_clu<-clusters_list[[j]]
                if(length(cur_clu)>1){
                    cur_beta=sapply(cur_clu, function(i){study_info[[i]]$Coeff})
                    cur_var=sapply(cur_clu, function(i){sqrt(study_info[[i]]$Covariance)})
                    cur_size=sapply(cur_clu, function(i){study_info[[i]]$Sample_size})
                    if(decor_method == "pip"){
                        maxid = which.max(cur_beta^2/cur_var)
                        Coeff=cur_beta[maxid]
                        Covariance=cur_var[maxid]
                        Sample_size=cur_size[maxid]
                        zpc=c(Coeff,Covariance,Sample_size,Z[,cur_clu[maxid]])
                    }else{
                        zpca=prcomp(Z[,cur_clu], rank. = 1)
                        rotate_=c(zpca$rotation)
                        Coeff=rotate_%*%cur_beta
                        Covariance=c(rotate_%*%prodv_rcpp(t(cur_var*corZ[cur_clu,cur_clu])*cur_var,rotate_))
                        #Covariance=c(rotate_%*%(t(cur_var*corZ[cur_clu,cur_clu])*cur_var)%*%rotate_)
                        Sample_size=rotate_%*%cur_size
                        zpc=c(Coeff,Covariance,Sample_size,zpca$x[,1])
                    }
                }else{
                    Coeff=study_info[[cur_clu]]$Coeff
                    Covariance=study_info[[cur_clu]]$Covariance
                    Sample_size=study_info[[cur_clu]]$Sample_size
                    zpc=c(Coeff,Covariance,Sample_size,Z[,cur_clu])
                }
                zpc
            })
            study_info_clu<-lapply(1:length(clusters_list), function(i){
                list("Coeff"=Z_clu[1,i],
                     "Covariance"=Z_clu[2,i],
                     "Sample_size"=Z_clu[3,i]
                )
            })
            Z_clu = Z_clu[-c(1,2,3),]

            beta_initial = NULL
            res_clu<-fm.htlgmm.default(y,Z_clu,W,study_info_clu,A,penalty_type,
                                       family,initial_with_type,beta_initial,
                                       hat_thetaA,V_thetaA,remove_penalty_Z,
                                       remove_penalty_W,inference,refine_C,
                                       sqrt_matrix,use_cv,type_measure,nfolds,
                                       fix_lambda,lambda_list,nlambda,lambda.min.ratio,
                                       tune_ratio,fix_ratio,ratio_list,
                                       gamma_adaptivelasso,
                                       use_sparseC,seed.use)
            c(res_clu,list("clusters_list"=clusters_list,"cor_cutoff"=cur_cor))
        })
}

        cur_cor = 0.9
        qval_cutoff_min = qchisq(0.05/sqrt(pZ),1,lower.tail = F)
        qval_cutoff_max = qchisq(0.05/pZ,1,lower.tail = F)
        qval_cutoff_seq = seq(qval_cutoff_min,qval_cutoff_max,length.out=ncor)
        res_clu_bycor=lapply(qval_cutoff_seq, function(qval_cutoff){
            #message(cur_cor)
            clusters_list=list()
            cur_id = 1
            cur_clu = 1
            remain_ct = pZ
            while(remain_ct>0 & cur_id < pZ){
                if(qval_cutoff>s[1,cur_id]){break}
                cur_cutoff = sqrt(1-qval_cutoff/s[1,cur_id])
                if(cur_cutoff<cur_cor){break}
                stmp = s[,s[3,]==1]
                cur_ids=which(abscorZ[s[2,cur_id],stmp[2,]]>cur_cutoff)
                clusters_list[[cur_clu]] = stmp[2,cur_ids]
                s[3,stmp[4,cur_ids]]=0
                cur_clu = cur_clu+1
                remain_ct=remain_ct-length(cur_ids)
                while (s[3,cur_id] == 0 & cur_id < pZ) {
                    cur_id = cur_id+1
                }
            }
            if(remain_ct>0 & cur_id < pZ){
                cur_clu = cur_clu-1
                s=s[,s[3,]==1]
                for(i in 1:ncol(s)){
                    clusters_list[[cur_clu+i]]=s[2,i]
                }
            }
            Z_clu<-sapply(1:length(clusters_list), function(j){
                cur_clu<-clusters_list[[j]]
                if(length(cur_clu)>1){
                    cur_beta=sapply(cur_clu, function(i){study_info[[i]]$Coeff})
                    cur_var=sapply(cur_clu, function(i){sqrt(study_info[[i]]$Covariance)})
                    cur_size=sapply(cur_clu, function(i){study_info[[i]]$Sample_size})
                    if(decor_method == "pip"){
                        maxid = which.max(cur_beta^2/cur_var)
                        Coeff=cur_beta[maxid]
                        Covariance=cur_var[maxid]
                        Sample_size=cur_size[maxid]
                        zpc=c(Coeff,Covariance,Sample_size,Z[,cur_clu[maxid]])
                    }else{
                        zpca=prcomp(Z[,cur_clu], rank. = 1)
                        rotate_=c(zpca$rotation)
                        Coeff=rotate_%*%cur_beta
                        Covariance=c(rotate_%*%prodv_rcpp(t(cur_var*corZ[cur_clu,cur_clu])*cur_var,rotate_))
                        #Covariance=c(rotate_%*%(t(cur_var*corZ[cur_clu,cur_clu])*cur_var)%*%rotate_)
                        Sample_size=rotate_%*%cur_size
                        zpc=c(Coeff,Covariance,Sample_size,zpca$x[,1])
                    }
                }else{
                    Coeff=study_info[[cur_clu]]$Coeff
                    Covariance=study_info[[cur_clu]]$Covariance
                    Sample_size=study_info[[cur_clu]]$Sample_size
                    zpc=c(Coeff,Covariance,Sample_size,Z[,cur_clu])
                }
                zpc
            })
            study_info_clu<-lapply(1:length(clusters_list), function(i){
                list("Coeff"=Z_clu[1,i],
                     "Covariance"=Z_clu[2,i],
                     "Sample_size"=Z_clu[3,i]
                )
            })
            Z_clu = Z_clu[-c(1,2,3),]

            beta_initial = NULL
            res_clu<-fm.htlgmm.default(y,Z_clu,W,study_info_clu,A,penalty_type,
                                       family,initial_with_type,beta_initial,
                                       hat_thetaA,V_thetaA,remove_penalty_Z,
                                       remove_penalty_W,inference,refine_C,
                                       sqrt_matrix,use_cv,type_measure,nfolds,
                                       fix_lambda,lambda_list,nlambda,lambda.min.ratio,
                                       tune_ratio,fix_ratio,ratio_list,
                                       gamma_adaptivelasso,
                                       use_sparseC,seed.use)
            c(res_clu,list("clusters_list"=clusters_list,"cor_cutoff"=cur_cor))
        })



        if(family == "gaussian"){
            best_cor_id = which.min(sapply(1:ncor, function(i){min(res_clu_bycor[[i]]$cv_mse)}))
        }else{
            if(type_measure%in%c("default","dev")){
                best_cor_id = which.min(sapply(1:ncor, function(i){min(res_clu_bycor[[i]]$cv_dev)}))
            }else{
                best_cor_id = which.max(sapply(1:ncor, function(i){max(res_clu_bycor[[i]]$cv_auc)}))
            }
        }
        output=res_clu_bycor[[best_cor_id]]
    }else{
        dissimilarity=as.dist(1-abscorZ)
        hc=hclust(dissimilarity, method = "complete")
        res_clu_bycor=lapply(cor_seq, function(cur_cor){
            #message(cur_cor)
            height_cutoff <- 1 - cur_cor
            clusters <- cutree(hc, h = height_cutoff)
            clusters_list <- split(1:length(clusters), clusters)
            if(decor_method == "pca"){
                Z_clu<-sapply(1:length(clusters_list), function(j){
                    cur_clu<-clusters_list[[j]]
                    if(length(cur_clu)>1){
                        zpca=prcomp(Z[,cur_clu], rank. = 1)
                        rotate_=c(zpca$rotation)
                        cur_beta=sapply(cur_clu, function(i){study_info[[i]]$Coeff})
                        cur_var=sapply(cur_clu, function(i){sqrt(study_info[[i]]$Covariance)})
                        cur_size=sapply(cur_clu, function(i){study_info[[i]]$Sample_size})
                        Coeff=rotate_%*%cur_beta
                        Covariance=c(rotate_%*%prodv_rcpp(t(cur_var*corZ[cur_clu,cur_clu])*cur_var,rotate_))
                        #Covariance=c(rotate_%*%(t(cur_var*corZ[cur_clu,cur_clu])*cur_var)%*%rotate_)
                        Sample_size=rotate_%*%cur_size
                        zpc=c(Coeff,Covariance,Sample_size,zpca$x[,1])
                    }else{
                        Coeff=study_info[[cur_clu]]$Coeff
                        Covariance=study_info[[cur_clu]]$Covariance
                        Sample_size=study_info[[cur_clu]]$Sample_size
                        zpc=c(Coeff,Covariance,Sample_size,Z[,cur_clu])
                    }
                    zpc
                })
            }else{
                Z_clu<-sapply(1:length(clusters_list), function(j){
                    cur_clu<-clusters_list[[j]]
                    if(length(cur_clu)>1){
                        cur_beta=sapply(cur_clu, function(i){study_info[[i]]$Coeff})
                        cur_var=sapply(cur_clu, function(i){sqrt(study_info[[i]]$Covariance)})
                        cur_size=sapply(cur_clu, function(i){study_info[[i]]$Sample_size})
                        maxid = which.max(abs(cur_beta^2/cur_var))
                        Coeff=cur_beta[maxid]
                        Covariance=cur_var[maxid]
                        Sample_size=cur_size[maxid]
                        zpc=c(Coeff,Covariance,Sample_size,Z[,cur_clu[maxid]])
                    }else{
                        Coeff=study_info[[cur_clu]]$Coeff
                        Covariance=study_info[[cur_clu]]$Covariance
                        Sample_size=study_info[[cur_clu]]$Sample_size
                        zpc=c(Coeff,Covariance,Sample_size,Z[,cur_clu])
                    }
                    zpc
                })
            }
            study_info_clu<-lapply(1:length(clusters_list), function(i){
                list("Coeff"=Z_clu[1,i],
                     "Covariance"=Z_clu[2,i],
                     "Sample_size"=Z_clu[3,i]
                )
            })
            Z_clu = Z_clu[-c(1,2,3),]

            beta_initial = NULL
            res_clu<-fm.htlgmm.default(y,Z_clu,W,study_info_clu,A,penalty_type,
                                       family,initial_with_type,beta_initial,
                                       hat_thetaA,V_thetaA,remove_penalty_Z,
                                       remove_penalty_W,inference,refine_C,
                                       sqrt_matrix,use_cv,type_measure,nfolds,
                                       fix_lambda,lambda_list,nlambda,lambda.min.ratio,
                                       tune_ratio,fix_ratio,ratio_list,
                                       gamma_adaptivelasso,
                                       use_sparseC,seed.use)
            c(res_clu,list("clusters_list"=clusters_list,"cor_cutoff"=cur_cor))
        })
        if(family == "gaussian"){
            best_cor_id = which.min(sapply(1:ncor, function(i){min(res_clu_bycor[[i]]$cv_mse)}))
        }else{
            if(type_measure%in%c("default","dev")){
                best_cor_id = which.min(sapply(1:ncor, function(i){min(res_clu_bycor[[i]]$cv_dev)}))
            }else{
                best_cor_id = which.max(sapply(1:ncor, function(i){max(res_clu_bycor[[i]]$cv_auc)}))
            }
        }
        output=res_clu_bycor[[best_cor_id]]
    }
    output
}


