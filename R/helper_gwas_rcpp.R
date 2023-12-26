Delta_opt_fast1<-function(y,Z,W,A,family,pA,
                         X_beta,A_thetaA,Z_thetaZ,V_thetaA,V_thetaZ,inv_GammaAA){
    n_main=length(y)
    X=cbind(A,W,Z)
    XR_theta=add_rcpp(A_thetaA,Z_thetaZ)
    if(family == "binomial"){
        X_beta=expit_rcpp(X_beta)
        XR_theta=expit_rcpp(XR_theta)
        mu_prime_XR_theta=timesv_rcpp(XR_theta,(1-XR_theta))
        A_thetaA=expit_rcpp(A_thetaA)
    }else{mu_prime_XR_theta=1}
    DX = X*c(X_beta-y)
    DZ = Z*c(X_beta-XR_theta)
    DZ2 = Z*c(mu_prime_XR_theta)
    V_U1=(1/n_main)*self_crossprod_rcpp(DX)
    V_U2=(1/n_main)*self_crossprod_rcpp(DZ)
    Cov_U1U2=(1/n_main)*crossprod_rcpp(DX,DZ)
    GammaZZ=(1/n_main)*crossprod_rcpp(DZ2,Z)
    Delta22=add_rcpp(V_U2,prod_rcpp(prod_rcpp(GammaZZ,(n_main*V_thetaZ)),t(GammaZZ)))
    Delta12=Cov_U1U2
    if(pA!=0){
        GammaAZ=(1/n_main)*crossprod_rcpp(A,DZ2)
        #inv_GammaAA=ginv((1/n_main)*crossprod(A*c(mu_prime_A_thetaA),A))
        Gammas=prod_rcpp(inv_GammaAA,GammaAZ)
        Cov_U1theta=(1/n_main)*crossprod_rcpp(DX,prod_rcpp(A*c(A_thetaA-y),Gammas))
        Cov_U2theta=(1/n_main)*crossprod_rcpp(DZ,prod_rcpp(A*c(A_thetaA-y),Gammas))
        Delta22 = add_rcpp(Delta22,prod_rcpp(crossprod_rcpp(GammaAZ,n_main*V_thetaA),GammaAZ))
        Delta22 = add_rcpp(Delta22,add_rcpp(Cov_U2theta,t(Cov_U2theta)))
        Delta12 = add_rcpp(Delta12, Cov_U1theta)
    }
    Delta = rbind(cbind(V_U1,Delta12),cbind(t(Delta12),Delta22))
    Delta
}


direct_fast1<-function(
        inv_C,y,Z,W,A,XR_theta,
        X_beta,family){
    Cn = spdinv(inv_C+diag(1e-15,nrow(inv_C)))
    X=cbind(A,W,Z)
    scale_factor = length(y)/sqrt(nrow(inv_C))
    if(family == "binomial"){
        XR_theta=expit_rcpp(c(XR_theta))
        expit_beta=expit_rcpp(c(X_beta))
        dexpit_beta=timesv_rcpp(expit_beta,(1-expit_beta))
        ps_y0=c(crossprod_rcpp(cbind(X,Z),c(dexpit_beta*X_beta-expit_beta)))
    }else{
        dexpit_beta=1
        ps_y0=0
    }
    pseudo_X0=crossprod_rcpp(cbind(X,Z),X*dexpit_beta)/scale_factor
    ps_XtX = prod_rcpp(crossprod_rcpp(pseudo_X0,Cn),pseudo_X0)
    ps_y1=crossprodv_rcpp(X,y)
    ps_y2=crossprodv_rcpp(Z,XR_theta)
    ps_y=add_rcpp(ps_y0,c(ps_y1,ps_y2))/scale_factor
    ps_Xty = prodv_rcpp(crossprod_rcpp(pseudo_X0,Cn),ps_y)
    beta = prodv_rcpp(spdinv(ps_XtX),ps_Xty)
    beta
}
direct_fast12<-function(inv_C,y,Z,W,A,XR_theta,
                       X_beta,family){
    Cn = spdinv(inv_C+diag(1e-15,nrow(inv_C)))
    X=cbind(A,W,Z)
    scale_factor = length(y)
    if(family == "binomial"){
        XR_theta=expit_rcpp(c(XR_theta))
        expit_beta=expit_rcpp(c(X_beta))
        dexpit_beta=timesv_rcpp(expit_beta,(1-expit_beta))
        ps_y0=c(crossprod_rcpp(cbind(X,Z),c(dexpit_beta*X_beta-expit_beta)))
    }else{
        dexpit_beta=1
        ps_y0=0
    }
    pseudo_X0=crossprod_rcpp(cbind(X,Z),X*dexpit_beta)/scale_factor
    ps_XtX = prod_rcpp(crossprod_rcpp(pseudo_X0,Cn),pseudo_X0)
    inv_ps_XtX=spdinv(ps_XtX)
    diag(inv_ps_XtX)/scale_factor
}

htlgmm.gwas.default3<-function(
        y,Z,W=NULL,
        study_info=NULL,
        A=NULL,
        family = "gaussian",
        AW_betaAW=NULL,
        A_thetaA=NULL,
        V_thetaA=NULL,
        inv_GammaAA=NULL,
        output_SNP_only=TRUE,
        seed.use = 97,
        verbose = FALSE,
        output_tmp = FALSE,
        stable=FALSE
){
    set.seed(seed.use)
    if (is.null(study_info)){stop("Please input study_info as trained model")}

    if(is.null(dim(Z)[1])){
        warning("Z is input as a vector, convert Z into matrix with size nZ*1")
        Z=matrix(Z,ncol=1)
    }
    nZ=nrow(Z)
    pZ=ncol(Z) # how many study_info
    if(is.null(W)){pW=0}else{pW=ncol(W)}
    if(is.null(A)){pA=0}else{
        if(is.null(dim(A)[1])){
            if(length(A)==1){
                if(A==1){A=matrix(1,nrow=nZ,ncol=1)}}}
        pA=ncol(A)
    }
    idZ=pA+pW+1
    if(length(study_info)!=pZ){
        stop("When using htlgmm.gwas, input Z as a matrix with size of sample*SNP, and study_info as a list of summary statistics. The columns of Z need to match the study_info.")
    }
    Acolnames=NULL
    if(pA>0){
        Acolnames=colnames(A)
        if(is.null(Acolnames[1])){
            Acolnames=paste0('A',1:pA)}
        if(length(unique(A[,1])) == 1){
            if(unique(A[,1]) == 1){
                Acolnames[1]='intercept'}}
        colnames(A)=Acolnames
    }
    Zcolnames=colnames(Z)
    Wcolnames=colnames(W)
    if(is.null(Zcolnames[1])){
        Zcolnames=paste0('Z',1:pZ)
        colnames(Z)=Zcolnames
    }
    if(is.null(Wcolnames[1])){
        Wcolnames=paste0('W',1:pW)
        colnames(W)=Wcolnames
    }
    Xcolnames<-c(Acolnames,Wcolnames,Zcolnames)

    # unique thetaA
    if(pA!=0){
        if(is.null(A_thetaA)){
            df=data.frame(y,A)
            if(family=="binomial"){
                hat_thetaA_glm=speedglm(y~0+.,data = df,family = binomial())
            }else if(family=="gaussian"){
                hat_thetaA_glm=speedlm(y~0+.,data = df)
            }
            hat_thetaA=c(hat_thetaA_glm$coefficients)
            V_thetaA=vcov(hat_thetaA_glm)
            A_thetaA = prodv_rcpp(A,hat_thetaA)
            if(family == "binomial"){
                expit_A_thetaA=expit_rcpp(A_thetaA)
                mu_prime_A_thetaA=timesv_rcpp(expit_A_thetaA,(1-expit_A_thetaA))
            }else{mu_prime_A_thetaA=1}
            if(stable){
                inv_GammaAA=ginv((1/nZ)*crossprod(A*c(mu_prime_A_thetaA),A))
            }else{
                inv_GammaAA=spdinv((1/nZ)*crossprod_rcpp(A*c(mu_prime_A_thetaA),A))
            }
        }else if(is.null(V_thetaA)){
            stop("When inputing A_thetaA, V_thetaA is needed.")
        }
    }else{
        V_thetaA = NULL
        A_thetaA = 0
        inv_GammaAA = NULL
    }

    # unique beta initial
    if(is.null(AW_betaAW)){
        df=data.frame(y,A,W)
        if(family=="binomial"){
            fit_initial=speedglm(y~0+.,data = df,family = binomial())
        }else if(family=="gaussian"){
            fit_initial=speedlm(y~0+.,data = df)
        }
        beta_initial_AW=c(fit_initial$coefficients)
        AW_betaAW = prod_rcpp(cbind(A,W),beta_initial_AW)
    }

    # Estimation of C
    beta_var_list<-lapply(1:pZ, function(id){
        if(verbose){message(id)}
        if(is.null(study_info[[id]]$Coeff[1])){
            return_list<-list("beta"=NULL,
                              "variance"=NULL)
        }else{
            Zid = Z[,id]
            Z_thetaZ = c(Zid*study_info[[id]]$Coeff)
            V_thetaZ = study_info[[id]]$Covariance
            X_beta = add_rcpp(AW_betaAW,Z_thetaZ)
            XR_theta = add_rcpp(A_thetaA,Z_thetaZ)
            # start model
            inv_C = Delta_opt_fast1(y=y,Z=Zid,W=W,A=A,
                                   family=family,pA=pA,
                                   X_beta=X_beta,
                                   A_thetaA=A_thetaA,
                                   Z_thetaZ=Z_thetaZ,
                                   V_thetaA=V_thetaA,
                                   V_thetaZ=V_thetaZ,
                                   inv_GammaAA=inv_GammaAA)

            if(stable){
                inv_C_svd<-fast.svd(inv_C+diag(1e-15,nrow(inv_C)))
                C_half<-inv_C_svd$v%*%diag(1/sqrt(inv_C_svd$d))%*%t(inv_C_svd$u)
                # Prepare for final model
                pseudo_Xy_list<-pseudo_Xy_fast(C_half=C_half,y=y,Z=Zid,
                                               W=W,A=A,
                                               XR_theta=XR_theta,
                                               X_beta=X_beta,
                                               family=family)
                initial_sf<-nZ/sqrt(nrow(pseudo_Xy_list$pseudo_X))
                pseudo_X<-pseudo_Xy_list$pseudo_X/initial_sf
                pseudo_y<-pseudo_Xy_list$pseudo_y/initial_sf
                fit_final_ols=lm(y~0+.,data = data.frame(y= pseudo_y,pseudo_X))
                beta=fit_final_ols$coefficients
            }else{
                beta=direct_fast1(inv_C,y,Z,W,A,XR_theta,
                                 X_beta,family)
            }

            return_list<-list("beta"=beta)


            # refine C
            X_beta = prodv_rcpp(cbind(A,W,Zid),beta)
            inv_C = Delta_opt_fast1(y=y,Z=Zid,W=W,A=A,
                                   family=family,pA=pA,
                                   X_beta=X_beta,
                                   A_thetaA=A_thetaA,
                                   Z_thetaZ=Z_thetaZ,
                                   V_thetaA=V_thetaA,
                                   V_thetaZ=V_thetaZ,
                                   inv_GammaAA=inv_GammaAA)
            if(stable){
                inv_C_svd<-fast.svd(inv_C+diag(1e-15,nrow(inv_C)))
                C_half<-inv_C_svd$v%*%diag(1/sqrt(inv_C_svd$d))%*%t(inv_C_svd$u)
                ## This need a new formulation
                pseudo_Xy_list<-pseudo_Xy_fast(C_half=C_half,y=y,Z=Zid,
                                               W=W,A=A,
                                               XR_theta=XR_theta,
                                               X_beta=X_beta,
                                               family=family)
                Sigsum_half<-pseudo_Xy_list$pseudo_X/nZ
                Sigsum_scaled<-self_crossprod_rcpp(Sigsum_half)
                inv_Sigsum_scaled<-spdinv(Sigsum_scaled)
                final_v<-diag(inv_Sigsum_scaled)/nZ
            }else{
                final_v<-direct_fast12(inv_C,
                                      y,Z,W,A,XR_theta,
                                      X_beta,family)
            }

            return_list<-c(return_list,list("variance"=final_v))
        }
        return_list
    })
    if(output_SNP_only){
        beta_var_mat<-sapply(1:pZ, function(id){
            if(is.null(beta_var_list[[id]]$beta[1])){
                res =c(NA,NA)
            }else{
                res=c(beta_var_list[[id]]$beta[idZ],beta_var_list[[id]]$variance[idZ])
            }
            res
        })
        beta_var_list<-list("beta"=beta_var_mat[1,],
                            "variance"=beta_var_mat[2,])
    }
    if(output_tmp){
        beta_var_list<-c(beta_var_list,list("AW_betaAW"=AW_betaAW,
                                            "A_thetaA"=A_thetaA,
                                            "V_thetaA"=V_thetaA,
                                            "inv_GammaAA"=inv_GammaAA
        ))
    }
    return(beta_var_list)
}
