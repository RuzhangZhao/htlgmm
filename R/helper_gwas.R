Delta_opt_fast<-function(y,Z,W,A,family,pA,
                    X_beta,XR_theta,V_thetaA,V_thetaZ){
    n_main=length(y)
    X=cbind(A,W,Z)
    if(family == "binomial"){
        X_beta=expit(X_beta)
        XR_theta=expit(XR_theta)
        mu_prime_XR_theta=XR_theta*(1-XR_theta)
    }else{mu_prime_XR_theta=1}
    V_U1=(1/n_main)*crossprod(X*c(X_beta-y))
    V_U2=(1/n_main)*crossprod(Z*c(X_beta-XR_theta))
    Cov_U1U2=(1/n_main)*crossprod(X*c(X_beta-y),Z*c(X_beta-XR_theta))
    GammaZZ=(1/n_main)*crossprod(Z*c(mu_prime_XR_theta),Z)
    Delta22=V_U2 +GammaZZ%*%(n_main*V_thetaZ)%*%t(GammaZZ)
    Delta12=Cov_U1U2
    if(pA!=0){
        GammaZA=(1/n_main)*crossprod(Z*c(mu_prime_XR_theta),A)
        inv_GammaAA=ginv((1/n_main)*crossprod(A*c(mu_prime_XR_theta),A))
        Cov_U1theta=(1/n_main)*crossprod(X*c(X_beta-y),A*c(XR_theta-y))%*%(inv_GammaAA%*%t(GammaZA))
        Cov_U2theta=(1/n_main)*crossprod(Z*c(X_beta-XR_theta),A*c(XR_theta-y))%*%(inv_GammaAA%*%t(GammaZA))
        Delta22 = Delta22 + GammaZA%*%(n_main*V_thetaA)%*%t(GammaZA)
        + Cov_U2theta+t(Cov_U2theta)
        Delta12 = Delta12 + Cov_U1theta
    }
    Delta = rbind(cbind(V_U1,Delta12),cbind(t(Delta12),Delta22))
    Delta
}

pseudo_Xy_fast<-function(
        C_half,y,Z,W,A,XR_theta,
        X_beta,family){
    X=cbind(A,W,Z)
    if(family == "binomial"){
        XR_theta=c(expit(XR_theta))
        expit_beta=c(expit(X_beta))
        dexpit_beta=c(expit_beta*(1-expit_beta))
        ps_y0=c(crossprod(cbind(X,Z),c(dexpit_beta*X_beta-expit_beta)))
    }else{
        dexpit_beta=1
        ps_y0=0
    }
    pseudo_X=C_half%*%crossprod(cbind(X,Z),X*dexpit_beta)

    ps_y1=c(crossprod(y,X))
    ps_y2=c(crossprod(XR_theta,Z))
    pseudo_y=c(c(ps_y0+c(ps_y1,ps_y2))%*%C_half)
    list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
}


htlgmm.gwas.default<-function(
        y,Z,W=NULL,
        study_info=NULL,
        A=NULL,
        family = "gaussian",
        output_SNP_only=TRUE,
        seed.use = 97,
        verbose = FALSE
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
        df=data.frame(y,A)
        if(family=="binomial"){
            hat_thetaA_glm=speedglm(y~0+.,data = df,family = binomial())
        }else if(family=="gaussian"){
            hat_thetaA_glm=speedlm(y~0+.,data = df)
        }
        hat_thetaA=hat_thetaA_glm$coefficients
        V_thetaA=vcov(hat_thetaA_glm)
        A_thetaA = c(A%*%hat_thetaA)
    }else{
        V_thetaA = NULL
        A_thetaA = 0
    }

    # unique beta initial
    df=data.frame(y,A,W)
    if(family=="binomial"){
        fit_initial=speedglm(y~0+.,data = df,family = binomial())
    }else if(family=="gaussian"){
        fit_initial=speedlm(y~0+.,data = df)
    }
    beta_initial_AW=fit_initial$coefficients
    AW_betaAW = c(cbind(A,W)%*%beta_initial_AW)

    # Estimation of C
    beta_var_list<-lapply(1:pZ, function(id){
        if(verbose){message(id)}
        if(is.null(study_info[[id]]$Coeff[1])){
            return_list<-list("beta"=NULL,
                                "variance"=NULL)
        }else{
            Zid = Z[,id]
            Z_thetaZ = c(Zid*study_info[[id]]$Coeff)
            X_beta = AW_betaAW+Z_thetaZ
            XR_theta = A_thetaA+Z_thetaZ

            # start model
            inv_C = Delta_opt_fast(y=y,Z=Zid,W=W,A=A,
                                   family=family,pA=pA,
                                   X_beta=X_beta,
                                   XR_theta=XR_theta,
                                   V_thetaA=V_thetaA,
                                   V_thetaZ=study_info[[id]]$Covariance)


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
            return_list<-list("beta"=beta)


            # refine C
            X_beta = cbind(A,W,Zid)%*%beta
            inv_C = Delta_opt_fast(y=y,Z=Zid,W=W,A=A,
                                   family=family,pA=pA,
                                   X_beta=X_beta,
                                   XR_theta=XR_theta,
                                   V_thetaA=V_thetaA,
                                   V_thetaZ=study_info[[id]]$Covariance)
            inv_C_svd<-fast.svd(inv_C+diag(1e-15,nrow(inv_C)))
            C_half<-inv_C_svd$v%*%diag(1/sqrt(inv_C_svd$d))%*%t(inv_C_svd$u)
            ## This need a new formulation
            pseudo_Xy_list<-pseudo_Xy_fast(C_half=C_half,y=y,Z=Zid,
                                           W=W,A=A,
                                           XR_theta=XR_theta,
                                           X_beta=X_beta,
                                           family=family)
            Sigsum_half<-pseudo_Xy_list$pseudo_X/nZ
            Sigsum_scaled<-crossprod(Sigsum_half)
            inv_Sigsum_scaled<-solve(Sigsum_scaled)
            final_v<-diag(inv_Sigsum_scaled)/nZ
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
    return(beta_var_list)
}

