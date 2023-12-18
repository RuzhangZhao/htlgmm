
dexpit1<-function(x){return (expit(x)*(1-expit(x)))} # rank3
dexpit2<-function (x) {1/(2+exp(x)+exp(-x))} # rank2
expit1<-function (x)
{
    y <- x
    ix <- (x < 0)
    y[ix] <- exp(x[ix])/(1 + exp(x[ix]))
    y[!ix] <- 1/(1 + exp(-x[!ix]))
    y
} # twice slower
pseudo_Xy_gaussian_finemapping<-function(
        C_half,Z,W,A,y,beta=NULL,hat_thetaA=NULL,study_info=NULL){
    X=cbind(A,Z,W)
    A_thetaA = c(A%*%hat_thetaA)
    pseudo_X=C_half%*%crossprod(cbind(X,Z),X)
    pseudo_y1=c(crossprod(y,X))
    pseudo_y2<-sapply(1:ncol(Z), function(id){
        c(c(A_thetaA+Z[,id]*study_info[[id]]$Coeff)%*%Z[,id])
    })
    pseudo_y=c(c(pseudo_y1,pseudo_y2)%*%C_half)
    list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
}
pseudo_Xy_binomial_finemapping<-function(
        C_half,Z,W,A,y,beta=NULL,hat_thetaA=NULL,study_info=NULL){
    X=cbind(A,Z,W)
    X_beta=X%*%beta
    A_thetaA=A%*%hat_thetaA
    expit_beta=c(expit(X_beta))
    dexpit_beta=c(expit_beta*(1-expit_beta))
    pseudo_X=C_half%*%crossprod(cbind(X,Z),X*dexpit_beta)
    ps_y0=c(crossprod(cbind(X,Z),c(dexpit_beta*X_beta-expit_beta)))
    ps_y1=c(crossprod(y,X))
    ps_y2<-sapply(1:ncol(Z), function(id){
        c(expit(c(A_thetaA+Z[,id]*study_info[[id]]$Coeff))%*%Z[,id])
    })
    pseudo_y=c(c(ps_y0+c(ps_y1,ps_y2))%*%C_half)
    list("pseudo_X"=pseudo_X,"pseudo_y"=pseudo_y)
}
Delta_opt_finemapping<-function(y,Z,W,family,
                       study_info,A=NULL,pA=NULL,
                       beta=NULL,hat_thetaA=NULL,
                       V_thetaA=NULL,use_offset=NULL){
    n_main=length(y)
    DX=cbind(A,W,Z)
    X_beta = DX%*%beta
    A_thetaA = A%*%hat_thetaA
    if(family == "binomial"){
        X_beta=expit(X_beta)
        DZ<-sapply(1:ncol(Z), function(id){
            XR_thetaid<-expit(c(A_thetaA+Z[,id]*study_info[[id]]$Coeff))
            X_beta-XR_thetaid})*Z #the same shape with Z
        DDZ<-sapply(1:ncol(Z), function(id){
            dexpit(A_thetaA+Z[,id]*study_info[[id]]$Coeff)})*Z
        GammaZZ_vec =colmeans(DDZ*Z)
        A_thetaA = expit(A_thetaA)
        mu_prime_A_thetaA = A_thetaA*(1-A_thetaA)
    }else{
        DZ<-sapply(1:ncol(Z), function(id){
            XR_thetaid<-c(A_thetaA+Z[,id]*study_info[[id]]$Coeff)
            X_beta-XR_thetaid
        }) #the same shape with Z
        DDZ = Z
        GammaZZ_vec = colmeans(Z^2)
        mu_prime_A_thetaA = 1
    }
    DX = DX*c(X_beta-y)
    V_U1=(1/n_main)*crossprod(DX)
    V_U2=(1/n_main)*crossprod(DZ)
    Cov_U1U2=(1/n_main)*crossprod(DX,DZ)

    V_thetaZ_vec<-sapply(1:ncol(Z), function(id){
        sqrt(study_info[[id]]$Covariance)})

    Delta22=V_U2+n_main*t((GammaZZ_vec*V_thetaZ_vec)*cor(Z))*(GammaZZ_vec*V_thetaZ_vec)
    #Delta22=V_U2+GammaZZ%*%(n_main*V_thetaZ)%*%t(GammaZZ)
    Delta12=Cov_U1U2
    if(pA!=0){
        GammaZA=(1/n_main)*crossprod(DDZ,A)
        inv_GammaAA=ginv((1/n_main)*crossprod(A*c(mu_prime_A_thetaA),A))
        Cov_U1theta=(1/n_main)*crossprod(DX,A*c(A_thetaA-y))%*%(inv_GammaAA%*%t(GammaZA))
        Cov_U2theta=(1/n_main)*crossprod(DZ,A*c(A_thetaA-y))%*%(inv_GammaAA%*%t(GammaZA))
        Delta22 = Delta22 + GammaZA%*%(n_main*V_thetaA)%*%t(GammaZA)
        + Cov_U2theta+t(Cov_U2theta)
        Delta12 = Delta12 + Cov_U1theta
    }
    Delta = rbind(cbind(V_U1,Delta12),cbind(t(Delta12),Delta22))
    Delta
}

