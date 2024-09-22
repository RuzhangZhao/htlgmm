#' htlgmm
#'
#' htlgmm fits a generalized linear model or cox proportional hazard model via penalized generalized method of moments,
#' i.e. Heterogeneous Transfer Learning via Generalized Method of Moments.
#' The input requires a main study and an external study.
#'
#'
#' @details htlgmm: Heterogeneous Transfer Learning via generalized method of moments(GMM).
#'
#' @param y The outcome variable, which can be continuous, binary or time-to-event data. For coxph model, y should be a list with two items including 'time' and 'event'.
#' @param Z The overlapping features in both main and external studies.
#' @param W The unmatched features only in main study.
#' @param ext_study_info The trained model from external study, including estimated coefficients, and estimated variance-covariance matrix.
#' The 'ext_study_info' is in the format of list. The first item is 'Coeff', the second item is 'Covariance'.
#' E.g. ext_study_info = list(list("Coeff"=coeff,"Covariance"=covariance))
#' @param A The covariates for study-specific adjustment. The default is 'default', which is 'NULL' for 'gaussian' family or 'cox' family, '1' for 'binomial' family.
#' Other than c('default',NULL,1), A must be a matrix whose dimension is the same as the sample dimension or Z and W.
#' For continuous variable, we suggest scaling the features Z, W to eliminate intercept term.  If 'A = NULL', there is no intercept term included.
#' For binary variable, we use intercept term by 'A=1' to adjust for different binary trait ratios in main and external studies.
#' If there is only intercept term in A, we use 'A=1'.
#' A are the features working for adjustment in reduced model, but A is not summarized in summary statistics(input:study_info).
#' @param penalty_type The penalty type for htlgmm, chosen from c("none","lasso","adaptivelasso","ridge","elasticnet"). The default is "none".
#' If 'penalty_type = 'none' ', we use without penalty. (For continuous y, we use ordinary least square, and for binary y, we use logistic regression without penalty.)
#' @param family The family is chosen from c("gaussian","binomial"). Linear regression for "gaussian" and logistic regression for "binomial".
#' @param initial_with_type Get initial estimation for beta using main study data only
#' by cross validation using (penalty) regression, chosen from c("ridge","lasso","glm"). The default is "glm". If penalty_type = 'glm',
#' for continuous y, we use ordinary least square, and for binary y, we use logistic regression without penalty.)
#' @param beta_initial The initial estimation for beta if a consistent estimator is available.
#' The default is NULL, and main study is used for initial estimation according to 'initial_with_type'.
#' @param weight_adaptivelasso The customized adaptive-weight adaptive lasso when using penalty_type = 'adaptivelasso'. The default is NULL. When using default weight_adaptivelasso, user may follow the instruction from 'Zou, H. (2006). The adaptive lasso and its oracle properties.'
#' The input of 'weight_adaptivelasso' can either be with all features (the same length as beta_initial) or without intercept.
#' @param alpha The elasticnet mixing parameter, between 0 and 1, which is corresponding to the alpha for glmnet. Only used when penalty_type = 'elasticnet'.
#' @param hat_thetaA If A is not NULL, one can provide hat_thetaA as the input. If 'hat_thetaA = NULL', we estimate hat_thetaA with glm by main study.
#' @param V_thetaA If A is not NULL, one can provide V_thetaA as the input. If 'V_thetaA = NULL', we estimate V_thetaA with glm by main study.
#' @param use_offset Whether to use offset regarding the external model estimated coefficient. The default is FALSE.
#' @param robust Whether to apply sandwich formula to compute the variance-covariance matrix of hat_thetaA.The default is TRUE.
#' For coxph model, robust is also about whether we apply the robust variance for the estimating equations.
#' @param remove_penalty_Z Do not penalize Z if it is TRUE. The default is FALSE.
#' @param remove_penalty_W Do not penalize W if it is TRUE. The default is FALSE.
#' @param inference Whether to do inference without penalty or post-selection inference with adaptive lasso penalty. The default is 'default', which is TRUE when penalty_type = 'none' or 'adaptivelasso', and FALSE otherwise.
#' @param fix_C When fix_C = NULL, the optimal C is computed. When user wants to customize the fix_C, please match its dimension as dim(A)+2*dim(Z)+dim(W) and make sure it is positive definite.
#' @param fix_inv_C When fix_inv_C = NULL, the optimal C is computed. When user wants to customize the fix_inv_C, please match its dimension as dim(A)+2*dim(Z)+dim(W) and make sure it is positive definite. When fix_C and fix_inv_C are both given, the fix_C will be used.
#' @param refine_C When computing the variance, whether recompute the weighting matrix C using final estimated beta. The default is FALSE, which means we do not update the C. Instead, keep it as what we use in the algorithm.
#' @param sqrt_matrix The method to split weighting matrix into square root matrix. Select from c('svd','cholesky'), where 'cholesky' generates faster computation.
#' @param fix_lambda Without cross validation, fix the lambda. The default is NULL.
#' @param lambda_list Customize the input lambda list for validation. The default is NULL to generate lambda list according to glmnet.
#' @param fix_ratio The fixed ratio for two-lambda strategy. The ratio is multiplied for Z features. The default is NULL. If it is NULL, select the best ratio via cross validation.
#' @param fix_weight The fixed weight for weighting matrix of external study.
#' @param gamma_adaptivelasso The gamma for adaptive lasso. Select from c(1/2,1,2). The default is 1/2.
#' @param use_sparseC Whether to use approximate version of weighting matrix C using only diagonal terms as nonzeros.
#' If approximation, use the diagonal of inverse of C(inv_C) to approximate the inv_C. The default is FALSE.
#' When main study sample size is limited, use_sparseC = TRUE is recommended.
#' When main study sample size is large enough, use_sparseC = FALSE is recommended.
#' @param seed.use The seed for  97.
#' @param output_all_betas Output all betas.
#'
#' @return \itemize{
#'  \item{beta:} The target coefficient estimation, the features will go in the order of (A,Z,W).
#'  \item{lambda_list:} The lambda list for cross validation.
#'  \item{ratio_list:} The ratio list for cross validation.
#'  \item{fix_lambda:} If the fix_lambda is not null, we output fix_lambda.
#'  \item{fix_ratio:} If the fix_ratio is not null, we output fix_ratio.
#'  \item{selected_vars:} For inference or post-selection inference, we output the inference results by a list. \itemize{
#'  \item{position:} The index of nonzero positions, the index comes from X = (A,Z,W).
#'  \item{name:} The feature name of nonzero positions. If there is no default name, we name it after Ai, Zi, Wi.
#'  \item{coef:} The coefficients of nonzero positions.
#'  \item{variance:} The variances for features with glm inference, for selected features with post-selection inference.
#'  \item{pval:} For p values for nonzero positions.
#'  \item{FDR_adjust_position:} The FDR adjusted positions passing significant level 0.05 after BH adjustment (Benjamini & Hochberg).
#'  \item{FDR_adjust_name:} The feature name based on FDR_adjust_position.
#'  }
#'  }
#'
#' @import glmnet
#' @import stats
#' @importFrom caret createFolds createDataPartition
#' @importFrom corpcor fast.svd
#' @importFrom magic adiag
#' @importFrom MASS ginv
#' @importFrom speedglm speedglm speedlm
#' @importFrom survival coxph Surv
#' @export
#'
#'
htlgmm<-function(
        y,Z,W=NULL,
        ext_study_info=NULL,
        A="default",
        penalty_type = "none",
        family = "gaussian",
        initial_with_type = "glm",
        beta_initial = NULL,
        weight_adaptivelasso = NULL,
        alpha = NULL,
        hat_thetaA = NULL,
        V_thetaA = NULL,
        use_offset = TRUE,
        robust = FALSE,
        remove_penalty_Z = FALSE,
        remove_penalty_W = FALSE,
        inference = 'default',
        fix_C = NULL,
        fix_inv_C = NULL,
        refine_C = FALSE,
        sqrt_matrix ="cholesky",
        fix_lambda = NULL,
        lambda_list = NULL,
        fix_ratio = NULL,
        fix_weight = NULL,
        gamma_adaptivelasso = 1/2,
        use_sparseC = FALSE,
        seed.use = 97,
        output_all_betas=FALSE
){

    if(!family %in% c("gaussian","binomial","cox")){
        stop("Select family from c('gaussian','binomial','cox')")
    }
    if(is.null(dim(A)[1])){
        if(length(A)==1){
            if(A=='default'){if(family == "binomial"){A=1}else{A=NULL}}else{
                if(A!=1){warning("If A is not a matrix, A should be selected from c('default',NULL,1).")}
                if(A==1&family=='cox'){warning("Coxph model usually does not include intercept term.")}
            }
        }else{
            warning("If A is not selected from c('default',NULL,1), A must be a matrix.")
        }}
    use_cv = FALSE
    nfolds = 10
    foldid = NULL
    tune_ratio = FALSE
    ratio_list = NULL
    tune_weight = FALSE
    weight_list = NULL
    tune_weight_method = 1
    type_measure = "default"
    nlambda = 100
    lambda.min.ratio = 0.0001

    if(family == 'cox'){
        res<-htlgmm.cox.default(y,Z,W,ext_study_info,A,penalty_type,
                                initial_with_type,beta_initial,
                                weight_adaptivelasso,alpha,
                                hat_thetaA,V_thetaA,robust,remove_penalty_Z,
                                remove_penalty_W,inference,fix_C,fix_inv_C,refine_C,
                                sqrt_matrix,use_cv,type_measure,nfolds,foldid,fix_lambda,
                                lambda_list,nlambda,lambda.min.ratio,
                                tune_weight,fix_weight,weight_list,tune_weight_method,
                                gamma_adaptivelasso,use_sparseC,seed.use)
    }else{
        V_thetaA_sandwich=robust
        res<-htlgmm.default(y,Z,W,ext_study_info,A,penalty_type,
                            family,initial_with_type,beta_initial,
                            weight_adaptivelasso,alpha,
                            hat_thetaA,V_thetaA,use_offset,
                            V_thetaA_sandwich,remove_penalty_Z,
                            remove_penalty_W,inference,fix_C,fix_inv_C,refine_C,
                            sqrt_matrix,use_cv,type_measure,nfolds,foldid,
                            fix_lambda,lambda_list,nlambda,
                            lambda.min.ratio,tune_ratio,fix_ratio,
                            ratio_list,tune_weight,fix_weight,
                            weight_list,tune_weight_method,gamma_adaptivelasso,
                            use_sparseC,seed.use,output_all_betas)
    }
    return(res)
}



#' Cross validation for htlgmm.
#'
#' htlgmm fits a generalized linear model via penalized generalized method of moments,
#' i.e. Heterogeneous Transfer Learning via Generalized Method of Moments.
#' The input requires main study and external study.
#' cv.htlgmm does k-fold cross validation for htlgmm.
#'
#'
#' @details Cross validation for htlgmm.
#'
#' @param y The variable of interest, which can be continuous or binary.
#' @param Z The overlapping features in both main and external studies.
#' @param W The unmatched features only in main study.
#' @param ext_study_info The trained model from external study, including estimated coefficients, and estimated variance-covariance matrix.
#' The 'ext_study_info' is in the format of list. The first item is 'Coeff', the second item is 'Covariance'.
#' E.g. ext_study_info = list(list("Coeff"=coeff,"Covariance"=covariance))
#' @param A The covariates for study-specific adjustment. The default is 'default', which is 'NULL' for 'gaussian' family, '1' for 'binomial' family.
#' Other than c('default',NULL,1), A must be a matrix whose dimension is the same as the sample dimension or Z and W.
#' For continuous variable, we suggest scaling the features Z, W to eliminate intercept term.  If 'A = NULL', there is no intercept term included.
#' For binary variable, we use intercept term by 'A=1' to adjust for different binary trait ratios in main and external studies.
#' If there is only intercept term in A, we use 'A=1'.
#' A are the features working for adjustment in reduced model, but A is not summarized in summary statistics(input:ext_study_info).
#' @param penalty_type The penalty type for htlgmm, chosen from c("none","lasso","adaptivelasso","ridge","elasticnet"). The default is "lasso".
#' If 'penalty_type = 'none' ', we use without penalty. (For continuous y, we use ordinary least square, and for binary y, we use logistic regression without penalty.)
#' @param family The family is chosen from c("gaussian","binomial"). Linear regression for "gaussian" and logistic regression for "binomial".
#' @param initial_with_type Get initial estimation for beta using main study data only
#' by cross validation using penalty regression, chosen from c("ridge","lasso","glm"). The default is "lasso". If penalty_type = 'glm',
#' for continuous y, we use ordinary least square, and for binary y, we use logistic regression without penalty.
#' @param beta_initial The initial estimation for beta if a consistent estimator is available.
#' The default is NULL, and main study is used for initial estimation according to 'initial_with_type'.
#' @param weight_adaptivelasso The customized adaptive-weight of adaptive lasso when using penalty_type = 'adaptivelasso'. The default is NULL. When using default weight_adaptivelasso, user may follow the instruction from 'Zou, H. (2006). The adaptive lasso and its oracle properties.'
#' The input of 'weight_adaptivelasso' can either be with all features (the same length as beta_initial) or without intercept.
#' @param alpha The elasticnet mixing parameter, between 0 and 1, which is corresponding to the alpha for glmnet. Only used when penalty_type = 'elasticnet'.
#' @param hat_thetaA If A is not NULL, one can provide hat_thetaA as the input. If 'hat_thetaA = NULL', we estimate hat_thetaA with glm by main study.
#' @param V_thetaA If A is not NULL, one can provide V_thetaA as the input. If 'V_thetaA = NULL', we estimate V_thetaA with glm by main study.
#' @param use_offset Whether to use offset regarding the external model estimated coefficient. The default is FALSE.
#' @param robust Whether to apply sandwich formula to compute the variance-covariance matrix of hat_thetaA.The default is 'default', which is TRUE for coxph model, and FALSE for others.
#' For coxph model, robust is also about whether we apply the robust variance for the estimating equations.
#' @param remove_penalty_Z Do not penalize Z if it is TRUE. The default is FALSE.
#' @param remove_penalty_W Do not penalize W if it is TRUE. The default is FALSE.
#' @param inference Whether to do inference without penalty or post-selection inference with adaptive lasso penalty. The default is 'default', which is TRUE when penalty_type = 'none' or 'adaptivelasso', and FALSE otherwise.
#' @param fix_C When fix_C = NULL, the optimal C is computed. When user wants to customize the fix_C, please match its dimension as dim(A)+2*dim(Z)+dim(W) and make sure it is positive definite.
#' @param fix_inv_C When fix_inv_C = NULL, the optimal C is computed. When user wants to customize the fix_inv_C, please match its dimension as dim(A)+2*dim(Z)+dim(W) and make sure it is positive definite. When fix_C and fix_inv_C are both given, the fix_C will be used.
#' @param refine_C When computing the variance, whether recompute the weighting matrix C using final estimated beta.
#' @param sqrt_matrix The method to split weighting matrix into square root matrix. Select from c('svd','cholesky'), where 'cholesky' generates faster computation.
#' @param use_cv Whether to use cross validation to determine the best lambda (or ratio).
#' @param type_measure Select from c("default", "mse", "deviance", "auc"). Default is mse(liner), auc(logistic). 'deviance' is another choice for binary y.
#' @param nfolds The fold number for cross validation. Only work for use_cv = TRUE.The default is 10.
#' @param foldid An optional vector of values with the length of sample size, identifying what fold each observation is in. If supplied, nfolds can be missing.
#' @param fix_lambda Without cross validation, fix the lambda. The default is NULL.
#' @param lambda_list Customize the input lambda list for validation. The default is NULL to generate lambda list according to glmnet.
#' @param nlambda The number of lambda values - default is 100.
#' @param lambda.min.ratio Smallest value for lambda, as a fraction of lambda.max. The default is 0.0001.
#' @param tune_ratio Whether to use two-lambda stratgey. The default is FALSE. This is not applied to coxph model.
#' @param fix_ratio The fixed ratio for two-lambda strategy. The ratio is multiplied for Z features. The default is NULL. If it is NULL, select the best ratio via cross validation.
#' @param ratio_list The ratio list if it is preset. The default is NULL and ratio list will be generated.
#' @param tune_weight Whether to assign tuning weight for the regularization of weighting matrix. The default is FALSE.
#' @param fix_weight The fixed weight for the regularization of weighting matrix. The default is NULL. If it is NULL, select the best weight via cross validation.
#' @param weight_list The weight list if it is preset. The default is NULL and weight list will be generated.
#' @param tune_weight_method Method for weighting matrix regularization. The default is 'multiplicative shrinkage', which can also be used as 'mshrink' or 'ms'. The other choice is 'ridge'.
#' @param gamma_adaptivelasso The gamma for adaptive lasso. Select from c(1/2,1,2). The default is 1/2.
#' @param use_sparseC Whether to use approximate version of weighting matrix C.
#' If approximation, use the diagonal of inverse of C(inv_C) to approximate the inv_C. The default is FALSE.
#' When main study sample size is limited, use_sparseC = TRUE is recommended.
#' When main study sample size is large enough, use_sparseC = FALSE is recommended.
#' @param seed.use The seed for  97.
#' @param output_all_betas output_all_betas
#'
#' @return \itemize{
#'  \item{beta:} The target coefficient estimation, the features will go in the order of (A,Z,W).
#'  \item{lambda_list:} The lambda list for cross validation.
#'  \item{ratio_list:} The ratio list for cross validation.
#'  \item{fix_lambda:} If the fix_lambda is not null, we output fix_lambda.
#'  \item{fix_ratio:} If the fix_ratio is not null, we output fix_ratio.
#'  \item{lambda_min:} The selected best lambda by cross validation.
#'  \item{ratio_min:} The selected best ratio by cross validation.
#'  \item{cv_mse:} The mean square error(mse) when family = "gaussian", and use_cv = TRUE.
#'  \item{cv_dev:} The deviance(dev) when family = "binomial", and use_cv = TRUE.
#'  \item{cv_auc:} The area under the curve of sensitivity specificity when family = "binomial", and use_cv = TRUE.
#'  \item{selected_vars:} For inference or post-selection inference, we output the inference results by a list. \itemize{
#'  \item{position:} The index of nonzero positions, the index comes from X = (A,Z,W).
#'  \item{name:} The feature name of nonzero positions. If there is no default name, we name it after Ai, Zi, Wi.
#'  \item{coef:} The coefficients of nonzero positions.
#'  \item{variance:} The variances for features with none penalty inference or for selected features with post-selection inference.
#'  \item{variance_covariance:} The variance-covariance matrix for features with none penalty inference or for selected features with post-selection inference.
#'  \item{pval:} For p values for nonzero positions.
#'  \item{FDR_adjust_position:} The FDR adjusted positions passing significant level 0.05 after BH adjustment (Benjamini & Hochberg).
#'  \item{FDR_adjust_name:} The feature name based on FDR_adjust_position.
#'  }
#'  }
#'
#'
#'
#' @import glmnet
#' @import stats
#' @importFrom caret createFolds createDataPartition
#' @importFrom corpcor fast.svd
#' @importFrom magic adiag
#' @importFrom MASS ginv
#' @importFrom pROC auc
#' @importFrom speedglm speedglm speedlm
#' @importFrom survival coxph Surv
#' @export
#'
#'
cv.htlgmm<-function(
        y,Z,W=NULL,
        ext_study_info=NULL,
        A="default",
        penalty_type = "lasso",
        family = "gaussian",
        initial_with_type = "lasso",
        beta_initial = NULL,
        weight_adaptivelasso=NULL,
        alpha = NULL,
        hat_thetaA = NULL,
        V_thetaA = NULL,
        use_offset = TRUE,
        robust = FALSE,
        remove_penalty_Z = FALSE,
        remove_penalty_W = FALSE,
        inference = 'default',
        fix_C = NULL,
        fix_inv_C=NULL,
        refine_C = FALSE,
        sqrt_matrix = 'cholesky',
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
        tune_weight_method = 'mshrink',
        gamma_adaptivelasso = 1/2,
        use_sparseC = FALSE,
        seed.use = 97,
        output_all_betas=FALSE
){
    if(tune_weight_method %in% c("mshrink","ms","multiplicative shrinkage") ){
        tune_weight_method = 1
    }else if(tune_weight_method == "ridge"){
        tune_weight_method = 2
    }
    if(!family %in% c("gaussian","binomial","cox")){
        stop("Select family from c('gaussian','binomial','cox')")
    }

    if(is.null(dim(A)[1])){
        if(length(A)==1){
            if(A=='default'){if(family == "binomial"){A=1}else{A=NULL}}else{
                if(A!=1){warning("If A is not a matrix, A should be selected from c('default',NULL,1).")}
                if(A==1&family=='cox'){warning("Coxph model usually does not include intercept term.")}
            }

        }else{
            warning("If A is not selected from c('default',NULL,1), A must be a matrix.")
        }}

    if(family == 'cox'){
        if(robust == 'default'){robust = TRUE}
        res<-htlgmm.cox.default(y,Z,W,ext_study_info,A,penalty_type,
                                initial_with_type,beta_initial,
                                weight_adaptivelasso,alpha,
                                hat_thetaA,V_thetaA,robust,remove_penalty_Z,
                                remove_penalty_W,inference,fix_C,fix_inv_C,refine_C,
                                sqrt_matrix,use_cv,type_measure,nfolds,foldid,fix_lambda,
                                lambda_list,nlambda,lambda.min.ratio,
                                tune_weight,fix_weight,weight_list,tune_weight_method,
                                gamma_adaptivelasso,use_sparseC,seed.use)
    }else{
        if(robust == 'default'){robust = FALSE}
        V_thetaA_sandwich=robust
        res<-htlgmm.default(y,Z,W,ext_study_info,A,penalty_type,
                            family,initial_with_type,beta_initial,
                            weight_adaptivelasso,alpha,hat_thetaA,V_thetaA,
                            use_offset,V_thetaA_sandwich,remove_penalty_Z,
                            remove_penalty_W,inference,fix_C,fix_inv_C,refine_C,
                            sqrt_matrix,use_cv,type_measure,nfolds,foldid,
                            fix_lambda,lambda_list,nlambda,
                            lambda.min.ratio,tune_ratio,fix_ratio,
                            ratio_list,tune_weight,fix_weight,
                            weight_list,tune_weight_method,gamma_adaptivelasso,
                            use_sparseC,seed.use,output_all_betas)
    }

    return(res)
}

#' GWAS version of htlgmm.
#'
#' GWAS version of htlgmm fits a generalized linear model via generalized method of moments,
#' The input requires main study and external study.
#'
#'
#' @details GWAS for htlgmm.
#'
#' @param y The variable of interest, which can be continuous or binary.
#' @param Z The overlapping features in both main and external studies.
#' @param W The unmatched features only in main study, the default is NULL.
#' @param ext_study_info The trained model from external GWAS study, including estimated coefficients, and estimated variance.
#' The 'ext_study_info' is in the format of list. The first item is 'Coeff', the second item is 'Covariance'.
#' E.g. ext_study_info = list(list("Coeff"=coeff,"Covariance"=variance))
#' @param A The covariates for study-specific adjustment. The default is 'default', which is 'NULL' for 'gaussian' family, '1' for 'binomial' family.
#' Other than c('default',NULL,1), A must be a matrix whose dimension is the same as the sample dimension or Z and W.
#' For continuous variable, we suggest scaling the features Z, W to eliminate intercept term.  If 'A = NULL', there is no intercept term included.
#' For binary variable, we use intercept term by 'A=1' to adjust for different binary trait ratios in main and external studies.
#' If there is only intercept term in A, we use 'A=1'.
#' A are the features working for adjustment in reduced model, but A is not summarized in summary statistics(input:ext_study_info).
#' @param family The family is chosen from c("gaussian","binomial"). Linear regression for "gaussian" and logistic regression for "binomial".
#' @param beta_initial The beta_initial list matching ext_study_info, which provides the initial value of beta, e.g. fitting full model using plink2. If there is only one SNP, the beta_initial can be a vector.
#' The default is NULL. The beta_initial should go in the order of (A,W,Z), where Z is for each SNP.
#' @param repeated_term Default is NULL, which includes A_hat_thetaA, V_thetaA, and inv_GammaAA.
#' @param refine_C When computing the variance, whether recompute the weighting matrix C using final estimated beta.
#' @param output_SNP_only Default is TRUE.
#' @param seed.use The seed for  97.
#' @param verbose Default is FALSE.
#' @param output_tmp Default is FALSE
#'
#' @return \itemize{
#'  \item{beta:} The target coefficient estimation, the features will go in the order of (A,Z,W).
#'  \item{variance:} The lambda list for cross validation.
#'  }
#'
#' @useDynLib htlgmm
#' @import glmnet
#' @import stats
#' @import Rcpp
#' @import RcppEigen
#' @importFrom corpcor fast.svd
#' @importFrom magic adiag
#' @importFrom MASS ginv
#' @importFrom speedglm speedglm speedlm
#' @export
#'
gwas.htlgmm<-function(
        y,Z,W=NULL,
        ext_study_info=NULL,
        A='default',
        family = "gaussian",
        beta_initial=NULL,
        repeated_term = NULL,
        refine_C = FALSE,
        output_SNP_only=TRUE,
        seed.use = 97,
        verbose = FALSE,
        output_tmp=FALSE
){
    if(!family %in% c("gaussian","binomial")){
        stop("Select family from c('gaussian','binomial')")
    }

    if(is.null(dim(A)[1])){
        if(length(A)==1){
            if(A=='default'){if(family == "gaussian"){A=NULL}else{A=1}}else{
                if(!is.null(A)){
                    if(A!=1){warning("If A is not a matrix, A should be selected from c('default',NULL,1).")}}
            }}else{
                warning("If A is not selected from c('default',NULL,1), A must be a matrix.")
            }}

    res<-htlgmm.gwas.default(y,Z,W,ext_study_info,A,family,beta_initial,
                             repeated_term,refine_C,
                             output_SNP_only,seed.use,
                             verbose,output_tmp)
    return(res)
}


#' Function to compute repeated term in gwas.htlgmm
#'
#' @details expit(s)=exp(s)/{1+exp(s)}
#'
#' @param y response variable
#' @param A covariates adjusted in external GWAS
#' @param family gaussian or binomial
#' @return the repeated term list.
#'
#'
#' @useDynLib htlgmm
#' @import Rcpp
#' @importFrom speedglm speedglm speedlm
#' @export
#'
#'
gwas.htlgmm.repeated.term<-function(y,A=NULL,family = "gaussian"){
    compute_repeated_term(y,A,family)
}



#' Expit Function from Rcpp
#'
#' @details expit(s)=exp(s)/{1+exp(s)}
#'
#' @param x x can be a single number or a vector or a matrix.
#' @return the element-wise expit value.
#'
#'
#' @useDynLib htlgmm
#' @import Rcpp
#' @export
#'
#'
expit<-function(x){
    if (is.null(dim(x)[1])){
        expitx=expit_rcpp(x)
    }else{
        expitx=expitm_rcpp(x)
    }
    expitx
}


#' Fast Matrix Inverse for Positive Definite Matrix by Cholesky Decomposition
#'
#' @details Implemented by RcppEigen
#'
#' @param X The positive definite matrix.
#' @return The inverse of X.
#'
#'
#' @useDynLib htlgmm
#' @import Rcpp
#' @export
#'
#'
#'
choinv<-function(X){
    if(!is.matrix(X)){
        stop("The input must be in matrix format")
    }
    choinv_rcpp(X)
}


#' Fast Matrix Square Root Inverse for Positive Definite Matrix by Cholesky Decomposition
#'
#' @details Implemented by RcppEigen
#'
#' @param X The positive definite matrix.
#' @return The inverse of X.
#'
#'
#' @useDynLib htlgmm
#' @import Rcpp
#' @export
#'
#'
#'
sqrtchoinv<-function(X){
    if(!is.matrix(X)){
        stop("The input must be in matrix format")
    }
    sqrtchoinv_rcpp(X)
}

#' Fast Matrix Dot Product with Transpose to Replace "crossprod".
#'
#' @details Implemented by RcppEigen
#'
#' @param x The input matrix.
#' @param y The input matrix or vector or NULL.
#' @return The dot product of the transpose of 'x' and 'y'.
#'
#'
#' @useDynLib htlgmm
#' @import Rcpp
#' @export
#'
#'
#'
crossprod_fast<-function(x,y=NULL){
    if(is.null(y)){
        if(!is.matrix(x)){
            x<-matrix(x,ncol = 1)
        }
        return(self_crossprod_rcpp(x))
    }
    if(is.matrix(x)){
        if(is.matrix(y)){
            if(nrow(x)!=nrow(y)){
                stop("dimensions do not match!")
            }else{
                return(crossprod_rcpp(x,y))
            }
        }else{
            if(nrow(x)!=length(y)){
                stop("dimensions do not match!")
            }else{
                return(crossprodv_rcpp(x,y))
            }
        }
    }else{
        if(is.matrix(y)){
            stop("dimensions do not match!")
        }else{
            return(sum(timesv_rcpp(x,y)))
        }
    }
}

#' Fast Matrix Dot Product with Transpose to Replace "\%*\%".
#'
#' @details Implemented by RcppEigen
#'
#' @param x The input matrix.
#' @param y The input matrix or vector or NULL.
#' @return The dot product of 'x' and 'y'.
#'
#'
#' @useDynLib htlgmm
#' @import Rcpp
#' @export
#'
#'
#'
prod_fast<-function(x,y){
    if(is.matrix(x)){
        if(is.matrix(y)){
            if(ncol(x)!=nrow(y)){
                stop("dimensions do not match!")
            }else{
                return(prod_rcpp(x,y))
            }
        }else{
            if(ncol(x)!=length(y)){
                stop("dimensions do not match!")
            }else{
                return(prodv_rcpp(x,y))
            }
        }
    }else{
        if(is.matrix(y)){
            stop("dimensions do not match!")
        }else{
            return(sum(timesv_rcpp(x,y)))
        }
    }
}




