#' htlgmm:
#'
#' htlgmm fits a generalized linear model via penalized generalized method of moments,
#' i.e. Heterogeneous Transfer Learning via Generalized Method of Moments.
#' The input requires main study and external study.
#'
#'
#' @details htlgmm: Heterogeneous Transfer Learning via generalized method of moments(GMM).
#'
#' @param y The variable of interest, which can be continuous or binary.
#' @param Z The overlapping features in both main and external studies.
#' @param W The unmatched features only in main study, the default is NULL.
#' @param study_info The trained model from external study, including estimate coefficients, estimated variance-covariance matrix and sample size.
#' The 'study_info' is in the format of list. The first item is 'Coeff', the second iterm is 'Covariance', and the third item is 'Sample_size'.
#' @param A The covariates for study-specific adjustment. The default is 'default', which is 'NULL' for 'gaussian' family, '1' for 'binomial' family.
#' Other than c('default',NULL,1), A must be a matrix whose dimension is the same as the sample dimension or Z and W.
#' For continuous variable, we suggest scaling the features Z, W to eliminate intercept term.  If 'A = NULL', there is no intercept term included.
#' For binary variable, we use intercept term by 'A=1' to adjust for different binary trait ratios in main and external studies.
#' If there is only intercept term in A, we use 'A=1'.
#' A are the features working for adjustment in reduced model, but A is not summarized in summary statistics(input:study_info).
#' @param penalty_type The penalty type for htlgmm, chosen from c("none","lasso","adaptivelasso","ridge"). The default is "none".
#' If 'penalty_type = 'none' ', we use without penalty. (For continuous y, we use ordinary least square, and for binary y, we use logistic regression without penalty.)
#' @param family The family is chosen from c("gaussian","binomial"). Linear regression for "gaussian" and logistic regression for "binomial".
#' @param initial_with_type Get initial estimation for beta using main study data only
#' by cross validation using (penalty) regression, chosen from c("ridge","lasso","glm"). The default is "ridge". If penalty_type = 'glm',
#' for continuous y, we use ordinary least square, and for binary y, we use logistic regression without penalty.)
#' @param beta_initial The initial estimation for beta if a consistent estimator is available.
#' E.g., one may input htlgmm result as beta_initial for more rounds to refine the final estimation.
#' The default is NULL, and main study is used for initial estimation according to 'initial_with_type'.
#' @param hat_thetaA If A is not NULL, one can provide hat_thetaA as the input. If 'hat_thetaA = NULL', we estimate hat_thetaA with glm by main study.
#' @param V_thetaA If A is not NULL, one can provide V_thetaA as the input. If 'V_thetaA = NULL', we estimate V_thetaA with glm by main study.
#' @param use_offset Whether to use offset regarding the external model estimated coefficient. The default is TRUE.
#' @param V_thetaA_sandwich Whether to apply sandwich formula to compute the variance-covariance matrix if hat_thetaA.The default is TRUE.
#' @param remove_penalty_Z Not penalize Z if it is TRUE. The default is FALSE.
#' @param remove_penalty_W Not penalize W if it is TRUE. The default is FALSE.
#' @param inference Whether to do inference without penalty or post-selection inference with adaptive lasso penalty. The default is TRUE.
#' @param fix_lambda Without cross validation, fix the lambda. The default is NULL.
#' @param lambda_list Customize the input lambda list for validation. The default is NULL to generate lambda list according to glmnet.
#' @param fix_ratio The fixed ratio for two-lambda strategy. The ratio is multiplied for Z features. The default is NULL. If it is NULL, select the best ratio via cross validation or holdout validation.
#' @param gamma_adaptivelasso The gamma for adaptive lasso. Select from c(1/2,1,2). The default is 1/2.
#' @param use_sparseC Whether to use approximate version of weighting matrix C.
#' If approximation, use the diagonal of inverse of C(inv_C) to approximate the inv_C. The default is TRUE.
#' When main study sample size is limited, use_sparseC = TRUE is recommended.
#' When main study sample size is large enough, use_sparseC = FALSE is recommended.
#' @param seed.use The seed for  97.
#'
#' @return \itemize{
#'  \item{beta:} The target coefficient estimation, the features will go in the order of (A,Z,W).
#'  \item{lambda_list:} The lambda list for cross validation.
#'  \item{ratio_list:} The ratio list for validation (cross validation or holdout validation).
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
#' @export
#'
#'

htlgmm<-function(
        y,Z,W=NULL,
        study_info=NULL,
        A="default",
        penalty_type = "none",
        family = "gaussian",
        initial_with_type = "glm",
        beta_initial = NULL,
        hat_thetaA = NULL,
        V_thetaA = NULL,
        use_offset = TRUE,
        V_thetaA_sandwich = TRUE,
        remove_penalty_Z = FALSE,
        remove_penalty_W = FALSE,
        inference = TRUE,
        fix_lambda = NULL,
        lambda_list = NULL,
        fix_ratio = NULL,
        gamma_adaptivelasso = 1/2,
        use_sparseC = TRUE,
        seed.use = 97
){

    if(!family %in% c("gaussian","binomial")){
        stop("Select family from c('gaussian','binomial')")
    }

    if(is.null(dim(A)[1])){
    if(length(A)==1){
    if(A=='default'){if(family == "gaussian"){A=NULL}else{A=1}}else{
        if(!is.null(A)){
            if(A!=1){warnings("If A is not a matrix, A should be selected from c('default',NULL,1).")}}
    }}else{
        warnings("If A is not selected from c('default',NULL,1), A must be a matrix.")
    }}
    use_cv = FALSE
    nfolds = 10
    tune_ratio = FALSE
    ratio_list = NULL
    type_measure = "default"
    res<-htlgmm.default(y,Z,W,study_info,A,penalty_type,
                        family,initial_with_type,beta_initial,
                        hat_thetaA,V_thetaA,use_offset,
                        V_thetaA_sandwich,remove_penalty_Z,
                        remove_penalty_W,inference,use_cv,
                        type_measure,nfolds,fix_lambda,
                        lambda_list,tune_ratio,fix_ratio,
                        ratio_list,gamma_adaptivelasso,
                        use_sparseC,seed.use)
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
#' @param W The unmatched features only in main study, the default is NULL.
#' @param study_info The trained model from external study, including estimate coefficients, estimated variance-covariance matrix and sample size.
#' The 'study_info' is in the format of list. The first item is 'Coeff', the second iterm is 'Covariance', and the third item is 'Sample_size'.
#' @param A The covariates for study-specific adjustment. The default is 'default', which is 'NULL' for 'gaussian' family, '1' for 'binomial' family.
#' Other than c('default',NULL,1), A must be a matrix whose dimension is the same as the sample dimension or Z and W.
#' For continuous variable, we suggest scaling the features Z, W to eliminate intercept term.  If 'A = NULL', there is no intercept term included.
#' For binary variable, we use intercept term by 'A=1' to adjust for different binary trait ratios in main and external studies.
#' If there is only intercept term in A, we use 'A=1'.
#' A are the features working for adjustment in reduced model, but A is not summarized in summary statistics(input:study_info).
#' @param penalty_type The penalty type for htlgmm, chosen from c("none","lasso","adaptivelasso","ridge"). The default is "adaptivelasso".
#' If 'penalty_type = 'none' ', we use without penalty. (For continuous y, we use ordinary least square, and for binary y, we use logistic regression without penalty.)
#' @param family The family is chosen from c("gaussian","binomial"). Linear regression for "gaussian" and logistic regression for "binomial".
#' @param initial_with_type Get initial estimation for beta using main study data only
#' by cross validation using penalty regression, chosen from c("ridge","lasso","glm"). The default is "ridge". If penalty_type = 'glm',
#' for continuous y, we use ordinary least square, and for binary y, we use logistic regression without penalty.)
#' @param beta_initial The initial estimation for beta if a consistent estimator is available.
#' E.g., one may input htlgmm result as beta_initial for more rounds to refine the final estimation.
#' The default is NULL, and main study is used for initial estimation according to 'initial_with_type'.
#' @param hat_thetaA If A is not NULL, one can provide hat_thetaA as the input. If 'hat_thetaA = NULL', we estimate hat_thetaA with glm by main study.
#' @param V_thetaA If A is not NULL, one can provide V_thetaA as the input. If 'V_thetaA = NULL', we estimate V_thetaA with glm by main study.
#' @param use_offset Whether to use offset regarding the external model estimated coefficient. The default is TRUE.
#' @param V_thetaA_sandwich Whether to apply sandwich formula to compute the variance-covariance matrix if hat_thetaA.The default is TRUE.
#' @param remove_penalty_Z Not penalize Z if it is TRUE. The default is FALSE.
#' @param remove_penalty_W Not penalize W if it is TRUE. The default is FALSE.
#' @param inference Whether to do inference without penalty or post-selection inference with adaptive lasso penalty. The default is TRUE.
#' @param use_cv Whether to use cross validation to determine the best lambda (or ratio).
#' @param type_measure Select from c("default", "mse", "deviance", "auc"). Default is mse(liner), deviance(logistic). 'auc' is another choice for binary y.
#' @param nfolds The fold number for cross validation. Only work for use_cv = TRUE.The default is 10.
#' @param fix_lambda Without cross validation, fix the lambda. The default is NULL.
#' @param lambda_list Customize the input lambda list for validation. The default is NULL to generate lambda list according to glmnet.
#' @param tune_ratio Whether to use two-lambda stratgey. The default is TRUE.
#' @param fix_ratio The fixed ratio for two-lambda strategy. The ratio is multiplied for Z features. The default is NULL. If it is NULL, select the best ratio via cross validation or holdout validation.
#' @param ratio_list The ratio list if it is preset. The default is NULL and ratio list will be generated.
#' @param gamma_adaptivelasso The gamma for adaptive lasso. Select from c(1/2,1,2). The default is 1/2.
#' @param use_sparseC Whether to use approximate version of weighting matrix C.
#' If approximation, use the diagonal of inverse of C(inv_C) to approximate the inv_C. The default is TRUE.
#' When main study sample size is limited, use_sparseC = TRUE is recommended.
#' When main study sample size is large enough, use_sparseC = FALSE is recommended.
#' @param seed.use The seed for  97.
#'
#' @return \itemize{
#'  \item{beta:} The target coefficient estimation, the features will go in the order of (A,Z,W).
#'  \item{lambda_list:} The lambda list for cross validation.
#'  \item{ratio_list:} The ratio list for validation (cross validation or holdout validation).
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
#'  \item{variance:} The variances for features with glm inference, for selected features with post-selection inference.
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
#' @export
#'
#'
cv.htlgmm<-function(
        y,Z,W=NULL,
        study_info=NULL,
        A="default",
        penalty_type = "adaptivelasso",
        family = "gaussian",
        initial_with_type = "ridge",
        beta_initial = NULL,
        hat_thetaA = NULL,
        V_thetaA = NULL,
        use_offset = TRUE,
        V_thetaA_sandwich = TRUE,
        remove_penalty_Z = FALSE,
        remove_penalty_W = FALSE,
        inference = TRUE,
        use_cv = TRUE,
        type_measure = "default",
        nfolds = 10,
        fix_lambda = NULL,
        lambda_list = NULL,
        tune_ratio = FALSE,
        fix_ratio = NULL,
        ratio_list = NULL,
        gamma_adaptivelasso = 1/2,
        use_sparseC = TRUE,
        seed.use = 97
){
    if(!family %in% c("gaussian","binomial")){
        stop("Select family from c('gaussian','binomial')")
    }

    if(is.null(dim(A)[1])){
    if(length(A)==1){
    if(A=='default'){if(family == "gaussian"){A=NULL}else{A=1}}else{
        if(!is.null(A)){
            if(A!=1){warnings("If A is not a matrix, A should be selected from c('default',NULL,1).")}}
    }}else{
        warnings("If A is not selected from c('default',NULL,1), A must be a matrix.")
    }}

    res<-htlgmm.default(y,Z,W,study_info,A,penalty_type,
                        family,initial_with_type,beta_initial,
                        hat_thetaA,V_thetaA,use_offset,
                        V_thetaA_sandwich,remove_penalty_Z,
                        remove_penalty_W,inference,use_cv,
                        type_measure,nfolds,fix_lambda,
                        lambda_list,tune_ratio,fix_ratio,
                        ratio_list,gamma_adaptivelasso,
                        use_sparseC,seed.use,pseudo_Xy_gaussian,
                        pseudo_Xy_binomial,Delta_opt)

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
#' @param study_info The trained model from external study, including estimate coefficients, estimated variance-covariance matrix and sample size.
#' The 'study_info' is in the format of list. The first item is 'Coeff', the second iterm is 'Covariance', and the third item is 'Sample_size'.
#' @param A The covariates for study-specific adjustment. The default is 'default', which is 'NULL' for 'gaussian' family, '1' for 'binomial' family.
#' Other than c('default',NULL,1), A must be a matrix whose dimension is the same as the sample dimension or Z and W.
#' For continuous variable, we suggest scaling the features Z, W to eliminate intercept term.  If 'A = NULL', there is no intercept term included.
#' For binary variable, we use intercept term by 'A=1' to adjust for different binary trait ratios in main and external studies.
#' If there is only intercept term in A, we use 'A=1'.
#' A are the features working for adjustment in reduced model, but A is not summarized in summary statistics(input:study_info).
#' @param family The family is chosen from c("gaussian","binomial"). Linear regression for "gaussian" and logistic regression for "binomial".
#' @param AW_betaAW Default is NULL.
#' @param A_thetaA Default is NULL.
#' @param V_thetaA Default is NULL.
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
#'
#'
#' @import glmnet
#' @import stats
#' @importFrom corpcor fast.svd
#' @importFrom magic adiag
#' @importFrom MASS ginv
#' @importFrom speedglm speedglm speedlm
#' @export
#'
gwas.htlgmm<-function(
        y,Z,W=NULL,
        study_info=NULL,
        A='default',
        family = "gaussian",
        AW_betaAW=NULL,
        A_thetaA=NULL,
        V_thetaA=NULL,
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
                    if(A!=1){warnings("If A is not a matrix, A should be selected from c('default',NULL,1).")}}
            }}else{
                warnings("If A is not selected from c('default',NULL,1), A must be a matrix.")
            }}

    res<-htlgmm.gwas.default(y,Z,W,study_info,A,family,
                             AW_betaAW,A_thetaA,V_thetaA,
                             output_SNP_only,seed.use,
                             verbose,output_tmp)
    return(res)
}

#' Fine-mapping for htlgmm.
#'
#' htlgmm.fm fits a generalized linear model via penalized generalized method of moments,
#' i.e. Heterogeneous Transfer Learning via Generalized Method of Moments,
#' with special focus on fine-mapping, where fm stands for fine-mapping.
#' The input requires main study and GWAS external studies.
#' The htlgmm.fm support cross-validation.
#'
#'
#' @details Fine-mapping for htlgmm.
#'
#' @param y The variable of interest, which can be continuous or binary.
#' @param Z The overlapping features in both main and external studies.
#' @param W The unmatched features only in main study, the default is NULL.
#' @param study_info The trained model from external study, including estimate coefficients, estimated variance-covariance matrix and sample size.
#' The 'study_info' is in the format of list. The first item is 'Coeff', the second iterm is 'Covariance', and the third item is 'Sample_size'.
#' @param A The covariates for study-specific adjustment. The default is 'default', which is 'NULL' for 'gaussian' family, '1' for 'binomial' family.
#' Other than c('default',NULL,1), A must be a matrix whose dimension is the same as the sample dimension or Z and W.
#' For continuous variable, we suggest scaling the features Z, W to eliminate intercept term.  If 'A = NULL', there is no intercept term included.
#' For binary variable, we use intercept term by 'A=1' to adjust for different binary trait ratios in main and external studies.
#' If there is only intercept term in A, we use 'A=1'.
#' A are the features working for adjustment in reduced model, but A is not summarized in summary statistics(input:study_info).
#' @param penalty_type The penalty type for htlgmm, chosen from c("none","lasso","adaptivelasso","ridge"). The default is "adaptivelasso".
#' If 'penalty_type = 'none' ', we use without penalty. (For continuous y, we use ordinary least square, and for binary y, we use logistic regression without penalty.)
#' @param family The family is chosen from c("gaussian","binomial"). Linear regression for "gaussian" and logistic regression for "binomial".
#' @param initial_with_type Get initial estimation for beta using main study data only
#' by cross validation using penalty regression, chosen from c("ridge","lasso","glm"). The default is "ridge". If penalty_type = 'glm',
#' for continuous y, we use ordinary least square, and for binary y, we use logistic regression without penalty.)
#' @param beta_initial The initial estimation for beta if a consistent estimator is available.
#' E.g., one may input htlgmm result as beta_initial for more rounds to refine the final estimation.
#' The default is NULL, and main study is used for initial estimation according to 'initial_with_type'.
#' @param hat_thetaA If A is not NULL, one can provide hat_thetaA as the input. If 'hat_thetaA = NULL', we estimate hat_thetaA with glm by main study.
#' @param V_thetaA If A is not NULL, one can provide V_thetaA as the input. If 'V_thetaA = NULL', we estimate V_thetaA with glm by main study.
#' @param use_offset Whether to use offset regarding the external model estimated coefficient. The default is TRUE.
#' @param V_thetaA_sandwich Whether to apply sandwich formula to compute the variance-covariance matrix if hat_thetaA.The default is TRUE.
#' @param remove_penalty_Z Not penalize Z if it is TRUE. The default is FALSE.
#' @param remove_penalty_W Not penalize W if it is TRUE. The default is FALSE.
#' @param inference Whether to do inference without penalty or post-selection inference with adaptive lasso penalty. The default is TRUE.
#' @param use_cv Whether to use cross validation to determine the best lambda (or ratio).
#' @param type_measure Select from c("default", "mse", "deviance", "auc"). Default is mse(liner), deviance(logistic). 'auc' is another choice for binary y.
#' @param nfolds The fold number for cross validation. Only work for use_cv = TRUE.The default is 10.
#' @param fix_lambda Without cross validation, fix the lambda. The default is NULL.
#' @param lambda_list Customize the input lambda list for validation. The default is NULL to generate lambda list according to glmnet.
#' @param tune_ratio Whether to use two-lambda stratgey. The default is TRUE.
#' @param fix_ratio The fixed ratio for two-lambda strategy. The ratio is multiplied for Z features. The default is NULL. If it is NULL, select the best ratio via cross validation or holdout validation.
#' @param ratio_list The ratio list if it is preset. The default is NULL and ratio list will be generated.
#' @param gamma_adaptivelasso The gamma for adaptive lasso. Select from c(1/2,1,2). The default is 1/2.
#' @param use_sparseC Whether to use approximate version of weighting matrix C.
#' If approximation, use the diagonal of inverse of C(inv_C) to approximate the inv_C. The default is TRUE.
#' When main study sample size is limited, use_sparseC = TRUE is recommended.
#' When main study sample size is large enough, use_sparseC = FALSE is recommended.
#' @param seed.use The seed for  97.
#'
#' @return \itemize{
#'  \item{beta:} The target coefficient estimation, the features will go in the order of (A,Z,W).
#'  \item{lambda_list:} The lambda list for cross validation.
#'  \item{ratio_list:} The ratio list for validation (cross validation or holdout validation).
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
#'  \item{variance:} The variances for features with glm inference, for selected features with post-selection inference.
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
#' @importFrom Rfast colmeans
#' @export
#'
#'
fm.htlgmm<-function(
        y,Z,W=NULL,
        study_info=NULL,
        A="default",
        penalty_type = "adaptivelasso",
        family = "gaussian",
        initial_with_type = "ridge",
        beta_initial = NULL,
        hat_thetaA = NULL,
        V_thetaA = NULL,
        use_offset = TRUE,
        V_thetaA_sandwich = TRUE,
        remove_penalty_Z = FALSE,
        remove_penalty_W = FALSE,
        inference = TRUE,
        use_cv = TRUE,
        type_measure = "default",
        nfolds = 10,
        fix_lambda = NULL,
        lambda_list = NULL,
        tune_ratio = FALSE,
        fix_ratio = NULL,
        ratio_list = NULL,
        gamma_adaptivelasso = 1/2,
        use_sparseC = TRUE,
        seed.use = 97
){
    if(!family %in% c("gaussian","binomial")){
        stop("Select family from c('gaussian','binomial')")
    }

    if(is.null(dim(A)[1])){
        if(length(A)==1){
            if(A=='default'){if(family == "gaussian"){A=NULL}else{A=1}}else{
                if(!is.null(A)){
                    if(A!=1){warnings("If A is not a matrix, A should be selected from c('default',NULL,1).")}}
            }}else{
                warnings("If A is not selected from c('default',NULL,1), A must be a matrix.")
            }}
    if(length(study_info)!=ncol(Z)){
        stop("When using htlgmm.finemapping, input Z as a matrix with size of sample*SNP, and study_info as a list of summary statistics. The columns of Z need to match the study_info.")
    }
    res<-htlgmm.default(y,Z,W,study_info,A,penalty_type,
                        family,initial_with_type,beta_initial,
                        hat_thetaA,V_thetaA,use_offset,
                        V_thetaA_sandwich,remove_penalty_Z,
                        remove_penalty_W,inference,use_cv,
                        type_measure,nfolds,fix_lambda,
                        lambda_list,tune_ratio,fix_ratio,
                        ratio_list,gamma_adaptivelasso,
                        use_sparseC,seed.use,pseudo_Xy_gaussian_finemapping,
                        pseudo_Xy_binomial_finemapping,Delta_opt_finemapping)
    return(res)
}
