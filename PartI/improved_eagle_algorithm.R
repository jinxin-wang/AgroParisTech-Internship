require('MM4LMM')
library("data.table")
library("dplyr")

# The name of functions follow such style : do.something(), do.some.vector(), do.some.Matrix()\
# The name of variables follow such style : some_data, some_vector, Some_Matrix

#### this function is only able to read dat and csv file
read.Geno <- function(geno_fname) {
  Geno <- as.matrix(fread(file = geno_fname, header = FALSE, sep = " ")) # - 1
  return(Geno)
}

read.GenoLable <- function(mapex_fname) {
  Mapex <- fread(mapex_fname, sep = " ", header = TRUE)
  return(Mapex)
}

read.Pheno <- function(pheno_fname) {
  Pheno <- read.csv(pheno_fname, sep = " ", header = TRUE)
  return(Pheno)
}

read.Zmat <- function(Zmat_fname) {
  Zmat <- as.matrix(fread(Zmat_fname, header = FALSE, sep = " "))
  return(Zmat)
}

write.Geno <- function(Geno_mat, geno_fname) {
  write.table(Geno_df, file = geno_fname, col.names = FALSE, row.names =  FALSE, sep = " ") # + 1
}

write_GenoLable <- function(Lable_df, lable_fname) {
  write.table(Lable_df, file = lable_fname, row.names =  FALSE, sep = " ", quote = FALSE)
}

write.Pheno <- function(Pheno_df, pheno_fname) {
  write.table(Pheno_df, file = pheno_fname, row.names = FALSE, sep = " ", quote = FALSE)
}

write.Pheno <- function(Zmat, zmat_fname) {
  write.table(Zmat, zmat_fname, row.names =  FALSE, col.names = FALSE, sep = " ", quote = FALSE)
}
### 

eagle.wo.Z <- function(Y, M, Cofactor, gamma = 0.5, max_itr_num = 5, epsilon = 1e-14) {
  # AA = 0, AB = 1, BB = 2
  # set -1, 0, 1 to AA, AB, BB 
  M_s <- as.matrix(M) - 1
  total_SNPs_num <- ncol(M)
  selected_locus <- c()
  selected_locus_names <- c()
  bic_history <- c()
  extBIC_history <- c() 
  vars_param_history <- list()
  loglik_history <- c()
  vg_history <- c()
  ve_history <- c()
  beta_hat_param_history <- list()
  snps_names <- colnames(M)
  score_names <- snps_names
  
  loci_name <- NULL
  Init <- NULL
  
  for (itr_num in 1:max_itr_num) {
    
    message(sprintf("\nIteration: %d", itr_num))
    
    if (length(selected_locus)) message(sprintf("new selected loci: %s, index: %d", loci_name, loci_idx))
    
    #### #### #### ####    STEP 1 model building    #### #### #### ####
    
    # rebuild X and M_s
    X <- as.matrix(cbind(Cofactor, M[, selected_locus] - 1)) 
    if (!is.null(selected_locus)) M_s <- as.matrix(M[, -selected_locus]) - 1
    
    # Kinship Matrix: 
    # Eagle Kinship Matrix : 
    # Attention, Eagle confuse MMt with kinship matrix
    MMt <- tcrossprod(M_s)
    K <- MMt/max(MMt) + diag(0.95, nrow(MMt))
    
    res <- MMEst(Y = Y, Cofactor = X, Method = "ML", 
                 VarList = list(Additive = K , Error = diag(1, length(Y))), 
                 Init = Init, MaxIter = 4000) #, CritLogLik = 0.0001)
    
    # Note: vg, ve, loglik are very slightly different from Eagle's
    vg = res$NullModel$Sigma2["Additive"]
    ve = res$NullModel$Sigma2["Error"]
    # To calculate log likelihood, MM4LMM use gradient descent, meanwhile Eagle use grid. 
    loglik <- res$NullModel$`LogLik (ML)` 
    
    message(sprintf("log likelihood: %f, var G: %f, var E: %f", loglik, vg, ve))
    
    L = ncol(M_s)
    Init = c(vg,ve)
    #### #### #### ####    STEP 2 model evaluation   #### #### #### #### 
    # Degree of freedom of Eagle's extend BIC equals to 1
    bic <- -2 * loglik + (ncol(X)+1) * log(length(Y))
    ext_bic <- bic + 2 * gamma * lchoose(total_SNPs_num, length(selected_locus))
    
    message(sprintf("BIC: %f, extend BIC: %f", bic, ext_bic))
    
    bic_history <- c(bic_history, bic)
    extBIC_history <- c(extBIC_history, ext_bic)
    vg_history <- c(vg_history, vg)
    ve_history <- c(ve_history, ve)
    loglik_history <- c(loglik_history, loglik)
    
    #### Termination Conditions: 
    if (!is.null(extBIC_history)) {
      if (ext_bic > min(extBIC_history)) {
        message("\n[Termination Message] The extend BIC of the new model cannot beat the previous one.")
        break
      }
    }
    
    #### #### #### ####    STEP 3 model selection   #### #### #### ####  
    
    # Refer to Tristan, Fabien, Alain 2022, page 12 
    Sigma_gamma <- diag(ve, nrow = nrow(K)) + vg * K
    Sigma_gamma_1 <- chol2inv(chol(Sigma_gamma))
    P_gamma <- Sigma_gamma_1 - crossprod(Sigma_gamma_1, X) %>% tcrossprod(.,solve(crossprod(X, Sigma_gamma_1) %>% tcrossprod(.,t(X)))) %>% tcrossprod(.,X) %>% tcrossprod(.,t(Sigma_gamma_1))
    beta_hat <- crossprod(M_s, crossprod(P_gamma, Y))
    
    var_Y <- Sigma_gamma
    P <- crossprod(t(P_gamma), tcrossprod(var_Y,P_gamma))
    diag_var_beta_hat <- rowSums(crossprod(M_s, P) * t(M_s))

    # calculate snps scores
    # bug: when an diagonal element of variance of the beta hat is much 
    #      too close to 0 and not equal to 0, the division becomes to na. 
    # indx <- which(diag_var_beta_hat != 0)
    # Fix: less than 1e-14 will be ignored as 0
    indx <- which(diag_var_beta_hat > epsilon)
    tsq <- beta_hat[indx]^2/diag_var_beta_hat[indx]
    names(tsq) <- seq(1, length(beta_hat))[indx]
    indx <- which(tsq == max(tsq, na.rm = TRUE))
    
    message(paste("the highest score snps local index : ", toString(indx)))
    
    # deal with the case of more than one snps reach the best score
    midpoint <- 1
    if (length(indx) > 2) midpoint <- trunc(length(indx)/2) + 1
    indx <- indx[midpoint]
    orig_indx <- seq(1, L)
    loci_idx <- orig_indx[as.numeric(names(tsq))[indx]]
    
    # snps name and global index
    if (!is.null(selected_locus)) score_names <- snps_names[-selected_locus]
    loci_name <- score_names[loci_idx]
    loci_idx <- which(snps_names==loci_name)
    
    selected_locus <- c(selected_locus, loci_idx)
    selected_locus_names <- c(selected_locus_names, loci_name)
    beta_hat_param_history[[itr_num]] <- list(beta_hat = beta_hat, diag_var_beta_hat = diag_var_beta_hat)
  }
  
  return(data.frame(
    "SNP" = c("NULL", selected_locus_names),
    "Index" = c(NA, selected_locus), 
    "loglik" = c(loglik_history),
    "vg" = c(vg_history), 
    "ve" = c(ve_history),
    "BIC" = c(bic_history),
    "extend BIC" = c(extBIC_history)))
}


eagle.w.Z <- function(Y, M, Cofactor, Z, gamma = 0.5, max_itr_num = 5, epsilon = 1e-14) {
  
  M_s <- as.matrix(M) - 1
  total_SNPs_num <- ncol(M)
  selected_locus <- c()
  selected_locus_names <- c()
  bic_history <- c()
  extBIC_history <- c() 
  vars_param_history <- list()
  loglik_history <- c()
  vg_history <- c()
  ve_history <- c()
  beta_hat_param_history <- list()
  snps_names <- colnames(M)
  score_names <- snps_names
  
  loci_name <- NULL
  Init <- NULL
  
  for (itr_num in 1:max_itr_num) {
    
    message(sprintf("\nIteration: %d", itr_num))
    
    if (length(selected_locus)) message(sprintf("new selected loci: %s, index: %d", loci_name, loci_idx))
    
    #### #### #### ####    STEP 1 model building    #### #### #### ####
    
    # rebuild X and M_s
    if (!is.null(selected_locus)) {
      X <- as.matrix(cbind(Cofactor, crossprod(t(Z), M[, selected_locus] - 1)))
    } else X <- Cofactor
    
    if (!is.null(selected_locus)) M_s <- as.matrix(M[, -selected_locus]) - 1
    
    # Eagle Kinship Matrix : 
    # Attention, Eagle confuse MMt with kinship matrix
    MMt <- tcrossprod(M_s)
    K <- MMt/max(MMt) + diag(0.95, nrow(MMt))
    
    res <- MMEst(Y = c(Y), Cofactor = X, Method = "ML",
                 VarList = list(Additive = K , Error = diag(1, length(Y))),
                 ZList = list(Additive=Z, Error = diag(1, length(Y))),
                 Init = Init, MaxIter = 400, CritLogLik = 0.0001)
    
    # Note: vg, ve, loglik are very slightly different from Eagle's
    vg <- res$NullModel$Sigma2["Additive"]
    ve <- res$NullModel$Sigma2["Error"]
    # To calculate log likelihood, MM4LMM use gradient descent, meanwhile Eagle use grid.   
    loglik <- res$NullModel$`LogLik (ML)` 
    
    message(sprintf("log likelihood: %f, var G: %f, var E: %f", loglik, vg, ve))
    
    L = ncol(M_s)
    Init = c(vg,ve)
    
    #### #### #### ####    STEP 2 model evaluation   #### #### #### #### 

    # Degree of freedom of Eagle's extend BIC equals to 1
    bic <- -2 * loglik + (ncol(X)+1) * log(length(Y))
    ext_bic <- bic + 2 * gamma * lchoose(total_SNPs_num, length(selected_locus))
    
    message(sprintf("BIC: %f, extend BIC: %f", bic, ext_bic))
    
    bic_history <- c(bic_history, bic)
    extBIC_history <- c(extBIC_history, ext_bic)
    vg_history <- c(vg_history, vg)
    ve_history <- c(ve_history, ve)
    loglik_history <- c(loglik_history, loglik)
    
    #### Termination Conditions: 
    if (!is.null(extBIC_history)) {
      if (ext_bic > min(extBIC_history)) {
        message("\n[Termination Message] The extend BIC of new model cannot beat the previous one.")
        break
      }
    }
    
    # Refer to Tristan, Fabien, Alain 2022, page 12 
    K <- crossprod(t(Z), tcrossprod(K,Z))
    Sigma_gamma <- diag(ve, nrow = nrow(K)) + vg * K
    Sigma_gamma_1 <- chol2inv(chol(Sigma_gamma))
    
    P_gamma <- Sigma_gamma_1 - crossprod(Sigma_gamma_1, X) %>% tcrossprod(.,solve(crossprod(X, Sigma_gamma_1) %>% tcrossprod(.,t(X)))) %>% tcrossprod(.,X) %>% tcrossprod(.,t(Sigma_gamma_1))
    beta_hat <- crossprod(M_s, crossprod(Z,crossprod(P_gamma, Y)))
    
    P <- crossprod(Z, P_gamma) %>% tcrossprod(., Sigma_gamma) %>% tcrossprod(., P_gamma) %>% tcrossprod(.,t(Z))
    diag_var_beta_hat <- rowSums(crossprod(M_s, P) * t(M_s))
    
    ####################################################################
    
    # calculate snps scores
    indx <- which(diag_var_beta_hat != 0)
    tsq <- beta_hat[indx]^2/diag_var_beta_hat[indx]
    names(tsq) <- seq(1, length(beta_hat))[indx]
    indx <- which(tsq == max(tsq, na.rm = TRUE))
    
    message(paste("the highest score snps local index : ", toString(indx)))
    
    # deal with the case of more than one snps reach the best score
    midpoint <- 1
    if (length(indx) > 2) midpoint <- trunc(length(indx)/2) + 1
    indx <- indx[midpoint]
    orig_indx <- seq(1, L)
    loci_idx <- orig_indx[as.numeric(names(tsq))[indx]]
    
    # snps name and global index
    if (!is.null(selected_locus)) score_names <- snps_names[-selected_locus]
    loci_name <- score_names[loci_idx]
    loci_idx <- which(snps_names==loci_name)
    
    selected_locus <- c(selected_locus, loci_idx)
    selected_locus_names <- c(selected_locus_names, loci_name)
    beta_hat_param_history[[itr_num]] <- list(beta_hat = beta_hat, diag_var_beta_hat = diag_var_beta_hat)
  }
  
  return(data.frame(
    "SNP" = c("NULL", selected_locus_names),
    "Index" = c(NA, selected_locus), 
    "loglik" = c(loglik_history),
    "vg" = c(vg_history), 
    "ve" = c(ve_history),
    "BIC" = c(bic_history),
    "extend BIC" = c(extBIC_history)))
}

do.eagle <- function(Y, M, Z = NULL, gamma = 0.5, max_itr_num = 5, epsilon = 1e-14) {
  
  if(is.null(dim(Y))) {
    intercept <- matrix(data = 1, nrow = length(Y), ncol = 1)
  } else intercept <- matrix(data = 1, nrow = nrow(Y), ncol = 1)
  colnames(intercept) <- c("intercept")
  
  if (is.null(Z)) {
    res <- eagle.wo.Z(Y = Y, M = M, Cofactor = intercept)
  } else {
    res <- eagle.w.Z(Y = Y, M = M, Cofactor = intercept, Z = Z)
  }
  return(res)
}

