library(MM4LMM)

EMG <- function(Y, G1, G2, Z1 = NULL, Z2 = NULL, genome_names = c("G1", "G2"), 
                max_itr_num = 30, gamma = 0.1, Method = "ML") {
  Pheno <- as.numeric(Y)
  Geno1 <- as.matrix(G1) - 1
  Geno2 <- as.matrix(G2) - 1
  
  X <- matrix(data = 1, nrow = length(Pheno), ncol = 1)
  colnames(X) <- c("intercept")
  
  Init  <- NULL
  total_SNPs_num <- ncol(Geno1) + ncol(Geno2)
  
  selected_locus <- list("G1" = c(), "G2" = c(), "order" = c())
  selected_locus_names <- list("G1" = c(), "G2" = c())
  bic_history <- c()
  extBIC_history <- c() 
  # vars_param_history <- list()
  loglik <- 0
  loglik_history <- c()
  vg_history <- list("G1"=c(), "G2"=c())
  ve_history <- c()
  
  G1_snps_names <- colnames(Geno1)
  G2_snps_names <- colnames(Geno2)
  
  G1_score_names <- G1_snps_names
  G2_score_names <- G2_snps_names
  
  loci_name <- NULL
  loci_genome <- NULL # equal to 1 or 2
  
  .build.ZMs <- function(Z,Geno,selected_locus) {
    ZMs <- crossprod(t(Z),Geno[, selected_locus])
    return(ZMs)
  }
  
  .build.K <- function(Geno, selected_locus) {
    if (is.null(selected_locus)) {
      M_s <- Geno
    } else {
      M_s <- Geno[, -selected_locus] 
    }
    MMt <- tcrossprod(M_s)
    # VERSION Eagle
    # K <- MMt/max(MMt) + diag(0.95, nrow(MMt))
    K <- MMt/mean(diag(MMt)) 
    return(list("M_s" = M_s, "K" = K))
  }
  
  .calc.a <- function(Z, M_s, Pheno, P_gamma, Sigma_gamma) {
    a_hat <- crossprod(M_s, crossprod(Z,crossprod(P_gamma, Pheno)))
    P <- crossprod(Z, P_gamma) %>% tcrossprod(., Sigma_gamma) %>% tcrossprod(., P_gamma) %>% tcrossprod(.,t(Z))
    diag_var_a_hat <- rowSums(crossprod(M_s, P) * t(M_s))
    return(list("a_hat" = a_hat, "var_a_hat" = diag_var_a_hat))
  }
  
  .calc.snps.scores <- function(a_hat, diag_var_a_hat, L, snps_names, selected_locus, genome_name) {
    # calculate snps scores
    indx <- which(diag_var_a_hat != 0)
    tsq <- a_hat[indx]^2/diag_var_a_hat[indx]
    names(tsq) <- seq(1, length(a_hat))[indx]
    best_score <- max(tsq, na.rm = TRUE)
    indx <- which(tsq == best_score)
    
    # deal with the case of more than one snps reach the best score
    midpoint <- 1
    if (length(indx) > 2) midpoint <- trunc(length(indx)/2) + 1
    indx <- indx[midpoint]
    orig_indx <- seq(1, L)
    loci_idx <- orig_indx[as.numeric(names(tsq))[indx]]
    
    # snps name and global index
    if (!is.null(selected_locus)) {
      score_names <- snps_names[-selected_locus]
    } else {
      score_names <- snps_names
    }
    
    loci_name <- score_names[loci_idx]
    loci_idx <- which(snps_names==loci_name)
    
    message(sprintf("[%s] the highest score snps index : %d, snp name: %s, snp score: %f", genome_name, loci_idx, loci_name, best_score))
    
    return(list("loci_name" = loci_name, "loci_idx" = loci_idx, "score" = best_score))
  }
  
  tryCatch({
    
    for (itr_num in 1:max_itr_num) {
      
      message(sprintf("\nIteration: %d", itr_num))
      
      if (!is.null(loci_name)) message(sprintf("new selected loci: %s, index: %d on %s", loci_name, loci_idx, genome_names[loci_genome]))
      
      #### #### #### ####    STEP 1 model building    #### #### #### ####
      if (!is.null(loci_genome)) {
        if (loci_genome == 1)  {
          X <- cbind(X, .build.ZMs(Z = Z1, Geno = Geno1, selected_locus = selected_locus$G1[length(selected_locus$G1)]))
        } else {
          X <- cbind(X, .build.ZMs(Z = Z2, Geno = Geno2, selected_locus = selected_locus$G2[length(selected_locus$G2)]))
        }
      }
  
      # Eagle Kinship Matrix : 
      # Attention, Eagle confuse MMt with kinship matrix
      
      if (is.null(loci_genome) || loci_genome == 1) {
        M1M1T <- .build.K(Geno1, selected_locus$G1)
      }
      
      if (is.null(loci_genome) || loci_genome == 2) {
        M2M2T <- .build.K(Geno2, selected_locus$G2)
      }
      
      VL <- list(G1 = M1M1T$K , G2 = M2M2T$K , Error = diag(1,length(Pheno)))
      ZL <- list(G1 = Z1, G2 = Z2, Error = diag(1,length(Pheno)))
      
      if (!is.null(genome_names)) {
        names(VL) <- c(genome_names, "Error")
        names(ZL) <- c(genome_names, "Error")
      }
      
      if (is.null(Init)) {
        Init <- rep(var(Pheno)/3,3)
      }
  
      res <- MMEst(Y = Pheno, Cofactor = X, Method = Method,
                   VarList = VL, ZList = ZL, Init = Init,
                   MaxIter = 400, CritLogLik = 0.0001)
      
      # Note: vg, ve, loglik are very slightly different from Eagle's
      vg1 <- res$NullModel$Sigma2[genome_names[1]]
      vg2 <- res$NullModel$Sigma2[genome_names[2]]
      ve <- res$NullModel$Sigma2["Error"]
      
      Init <- c(vg1, vg2, ve)
      
      # To calculate log likelihood, MM4LMM use gradient descent, meanwhile Eagle use grid.   
      if (!is.null(res$NullModel$`LogLik (ML)`)) {
        loglik <- res$NullModel$`LogLik (ML)`   
      } else {
        loglik <- res$NullModel$`LogLik (Reml)`
      }
      
      
      message(sprintf("log likelihood: %f, var G1: %f, var G2: %f, var E: %f", loglik, vg1, vg2, ve))
      
      #### #### #### ####    STEP 2 model evaluation   #### #### #### #### 
      # Degree of freedom of Eagle's extend BIC equals to 1
      bic <- -2 * loglik + (ncol(X)+1) * log(length(Pheno))
      ext_bic <- bic + 2 * gamma * lchoose(total_SNPs_num, length(selected_locus$order))
      
      message(sprintf("BIC: %f, extend BIC: %f", bic, ext_bic))
      
      bic_history <- c(bic_history, bic)
      extBIC_history <- c(extBIC_history, ext_bic)
      vg1_history <- c(vg_history, vg1)
      vg2_history <- c(vg_history, vg2)
      ve_history <- c(ve_history, ve)
      loglik_history <- c(loglik_history, loglik)
      
      #### Termination Conditions: 
  
      # if (!is.null(extBIC_history)) {
      #   if (ext_bic > min(extBIC_history)) {
      #     message("\n[Termination Message] The score of the new model cannot beat the previous one.")
      # 
      #     itr_num <- itr_num - 1
      #     selected_locus$order <- selected_locus$order[-itr_num]
      #     if (loci_genome == 1) {
      #       selected_locus$G1 <- selected_locus$G1[-itr_num]
      #       selected_locus_names$G1 <- selected_locus_names$G1[-itr_num]
      #     }
      # 
      #     if (loci_genome == 2) {
      #       selected_locus$G2 <- selected_locus$G1[-itr_num]
      #       selected_locus_names$G2 <- selected_locus_names$G2[-itr_num]
      #     }
      # 
      #     break
      #   }
      # }
  
      #### #### #### ####    STEP 3 model selection   #### #### #### ####  
      
      K1 <- crossprod(t(Z1), tcrossprod(M1M1T$K,Z1))
      K2 <- crossprod(t(Z2), tcrossprod(M2M2T$K,Z2))
      
      Sigma_gamma <- diag(ve, nrow = nrow(K1)) + vg1 * K1 + vg2 * K2
      Sigma_gamma_1 <- chol2inv(chol(Sigma_gamma))
      
      # P_gamma <- Sigma_gamma_1 - crossprod(Sigma_gamma_1, X) %>% 
      #                               tcrossprod(.,solve(crossprod(X, Sigma_gamma_1) %>% 
      #                                   tcrossprod(.,t(X)))) %>% 
      #                                       tcrossprod(.,X) %>% 
      #                                         tcrossprod(.,Sigma_gamma_1)
      
      XtS_1 <- crossprod(X, Sigma_gamma_1)
      P_gamma <-  Sigma_gamma_1 - solve(crossprod(X, t(XtS_1))) %>% crossprod(., XtS_1) %>% crossprod(XtS_1, .)
      
      a_vara1 <- .calc.a(Z = Z1, M_s = M1M1T$M_s, Pheno = Pheno, P_gamma = P_gamma, Sigma_gamma = Sigma_gamma)
      G1_score <- .calc.snps.scores(a_hat = a_vara1$a_hat, diag_var_a_hat = a_vara1$var_a_hat, 
                                      L = ncol(Geno1), snps_names = G1_snps_names, 
                                      selected_locus = selected_locus$G1, genome_name = genome_names[1]) 
      
      a_vara2 <- .calc.a(Z = Z2, M_s = M2M2T$M_s, Pheno = Pheno, P_gamma = P_gamma, Sigma_gamma = Sigma_gamma)
      G2_score <- .calc.snps.scores(a_hat = a_vara2$a_hat, diag_var_a_hat = a_vara2$var_a_hat, 
                                      L = ncol(Geno2), snps_names = G2_snps_names, 
                                      selected_locus = selected_locus$G2, genome_name = genome_names[2]) 
      
      if ( G1_score$score > G2_score$score) {
        loci_genome <- 1
        selected_locus$order <- c(selected_locus$order, 1)
        selected_locus$G1 <- c(selected_locus$G1, G1_score$loci_idx)
        selected_locus_names$G1 <- c(selected_locus_names$G1, G1_score$loci_name)
        loci_name <- G1_score$loci_name
        loci_idx  <- G1_score$loci_idx
        message(sprintf("[%s] new selected loci: %s, index: %d", genome_names[1], G1_score$loci_name, G1_score$loci_idx))
      } else {
        loci_genome <- 2
        selected_locus$order <- c(selected_locus$order, 2)
        selected_locus$G2 <- c(selected_locus$G2, G2_score$loci_idx)
        selected_locus_names$G2 <- c(selected_locus_names$G2, G2_score$loci_name)
        loci_name <- G2_score$loci_name
        loci_idx  <- G2_score$loci_idx
        message(sprintf("[%s] new selected loci: %s, index: %d", genome_names[2], G2_score$loci_name, G2_score$loci_idx))
      }
    }
  }, warning = function(warn) {
    print(paste("EMG WARNING MESSAGE: ", warn))
  }, error = function(err) {
    print(paste("EMG ERROR MESSAGE: ", err))
  }, finally = function(f) {
    selected_locus$extBIC <- extBIC_history
    selected_locus$loglik <- loglik_history
    return(selected_locus)
  })
  
  selected_locus$extBIC <- extBIC_history
  selected_locus$loglik <- loglik_history
  return(selected_locus)
}
