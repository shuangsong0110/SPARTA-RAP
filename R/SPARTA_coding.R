###########################################################
#  SPARTA - Step 4.1: Gene-centric CODING analysis  (UKB RAP / DNANexus)
###########################################################
# Run as a DNANexus array job: arrayid (1..195) selects a block of genes on one
# chromosome. Inputs are pulled from the RAP project with `dx download` and
# results pushed back with `dx upload` (project paths set in the User Input block).
# Authenticate by exporting your token first:  export DX_TOKEN=xxxxxxxx
###########################################################
print(commandArgs(TRUE))
args<- commandArgs(TRUE)
start_time<-Sys.time()
arrayid<- as.numeric(args[1])
print('arrayid:')
print(arrayid)
system("pip3 install setuptools-rust")
system("pip3 install --upgrade --ignore-installed pip")
system("pip3 install dxpy")
system(paste0("dx login --token ", Sys.getenv("DX_TOKEN")))
if(T){
  
  #------------------------------------------------------------
  find_root_spa <- function(gfun, t0 = 10, max_expand = 8, ...) {
    for (k in 0:max_expand) {
      lo <- -(t0 * 2^k)
      hi <-  (t0 * 2^k)
      glo <- gfun(lo, ...)
      ghi <- gfun(hi, ...)
      if (is.finite(glo) && is.finite(ghi) && sign(glo) != sign(ghi)) {
        return(uniroot(gfun, lower = lo, upper = hi, ...)$root)
      }
    }
    stop("Root not bracketed even after expanding interval.")
  }
  
  #------------------------------------------------------------
  # 主函数：双正态混合 SPA 校正 p 值
  #------------------------------------------------------------
  
  pv_normmix <- function(mu1, sigma21, mu2, sigma22, pi1,
                         cut, pv.unadj = 0.5,
                         t0 = 4, max_expand = 4) {
    pv.unadj <- pnorm(-abs(cut), sd=sqrt(min(sigma21, sigma22)))
    pv.unadj <- pnorm(-abs(cut), sd=sqrt(sum(sigma21*(1-pi1)+ sigma22*pi1)))
    print(pv.unadj)
    ## 如果未调 SPA（pv.unadj ≥ 0.001）就直接返回
    if (pv.unadj >= 0.01) {
      pv <- pnorm(-abs(cut), sd=sqrt((1-pi1)*sigma21+pi1*sigma22))
      return(list(pv = pv.unadj, pv_unadj = pv.unadj, SPA = FALSE))
    }
    
    ##---- 定义 g(t) = K'(t) - cut
    gfun <- function(tt, ...) {
      K_normmix_prime1(...) - cut
    }
    
    ##---- 求鞍点根（自适应区间）
    t_hat <- find_root_spa(
      gfun,
      t0 = t0,
      max_expand = max_expand,
      mu1 = mu1, sigma21 = sigma21,
      mu2 = mu2, sigma22 = sigma22,
      pi1 = pi1, tt = tt    # 注意: gfun 里 K' 会忽略多余参数
    )
    #print(t_hat)
    ##---- 计算 w, v, 校正 p 值
    K_val  <- K_normmix(mu1, sigma21, mu2, sigma22, pi1, tt = t_hat)
    K2_val <- K_normmix_prime2(mu1, sigma21, mu2, sigma22, pi1, tt = t_hat)
    
    w <- sign(t_hat) * sqrt(2 * (t_hat * cut - K_val))
    v <- t_hat * sqrt(K2_val)
    # phi_w <- dnorm(w)
    # pv <- 1 - pnorm(w) + phi_w * (1/w - 1/v)  
    # pv <- 1 - pnorm(w + (log(v / w)) / w)
    if (abs(w) < 1e-2) {  # 处理 w≈0 的情况
      print(w)
      pv <- 0.49
      spa_flag <- F
    } else {
      adj <- log(v / w) / w
      pv <- pnorm(-w - adj)
    }
    # print(w)
    # print(v)
    ## 如遇数值下溢回退
    if (is.na(pv) || pv == 0) {
      pv <- pv.unadj
      spa_flag <- FALSE
    } else {
      spa_flag <- TRUE
    }
    
    list(pv = pv, pv_unadj = pv.unadj, SPA = spa_flag)
  }
  
  K_normmix <- function(mu1, sigma21, mu2, sigma22,pi1, tt){
    A <- exp(mu1*tt+0.5*sigma21*tt^2)
    B <- exp(mu2*tt+0.5*sigma22*tt^2)
    D <- (1-pi1)*A+pi1*B
    return(sum(log(D)))
  }           
  K_normmix_prime2 <- function(mu1, sigma21,   # μ1, σ1²
                               mu2, sigma22,   # μ2, σ2²
                               pi1,            # π1 ∈ (0,1)
                               tt) {           # t  (可向量化)
    # 指数项
    eta1 <- mu1 * tt + 0.5 * sigma21 * tt^2
    eta2 <- mu2 * tt + 0.5 * sigma22 * tt^2
    A    <- exp(eta1)
    B    <- exp(eta2)
    
    # 公分母
    D <- (1 - pi1) * A + pi1 * B
    
    # 一阶导中用到的 dη/dt
    d1 <- mu1 + sigma21 * tt      # η1'
    d2 <- mu2 + sigma22 * tt      # η2'
    
    # 构造二阶导公式
    num1 <- (1 - pi1) * A * (d1^2 + sigma21)
    num2 <- pi1        * B * (d2^2 + sigma22)
    K1   <- ((1 - pi1) * A * d1 + pi1 * B * d2) / D   # 一阶导数
    K2   <- (num1 + num2) / D - K1^2                  # 二阶导数
    
    return(sum(K2))
  }
  K_normmix_prime1 <- function(mu1, sigma21,   # μ1, σ1²
                               mu2, sigma22,   # μ2, σ2²
                               pi1,            # π1 ∈ (0,1)
                               tt) {           # t  (vector OK)
    ## 1. 计算指数项
    eta1 <- mu1 * tt + 0.5 * sigma21 * tt^2
    eta2 <- mu2 * tt + 0.5 * sigma22 * tt^2
    A    <- exp(eta1)
    B    <- exp(eta2)
    
    ## 2. 公分母
    D <- (1 - pi1) * A + pi1 * B
    
    ## 3. dη/dt
    d1 <- mu1 + sigma21 * tt   # η1'
    d2 <- mu2 + sigma22 * tt   # η2'
    
    ## 4. 一阶导
    K1 <- ((1 - pi1) * A * d1 + pi1 * B * d2) / D
    return(sum(K1))
  }
  K_normmix_prime_adj <- function(mu1, sigma21, mu2, sigma22, pi1, tt, x){
    return(K_normmix_prime1(mu1, sigma21, mu2, sigma22, pi1, tt)-x)
  }
  
  
}

if(T){
  
  STAAR_sp_2norm_SPA <-   function (genotype_sp, MAF = NULL, obj_nullmodel, annotation_phred = NULL,
                                    rare_maf_cutoff = 0.01, rv_num_cutoff = 2, rv_num_cutoff_max = 1e+09,
                                    pi1, resid_ratio)
  {
    pp <- resid_ratio*pi1+1-pi1 ## 117
    
    if (!inherits(genotype_sp, "matrix") && !inherits(genotype_sp,
                                                      "Matrix")) {
      print(class(genotype_sp))
      stop("genotype is not a matrix!")
    }
    if (!inherits(genotype_sp, "sparseMatrix")) {
      genotype_sp <- as(genotype_sp, "dgCMatrix")
    }
    if (dim(genotype_sp)[2] == 1) {
      stop(paste0("Number of rare variant in the set is less than 2!"))
    }
    if (is.null(MAF)) {
      MAF <- Matrix::colMeans(genotype_sp, na.rm = TRUE)/2
    }
    annotation_phred <- as.data.frame(annotation_phred)
    if (dim(annotation_phred)[1] != 0 & (dim(genotype_sp)[2] !=
                                         dim(annotation_phred)[1] | dim(genotype_sp)[2] != length(MAF))) {
      stop(paste0("Dimensions don't match for genotype, MAF and annotation!"))
    }
    RV_label <- as.vector((MAF < rare_maf_cutoff) & (MAF > 0) &
                            (!is.na(MAF)))
    G <- genotype_sp[, RV_label]
    #print(dim(genotype_sp))
    #print(class(genotype_sp))
    #MAC0 <-  colSums(genotype_sp)
    MAC0 <-  colSums(genotype_sp[which(obj_nullmodel$R==1),])
    # print('MAC0')
    # print(summary(MAC0))
    rm(genotype_sp)
    gc()
    annotation_phred <- annotation_phred[RV_label, , drop = FALSE]
    if (sum(RV_label) >= rv_num_cutoff_max) {
      stop(paste0("Number of rare variant in the set is more than ",
                  rv_num_cutoff_max, "!"))
    }
    if (sum(RV_label) >= rv_num_cutoff) {
      MAF <- MAF[RV_label]
      MAC0 <- MAC0[RV_label]
      if (any(is.na(G@x))) {
        G <- na.replace.sp(G, m = 2 * MAF)
      }
      gc()
      annotation_rank <- 1 - 10^(-annotation_phred/10)
      w_1 <- dbeta(MAF, 1, 25)
      w_2 <- dbeta(MAF, 1, 1)
      if (dim(annotation_phred)[2] == 0) {
        w_B <- w_S <- as.matrix(cbind(w_1, w_2))
        w_A <- as.matrix(cbind(w_1^2/dbeta(MAF, 0.5, 0.5)^2,
                               w_2^2/dbeta(MAF, 0.5, 0.5)^2))
      }
      else {
        w_B_1 <- annotation_rank * w_1
        w_B_1 <- cbind(w_1, w_B_1)
        w_B_2 <- annotation_rank * w_2
        w_B_2 <- cbind(w_2, w_B_2)
        w_B <- cbind(w_B_1, w_B_2)
        w_B <- as.matrix(w_B)
        w_S_1 <- sqrt(annotation_rank) * w_1
        w_S_1 <- cbind(w_1, w_S_1)
        w_S_2 <- sqrt(annotation_rank) * w_2
        w_S_2 <- cbind(w_2, w_S_2)
        w_S <- cbind(w_S_1, w_S_2)
        w_S <- as.matrix(w_S)
        w_A_1 <- annotation_rank * w_1^2/dbeta(MAF, 0.5,
                                               0.5)^2
        w_A_1 <- cbind(w_1^2/dbeta(MAF, 0.5, 0.5)^2, w_A_1)
        w_A_2 <- annotation_rank * w_2^2/dbeta(MAF, 0.5,
                                               0.5)^2
        w_A_2 <- cbind(w_2^2/dbeta(MAF, 0.5, 0.5)^2, w_A_2)
        w_A <- cbind(w_A_1, w_A_2)
        w_A <- as.matrix(w_A)
      }
      #newadd <- T
      if(T){
        
        
        if (obj_nullmodel$relatedness) {
          if (!obj_nullmodel$sparse_kins) {
            P <- obj_nullmodel$P
            residuals.phenotype <- obj_nullmodel$scaled.residuals
            pvalues <- STAAR_O_SMMAT(G, P, residuals.phenotype,
                                     weights_B = w_B, weights_S = w_S, weights_A = w_A,
                                     mac = as.integer(round(MAF * 2 * dim(G)[1])))
          }
          else {
            Sigma_i <- obj_nullmodel$Sigma_i
            Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
            cov <- obj_nullmodel$cov
            residuals.phenotype <- obj_nullmodel$scaled.residuals
            pvalues <- STAAR_O_SMMAT_sparse(G, Sigma_i, Sigma_iX,
                                            cov, residuals.phenotype, weights_B = w_B,
                                            weights_S = w_S, weights_A = w_A, mac = as.integer(round(MAF *
                                                                                                       2 * dim(G)[1])))
          }
        }
        else {
          X <- model.matrix(obj_nullmodel)
          working <- obj_nullmodel$weights
          sigma <- sqrt(summary(obj_nullmodel)$dispersion)
          if (obj_nullmodel$family[1] == "binomial") {
            fam <- 1
          }
          else if (obj_nullmodel$family[1] == "gaussian") {
            fam <- 0
          }
          residuals.phenotype <- obj_nullmodel$y - obj_nullmodel$fitted.values
          pvalues <- STAAR_O(G, X, working, sigma, fam, residuals.phenotype,
                             weights_B = w_B, weights_S = w_S, weights_A = w_A,
                             mac = as.integer(round(MAF * 2 * dim(G)[1])))
        }
        
      }
      if(min(pvalues)<1e-3){
        newadd <- T
        print('use SPA')
      }else{
        newadd <- F
        print('no SPA')
      }
      if(newadd){
        Sigma_i <- obj_nullmodel$Sigma_i
        Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
        cov <- obj_nullmodel$cov
        
        residuals.phenotype <- obj_nullmodel$scaled.residuals
        residuals <- residuals.phenotype
        # pvalues <- STAAR_O_SMMAT_sparse(G,Sigma_i,Sigma_iX,cov,residuals.phenotype,
        #                                 weights_B=w_B,weights_S=w_S,weights_A=w_A,
        #                                 mac=as.integer(round(MAF*2*dim(G)[1])))
        weights_B=w_B
        weights_S=w_S
        weights_A=w_A
        mac=as.integer(round(MAF*2*dim(G)[1]))
        mac_thres <- 10
        ##########################################################################
        ## ---------------- 基本尺寸 ----------------
        un <- ncol(G)               # m: variant 数
        wn <- ncol(weights_B)       # w: 权重组数
        res <- numeric(3 * wn)      # 输出
        
        lower <- FALSE; logp <- FALSE
        
        ## ---------------- 协方差矩阵 Cov ----------------
        ## Cov = (G' Σ⁻¹ G) - (G' Σ⁻¹ X)(X' Σ⁻¹ X)⁻¹ (X' Σ⁻¹ G)
        tSigma_iX_G <- Matrix::t(Sigma_iX) %*% G
        Cov <- Matrix::t(Sigma_i %*% G) %*% G -
          Matrix::t(tSigma_iX_G) %*% cov %*% tSigma_iX_G
        ## 取对角向量可加速
        Cov_diag <- Matrix::diag(Cov)
        
        ## ---------------- 预计算 ----------------
        x <- as.numeric(Matrix::t(residuals) %*% G)  # 1 × m
        
        id_veryrare <- which(MAC0 <= mac_thres)
        id_common   <- which(MAC0 >  mac_thres)
        n0 <- length(id_veryrare)
        n1 <- length(id_common)
        
        ## p 值序列（常见等位基因）
        pseq0 <- numeric(un)
        pseq <- numeric(un)
        chi_vals0 <- x^2 / Cov_diag
        for(j in 1:un){
          G_sum <- G[,j]
          kk <- length(which(G_sum!=0))
          G_sum_use <- G_sum[which(G_sum!=0)]
          sigma21 <- 1/pp
          sigma22 <- 1/pp*resid_ratio
          aa <- try(pseq0[j] <- pv_normmix(rep(0,kk),
                                           G_sum_use^2*sigma21,
                                           rep(0,kk),
                                           G_sum_use^2*sigma22,
                                           pi1=rep(pi1,kk),
                                           cut=sqrt(chi_vals0[j])*sqrt(1*sum(G_sum_use^2)),
                                           pv.unadj=1e-5,max_expand = 4)$pv*2)
          if(inherits(aa, "try-error")){
            pseq0[j] <- 0.99
          }
          
        }
        pseq0 <- pmin(pseq0,0.99)
        if (n1 > 0) {
          pseq[1:n1] <- pseq0[id_common]
          
        }
        
        ## ---------------- 主循环：每组权重 ----------------
        for (i in seq_len(wn)) {
          
          ## -------- SKAT 部分 --------
          wS <- weights_S[, i]/sum(weights_S[, i])
          if(newadd){
            ## 这里assume x~N(0,diag(Cov))
            xnew <- qnorm(pseq0/2)*sqrt(Matrix::diag(Cov))
            
            Covw <- (wS %o% wS) * Cov              # 外积再逐元素乘
            Q_skat <- sum((xnew^2) * (wS^2))
            eigenvals <- eigen(Covw, symmetric = TRUE,
                               only.values = TRUE)$values
            eigenvals[eigenvals < 1e-8] <- 0
            
            p_skat <- Saddle(Q_skat, eigenvals)
            p_skat
            ## Davies 可能返回 2（数值不稳定的标记），
            ## 这里跟原代码相同的补救:
            if (identical(p_skat, 2)) {
              c1  <- sum(Matrix::diag(Covw))
              c2  <- sum(Matrix::diag(Covw %*% Covw))
              c4  <- sum(Matrix::diag(Covw %*% Covw %*% Covw %*% Covw))
              Qn  <- (Q_skat - c1) / sqrt(2 * c2)
              l   <- (c2^2) / c4
              p_skat <- pchisq(Qn * sqrt(2 * l) + l, df = l,
                               lower.tail = lower, log.p = logp)
            }
            
          }
          res[i] <- p_skat
          
          ## -------- Burden 部分 --------
          wB <- weights_B[, i]/sum(weights_B[, i])
          Covw_b <- (wB %o% wB) * Cov
          sumw_b <- sum(Covw_b)
          T_burden <- (sum(x * wB))^2 / sumw_b
          res[wn + i] <- pchisq(T_burden, df = 1,lower.tail = lower, log.p = logp)
          if(newadd){
            G_sum <- G%*%wB
            kk <- length(which(G_sum[,1]!=0))
            G_sum_use <- G_sum[which(G_sum[,1]!=0),1]
            #pi1 <- mean(R)
            sigma21 <- 1/pp
            sigma22 <- 1/pp*resid_ratio
            aa <- try(res[wn + i] <- pv_normmix(rep(0,kk),
                                                G_sum_use^2*sigma21,
                                                rep(0,kk),
                                                G_sum_use^2*sigma22,
                                                pi1=rep(pi1,kk),
                                                cut=sqrt(T_burden)*sqrt(1*sum(G_sum_use^2)),
                                                pv.unadj=1e-5,max_expand = 4)$pv*2)
            if(inherits(aa, "try-error")){
              res[wn + i] <- 0.99
            }else{
              res[wn + i] <- pmin(0.99, res[wn + i])
            }
            
          }
          ## -------- ACAT (聚合) --------
          wseq_common <- weights_A[id_common, i]
          if (n0 == 0) {
            ## 没有 very rare
            res[2 * wn + i] <- CCT_pval(pseq[id_common], wseq_common)
          } else {
            ## 先组合 very rare burden p 值
            wA_vr <- weights_A[id_veryrare, i]
            sum0_a <- sum(x[id_veryrare] * wB[id_veryrare])
            denom  <- sum(Covw_b[id_veryrare, id_veryrare])
            T_vr   <- (sum0_a)^2 / denom
            p_vr   <- pchisq(T_vr, df = 1,
                             lower.tail = lower, log.p = logp)
            if(newadd){
              if(length(id_veryrare)==1){
                G_sum <- G[,id_veryrare]*wB[id_veryrare]  
                kk <- length(which(G_sum!=0))
                G_sum_use <- G_sum[which(G_sum!=0)]
              }else{
                G_sum <- G[,id_veryrare]%*%wB[id_veryrare]  
                kk <- length(which(G_sum[,1]!=0))
                G_sum_use <- G_sum[which(G_sum[,1]!=0),1]
              }
              
              
              sigma21 <- 1/pp
              sigma22 <- 1/pp*resid_ratio
              
              aa <- try(p_vr <- pv_normmix(rep(0,kk),
                                           G_sum_use^2*sigma21,
                                           rep(0,kk),
                                           G_sum_use^2*sigma22,
                                           pi1=rep(pi1,kk),
                                           cut=sqrt(T_vr)*sqrt(1*sum(G_sum_use^2)),
                                           pv.unadj=1e-5,max_expand = 4)$pv*2)
              if(inherits(aa, "try-error")){
                p_vr <- 0.99
              }else{
                p_vr <- pmin(0.99, res[wn + i])
              }
              
            }
            p_all  <- c(pseq0[id_common], p_vr)
            w_all  <- c(wseq_common, sum(wA_vr) / n0)
            
            res[2 * wn + i] <- CCT_pval(p_all, w_all)
          }
        }
        
        pvalues <- res
      }
      
      
      num_variant <- sum(RV_label)
      cMAC <- sum(G@x)
      num_annotation <- dim(annotation_phred)[2] + 1
      results_STAAR_O <- CCT(pvalues)
      results_ACAT_O <- CCT(pvalues[c(1, num_annotation + 1,
                                      2 * num_annotation + 1, 3 * num_annotation + 1, 4 *
                                        num_annotation + 1, 5 * num_annotation + 1)])
      pvalues_STAAR_S_1_25 <- CCT(pvalues[1:num_annotation])
      pvalues_STAAR_S_1_1 <- CCT(pvalues[(num_annotation +
                                            1):(2 * num_annotation)])
      pvalues_STAAR_B_1_25 <- CCT(pvalues[(2 * num_annotation +
                                             1):(3 * num_annotation)])
      pvalues_STAAR_B_1_1 <- CCT(pvalues[(3 * num_annotation +
                                            1):(4 * num_annotation)])
      pvalues_STAAR_A_1_25 <- CCT(pvalues[(4 * num_annotation +
                                             1):(5 * num_annotation)])
      pvalues_STAAR_A_1_1 <- CCT(pvalues[(5 * num_annotation +
                                            1):(6 * num_annotation)])
      results_STAAR_S_1_25 <- c(pvalues[1:num_annotation],
                                pvalues_STAAR_S_1_25)
      results_STAAR_S_1_25 <- data.frame(t(results_STAAR_S_1_25))
      results_STAAR_S_1_1 <- c(pvalues[(num_annotation + 1):(2 *
                                                               num_annotation)], pvalues_STAAR_S_1_1)
      results_STAAR_S_1_1 <- data.frame(t(results_STAAR_S_1_1))
      results_STAAR_B_1_25 <- c(pvalues[(2 * num_annotation +
                                           1):(3 * num_annotation)], pvalues_STAAR_B_1_25)
      results_STAAR_B_1_25 <- data.frame(t(results_STAAR_B_1_25))
      results_STAAR_B_1_1 <- c(pvalues[(3 * num_annotation +
                                          1):(4 * num_annotation)], pvalues_STAAR_B_1_1)
      results_STAAR_B_1_1 <- data.frame(t(results_STAAR_B_1_1))
      results_STAAR_A_1_25 <- c(pvalues[(4 * num_annotation +
                                           1):(5 * num_annotation)], pvalues_STAAR_A_1_25)
      results_STAAR_A_1_25 <- data.frame(t(results_STAAR_A_1_25))
      results_STAAR_A_1_1 <- c(pvalues[(5 * num_annotation +
                                          1):(6 * num_annotation)], pvalues_STAAR_A_1_1)
      results_STAAR_A_1_1 <- data.frame(t(results_STAAR_A_1_1))
      if (dim(annotation_phred)[2] == 0) {
        colnames(results_STAAR_S_1_25) <- c("SKAT(1,25)",
                                            "STAAR-S(1,25)")
        colnames(results_STAAR_S_1_1) <- c("SKAT(1,1)", "STAAR-S(1,1)")
        colnames(results_STAAR_B_1_25) <- c("Burden(1,25)",
                                            "STAAR-B(1,25)")
        colnames(results_STAAR_B_1_1) <- c("Burden(1,1)",
                                           "STAAR-B(1,1)")
        colnames(results_STAAR_A_1_25) <- c("ACAT-V(1,25)",
                                            "STAAR-A(1,25)")
        colnames(results_STAAR_A_1_1) <- c("ACAT-V(1,1)",
                                           "STAAR-A(1,1)")
      }
      else {
        colnames(results_STAAR_S_1_25) <- c("SKAT(1,25)",
                                            paste0("SKAT(1,25)-", colnames(annotation_phred)),
                                            "STAAR-S(1,25)")
        colnames(results_STAAR_S_1_1) <- c("SKAT(1,1)", paste0("SKAT(1,1)-",
                                                               colnames(annotation_phred)), "STAAR-S(1,1)")
        colnames(results_STAAR_B_1_25) <- c("Burden(1,25)",
                                            paste0("Burden(1,25)-", colnames(annotation_phred)),
                                            "STAAR-B(1,25)")
        colnames(results_STAAR_B_1_1) <- c("Burden(1,1)",
                                           paste0("Burden(1,1)-", colnames(annotation_phred)),
                                           "STAAR-B(1,1)")
        colnames(results_STAAR_A_1_25) <- c("ACAT-V(1,25)",
                                            paste0("ACAT-V(1,25)-", colnames(annotation_phred)),
                                            "STAAR-A(1,25)")
        colnames(results_STAAR_A_1_1) <- c("ACAT-V(1,1)",
                                           paste0("ACAT-V(1,1)-", colnames(annotation_phred)),
                                           "STAAR-A(1,1)")
      }
      return(list(num_variant = num_variant, cMAC = cMAC, RV_label = RV_label,
                  results_STAAR_O = results_STAAR_O, results_ACAT_O = results_ACAT_O,
                  results_STAAR_S_1_25 = results_STAAR_S_1_25, results_STAAR_S_1_1 = results_STAAR_S_1_1,
                  results_STAAR_B_1_25 = results_STAAR_B_1_25, results_STAAR_B_1_1 = results_STAAR_B_1_1,
                  results_STAAR_A_1_25 = results_STAAR_A_1_25, results_STAAR_A_1_1 = results_STAAR_A_1_1))
    }
    else {
      stop(paste0("Number of rare variant in the set is less than ",
                  rv_num_cutoff, "!"))
    }
  }
  
  coding_PheWAS_hwe_2norm_SPA <- function(chr,gene_name,genofile,obj_nullmodel_list,genes,
                                          rare_maf_cutoff=0.01,rv_num_cutoff=2,rv_num_cutoff_max=1e9,rv_num_cutoff_max_prefilter=1e9,
                                          QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                          Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                          Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                          SPA_p_filter=FALSE,p_filter_cutoff=0.05,silent=FALSE,
                                          pi1_all, ratio_all){
    
    ## evaluate choices
    variant_type <- match.arg(variant_type)
    geno_missing_imputation <- match.arg(geno_missing_imputation)
    
    phenotype.id <- Reduce(union, lapply(obj_nullmodel_list, function(x){
      as.character(x$id_include)
    }))
    
    ## get SNV id, position, REF, ALT (whole genome)
    filter <- seqGetData(genofile, QC_label)
    if(variant_type=="variant")
    {
      SNVlist <- filter == "PASS"
    }
    
    if(variant_type=="SNV")
    {
      SNVlist <- (filter == "PASS") & isSNV(genofile)
    }
    
    if(variant_type=="Indel")
    {
      SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    }
    
    position <- as.numeric(seqGetData(genofile, "position"))
    variant.id <- seqGetData(genofile, "variant.id")
    
    rm(filter)
    gc()
    
    ### Gene
    kk <- which(genes[,1]==gene_name)
    
    sub_start_loc <- genes[kk,3]
    sub_end_loc <- genes[kk,4]
    
    is.in <- (SNVlist)&(position>=sub_start_loc)&(position<=sub_end_loc)
    variant.id.gene <- variant.id[is.in]
    
    rm(position)
    gc()
    
    seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)
    
    ## HWE is assumed to be handled upstream in variant QC (QC_label); no HWE filtering here.
    seqSetFilter(genofile,variant.id=variant.id.gene)
    
    ## Gencode_Exonic
    GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
    ## Gencode
    GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
    ## Meta.SVM.Pred
    MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))
    
    ################################################
    #           Coding
    ################################################
    variant.id.gene <- seqGetData(genofile, "variant.id")
    lof.in.coding <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|(GENCODE.EXONIC.Category=="nonsynonymous SNV")|(GENCODE.EXONIC.Category=="synonymous SNV")
    variant.id.gene <- variant.id.gene[lof.in.coding]
    
    if(length(variant.id.gene)==0)
    {
      results_list_sub <- rep(list(c()), length(obj_nullmodel_list))
      results_coding_list <- list(plof = results_list_sub,
                                  plof_ds = results_list_sub,
                                  missense = results_list_sub,
                                  disruptive_missense = results_list_sub,
                                  synonymous = results_list_sub)
      seqResetFilter(genofile)
      return(results_coding_list)
    }
    
    seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)
    
    ## Gencode_Exonic
    GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
    ## Gencode
    GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
    ## Meta.SVM.Pred
    MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))
    
    ## Annotation
    Anno.Int.PHRED.sub <- NULL
    Anno.Int.PHRED.sub.name <- NULL
    
    if(variant_type=="SNV")
    {
      if(Use_annotation_weights)
      {
        for(k in 1:length(Annotation_name))
        {
          if(Annotation_name[k]%in%Annotation_name_catalog$name)
          {
            Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
            Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))
            
            if(Annotation_name[k]=="CADD")
            {
              Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
            }
            
            if(Annotation_name[k]=="aPC.LocalDiversity")
            {
              Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
              Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
              Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
            }
            Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
          }
        }
        
        Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
        colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
      }
    }
    
    ## get AF, Missing rate
    AF_AC_Missing <- seqGetAF_AC_Missing(genofile,minor=FALSE,parallel=FALSE)
    REF_AF <- AF_AC_Missing$af
    Missing_rate <- AF_AC_Missing$miss
    variant_maf_cutoff_filter <- ifelse(rare_maf_cutoff<=0.01,0.05,1)
    
    ################################################
    #                  plof_ds
    ################################################
    print('plof_ds')
    variant.id.gene <- seqGetData(genofile, "variant.id")
    lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|((GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D"))
    variant.id.gene.category <- variant.id.gene[lof.in.plof]
    
    ## Annotation
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.plof,]
    
    Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.id.gene.category,
                                          sample.id=phenotype.id,
                                          REF_AF=REF_AF[lof.in.plof],variant_maf_cutoff_filter=variant_maf_cutoff_filter,
                                          Missing_rate=Missing_rate[lof.in.plof],
                                          rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                          annotation_phred=Anno.Int.PHRED.sub.category)
    Geno <- Genotype_sp$Geno
    Anno.Int.PHRED.sub.category <- Genotype_sp$annotation_phred
    rm(Genotype_sp)
    gc()
    
    results_list_plof_ds <- rep(list(c()), length(obj_nullmodel_list))
    phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))
    
    if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
    {
      for(num in 1:length(obj_nullmodel_list))
      {
        id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
        id.phenotype.num.merge <- data.frame(id.phenotype.num)
        id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
        phenotype.id.match <- id.phenotype.num.merge$index
        samplesize.num <- length(phenotype.id.match)
        
        ## Extract the genotype matrix and variant information of the current trait
        Geno.num <- Geno[phenotype.id.match, , drop = FALSE]
        MAC.in <- Matrix::colSums(Geno.num,na.rm = TRUE)
        #MAC0.in <- Matrix::colSums(Geno[phenotype.id.match[which(obj_nullmodel$R==1)], , drop = FALSE],na.rm = TRUE)
        Missing_num.in <- Missing_num.sp(Geno.num)
        MAF.in <- MAC.in/(2*(samplesize.num-Missing_num.in))
        
        if (geno_missing_imputation == "mean")
        {
          Geno.num <- na.replace.sp(Geno.num,m=2*MAF.in)
        }
        if (geno_missing_imputation == "minor")
        {
          Geno.num <- na.replace.sp(Geno.num,is_NA_to_Zero=TRUE)
          MAF.in <- MAC.in/(2*samplesize.num)
        }
        
        ## number of traits in analysis
        n_pheno <- obj_nullmodel_list[[num]]$n.pheno
        
        ## SPA status
        if(!is.null(obj_nullmodel_list[[num]]$use_SPA))
        {
          use_SPA <- obj_nullmodel_list[[num]]$use_SPA
        }else
        {
          use_SPA <- FALSE
        }
        
        pvalues <- 0
        if(n_pheno == 1)
        {
          if(!use_SPA)
          {
            try(pvalues <- STAAR_sp_2norm_SPA(Geno.num,MAF.in,obj_nullmodel_list[[num]],
                                              Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,
                                              rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,pi1=pi1_all[num], resid_ratio=ratio_all[num]),silent=silent)
          }else{
            try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
          }
        }else
        {
          try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
        }
        
        results_plof_ds <- c()
        if(inherits(pvalues, "list"))
        {
          results_temp <- as.vector(genes[kk,])
          results_temp[3] <- "plof_ds"
          results_temp[2] <- chr
          results_temp[1] <- as.character(genes[kk,1])
          results_temp[4] <- pvalues$num_variant
          
          if(!use_SPA)
          {
            results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                              pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                              pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
          }else
          {
            results_temp <- c(results_temp,pvalues$cMAC,
                              pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
          }
          
          results_plof_ds <- rbind(results_plof_ds,results_temp)
        }
        
        if(!is.null(results_plof_ds))
        {
          if(!use_SPA)
          {
            colnames(results_plof_ds) <- colnames(results_plof_ds, do.NULL = FALSE, prefix = "col")
            colnames(results_plof_ds)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
            colnames(results_plof_ds)[(dim(results_plof_ds)[2]-1):dim(results_plof_ds)[2]] <- c("ACAT-O","STAAR-O")
          }else
          {
            colnames(results_plof_ds) <- colnames(results_plof_ds, do.NULL = FALSE, prefix = "col")
            colnames(results_plof_ds)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
            colnames(results_plof_ds)[dim(results_plof_ds)[2]] <- c("STAAR-B")
          }
          
          results_list_plof_ds[[num]] <- results_plof_ds
        }
      }
    }
    
    #####################################################
    #                      plof
    #####################################################
    print('plof')
    
    lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")
    variant.id.gene.category <- variant.id.gene[lof.in.plof]
    
    ## Annotation
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.plof,]
    
    Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.id.gene.category,
                                          sample.id=phenotype.id,
                                          REF_AF=REF_AF[lof.in.plof],variant_maf_cutoff_filter=variant_maf_cutoff_filter,
                                          Missing_rate=Missing_rate[lof.in.plof],
                                          rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                          annotation_phred=Anno.Int.PHRED.sub.category)
    Geno <- Genotype_sp$Geno
    Anno.Int.PHRED.sub.category <- Genotype_sp$annotation_phred
    rm(Genotype_sp)
    gc()
    
    results_list_plof <- rep(list(c()), length(obj_nullmodel_list))
    phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))
    
    if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
    {
      for(num in 1:length(obj_nullmodel_list))
      {
        id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
        id.phenotype.num.merge <- data.frame(id.phenotype.num)
        id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
        phenotype.id.match <- id.phenotype.num.merge$index
        samplesize.num <- length(phenotype.id.match)
        
        ## Extract the genotype matrix and variant information of the current trait
        Geno.num <- Geno[phenotype.id.match, , drop = FALSE]
        MAC.in <- Matrix::colSums(Geno.num,na.rm = TRUE)
        Missing_num.in <- Missing_num.sp(Geno.num)
        MAF.in <- MAC.in/(2*(samplesize.num-Missing_num.in))
        
        if (geno_missing_imputation == "mean")
        {
          Geno.num <- na.replace.sp(Geno.num,m=2*MAF.in)
        }
        if (geno_missing_imputation == "minor")
        {
          Geno.num <- na.replace.sp(Geno.num,is_NA_to_Zero=TRUE)
          MAF.in <- MAC.in/(2*samplesize.num)
        }
        
        ## number of traits in analysis
        n_pheno <- obj_nullmodel_list[[num]]$n.pheno
        
        ## SPA status
        if(!is.null(obj_nullmodel_list[[num]]$use_SPA))
        {
          use_SPA <- obj_nullmodel_list[[num]]$use_SPA
        }else
        {
          use_SPA <- FALSE
        }
        
        pvalues <- 0
        if(n_pheno == 1)
        {
          if(!use_SPA)
          {
            try(pvalues <- STAAR_sp_2norm_SPA(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,
                                              rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,pi1=pi1_all[num], resid_ratio=ratio_all[num]),silent=silent)
          }else{
            try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
          }
        }else
        {
          try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
        }
        
        results_plof <- c()
        if(inherits(pvalues, "list"))
        {
          results_temp <- as.vector(genes[kk,])
          results_temp[3] <- "plof"
          results_temp[2] <- chr
          results_temp[1] <- as.character(genes[kk,1])
          results_temp[4] <- pvalues$num_variant
          
          if(!use_SPA)
          {
            results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                              pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                              pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
          }else
          {
            results_temp <- c(results_temp,pvalues$cMAC,
                              pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
          }
          
          results_plof <- rbind(results_plof,results_temp)
        }
        
        if(!is.null(results_plof))
        {
          if(!use_SPA)
          {
            colnames(results_plof) <- colnames(results_plof, do.NULL = FALSE, prefix = "col")
            colnames(results_plof)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
            colnames(results_plof)[(dim(results_plof)[2]-1):dim(results_plof)[2]] <- c("ACAT-O","STAAR-O")
          }else
          {
            colnames(results_plof) <- colnames(results_plof, do.NULL = FALSE, prefix = "col")
            colnames(results_plof)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
            colnames(results_plof)[dim(results_plof)[2]] <- c("STAAR-B")
          }
          
          results_list_plof[[num]] <- results_plof
        }
      }
    }
    
    #############################################
    #             synonymous
    #############################################
    if(T){
      print('synonymous')
      lof.in.synonymous <- (GENCODE.EXONIC.Category=="synonymous SNV")
      variant.id.gene.category <- variant.id.gene[lof.in.synonymous]
      
      ## Annotation
      Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.synonymous,]
      
      Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.id.gene.category,
                                            sample.id=phenotype.id,
                                            REF_AF=REF_AF[lof.in.synonymous],variant_maf_cutoff_filter=variant_maf_cutoff_filter,
                                            Missing_rate=Missing_rate[lof.in.synonymous],
                                            rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                            annotation_phred=Anno.Int.PHRED.sub.category)
      Geno <- Genotype_sp$Geno
      Anno.Int.PHRED.sub.category <- Genotype_sp$annotation_phred
      rm(Genotype_sp)
      gc()
      
      results_list_synonymous <- rep(list(c()), length(obj_nullmodel_list))
      phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))
      
      if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
      {
        for(num in 1:length(obj_nullmodel_list))
        {
          id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
          id.phenotype.num.merge <- data.frame(id.phenotype.num)
          id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
          phenotype.id.match <- id.phenotype.num.merge$index
          samplesize.num <- length(phenotype.id.match)
          
          ## Extract the genotype matrix and variant information of the current trait
          Geno.num <- Geno[phenotype.id.match, , drop = FALSE]
          MAC.in <- Matrix::colSums(Geno.num,na.rm = TRUE)
          Missing_num.in <- Missing_num.sp(Geno.num)
          MAF.in <- MAC.in/(2*(samplesize.num-Missing_num.in))
          
          if (geno_missing_imputation == "mean")
          {
            Geno.num <- na.replace.sp(Geno.num,m=2*MAF.in)
          }
          if (geno_missing_imputation == "minor")
          {
            Geno.num <- na.replace.sp(Geno.num,is_NA_to_Zero=TRUE)
            MAF.in <- MAC.in/(2*samplesize.num)
          }
          
          ## number of traits in analysis
          n_pheno <- obj_nullmodel_list[[num]]$n.pheno
          
          ## SPA status
          if(!is.null(obj_nullmodel_list[[num]]$use_SPA))
          {
            use_SPA <- obj_nullmodel_list[[num]]$use_SPA
          }else
          {
            use_SPA <- FALSE
          }
          
          pvalues <- 0
          if(n_pheno == 1)
          {
            if(!use_SPA)
            {
              try(pvalues <- STAAR_sp_2norm_SPA(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,
                                                rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max, pi1=pi1_all[num], resid_ratio=ratio_all[num]),silent=silent)
            }else{
              try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
            }
          }else
          {
            try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
          }
          
          results_synonymous <- c()
          if(inherits(pvalues, "list"))
          {
            results_temp <- as.vector(genes[kk,])
            results_temp[3] <- "synonymous"
            results_temp[2] <- chr
            results_temp[1] <- as.character(genes[kk,1])
            results_temp[4] <- pvalues$num_variant
            
            if(!use_SPA)
            {
              results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                                pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                                pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
            }else
            {
              results_temp <- c(results_temp,pvalues$cMAC,
                                pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
            }
            
            results_synonymous <- rbind(results_synonymous,results_temp)
          }
          
          if(!is.null(results_synonymous))
          {
            if(!use_SPA)
            {
              colnames(results_synonymous) <- colnames(results_synonymous, do.NULL = FALSE, prefix = "col")
              colnames(results_synonymous)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
              colnames(results_synonymous)[(dim(results_synonymous)[2]-1):dim(results_synonymous)[2]] <- c("ACAT-O","STAAR-O")
            }else
            {
              colnames(results_synonymous) <- colnames(results_synonymous, do.NULL = FALSE, prefix = "col")
              colnames(results_synonymous)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
              colnames(results_synonymous)[dim(results_synonymous)[2]] <- c("STAAR-B")
            }
            
            results_list_synonymous[[num]] <- results_synonymous
          }
        }
      }
    }
    #################################################
    #        missense
    #################################################
    print('missense')
    
    lof.in.missense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")
    variant.id.gene.category <- variant.id.gene[lof.in.missense]
    
    ## Annotation
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.missense,]
    
    Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.id.gene.category,
                                          sample.id=phenotype.id,
                                          REF_AF=REF_AF[lof.in.missense],variant_maf_cutoff_filter=variant_maf_cutoff_filter,
                                          Missing_rate=Missing_rate[lof.in.missense],
                                          rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                          annotation_phred=Anno.Int.PHRED.sub.category)
    Geno <- Genotype_sp$Geno
    Anno.Int.PHRED.sub.category <- Genotype_sp$annotation_phred
    rm(Genotype_sp)
    gc()
    
    results_list_missense <- rep(list(c()), length(obj_nullmodel_list))
    phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))
    
    if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
    {
      for(num in 1:length(obj_nullmodel_list))
      {
        id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
        id.phenotype.num.merge <- data.frame(id.phenotype.num)
        id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
        phenotype.id.match <- id.phenotype.num.merge$index
        samplesize.num <- length(phenotype.id.match)
        
        ## Extract the genotype matrix and variant information of the current trait
        Geno.num <- Geno[phenotype.id.match, , drop = FALSE]
        MAC.in <- Matrix::colSums(Geno.num,na.rm = TRUE)
        Missing_num.in <- Missing_num.sp(Geno.num)
        MAF.in <- MAC.in/(2*(samplesize.num-Missing_num.in))
        
        if (geno_missing_imputation == "mean")
        {
          Geno.num <- na.replace.sp(Geno.num,m=2*MAF.in)
        }
        if (geno_missing_imputation == "minor")
        {
          Geno.num <- na.replace.sp(Geno.num,is_NA_to_Zero=TRUE)
          MAF.in <- MAC.in/(2*samplesize.num)
        }
        
        ## number of traits in analysis
        n_pheno <- obj_nullmodel_list[[num]]$n.pheno
        
        ## SPA status
        if(!is.null(obj_nullmodel_list[[num]]$use_SPA))
        {
          use_SPA <- obj_nullmodel_list[[num]]$use_SPA
        }else
        {
          use_SPA <- FALSE
        }
        
        pvalues <- 0
        if(n_pheno == 1)
        {
          if(!use_SPA)
          {
            try(pvalues <- STAAR_sp_2norm_SPA(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,
                                              rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,pi1=pi1_all[num], resid_ratio=ratio_all[num]),silent=silent)
          }else{
            try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
          }
        }else
        {
          try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
        }
        
        results <- c()
        if(inherits(pvalues, "list"))
        {
          results_temp <- as.vector(genes[kk,])
          results_temp[3] <- "missense"
          results_temp[2] <- chr
          results_temp[1] <- as.character(genes[kk,1])
          results_temp[4] <- pvalues$num_variant
          
          if(!use_SPA)
          {
            results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                              pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                              pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
          }else
          {
            results_temp <- c(results_temp,pvalues$cMAC,
                              pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
          }
          
          results <- rbind(results,results_temp)
          results_list_missense[[num]] <- results
        }
      }
    }
    
    #################################################
    #         disruptive missense
    #################################################
    print('d-missense')
    lof.in.dmissense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D")
    variant.id.gene.category <- variant.id.gene[lof.in.dmissense]
    
    ## Annotation
    Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.dmissense,]
    
    Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.id.gene.category,
                                          sample.id=phenotype.id,
                                          REF_AF=REF_AF[lof.in.dmissense],variant_maf_cutoff_filter=variant_maf_cutoff_filter,
                                          Missing_rate=Missing_rate[lof.in.dmissense],
                                          rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                          annotation_phred=Anno.Int.PHRED.sub.category)
    Geno <- Genotype_sp$Geno
    Anno.Int.PHRED.sub.category <- Genotype_sp$annotation_phred
    rm(Genotype_sp)
    gc()
    
    results_list_disruptive_missense <- rep(list(c()), length(obj_nullmodel_list))
    phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))
    
    for (num in 1:length(obj_nullmodel_list))
    {
      id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
      id.phenotype.num.merge <- data.frame(id.phenotype.num)
      id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
      phenotype.id.match <- id.phenotype.num.merge$index
      samplesize.num <- length(phenotype.id.match)
      
      if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
      {
        ## Extract the genotype matrix and variant information of the current trait
        Geno.num <- Geno[phenotype.id.match, , drop = FALSE]
        MAC.in <- Matrix::colSums(Geno.num,na.rm = TRUE)
        Missing_num.in <- Missing_num.sp(Geno.num)
        MAF.in <- MAC.in/(2*(samplesize.num-Missing_num.in))
        
        if (geno_missing_imputation == "mean")
        {
          Geno.num <- na.replace.sp(Geno.num,m=2*MAF.in)
        }
        if (geno_missing_imputation == "minor")
        {
          Geno.num <- na.replace.sp(Geno.num,is_NA_to_Zero=TRUE)
          MAF.in <- MAC.in/(2*samplesize.num)
        }
      } else
      {
        Geno.num <- NULL
        MAF.in <- numeric()
      }
      
      ## number of traits in analysis
      n_pheno <- obj_nullmodel_list[[num]]$n.pheno
      
      ## SPA status
      if(!is.null(obj_nullmodel_list[[num]]$use_SPA))
      {
        use_SPA <- obj_nullmodel_list[[num]]$use_SPA
      }else
      {
        use_SPA <- FALSE
      }
      
      pvalues <- 0
      if(n_pheno == 1)
      {
        if(!use_SPA)
        {
          try(pvalues <- STAAR_sp_2norm_SPA(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,
                                            rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max, pi1=pi1_all[num], resid_ratio=ratio_all[num]),silent=silent)
        }else{
          try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
        }
      }else
      {
        try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
      }
      
      if(inherits(pvalues, "list"))
      {
        results_temp <- as.vector(genes[kk,])
        results_temp[3] <- "disruptive_missense"
        results_temp[2] <- chr
        results_temp[1] <- as.character(genes[kk,1])
        results_temp[4] <- pvalues$num_variant
        
        if(!use_SPA)
        {
          results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
                            pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
        }else
        {
          results_temp <- c(results_temp,pvalues$cMAC,
                            pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
        }
        
        results <- rbind(results_list_missense[[num]],results_temp)
      }else
      {
        results <- results_list_missense[[num]]
      }
      
      if(!is.null(results))
      {
        if(!use_SPA)
        {
          colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
          colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("ACAT-O","STAAR-O")
        }else
        {
          colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
          colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
          colnames(results)[dim(results)[2]] <- c("STAAR-B")
        }
        
        if(dim(results)[1]==1)
        {
          if(results[3]!="disruptive_missense")
          {
            if(!use_SPA)
            {
              results <- cbind(results,matrix(1,1,6))
              colnames(results)[(dim(results)[2]-5):dim(results)[2]] <- c("SKAT(1,25)-Disruptive","SKAT(1,1)-Disruptive","Burden(1,25)-Disruptive","Burden(1,1)-Disruptive","ACAT-V(1,25)-Disruptive","ACAT-V(1,1)-Disruptive")
              results_missense <- results
              results_ds <- c()
            }else{
              results <- cbind(results,matrix(1,1,2))
              colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("Burden(1,25)-Disruptive","Burden(1,1)-Disruptive")
              results_missense <- results
              results_ds <- c()
            }
            
            results_list_missense[[num]] <- results
          }else
          {
            results_missense <- c()
            results_ds <- results
            results <- c()
            results_list_disruptive_missense[[num]] <- results_ds
          }
        }
        
        if(!is.null(results))
        {
          if(dim(results)[1]==2)
          {
            if(!use_SPA)
            {
              results_m <- c(results[1,],rep(0,6))
              names(results_m)[(length(results_m)-5):length(results_m)] <- c("SKAT(1,25)-Disruptive","SKAT(1,1)-Disruptive","Burden(1,25)-Disruptive","Burden(1,1)-Disruptive","ACAT-V(1,25)-Disruptive","ACAT-V(1,1)-Disruptive")
              results_m[(length(results_m)-5):length(results_m)] <- results[2,c("SKAT(1,25)","SKAT(1,1)","Burden(1,25)","Burden(1,1)","ACAT-V(1,25)","ACAT-V(1,1)")]
              apc_num <- (length(results_m)-19)/6
              p_seq <- c(1:apc_num,1:apc_num+(apc_num+1),1:apc_num+2*(apc_num+1),1:apc_num+3*(apc_num+1),1:apc_num+4*(apc_num+1),1:apc_num+5*(apc_num+1),(6*apc_num+9):(6*apc_num+14))
              results_m["STAAR-O"] <- CCT(as.numeric(results_m[6:length(results_m)][p_seq]))
              results_m["STAAR-S(1,25)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num,6*apc_num+9)]))
              results_m["STAAR-S(1,1)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+(apc_num+1),6*apc_num+10)]))
              results_m["STAAR-B(1,25)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+2*(apc_num+1),6*apc_num+11)]))
              results_m["STAAR-B(1,1)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+3*(apc_num+1),6*apc_num+12)]))
              results_m["STAAR-A(1,25)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+4*(apc_num+1),6*apc_num+13)]))
              results_m["STAAR-A(1,1)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+5*(apc_num+1),6*apc_num+14)]))
              
              results_ds <- c()
              results_ds <- rbind(results_ds,results[2,])
              
              results <- c()
              results <- rbind(results,results_m)
              
              results_list_missense[[num]] <- results
              results_list_disruptive_missense[[num]] <- results_ds
            }else
            {
              results_m <- c(results[1,],rep(0,2))
              names(results_m)[(length(results_m)-1):length(results_m)] <- c("Burden(1,25)-Disruptive","Burden(1,1)-Disruptive")
              results_m[(length(results_m)-1):length(results_m)] <- results[2,c("Burden(1,25)","Burden(1,1)")]
              
              ## check whether the p-values is NA. If so, set NA equals 1.
              if(is.na(results_m[(length(results_m)-1)]))
              {
                results_m[(length(results_m)-1)] <- 1
              }
              
              if(is.na(results_m[length(results_m)]))
              {
                results_m[length(results_m)] <- 1
              }
              
              apc_num <- (length(results_m)-10)/2
              p_seq <- c(1:apc_num,1:apc_num+(apc_num+1),(length(results_m)-6):(length(results_m)-5))
              
              ## calculate STAAR-B
              pvalues_sub <- as.numeric(results_m[6:length(results_m)][p_seq])
              if(sum(is.na(pvalues_sub))>0)
              {
                if(sum(is.na(pvalues_sub))==length(pvalues_sub))
                {
                  results_m["STAAR-B"] <- 1
                }else
                {
                  ## not all NAs
                  pvalues_sub <- pvalues_sub[!is.na(pvalues_sub)]
                  if(sum(pvalues_sub[pvalues_sub<1])>0)
                  {
                    ## not all ones
                    results_m["STAAR-B"] <- CCT(pvalues_sub[pvalues_sub<1])
                    
                  }else
                  {
                    results_m["STAAR-B"] <- 1
                    
                  }
                }
              }else
              {
                if(sum(pvalues_sub[pvalues_sub<1])>0)
                {
                  results_m["STAAR-B"] <- CCT(pvalues_sub[pvalues_sub<1])
                }else
                {
                  results_m["STAAR-B"] <- 1
                }
              }
              
              ## calculate STAAR-B(1,25)
              pvalues_sub <- as.numeric(results_m[6:length(results_m)][c(1:apc_num,(length(results_m)-6))])
              if(sum(is.na(pvalues_sub))>0)
              {
                if(sum(is.na(pvalues_sub))==length(pvalues_sub))
                {
                  results_m["STAAR-B(1,25)"] <- 1
                }else
                {
                  ## not all NAs
                  pvalues_sub <- pvalues_sub[!is.na(pvalues_sub)]
                  if(sum(pvalues_sub[pvalues_sub<1])>0)
                  {
                    ## not all ones
                    results_m["STAAR-B(1,25)"] <- CCT(pvalues_sub[pvalues_sub<1])
                    
                  }else
                  {
                    results_m["STAAR-B(1,25)"] <- 1
                    
                  }
                }
              }else
              {
                if(sum(pvalues_sub[pvalues_sub<1])>0)
                {
                  results_m["STAAR-B(1,25)"] <- CCT(pvalues_sub[pvalues_sub<1])
                }else
                {
                  results_m["STAAR-B(1,25)"] <- 1
                }
              }
              
              ## calculate STAAR-B(1,1)
              pvalues_sub <- as.numeric(results_m[6:length(results_m)][c(1:apc_num+(apc_num+1),(length(results_m)-5))])
              if(sum(is.na(pvalues_sub))>0)
              {
                if(sum(is.na(pvalues_sub))==length(pvalues_sub))
                {
                  results_m["STAAR-B(1,1)"] <- 1
                }else
                {
                  ## not all NAs
                  pvalues_sub <- pvalues_sub[!is.na(pvalues_sub)]
                  if(sum(pvalues_sub[pvalues_sub<1])>0)
                  {
                    ## not all ones
                    results_m["STAAR-B(1,1)"] <- CCT(pvalues_sub[pvalues_sub<1])
                    
                  }else
                  {
                    results_m["STAAR-B(1,1)"] <- 1
                    
                  }
                }
              }else
              {
                if(sum(pvalues_sub[pvalues_sub<1])>0)
                {
                  results_m["STAAR-B(1,1)"] <- CCT(pvalues_sub[pvalues_sub<1])
                }else
                {
                  results_m["STAAR-B(1,1)"] <- 1
                }
              }
              
              results_ds <- c()
              results_ds <- rbind(results_ds,results[2,])
              
              results <- c()
              results <- rbind(results,results_m)
              
              results_list_missense[[num]] <- results
              results_list_disruptive_missense[[num]] <- results_ds
            }
          }
        }
      }else
      {
        results <- c()
        results_ds <- c()
      }
    }
    
    results_coding_list <- list(plof = results_list_plof,
                                plof_ds = results_list_plof_ds,
                                missense = results_list_missense,
                                disruptive_missense = results_list_disruptive_missense,
                                synonymous = results_list_synonymous)
    
    seqResetFilter(genofile)
    
    return(results_coding_list)
  }
  
  
  
}

library(devtools)
library(STAARpipelinePheWAS)
library(data.table)
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
#library(MetaSTAAR)
#devtools::install_github("li-lab-genetics/MetaSTAARlite",ref="main")
# devtools::install_github(
#   "li-lab-genetics/MetaSTAARlite",
#   ref       = "main",
#   upgrade   = "never",        
# )
#library(MetaSTAARlite)

###########################################################
#           User Input
###########################################################
## Number of jobs for each chromosome

## --- input locations in your DNANexus (RAP) project ---
## Directory holding the per-phenotype null models from Step 3.2, named
##   obj.STAAR_500K_<field>_norm_pseudo_rf100_v5.rda  (each stores obj_nullmodel/resid_ratio/pi1)
dx_null_model_dir <- "/null_model/"
## Directory of annotated GDS (aGDS) files, one per chromosome:
dx_agds_dir <- "/aGDS/"
## Project directory to upload per-job results into (<dx_results_dir>/<arrayid>/):
dx_results_dir <- "/results/coding/"
## STAARpipeline annotation channel catalog (.RData providing `Annotation_name_catalog`).
annotation_catalog_file <- "Annotation_name_catalog.RData"
system(paste0("dx download -f /", annotation_catalog_file))
if (grepl("\\.csv$", annotation_catalog_file)) {
  Annotation_name_catalog <- read.csv(annotation_catalog_file)
} else {
  load(annotation_catalog_file)
}
#obj_nullmodel <- get(load(paste0('./CAD_250718.Rdata')))
rare_maf_cutoff=0.01;
rv_num_cutoff=2;
#varRatio=ratio;
QC_label= 'annotation/info/QC_label2'
#QC_label= 'annotation/filter'
variant_type=c("SNV");
geno_missing_imputation=c("mean");
Annotation_dir <- "annotation/info/FunctionalAnnotation"

## output file name
output_path <- ".//step3_output//"
output_file_name <- "Coding_"
dir.create(output_path, recursive = T)


#jobs_num <- get(load("./step0_output//jobs_num.Rdata"))
## aGDS directory
## Null model
#obj_nullmodel <- get(load(paste0("./step1_output/obj.STAAR_washU_callset_notopmed_v1.Rdata")))

## QC_label
#QC_label= 'annotation/info/QC_label2'



## Annotation channel
#Annotation_name_catalog <- get(load("./step0_output//Annotation_name_catalog.Rdata"))
# Or equivalently
# Annotation_name_catalog <- read.csv("/path_to_the_file/Annotation_name_catalog.csv")
## Use_annotation_weights
Use_annotation_weights <- T
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

#dir.create(output_path, recursive = T)
#if(!file.exists(paste0(output_path,output_file_name,"_sumstat_",arrayid,".Rdata"))){
### jobs_num
# load('./jobs_num.Rdata')
gene_num_in_array <- 100
## 1:81 : chr 1-7
## 82: 195: chr 8-22
group.num.allchr <- ceiling(table(genes_info[,2])/gene_num_in_array)
sum(group.num.allchr)

chr <- which.max(arrayid <= cumsum(group.num.allchr))

group.num <- group.num.allchr[chr]
if (chr == 1){
  groupid <- arrayid
}else{
  groupid <- arrayid - cumsum(group.num.allchr)[chr-1]
}

genes_info_chr <- genes_info[genes_info[,2]==chr,]
sub_seq_num <- dim(genes_info_chr)[1]

if(groupid < group.num)
{
  sub_seq_id <- ((groupid - 1)*gene_num_in_array + 1):(groupid*gene_num_in_array)
}else
{
  sub_seq_id <- ((groupid - 1)*gene_num_in_array + 1):sub_seq_num
}
print(paste("Chromosome:",chr))


## Resume: pull any already-computed results for this array job into ./<arrayid>/, so that
## genes finished in a previous (e.g. preempted) run are skipped by the check below.
system(paste0("dx mkdir -p ", dx_results_dir, arrayid, "/"))
system(paste0("dx download -r ", dx_results_dir, arrayid, "/"))





## Phenotypes to analyze (PheWAS). Each element is a phenotype field id; a pair
## (e.g. left/right IDPs) shares the null model keyed by its first id. For a
## single trait use  pheno_list <- list("<field>").
if(T){
  pheno_list <- list()
  for(k in 1:6){
    pheno_list[[k]] <- c(25055+k)
  }
  l <- 7
  for(k in c(1:21)){
    pheno_list[[l]] <- c(25062+2*k-2,25062+2*k-1)
    l <- l+1
  }
}


nn <- length(pheno_list)

ll <- 1
res <- list()
field_all <- c(#'25001','25003','25005','25007','25009','25010',
               unlist(lapply(pheno_list,function(x){x[1]})))


run_id <- c()
for(k in sub_seq_id){
  if(!file.exists(paste0(arrayid,"/",output_file_name,k,".Rdata"))){
    run_id <- c(run_id, k)
  }
}

sub_seq_id <- run_id

print(sub_seq_id)
if(length(sub_seq_id)>0){
  
  testlist <- list()
  pi1_all <- c()
  ratio_all <- c()
  l <- 1
  # field_all <- c('25001',
  #                '25003','25005','25007','25009','25010')
  
  for(ll in field_all){
    
    #for(l in 1:24){
    system(paste0("dx download -f ", dx_null_model_dir, "obj.STAAR_500K_",ll,"_norm_pseudo_rf100_v5.rda"))
    load(paste0("obj.STAAR_500K_",ll,"_norm_pseudo_rf100_v5.rda"))
    # obj_nullmodel$X                 <- NULL
    #   obj_nullmodel$Y                 <- NULL
    #   obj_nullmodel$residuals         <- NULL
    #   obj_nullmodel$linear.predictors <- NULL
    testlist[[l]] <-  obj_nullmodel
    ratio_all[l] <- resid_ratio
    pi1_all[l] <- pi1
    l <- l+1
    #testlist[[2]] <-  get(load(paste0("obj.STAAR_500K_26900_i2_norm.rda")))
  }
  
  
  ######################################################################################
  gds.path <- paste0("ukb.500k.wgs.chr",chr,".pass.annotated.gds")
  system(paste0("dx download -f ", dx_agds_dir, gds.path))
  genofile <- seqOpen(gds.path)
  #########################################################################################
  do_kk <- function(kk){
    print(kk)
    gene_name <- genes_info_chr[kk,1]
    start_time <- proc.time()
    results <- coding_PheWAS_hwe_2norm_SPA(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel_list=testlist,genes=genes_info,
                                           rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                           QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                           Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                           Use_annotation_weights=Use_annotation_weights,
                                           Annotation_name=Annotation_name[-2],
                                           pi1_all=pi1_all, ratio_all=ratio_all)
    end_time <- proc.time()
    print(end_time-start_time)
    save(results,file=paste0(output_path,output_file_name,kk,".Rdata"))
    system(paste0("dx mkdir -p ", dx_results_dir, arrayid, "/"))
    system(paste0("dx upload ", output_path, output_file_name, kk, ".Rdata --path ", dx_results_dir, arrayid, "/"))
  }
  library(parallel)
  res <- mclapply(sub_seq_id, do_kk, mc.cores=2)
  seqClose(genofile)
}










