###########################################################
#  SPARTA - Step 4.2: Gene-centric NON-CODING analysis  (UKB RAP / DNANexus)
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
   # print(pv.unadj)
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
      #print(w)
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
  # STAAR_sp_test <- function(genotype_sp,MAF=NULL,obj_nullmodel,annotation_phred=NULL,
  #                      rare_maf_cutoff=0.01,rv_num_cutoff=2,
  #                      rv_num_cutoff_max=1e9){
  #   print('check anno')
  #   print(dim(genotype_sp))
  #   #print(summary(annotation_phred))
  #   annotation_phred <- as.data.frame(annotation_phred)
  #   annotation_phred[is.na(annotation_phred)] <- 0 ## 防止报错sym()
  #   #print(summary(annotation_phred))
  #   return(STAAR_Binary_SPA_sp(genotype_sp=genotype_sp,MAF=MAF,obj_nullmodel=obj_nullmodel,
  #                              annotation_phred=annotation_phred,
  #                              rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
  #                              rv_num_cutoff_max=rv_num_cutoff_max))
  #   
  #   
  #   
  # }
  STAAR_sp_2norm_SPA_old <-   function (genotype_sp, MAF = NULL, obj_nullmodel, annotation_phred = NULL,
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
    MAC0 <-  colSums(genotype_sp)
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
        newadd <- T
        #print('no SPA')
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
  
  ### 调整了ACAT
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
       # print('no SPA')
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
        mac_thres <- 20
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
  
  
  noncoding_PheWAS_hwe_2norm_SPA <- function(chr,gene_name,genofile,obj_nullmodel_list,genes,
                                   rare_maf_cutoff=0.01,rv_num_cutoff=2,rv_num_cutoff_max=1e9,rv_num_cutoff_max_prefilter=1e9,
                                   QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                   Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                   Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                   SPA_p_filter=FALSE,p_filter_cutoff=0.05,silent=FALSE,
                                   id_all,pi1_all, ratio_all){
    
    ## evaluate choices
    variant_type <- match.arg(variant_type)
    geno_missing_imputation <- match.arg(geno_missing_imputation)
    
    phenotype.id <- Reduce(union, lapply(obj_nullmodel_list, function(x){
      as.character(x$id_include)
    }))
    
    ww <- length(id_all)
    if(ww<7){
      for(ll in (ww+1):7){
        id_all[[ll]] <- integer(0)
      }
    }
    
    ## get SNV id, position, REF, ALT (whole genome)
    # filter <- seqGetData(genofile, QC_label)
    # if(variant_type=="variant")
    # {
    #   SNVlist <- filter == "PASS"
    # }
    # 
    # if(variant_type=="SNV")
    # {
    #   SNVlist <- (filter == "PASS") & isSNV(genofile)
    # }
    # 
    # if(variant_type=="Indel")
    # {
    #   SNVlist <- (filter == "PASS") & (!isSNV(genofile))
    # }
    # 
    # position <- as.numeric(seqGetData(genofile, "position"))
    # variant.id <- seqGetData(genofile, "variant.id")
    # 
    # rm(filter)
    # gc()
    
    ### Gene
    # kk <- which(genes[,1]==gene_name)
    # 
    # sub_start_loc <- genes[kk,3]
    # sub_end_loc <- genes[kk,4]
    # 
    # is.in <- (SNVlist)&(position>=sub_start_loc)&(position<=sub_end_loc)
    # variant.id.gene <- variant.id[is.in]
    # 
    # rm(position)
    # gc()
    # 
    # seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)
    
    rv_num_cutoff_min_prefilter <- 1
    if(T){
      ########################################
      #   Downstream
      print('downstream')
      ### id list
      if(T){
        phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))
        
        id.phenotype.num.list <- list()
        for (num in 1:length(obj_nullmodel_list)){
          id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
          id.phenotype.num.merge <- data.frame(id.phenotype.num)
          id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
          phenotype.id.match <- id.phenotype.num.merge$index
          id.phenotype.num.list[[num]] <- phenotype.id.match
        }
      }
      results_list_downstream <- rep(list(c()), length(obj_nullmodel_list))
      
      if(length(id_all[[1]])>=rv_num_cutoff_min_prefilter & length(id_all[[1]])<rv_num_cutoff_max_prefilter){
        variant.is.in <- id_all[[1]]
        seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)
        
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
        
        if(length(variant.is.in)>=1)
        {
          ## get AF, Missing rate
          AF_AC_Missing <- seqGetAF_AC_Missing(genofile,minor=FALSE,parallel=FALSE)
          REF_AF <- AF_AC_Missing$af
          Missing_rate <- AF_AC_Missing$miss
          ## variant id
          variant.is.in <- seqGetData(genofile, "variant.id")
          variant_maf_cutoff_filter <- ifelse(rare_maf_cutoff<=0.01,0.05,1)
          
          #seqResetFilter(genofile)
          
          Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.is.in,
                                                sample.id=phenotype.id,
                                                REF_AF=REF_AF,variant_maf_cutoff_filter=variant_maf_cutoff_filter,
                                                Missing_rate=Missing_rate,
                                                rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                                annotation_phred=Anno.Int.PHRED.sub)
          Geno <- Genotype_sp$Geno
          Anno.Int.PHRED.sub <- Genotype_sp$annotation_phred
          rm(Genotype_sp)
          gc()
        } else
        {
          Geno <- NULL
        }
        
        
        if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
        {
          for (num in 1:length(obj_nullmodel_list))
          {
            # id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
            # id.phenotype.num.merge <- data.frame(id.phenotype.num)
            # id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
            # phenotype.id.match <- id.phenotype.num.merge$index
            phenotype.id.match <- id.phenotype.num.list[[num]]
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
                try(pvalues <- STAAR_sp_2norm_SPA(Geno.num,MAF.in,obj_nullmodel_list[[num]],
                                                  Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                                                  rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,pi1=pi1_all[num], resid_ratio=ratio_all[num]),silent=silent)
              }else
              {
                try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
              }
            }else
            {
              try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
            }
            
            results_downstream <- c()
            if(inherits(pvalues, "list"))
            {
              results_temp <- rep(NA,4)
              results_temp[3] <- "downstream"
              results_temp[2] <- chr
              results_temp[1] <- as.character(gene_name)
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
              
              results_downstream <- rbind(results_downstream,results_temp)
            }
            
            if(!is.null(results_downstream))
            {
              if(!use_SPA)
              {
                colnames(results_downstream) <- colnames(results_downstream, do.NULL = FALSE, prefix = "col")
                colnames(results_downstream)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
                colnames(results_downstream)[(dim(results_downstream)[2]-1):dim(results_downstream)[2]] <- c("ACAT-O","STAAR-O")
              }else
              {
                colnames(results_downstream) <- colnames(results_downstream, do.NULL = FALSE, prefix = "col")
                colnames(results_downstream)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
                colnames(results_downstream)[dim(results_downstream)[2]] <- c("STAAR-B")
              }
              
              results_list_downstream[[num]] <- results_downstream
            }
          }
        }
        #seqResetFilter(genofile)
      }else{
        seqSetFilter(genofile,sample.id=phenotype.id)
        
      }
      
      ########################################
      #   Upstream
      print('upstream')
      results_list_upstream <- rep(list(c()), length(obj_nullmodel_list))
      if(length(id_all[[2]])>=rv_num_cutoff_min_prefilter & length(id_all[[2]])<rv_num_cutoff_max_prefilter){
        variant.is.in <- id_all[[2]]
        seqSetFilter(genofile,variant.id=variant.is.in)
        
        
        
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
        
        if(length(variant.is.in)>=1)
        {
          ## get AF, Missing rate
          AF_AC_Missing <- seqGetAF_AC_Missing(genofile,minor=FALSE,parallel=FALSE)
          REF_AF <- AF_AC_Missing$af
          Missing_rate <- AF_AC_Missing$miss
          ## variant id
          variant.is.in <- seqGetData(genofile, "variant.id")
          variant_maf_cutoff_filter <- ifelse(rare_maf_cutoff<=0.01,0.05,1)
          #seqResetFilter(genofile)
          
          Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.is.in,
                                                sample.id=phenotype.id,
                                                REF_AF=REF_AF,variant_maf_cutoff_filter=variant_maf_cutoff_filter,
                                                Missing_rate=Missing_rate,
                                                rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                                annotation_phred=Anno.Int.PHRED.sub)
          Geno <- Genotype_sp$Geno
          Anno.Int.PHRED.sub <- Genotype_sp$annotation_phred
          rm(Genotype_sp)
          gc()
        } else
        {
          Geno <- NULL
        }
        
        #results_list_upstream <- rep(list(c()), length(obj_nullmodel_list))
        #phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))

        if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
        {
          for (num in 1:length(obj_nullmodel_list))
          {
            # id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
            # id.phenotype.num.merge <- data.frame(id.phenotype.num)
            # id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
            # phenotype.id.match <- id.phenotype.num.merge$index
            phenotype.id.match <- id.phenotype.num.list[[num]]
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
                try(pvalues <- STAAR_sp_2norm_SPA(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                                  rv_num_cutoff_max=rv_num_cutoff_max,pi1=pi1_all[num], resid_ratio=ratio_all[num]),silent=silent)
              }else
              {
                try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
              }
            }else
            {
              try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
            }
            
            results_upstream <- c()
            if(inherits(pvalues, "list"))
            {
              results_temp <- rep(NA,4)
              results_temp[3] <- "upstream"
              results_temp[2] <- chr
              results_temp[1] <- as.character(gene_name)
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
              
              results_upstream <- rbind(results_upstream,results_temp)
            }
            
            if(!is.null(results_upstream))
            {
              if(!use_SPA)
              {
                colnames(results_upstream) <- colnames(results_upstream, do.NULL = FALSE, prefix = "col")
                colnames(results_upstream)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
                colnames(results_upstream)[(dim(results_upstream)[2]-1):dim(results_upstream)[2]] <- c("ACAT-O","STAAR-O")
              }else
              {
                colnames(results_upstream) <- colnames(results_upstream, do.NULL = FALSE, prefix = "col")
                colnames(results_upstream)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
                colnames(results_upstream)[dim(results_upstream)[2]] <- c("STAAR-B")
              }
              
              results_list_upstream[[num]] <- results_upstream
            }
          }
        }
        
        #seqResetFilter(genofile)
      }
      ########################################################
      #                UTR
      print('UTR')
      results_list_UTR <- rep(list(c()), length(obj_nullmodel_list))
      if(length(id_all[[3]])>=rv_num_cutoff_min_prefilter & length(id_all[[3]])<rv_num_cutoff_max_prefilter){
        variant.is.in <- id_all[[3]]
        seqSetFilter(genofile,variant.id=variant.is.in)
        
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
        
        if(length(variant.is.in)>=1)
        {
          ## get AF, Missing rate
          AF_AC_Missing <- seqGetAF_AC_Missing(genofile,minor=FALSE,parallel=FALSE)
          REF_AF <- AF_AC_Missing$af
          Missing_rate <- AF_AC_Missing$miss
          ## variant id
          variant.is.in <- seqGetData(genofile, "variant.id")
          variant_maf_cutoff_filter <- ifelse(rare_maf_cutoff<=0.01,0.05,1)
          #seqResetFilter(genofile)
          
          Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.is.in,
                                                sample.id=phenotype.id,
                                                REF_AF=REF_AF,variant_maf_cutoff_filter=variant_maf_cutoff_filter,
                                                Missing_rate=Missing_rate,
                                                rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                                annotation_phred=Anno.Int.PHRED.sub)
          Geno <- Genotype_sp$Geno
          Anno.Int.PHRED.sub <- Genotype_sp$annotation_phred
          rm(Genotype_sp)
          gc()
        } else
        {
          Geno <- NULL
        }
        
        #results_list_UTR <- rep(list(c()), length(obj_nullmodel_list))
        #phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))

        if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
        {
          for (num in 1:length(obj_nullmodel_list))
          {
            # id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
            # id.phenotype.num.merge <- data.frame(id.phenotype.num)
            # id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
            # phenotype.id.match <- id.phenotype.num.merge$index
            phenotype.id.match <- id.phenotype.num.list[[num]]
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
                try(pvalues <- STAAR_sp_2norm_SPA(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                                  rv_num_cutoff_max=rv_num_cutoff_max,pi1=pi1_all[num], resid_ratio=ratio_all[num]),silent=silent)
              }else
              {
                try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
              }
            }else
            {
              try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
            }
            
            results_UTR <- c()
            if(inherits(pvalues, "list"))
            {
              results_temp <- rep(NA,4)
              results_temp[3] <- "UTR"
              results_temp[2] <- chr
              results_temp[1] <- as.character(gene_name)
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
              
              results_UTR <- rbind(results_UTR,results_temp)
            }
            
            if(!is.null(results_UTR))
            {
              if(!use_SPA)
              {
                colnames(results_UTR) <- colnames(results_UTR, do.NULL = FALSE, prefix = "col")
                colnames(results_UTR)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
                colnames(results_UTR)[(dim(results_UTR)[2]-1):dim(results_UTR)[2]] <- c("ACAT-O","STAAR-O")
              }else
              {
                colnames(results_UTR) <- colnames(results_UTR, do.NULL = FALSE, prefix = "col")
                colnames(results_UTR)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
                colnames(results_UTR)[dim(results_UTR)[2]] <- c("STAAR-B")
              }
              
              results_list_UTR[[num]] <- results_UTR
            }
          }
        }
        
        #seqResetFilter(genofile)
      }
      
      #############################################
      #   Promoter-CAGE
      print('Promoter-CAGE')
      
      ## Promoter
      results_list_promoter_CAGE <- rep(list(c()), length(obj_nullmodel_list))
      
      if(length(id_all[[4]])>=rv_num_cutoff_min_prefilter & length(id_all[[4]])<rv_num_cutoff_max_prefilter){
        variant.is.in <- id_all[[4]]
        seqSetFilter(genofile,variant.id=variant.is.in)
        
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
        
        if(length(variant.is.in)>=1)
        {
          ## get AF, Missing rate
          AF_AC_Missing <- seqGetAF_AC_Missing(genofile,minor=FALSE,parallel=FALSE)
          REF_AF <- AF_AC_Missing$af
          Missing_rate <- AF_AC_Missing$miss
          ## variant id
          variant.is.in <- seqGetData(genofile, "variant.id")
          variant_maf_cutoff_filter <- ifelse(rare_maf_cutoff<=0.01,0.05,1)
          #seqResetFilter(genofile)
          
          Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.is.in,
                                                sample.id=phenotype.id,
                                                REF_AF=REF_AF,variant_maf_cutoff_filter=variant_maf_cutoff_filter,
                                                Missing_rate=Missing_rate,
                                                rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                                annotation_phred=Anno.Int.PHRED.sub)
          Geno <- Genotype_sp$Geno
          Anno.Int.PHRED.sub <- Genotype_sp$annotation_phred
          rm(Genotype_sp)
          gc()
        } else
        {
          Geno <- NULL
        }
        
        #results_list_promoter_CAGE <- rep(list(c()), length(obj_nullmodel_list))
        #phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))
        
        if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
        {
          for (num in 1:length(obj_nullmodel_list))
          {
            # id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
            # id.phenotype.num.merge <- data.frame(id.phenotype.num)
            # id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
            # phenotype.id.match <- id.phenotype.num.merge$index
            phenotype.id.match <- id.phenotype.num.list[[num]]
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
                try(pvalues <- STAAR_sp_2norm_SPA(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                                  rv_num_cutoff_max=rv_num_cutoff_max,pi1=pi1_all[num], resid_ratio=ratio_all[num]),silent=silent)
              }else
              {
                try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
              }
            }else
            {
              try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
            }
            
            results_promoter_CAGE <- c()
            if(inherits(pvalues, "list"))
            {
              results_temp <- c()
              results_temp[3] <- "promoter_CAGE"
              results_temp[2] <- chr
              results_temp[1] <- as.character(gene_name)
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
              
              results_promoter_CAGE <- rbind(results_promoter_CAGE,results_temp)
            }
            
            if(!is.null(results_promoter_CAGE))
            {
              if(!use_SPA)
              {
                colnames(results_promoter_CAGE) <- colnames(results_promoter_CAGE, do.NULL = FALSE, prefix = "col")
                colnames(results_promoter_CAGE)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
                colnames(results_promoter_CAGE)[(dim(results_promoter_CAGE)[2]-1):dim(results_promoter_CAGE)[2]] <- c("ACAT-O","STAAR-O")
              }else
              {
                colnames(results_promoter_CAGE) <- colnames(results_promoter_CAGE, do.NULL = FALSE, prefix = "col")
                colnames(results_promoter_CAGE)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
                colnames(results_promoter_CAGE)[dim(results_promoter_CAGE)[2]] <- c("STAAR-B")
              }
              
              results_list_promoter_CAGE[[num]] <- results_promoter_CAGE
            }
          }
        }
        
        # seqResetFilter(genofile)
      }
      ##################################################
      #       Promoter-DHS
      print('Promoter-DHS')
      
      # Subsetting Promoters that within +/-3kb of TSS and have rOCRs signals
      results_list_promoter_DHS <- rep(list(c()), length(obj_nullmodel_list))
      if(length(id_all[[5]])>=rv_num_cutoff_min_prefilter & length(id_all[[5]])<rv_num_cutoff_max_prefilter){
        #seqSetFilter(genofile,variant.id=id_all[[5]])
        variant.is.in <- id_all[[5]]
        print('Promoter-DHS')
        seqSetFilter(genofile,variant.id=variant.is.in)
        
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
        
        if(length(variant.is.in)>=1)
        {
          ## get AF, Missing rate
          AF_AC_Missing <- seqGetAF_AC_Missing(genofile,minor=FALSE,parallel=FALSE)
          REF_AF <- AF_AC_Missing$af
          Missing_rate <- AF_AC_Missing$miss
          ## variant id
          variant.is.in <- seqGetData(genofile, "variant.id")
          variant_maf_cutoff_filter <- ifelse(rare_maf_cutoff<=0.01,0.05,1)
          #seqResetFilter(genofile)
          
          Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.is.in,
                                                sample.id=phenotype.id,
                                                REF_AF=REF_AF,variant_maf_cutoff_filter=variant_maf_cutoff_filter,
                                                Missing_rate=Missing_rate,
                                                rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                                annotation_phred=Anno.Int.PHRED.sub)
          Geno <- Genotype_sp$Geno
          Anno.Int.PHRED.sub <- Genotype_sp$annotation_phred
          rm(Genotype_sp)
          gc()
        } else
        {
          Geno <- NULL
        }
        
        #results_list_promoter_DHS <- rep(list(c()), length(obj_nullmodel_list))
        #phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))
        
        if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
        {
          for (num in 1:length(obj_nullmodel_list))
          {
            # id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
            # id.phenotype.num.merge <- data.frame(id.phenotype.num)
            # id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
            # phenotype.id.match <- id.phenotype.num.merge$index
            # id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
            # id.phenotype.num.merge <- data.frame(id.phenotype.num)
            # id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
            # phenotype.id.match <- id.phenotype.num.merge$index
            phenotype.id.match <- id.phenotype.num.list[[num]]
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
                try(pvalues <- STAAR_sp_2norm_SPA(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,
                                                  rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,pi1=pi1_all[num], resid_ratio=ratio_all[num]),silent=silent)
              }else
              {
                try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
              }
            }else
            {
              try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
            }
            
            results_promoter_DHS <- c()
            if(inherits(pvalues, "list"))
            {
              results_temp <- c()
              results_temp[3] <- "promoter_DHS"
              results_temp[2] <- chr
              results_temp[1] <- as.character(gene_name)
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
              
              results_promoter_DHS <- rbind(results_promoter_DHS,results_temp)
            }
            
            if(!is.null(results_promoter_DHS))
            {
              if(!use_SPA)
              {
                colnames(results_promoter_DHS) <- colnames(results_promoter_DHS, do.NULL = FALSE, prefix = "col")
                colnames(results_promoter_DHS)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
                colnames(results_promoter_DHS)[(dim(results_promoter_DHS)[2]-1):dim(results_promoter_DHS)[2]] <- c("ACAT-O","STAAR-O")
              }else
              {
                colnames(results_promoter_DHS) <- colnames(results_promoter_DHS, do.NULL = FALSE, prefix = "col")
                colnames(results_promoter_DHS)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
                colnames(results_promoter_DHS)[dim(results_promoter_DHS)[2]] <- c("STAAR-B")
              }
              
              results_list_promoter_DHS[[num]] <- results_promoter_DHS
            }
          }
        }
        
        #seqResetFilter(genofile)
      }
      ###########################################
      #        Enhancer-CAGE
      print('Enhancer-CAGE')
      
      #Now extract the GeneHancer with CAGE Signal Overlay
      results_list_enhancer_CAGE <- rep(list(c()), length(obj_nullmodel_list))
      if(length(id_all[[6]])>=rv_num_cutoff_min_prefilter & length(id_all[[6]])<rv_num_cutoff_max_prefilter){
        #seqSetFilter(genofile,variant.id=id_all[[6]])
        variant.is.in <- id_all[[6]]
        seqSetFilter(genofile,variant.id=variant.is.in)
        
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
        
        if(length(variant.is.in)>=1)
        {
          ## get AF, Missing rate
          AF_AC_Missing <- seqGetAF_AC_Missing(genofile,minor=FALSE,parallel=FALSE)
          REF_AF <- AF_AC_Missing$af
          Missing_rate <- AF_AC_Missing$miss
          ## variant id
          variant.is.in <- seqGetData(genofile, "variant.id")
          variant_maf_cutoff_filter <- ifelse(rare_maf_cutoff<=0.01,0.05,1)
          #seqResetFilter(genofile)
          
          Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.is.in,
                                                sample.id=phenotype.id,
                                                REF_AF=REF_AF,variant_maf_cutoff_filter=variant_maf_cutoff_filter,
                                                Missing_rate=Missing_rate,
                                                rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                                annotation_phred=Anno.Int.PHRED.sub)
          Geno <- Genotype_sp$Geno
          Anno.Int.PHRED.sub <- Genotype_sp$annotation_phred
          rm(Genotype_sp)
          gc()
        } else
        {
          Geno <- NULL
        }
        
        #results_list_enhancer_CAGE <- rep(list(c()), length(obj_nullmodel_list))
        #phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))
        
        if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
        {
          for (num in 1:length(obj_nullmodel_list))
          {
            # id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
            # id.phenotype.num.merge <- data.frame(id.phenotype.num)
            # id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
            # phenotype.id.match <- id.phenotype.num.merge$index
            phenotype.id.match <- id.phenotype.num.list[[num]]
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
                try(pvalues <- STAAR_sp_2norm_SPA(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,
                                                  rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,pi1=pi1_all[num], resid_ratio=ratio_all[num]),silent=silent)
              }else
              {
                try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
              }
            }else
            {
              try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
            }
            
            results_enhancer_CAGE <- c()
            if(inherits(pvalues, "list"))
            {
              results_temp <- c()
              results_temp[3] <- "enhancer_CAGE"
              results_temp[2] <- chr
              results_temp[1] <- as.character(gene_name)
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
              
              results_enhancer_CAGE <- rbind(results_enhancer_CAGE,results_temp)
            }
            
            if(!is.null(results_enhancer_CAGE))
            {
              if(!use_SPA)
              {
                colnames(results_enhancer_CAGE) <- colnames(results_enhancer_CAGE, do.NULL = FALSE, prefix = "col")
                colnames(results_enhancer_CAGE)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
                colnames(results_enhancer_CAGE)[(dim(results_enhancer_CAGE)[2]-1):dim(results_enhancer_CAGE)[2]] <- c("ACAT-O","STAAR-O")
              }else
              {
                colnames(results_enhancer_CAGE) <- colnames(results_enhancer_CAGE, do.NULL = FALSE, prefix = "col")
                colnames(results_enhancer_CAGE)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
                colnames(results_enhancer_CAGE)[dim(results_enhancer_CAGE)[2]] <- c("STAAR-B")
              }
              
              results_list_enhancer_CAGE[[num]] <- results_enhancer_CAGE
            }
          }
        }
        
        #seqResetFilter(genofile)
      }
      ##################################################
      #       Enhancer-DHS
      print('Enhancer-DHS')
      
      results_list_enhancer_DHS <- rep(list(c()), length(obj_nullmodel_list))
      if(length(id_all[[7]])>=rv_num_cutoff_min_prefilter & length(id_all[[7]])<rv_num_cutoff_max_prefilter){
        
        variant.is.in <- id_all[[7]]
        seqSetFilter(genofile,variant.id=variant.is.in)
        
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
        
        if(length(variant.is.in)>=1)
        {
          ## get AF, Missing rate
          AF_AC_Missing <- seqGetAF_AC_Missing(genofile,minor=FALSE,parallel=FALSE)
          REF_AF <- AF_AC_Missing$af
          Missing_rate <- AF_AC_Missing$miss
          ## variant id
          variant.is.in <- seqGetData(genofile, "variant.id")
          variant_maf_cutoff_filter <- ifelse(rare_maf_cutoff<=0.01,0.05,1)
          #seqResetFilter(genofile)
          
          Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.is.in,
                                                sample.id=phenotype.id,
                                                REF_AF=REF_AF,variant_maf_cutoff_filter=variant_maf_cutoff_filter,
                                                Missing_rate=Missing_rate,
                                                rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
                                                annotation_phred=Anno.Int.PHRED.sub)
          Geno <- Genotype_sp$Geno
          Anno.Int.PHRED.sub <- Genotype_sp$annotation_phred
          rm(Genotype_sp)
          gc()
        } else
        {
          Geno <- NULL
        }
        
        #results_list_enhancer_DHS <- rep(list(c()), length(obj_nullmodel_list))
        #phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))
        
        if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
        {
          for (num in 1:length(obj_nullmodel_list))
          {
            # id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
            # id.phenotype.num.merge <- data.frame(id.phenotype.num)
            # id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
            # phenotype.id.match <- id.phenotype.num.merge$index
            phenotype.id.match <- id.phenotype.num.list[[num]]
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
                try(pvalues <- STAAR_sp_2norm_SPA(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,pi1=pi1_all[num], resid_ratio=ratio_all[num]),silent=silent)
              }else
              {
                try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
              }
            }else
            {
              try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
            }
            
            results_enhancer_DHS <- c()
            if(inherits(pvalues, "list"))
            {
              results_temp <- c()
              results_temp[3] <- "enhancer_DHS"
              results_temp[2] <- chr
              results_temp[1] <- as.character(gene_name)
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
              
              results_enhancer_DHS <- rbind(results_enhancer_DHS,results_temp)
            }
            
            if(!is.null(results_enhancer_DHS))
            {
              if(!use_SPA)
              {
                colnames(results_enhancer_DHS) <- colnames(results_enhancer_DHS, do.NULL = FALSE, prefix = "col")
                colnames(results_enhancer_DHS)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
                colnames(results_enhancer_DHS)[(dim(results_enhancer_DHS)[2]-1):dim(results_enhancer_DHS)[2]] <- c("ACAT-O","STAAR-O")
              }else
              {
                colnames(results_enhancer_DHS) <- colnames(results_enhancer_DHS, do.NULL = FALSE, prefix = "col")
                colnames(results_enhancer_DHS)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
                colnames(results_enhancer_DHS)[dim(results_enhancer_DHS)[2]] <- c("STAAR-B")
              }
              
              results_list_enhancer_DHS[[num]] <- results_enhancer_DHS
            }
          }
        }
        
      }
      
      #seqResetFilter(genofile)
      
      ############################################
      #           results
      
      results_noncoding_list <- list(upstream = results_list_upstream,
                                     downstream = results_list_downstream,
                                     UTR = results_list_UTR,
                                     promoter_CAGE = results_list_promoter_CAGE,
                                     promoter_DHS = results_list_promoter_DHS,
                                     enhancer_CAGE = results_list_enhancer_CAGE,
                                     enhancer_DHS = results_list_enhancer_DHS)
      
      seqResetFilter(genofile)
      
      return(results_noncoding_list)
      
      
    }
    
    
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
if(T){
  
  ### modifications: remove resetfiltering
  Genotype_sp_extraction <- function(genofile,variant.id,sample.id,
                                     REF_AF,rare_maf_cutoff=0.01,
                                     variant_maf_cutoff_filter=1,
                                     Missing_rate,Missing_cutoff=0.01,
                                     subset_variants_num=1e3,
                                     rv_num_cutoff_max_prefilter=1e9,
                                     annotation_phred=NULL)
  {
    
    if((length(variant.id) != length(REF_AF)) | (length(Missing_rate) != length(REF_AF)))
    {
      stop(paste0("Lengths don't match for variant.id, REF_AF and Missing_rate!"))
    }
    
    annotation_phred <- as.data.frame(annotation_phred)
    if(dim(annotation_phred)[1] != 0 & length(variant.id) != dim(annotation_phred)[1])
    {
      stop(paste0("Dimensions don't match for genotype and annotation!"))
    }
    
    ## get AF
    ALT_AF <- 1-REF_AF
    MAF <- ifelse(REF_AF>=ALT_AF,ALT_AF,REF_AF)
    is.include <- !((MAF>=variant_maf_cutoff_filter) | (MAF==0) | is.na(MAF) | is.na(Missing_rate))
    
    variant.id.in <- variant.id[is.include]
    ALT_AF.in <- ALT_AF[is.include]
    MAF.in <- MAF[is.include]
    Missing_rate.in <- Missing_rate[is.include]
    
    annotation_phred.sub <- annotation_phred[is.include,,drop=FALSE]
    
    results <- list(Geno=NULL,results_information=NULL,annotation_phred=NULL)
    
    if(length(variant.id.in) == 0 | length(variant.id.in)>=rv_num_cutoff_max_prefilter)
    {
      return(results)
    }
    
    is.dosge <- (ALT_AF.in>0.5)
    is.dosge_alt <- (ALT_AF.in<=0.5) & (MAF.in>=rare_maf_cutoff | Missing_rate.in>=Missing_cutoff)
    is.dosge_sp <- (ALT_AF.in<=0.5) & (MAF.in<rare_maf_cutoff) & (Missing_rate.in<Missing_cutoff)
    
    ## dosage
    if(sum(is.dosge)>=1)
    {
      variant.id.dosage <- variant.id.in[is.dosge]
      subset.num.dosage <- ceiling(length(variant.id.dosage)/subset_variants_num)
      
      MAF.in.dosage <- MAF.in[is.dosge]
      ALT_AF.in.dosage <- ALT_AF.in[is.dosge]
      Missing_rate.in.dosage <- Missing_rate.in[is.dosge]
      
      Geno_dosage <- NULL
      results_dosage <- NULL
      
      for(kk in 1:subset.num.dosage)
      {
        if(kk < subset.num.dosage)
        {
          is.in.dosage <- ((kk-1)*subset_variants_num+1):(kk*subset_variants_num)
        }
        if(kk == subset.num.dosage)
        {
          is.in.dosage <- ((kk-1)*subset_variants_num+1):length(variant.id.dosage)
        }
        
        seqSetFilter(genofile,variant.id=variant.id.dosage[is.in.dosage],sample.id=sample.id)
        
        ## genotype id
        id.genotype <- as.character(seqGetData(genofile,"sample.id"))
        
        id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
        phenotype.id.merge <- data.frame(sample.id)
        phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("sample.id"="id.genotype"))
        id.genotype.match <- phenotype.id.merge$index
        
        Geno_dosage_temp <- seqGetData(genofile, "$dosage")
        Geno_dosage_temp <- as(Geno_dosage_temp,"dgCMatrix")
        Geno_dosage_temp <- Geno_dosage_temp[id.genotype.match,,drop=FALSE]
        Geno_dosage <- cbind(Geno_dosage,Geno_dosage_temp)
        gc()
        
        CHR.dosage <- as.numeric(seqGetData(genofile, "chromosome"))
        position.dosage <- as.numeric(seqGetData(genofile, "position"))
        REF.dosage <- as.character(seqGetData(genofile, "$ref"))
        ALT.dosage <- as.character(seqGetData(genofile, "$alt"))
        variant.id.in.dosage <- seqGetData(genofile, "variant.id")
        
        results_dosage_temp <- data.frame(CHR=CHR.dosage,position=position.dosage,
                                          REF=REF.dosage,ALT=ALT.dosage,
                                          variant.id=variant.id.in.dosage,
                                          stringsAsFactors = FALSE)
        results_dosage <- rbind(results_dosage,results_dosage_temp)
        
        # seqResetFilter(genofile)
      }
      
      annotation_phred.sub_dosage <- annotation_phred.sub[is.dosge,,drop=FALSE]
      ## Reorder
      variant_order.dosage <- match(results_dosage$variant.id, variant.id.dosage)
      annotation_phred.sub_dosage <- annotation_phred.sub_dosage[variant_order.dosage,,drop=FALSE]
      results_dosage <- dplyr::left_join(results_dosage,
                                         data.frame(MAF=MAF.in.dosage,
                                                    ALT_AF=ALT_AF.in.dosage,
                                                    Missing_rate=Missing_rate.in.dosage,
                                                    variant.id=variant.id.dosage),
                                         by="variant.id")
      Geno <- Geno_dosage
      results_info <- results_dosage
      annotation_phred.sub_sort <- annotation_phred.sub_dosage
    }
    
    ## dosage_alt
    if(sum(is.dosge_alt)>=1)
    {
      variant.id.dosage_alt <- variant.id.in[is.dosge_alt]
      subset.num.dosage_alt <- ceiling(length(variant.id.dosage_alt)/subset_variants_num)
      
      MAF.in.dosage_alt <- MAF.in[is.dosge_alt]
      ALT_AF.in.dosage_alt <- ALT_AF.in[is.dosge_alt]
      Missing_rate.in.dosage_alt <- Missing_rate.in[is.dosge_alt]
      
      Geno_dosage_alt <- NULL
      results_dosage_alt <- NULL
      
      for(kk in 1:subset.num.dosage_alt)
      {
        if(kk < subset.num.dosage_alt)
        {
          is.in.dosage_alt <- ((kk-1)*subset_variants_num+1):(kk*subset_variants_num)
        }
        if(kk == subset.num.dosage_alt)
        {
          is.in.dosage_alt <- ((kk-1)*subset_variants_num+1):length(variant.id.dosage_alt)
        }
        
        seqSetFilter(genofile,variant.id=variant.id.dosage_alt[is.in.dosage_alt],sample.id=sample.id)
        
        ## genotype id
        id.genotype <- as.character(seqGetData(genofile,"sample.id"))
        
        id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
        phenotype.id.merge <- data.frame(sample.id)
        phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("sample.id"="id.genotype"))
        id.genotype.match <- phenotype.id.merge$index
        
        Geno_dosage_alt_temp <- seqGetData(genofile, "$dosage_alt")
        Geno_dosage_alt_temp <- as(Geno_dosage_alt_temp,"dgCMatrix")
        Geno_dosage_alt_temp <- Geno_dosage_alt_temp[id.genotype.match,,drop=FALSE]
        Geno_dosage_alt <- cbind(Geno_dosage_alt,Geno_dosage_alt_temp)
        gc()
        
        CHR.dosage_alt <- as.numeric(seqGetData(genofile, "chromosome"))
        position.dosage_alt <- as.numeric(seqGetData(genofile, "position"))
        REF.dosage_alt <- as.character(seqGetData(genofile, "$ref"))
        ALT.dosage_alt <- as.character(seqGetData(genofile, "$alt"))
        variant.id.in.dosage_alt <- seqGetData(genofile, "variant.id")
        
        results_dosage_alt_temp <- data.frame(CHR=CHR.dosage_alt,position=position.dosage_alt,
                                              REF=REF.dosage_alt,ALT=ALT.dosage_alt,
                                              variant.id=variant.id.in.dosage_alt,
                                              stringsAsFactors = FALSE)
        results_dosage_alt <- rbind(results_dosage_alt,results_dosage_alt_temp)
        
        # seqResetFilter(genofile)
      }
      
      annotation_phred.sub_alt <- annotation_phred.sub[is.dosge_alt,,drop=FALSE]
      ## Reorder
      variant_order.dosage_alt <- match(results_dosage_alt$variant.id, variant.id.dosage_alt)
      annotation_phred.sub_alt <- annotation_phred.sub_alt[variant_order.dosage_alt,,drop=FALSE]
      results_dosage_alt <- dplyr::left_join(results_dosage_alt,
                                             data.frame(MAF=MAF.in.dosage_alt,
                                                        ALT_AF=ALT_AF.in.dosage_alt,
                                                        Missing_rate=Missing_rate.in.dosage_alt,
                                                        variant.id=variant.id.dosage_alt),
                                             by="variant.id")
      
      if(sum(is.dosge)>=1)
      {
        Geno <- cbind(Geno_dosage_alt,Geno_dosage)
        results_info <- rbind(results_dosage_alt,results_dosage)
        annotation_phred.sub_sort <- rbind(annotation_phred.sub_alt,annotation_phred.sub_dosage)
      } else
      {
        Geno <- Geno_dosage_alt
        results_info <- results_dosage_alt
        annotation_phred.sub_sort <- annotation_phred.sub_alt
      }
    }
    
    ## dosage_sp
    if(sum(is.dosge_sp)>=1)
    {
      variant.id.dosage_sp <- variant.id.in[is.dosge_sp]
      
      MAF.in.dosage_sp <- MAF.in[is.dosge_sp]
      ALT_AF.in.dosage_sp <- ALT_AF.in[is.dosge_sp]
      Missing_rate.in.dosage_sp <- Missing_rate.in[is.dosge_sp]
      
      Geno_dosage_sp <- NULL
      results_dosage_sp <- NULL
      
      seqSetFilter(genofile,variant.id=variant.id.dosage_sp,sample.id=sample.id)
      
      ## genotype id
      id.genotype <- as.character(seqGetData(genofile,"sample.id"))
      
      id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
      phenotype.id.merge <- data.frame(sample.id)
      phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("sample.id"="id.genotype"))
      id.genotype.match <- phenotype.id.merge$index
      
      Geno_dosage_sp <- seqGetData(genofile, "$dosage_sp")
      Geno_dosage_sp <- Geno_dosage_sp[id.genotype.match,,drop=FALSE]
      
      CHR.dosage_sp <- as.numeric(seqGetData(genofile, "chromosome"))
      position.dosage_sp <- as.numeric(seqGetData(genofile, "position"))
      REF.dosage_sp <- as.character(seqGetData(genofile, "$ref"))
      ALT.dosage_sp <- as.character(seqGetData(genofile, "$alt"))
      variant.id.in.dosage_sp <- seqGetData(genofile, "variant.id")
      
      annotation_phred.sub_sp <- annotation_phred.sub[is.dosge_sp,,drop=FALSE]
      results_dosage_sp <- data.frame(CHR=CHR.dosage_sp,position=position.dosage_sp,
                                      REF=REF.dosage_sp,ALT=ALT.dosage_sp,
                                      variant.id=variant.id.in.dosage_sp,
                                      stringsAsFactors = FALSE)
      
      ## Reorder
      variant_order.dosage_sp <- match(results_dosage_sp$variant.id, variant.id.dosage_sp)
      annotation_phred.sub_sp <- annotation_phred.sub_sp[variant_order.dosage_sp,,drop=FALSE]
      results_dosage_sp <- dplyr::left_join(results_dosage_sp,
                                            data.frame(MAF=MAF.in.dosage_sp,
                                                       ALT_AF=ALT_AF.in.dosage_sp,
                                                       Missing_rate=Missing_rate.in.dosage_sp,
                                                       variant.id=variant.id.dosage_sp),
                                            by="variant.id")
      
      # seqResetFilter(genofile)
      
      if(sum(is.dosge)>=1 | sum(is.dosge_alt)>=1)
      {
        Geno <- cbind(Geno_dosage_sp,Geno)
        results_info <- rbind(results_dosage_sp,results_info)
        annotation_phred.sub_sort <- rbind(annotation_phred.sub_sp,annotation_phred.sub_sort)
      } else
      {
        Geno <- Geno_dosage_sp
        results_info <- results_dosage_sp
        annotation_phred.sub_sort <- annotation_phred.sub_sp
      }
    }
    
    if(dim(annotation_phred.sub_sort)[2]==0)
    {
      annotation_phred.sub_sort <- NULL
    }
    
    results <- list(Geno=Geno,results_information=results_info,annotation_phred=annotation_phred.sub_sort)
    return(results)
  }
  
}

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
## Directory of precomputed per-gene noncoding variant-id sets
##   (<dx_noncoding_id_dir>/<arrayid>/Noncoding_<kk>.Rdata provides the list `id_all00`):
dx_noncoding_id_dir <- "/noncoding_id/"
## Project directory to upload per-job results into (<dx_results_dir>/<arrayid>/):
dx_results_dir <- "/results/noncoding/"
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
output_path <- ".//step4_output//"
output_file_name <- "Noncoding_"
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
  #   '25003','25005','25007','25009','25010')
  # 'postcentral',
  # 'transversetemporal')
  for(ll in field_all){
    #obj.STAAR_500K_26566
    
    #for(l in 1:24){
    system(paste0("dx download -f ", dx_null_model_dir, "obj.STAAR_500K_",ll,"_norm_pseudo_rf100_v5.rda"))
    load(paste0("obj.STAAR_500K_",ll,"_norm_pseudo_rf100_v5.rda"))
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
    ## precomputed per-gene noncoding variant-id sets (provides the list `id_all00`)
    system(paste0("dx download -f ", dx_noncoding_id_dir, arrayid, "/Noncoding_", kk, ".Rdata --overwrite"))
    load(paste0("Noncoding_",kk,".Rdata"))
    start_time <- proc.time()
    results <- noncoding_PheWAS_hwe_2norm_SPA(chr=chr,gene_name=gene_name,genofile=genofile,
                                              obj_nullmodel_list=testlist, genes=genes_info,
                                              rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
                                              rv_num_cutoff_max_prefilter=2000,
                                              QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                              Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                              Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
                                              id_all=id_all00, pi1_all=pi1_all, ratio_all=ratio_all)
    end_time <- proc.time()
    print(end_time-start_time)
    save(results,field_all,file=paste0(output_path,output_file_name,kk,".Rdata"))
    system(paste0("dx mkdir -p ", dx_results_dir, arrayid, "/"))
    system(paste0("dx upload ", output_path, output_file_name, kk, ".Rdata --path ", dx_results_dir, arrayid, "/"))
  }
  library(parallel)
  res <- mclapply(sub_seq_id, do_kk, mc.cores=2)
  seqClose(genofile)
}










