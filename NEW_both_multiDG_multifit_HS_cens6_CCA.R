# Available-case benchmark for data generated with 60% censoring.
# Usage: Rscript NEW_both_multiDG_multifit_HS_cens6_CCA.R <case_id> <replicate_id>

suppressPackageStartupMessages({
  library(R2jags)
  library(LaplacesDemon)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript NEW_both_multiDG_multifit_HS_cens6_CCA.R <case_id: 1|2> <replicate_id: positive integer>")
}
case_ids <- suppressWarnings(as.integer(args[1]))
replicate_ids <- suppressWarnings(as.integer(args[2]))
if (is.na(case_ids) || !case_ids %in% 1:2 || is.na(replicate_ids) || replicate_ids < 1L) {
  stop("case_id must be 1 or 2, and replicate_id must be a positive integer.")
}

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
results_dir <- file.path(script_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

for (replicate_id in replicate_ids){
  for (caseid in case_ids){
    sp <- caseid
    
    cr15 <- cr35 <- 0.6
    N <- 100
    n_AE <- 85
    
    ## AE-specific baseline probabilities
    p15A_AE_true <- c(rep(0.4, 20), rep(0.2, 30), rep(0.05, 20), rep(0.4, 15))
    p35A_AE_true <- c(rep(0.1, 20), rep(0.08, 30), rep(0.02, 20), rep(0.02, 15))
    stopifnot(length(p15A_AE_true) == n_AE,
              length(p35A_AE_true) == n_AE,
              all(p35A_AE_true <= p15A_AE_true))
    
    OR15_pos <- 2
    OR15_neg <- 1/2
    
    OR35_pos <- 3
    tau_study1 <- tau_study2 <- 0.5
    tau_eta1 <- tau_eta2 <- 0.2
    
    ## Sparse truth: four G15 ORs of 2, three G15 ORs of 1/2, and one
    ## G35 OR of 3 nested within a G15 OR of 1/2.
    signals_15 <- c(rep(OR15_neg, 3), rep(OR15_pos, 4), rep(1, n_AE - 7))
    signals_35 <- c(OR35_pos, rep(1, n_AE - 1))
    stopifnot(length(signals_15) == n_AE, length(signals_35) == n_AE)
    if (sp==1){ 
      study_id_A <- c(1:2, 3:(3+29-1))
      study_id_B <- c(1:2, 32:(32+13-1))
      n_study <- 2 + 29 + 13 #44
      
    } else if (sp==2){ 
      study_id_A <- c(1:2, 3:(3+21-1))
      study_id_B <- c(1:2, 23:(23+21-1))
      n_study <- 2 + 21 + 21 #44
    }
    
    dt_A <- data.frame(study=rep(study_id_A, each=n_AE),
                       AE=rep(1:n_AE, length(study_id_A)))
    dt_B <- data.frame(study=rep(study_id_B, each=n_AE),
                       AE=rep(1:n_AE, length(study_id_B)))
    
    dt_A$N <- N
    dt_B$N <- N
    
    set.seed(2025 + 123 * replicate_id)
    shuffle_index <- sample.int(n_AE)
    OR15_true <- signals_15[shuffle_index]
    OR35_true <- signals_35[shuffle_index]
    u_study1 <- rnorm(n_study, mean=0, sd=tau_study1)
    u_study2 <- rnorm(n_study, mean=0, sd=tau_study2)
    u_eta1 <- matrix(NA, nrow=n_study, ncol=n_AE)
    u_eta2 <- matrix(NA, nrow=n_study, ncol=n_AE)
    
    for (i in 1:n_study){
      for (j in 1:n_AE){
        u_eta1[i,j] <- u_study1[i] + rnorm(1, mean=0, sd=tau_eta1)
        u_eta2[i,j] <- u_study2[i] + rnorm(1, mean=0, sd=tau_eta2)
      }
    }
    
    mu1A_AE <- log((p15A_AE_true-p35A_AE_true)/(1-p15A_AE_true))
    mu2A_AE <- log(p35A_AE_true/(1-p15A_AE_true))
    
    p15A <- NULL
    p35A <- NULL
    p15B <- NULL
    p35B <- NULL
    count_A <- NULL
    count_B <- NULL
    for (i in 1:nrow(dt_A)){
      r1A <- exp(mu1A_AE[dt_A$AE[i]] + u_eta1[dt_A$study[i], dt_A$AE[i]])
      r2A <- exp(mu2A_AE[dt_A$AE[i]] + u_eta2[dt_A$study[i], dt_A$AE[i]])
      p15A[i] <- (r1A+r2A)/(1+r1A+r2A)
      p35A[i] <- r2A/(1+r1A+r2A)
      count_A <- rbind(count_A, as.numeric(rmultinom(n=1, size=N, prob=c(1-p15A[i], p15A[i]-p35A[i], p35A[i]))))
    }
    for (i in 1:nrow(dt_B)){
      r1A_mean <- exp(mu1A_AE[dt_B$AE[i]] + u_eta1[dt_B$study[i], dt_B$AE[i]])
      r2A_mean <- exp(mu2A_AE[dt_B$AE[i]] + u_eta2[dt_B$study[i], dt_B$AE[i]])
      p15A_mean <- (r1A_mean+r2A_mean)/(1+r1A_mean+r2A_mean)
      p35A_mean <- r2A_mean/(1+r1A_mean+r2A_mean)
      
      p15B[i] <- invlogit( logit(p15A_mean) + log(OR15_true[dt_B$AE[i]]) )
      p35B[i] <- min( invlogit( logit(p35A_mean) + log(OR35_true[dt_B$AE[i]]) ), p15B[i] )
      count_B <- rbind(count_B, as.numeric(rmultinom(n=1, size=N, prob=c(1-p15B[i], p15B[i]-p35B[i], p35B[i]))))
    }
    
    p15B_AE_true <- NULL
    p35B_AE_true <- NULL
    for (k in 1:n_AE){
      p15B_AE_true[k] <- invlogit( logit(p15A_AE_true[k]) + log(OR15_true[k]) )
      p35B_AE_true[k] <- min( invlogit( logit(p35A_AE_true[k]) + log(OR35_true[k]) ), p15B_AE_true[k] )
    }
    
    dt_A$G0 <- count_A[,1]
    dt_A$G12 <- count_A[,2]
    dt_A$G35 <- count_A[,3]
    dt_A$G15 <- N-dt_A$G0
    
    dt_B$G0 <- count_B[,1]
    dt_B$G12 <- count_B[,2]
    dt_B$G35 <- count_B[,3]
    dt_B$G15 <- N-dt_B$G0
    
    dt <- rbind(dt_A, dt_B)
    dt$drug <- rep(c(1,2), times=c(nrow(dt_A), nrow(dt_B)))
    
    compute_auc_from_credible <- function(posterior_draws, truth,
                                          c_grid = seq(0.50, 0.995, by = 0.005)) {
      
      stopifnot(is.matrix(posterior_draws), length(truth) == ncol(posterior_draws))
      
      ci_equal_tailed <- function(x, level) {
        a <- (1 - level) / 2
        c(quantile(x, a, names = FALSE), quantile(x, 1 - a, names = FALSE))
      }
      ci_fun <- ci_equal_tailed
      
      is_signal <- as.logical(truth!=1)
      
      stat_counts <- function(sel) {
        TP <- sum(sel & is_signal)
        FP <- sum(sel & !is_signal)
        FN <- sum(!sel & is_signal)
        TN <- sum(!sel & !is_signal)
        list(TP = TP, FP = FP, FN = FN, TN = TN)
      }
      
      sensitivity <- function(TP, FN) if ((TP + FN) == 0) NA_real_ else TP / (TP + FN)
      specificity <- function(TN, FP) if ((TN + FP) == 0) NA_real_ else TN / (TN + FP)
      fdr        <- function(TP, FP) if ((TP + FP) == 0) 0 else FP / (TP + FP)
      fwer       <- function(FP) as.numeric(FP > 0)
      
      curve_df <- do.call(rbind, lapply(c_grid, function(clev) {
        ci_mat <- t(apply(posterior_draws, 2, ci_fun, level = clev))
        sel <- (ci_mat[, 1] > 1) | (ci_mat[, 2] < 1)
        
        cnt <- stat_counts(sel)
        sens <- sensitivity(cnt$TP, cnt$FN)
        spec <- specificity(cnt$TN, cnt$FP)
        out <- data.frame(
          c_level = clev,
          sens = sens,
          FDR = fdr(cnt$TP, cnt$FP),
          FWER = fwer(cnt$FP),
          one_minus_spec = 1 - spec,
          TP = cnt$TP, FP = cnt$FP, FN = cnt$FN, TN = cnt$TN,
          selected = cnt$TP + cnt$FP
        )
        out
      }))
      
      trapz_auc <- function(x, y, pad_to01 = TRUE) {
        ok <- is.finite(x) & is.finite(y)
        x <- x[ok]; y <- y[ok]
        if (length(x) < 2) return(NA_real_)
        o <- order(x); x <- x[o]; y <- y[o]
        
        if (pad_to01) {
          if (min(x) > 0) { x <- c(0, x); y <- c(y[1], y) }
          if (max(x) < 1) { x <- c(x, 1); y <- c(y, y[length(y)]) }
        }
        sum(diff(x) * (head(y, -1) + tail(y, 1)) / 2)
      }
      
      auc_sens_fdr  <- trapz_auc(curve_df$FDR,  curve_df$sens,  pad_to01 = TRUE)   # Sens vs FDR
      auc_sens_fwer <- trapz_auc(curve_df$FWER, curve_df$sens,  pad_to01 = TRUE)   # Sens vs FWER
      auc_roc       <- trapz_auc(curve_df$one_minus_spec, curve_df$sens, pad_to01 = TRUE)  #ROC AUC
      
      list(
        curves = curve_df,
        auc = list(
          AUC_Sensitivity_vs_FDR  = auc_sens_fdr,
          AUC_Sensitivity_vs_FWER = auc_sens_fwer,
          AUC_ROC_Sens_vs_1MinusSpec = auc_roc
        )
      )
    }
    
    
    ############ Simulate censoring within each study ########################################################
    dt2 <- NULL
    for (i in 3:n_study){ ### Single-arm trials
      dt_temp <- dt[dt$study==i,]
      n_censored15 <- floor(cr15 * nrow(dt_temp))
      n_censored35 <- floor(cr35 * nrow(dt_temp))
      order15 <- order(dt_temp$G15, decreasing=FALSE)
      order35 <- order(dt_temp$G35, decreasing=FALSE)
      cut15 <- dt_temp$G15[order15[n_censored15]] # Censored counts are <= this cutoff.
      cut35 <- dt_temp$G35[order35[n_censored35]]
      dt2_temp <- dt_temp
      dt2_temp$G15[order15[seq_len(n_censored15)]] <- NA
      dt2_temp$cut15 <- cut15
      dt2_temp$G35[order35[seq_len(n_censored35)]] <- NA
      dt2_temp$cut35 <- cut35
      
      dt2 <- rbind(dt2, dt2_temp)
    }
    
    for (i in 1:2){ ### Two-arm randomized trials
      dt_temp1 <- dt[dt$study==i & dt$drug==1,]
      dt_temp2 <- dt[dt$study==i & dt$drug==2,]
      max_G15 <- apply(cbind(dt_temp1$G15, dt_temp2$G15), 1, max)
      max_G35 <- apply(cbind(dt_temp1$G35, dt_temp2$G35), 1, max)
      n_censored15 <- floor(cr15 * length(max_G15))
      n_censored35 <- floor(cr35 * length(max_G35))
      cut15 <- sort(max_G15)[n_censored15]
      cut35 <- sort(max_G35)[n_censored35]
      dt2_temp1 <- dt_temp1
      dt2_temp2 <- dt_temp2
      
      dt2_temp1$G15[order(dt2_temp1$G15)[seq_len(n_censored15)]] <- NA
      dt2_temp1$cut15 <- cut15
      dt2_temp1$G35[order(dt2_temp1$G35)[seq_len(n_censored35)]] <- NA
      dt2_temp1$cut35 <- cut35
      
      dt2_temp2$G15[order(dt2_temp2$G15)[seq_len(n_censored15)]] <- NA
      dt2_temp2$cut15 <- cut15
      dt2_temp2$G35[order(dt2_temp2$G35)[seq_len(n_censored35)]] <- NA
      dt2_temp2$cut35 <- cut35
      
      dt2 <- rbind(dt2, dt2_temp1, dt2_temp2)
    }
    
    dt2_11 <- dt2[!is.na(dt2$G15) & !is.na(dt2$G35),]
    dt2_10 <- dt2[!is.na(dt2$G15) & is.na(dt2$G35),]
    dt2_01 <- dt2[is.na(dt2$G15) & !is.na(dt2$G35),]
    dt2_00 <- dt2[is.na(dt2$G15) & is.na(dt2$G35),]
    
    ## When only G35 is censored, its cutoff cannot exceed the observed G15 count.
    dt2_10$cut35 <- ifelse(dt2_10$cut35>dt2_10$G15, dt2_10$G15, dt2_10$cut35)
    
    dt2_11_A <- dt2_11[dt2_11$drug==1,]
    dt2_10_A <- dt2_10[dt2_10$drug==1,]
    dt2_01_A <- dt2_01[dt2_01$drug==1,]
    dt2_00_A <- dt2_00[dt2_00$drug==1,]
    dt2_11_B <- dt2_11[dt2_11$drug==2,]
    dt2_10_B <- dt2_10[dt2_10$drug==2,]
    dt2_01_B <- dt2_01[dt2_01$drug==2,]
    dt2_00_B <- dt2_00[dt2_00$drug==2,]
    
    ########################################################################################
    
    ## Assemble data by censoring pattern and treatment group.
    
    res <- list()
    res$OR15_true <- OR15_true
    res$OR35_true <- OR35_true
    res$p15A_AE <- p15A_AE_true
    res$p35A_AE <- p35A_AE_true
    res$p15B_AE <- p15B_AE_true
    res$p35B_AE <- p35B_AE_true
    
    data <- list(J1_A=nrow(dt2_11_A), J2_A=nrow(dt2_10_A), J3_A=nrow(dt2_01_A), 
                 J1_B=nrow(dt2_11_B), J2_B=nrow(dt2_10_B), J3_B=nrow(dt2_01_B), 
                 n_study=n_study, n_AE=n_AE,
                 N_A=c(dt2_11_A$N, dt2_10_A$N, dt2_01_A$N), 
                 Y_A=rbind(dt2_11_A[,c("G0", "G12", "G35")], dt2_10_A[,c("G0", "G12", "G35")], dt2_01_A[,c("G0", "G12", "G35")]), 
                 N_B=c(dt2_11_B$N, dt2_10_B$N, dt2_01_B$N), 
                 Y_B=rbind(dt2_11_B[,c("G0", "G12", "G35")], dt2_10_B[,c("G0", "G12", "G35")], dt2_01_B[,c("G0", "G12", "G35")]),
                 AE_A=c(dt2_11_A$AE, dt2_10_A$AE, dt2_01_A$AE),
                 study_A=c(dt2_11_A$study, dt2_10_A$study, dt2_01_A$study),
                 AE_B=c(dt2_11_B$AE, dt2_10_B$AE, dt2_01_B$AE),
                 study_B=c(dt2_11_B$study, dt2_10_B$study, dt2_01_B$study) )
    
    
    inits <- function(){list(mu_study1 = -0.5, mu_study2 = -1,
                             u_AE1 = rep(-0.5, 85), u_AE2 = rep(-1, 85),
                             tau_INT1 = 0.5, tau_INT2 = 0.5,
                             tau_study1 = 0.5, tau_study2 = 0.5, 
                             tau15 = 1, tau35= 1, 
                             lambda15=rep(1, 85), lambda35=rep(1, 85))}
    system.time(fit1 <- try(jags.parallel(data=data, 
                                          inits = inits, 
                                          parameters.to.save=c("OR15", "OR35",
                                                               "p15A_AE", "p35A_AE", "p15B_AE", "p35B_AE"),
                                          model.file=file.path(script_dir, "multinomf_ORsel_HS_cens_CCA.txt"),
                                          n.chains = 2, n.iter = 100000, n.burnin = 50000, n.thin = 5,
                                          DIC=FALSE, jags.seed = 2025), silent=TRUE))
    if (inherits(fit1, "try-error")) stop("JAGS fitting failed: ", as.character(fit1))
    if (!inherits(fit1, "try-error")){
      res15 <- compute_auc_from_credible(posterior_draws = as.matrix(fit1$BUGSoutput$sims.list$OR15), truth = OR15_true)
      res$dt_curve15 <- res15$curves
      res$AUC15 <- res15$auc
      
      res35 <- compute_auc_from_credible(posterior_draws = as.matrix(fit1$BUGSoutput$sims.list$OR35), truth = OR35_true)
      res$dt_curve35 <- res35$curves
      res$AUC35 <- res35$auc
      
      res$OR15 <- fit1$BUGSoutput$summary[which(rownames(fit1$BUGSoutput$summary)=="OR15[1]"):which(rownames(fit1$BUGSoutput$summary)=="OR15[85]"),
                                          c(5,3,7,8,9)]
      res$OR35 <- fit1$BUGSoutput$summary[which(rownames(fit1$BUGSoutput$summary)=="OR35[1]"):which(rownames(fit1$BUGSoutput$summary)=="OR35[85]"),
                                          c(5,3,7,8,9)]
      
      res$p15A <- fit1$BUGSoutput$summary[which(rownames(fit1$BUGSoutput$summary)=="p15A_AE[1]"):which(rownames(fit1$BUGSoutput$summary)=="p15A_AE[85]"),
                                          c(5,3,7,8,9)]
      res$p15B <- fit1$BUGSoutput$summary[which(rownames(fit1$BUGSoutput$summary)=="p15B_AE[1]"):which(rownames(fit1$BUGSoutput$summary)=="p15B_AE[85]"),
                                          c(5,3,7,8,9)]
      res$p35A <- fit1$BUGSoutput$summary[which(rownames(fit1$BUGSoutput$summary)=="p35A_AE[1]"):which(rownames(fit1$BUGSoutput$summary)=="p35A_AE[85]"),
                                          c(5,3,7,8,9)]
      res$p35B <- fit1$BUGSoutput$summary[which(rownames(fit1$BUGSoutput$summary)=="p35B_AE[1]"):which(rownames(fit1$BUGSoutput$summary)=="p35B_AE[85]"),
                                          c(5,3,7,8,9)]
      
    }
    
    output_file <- file.path(results_dir, paste0(
      "NEW_both_multiDG_multifit_HS_c6_CCA_case", caseid,
      "_round", replicate_id, ".RData"
    ))
    save(res, file=output_file)
    
    
    message("Completed case ", caseid, ", replicate ", replicate_id,
            "; saved to ", output_file)
  }
}

