# Complete-data multinomial simulation with fixed-effect priors.
# Usage: Rscript NEW_multiDG_multifit_FX_comp.R <case_id> <replicate_id>

suppressPackageStartupMessages({
  library(R2jags)
  library(LaplacesDemon)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript NEW_multiDG_multifit_FX_comp.R <case_id: 1|2> <replicate_id: positive integer>")
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

for (caseid in case_ids){
  for (replicate_id in replicate_ids){
    sp <- caseid
    
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
    
    ########################################################################################
    
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
    
    
    ## Assemble the JAGS data list for both treatment groups.
    
    res <- list()
    res$OR15_true <- OR15_true
    res$OR35_true <- OR35_true
    res$p15A_AE <- p15A_AE_true
    res$p35A_AE <- p35A_AE_true
    res$p15B_AE <- p15B_AE_true
    res$p35B_AE <- p35B_AE_true
    
    data <- list(J1_A=nrow(dt_A), 
                 J1_B=nrow(dt_B), 
                 n_study=n_study, n_AE=n_AE,
                 N_A=dt_A$N, 
                 Y_A=dt_A[,c("G0", "G12", "G35")],
                 N_B=dt_B$N, 
                 Y_B=dt_B[,c("G0", "G12", "G35")],
                 AE_A=dt_A$AE, 
                 study_A=dt_A$study,
                 AE_B=dt_B$AE, 
                 study_B=dt_B$study )
    
    system.time(fit1 <- try(jags.parallel(data=data, 
                                          parameters.to.save=c("OR15", 
                                                               "OR35",
                                                               "p15A_AE", "p35A_AE", "p15B_AE", "p35B_AE"),
                                          model.file=file.path(script_dir, "multinomf_ORsel_FX_comp.txt"),
                                          n.chains = 2, n.iter = 100000, n.burnin = 50000, n.thin = 5,
                                          DIC=TRUE, jags.seed = 2025), silent=TRUE))
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
      "NEW_multiDG_multifit_FX_comp_case", caseid,
      "_round", replicate_id, ".RData"
    ))
    save(res, file=output_file)
    
    
  
    message("Completed case ", caseid, ", replicate ", replicate_id,
            "; saved to ", output_file)
  }
}

