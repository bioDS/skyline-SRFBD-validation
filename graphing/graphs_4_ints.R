######################################################
# Coverage / accuracy validation for skyline SRFBD sims
######################################################
library(ggplot2)
library(coda)
library(data.table)
library(Hmisc)   # needed for mean_cl_boot
options(digits = 7)

rm(list = ls())

setwd("/home/ket581/skyline/SUMMER/skyline-SRFBD-validation/graphing")
fig_dir <- "figures_sims_with_morph_and_dna/"
log_dir <- "../logs_with_morph_and_dna/"
dir.create(fig_dir, showWarnings = FALSE)

system("bash ../move_logs.sh")

true.rates <- fread("/home/ket581/skyline/SUMMER/skyline-SRFBD-validation/true_rates_morph_and_dna_sims.csv",
                    sep = ",", fill = TRUE)

n       <- 200
ints    <- 4
burnin  <- 0.1
min_ess <- 100
probs   <- seq(0, 1, 0.05)

######################################################
# Parameter definitions: one row per estimated quantity
######################################################
params <- data.frame(
  name     = c(paste0("diversificationRate", 1:ints),
               paste0("turnover", 1:ints),
               paste0("samplingProportion", 1:ints),
               "samplingAtPresentProb", "mrca", "origin"),
  log_col  = c(paste0("netDiversification.", 1:ints),
               paste0("turnOver.", 1:ints),
               paste0("samplingProportion.", 1:ints),
               "rho", "TreeHeight", "origin"),
  true_col = c(paste0("div_rate.", 1:ints, "."),
               paste0("turnover.", 1:ints, "."),
               paste0("sampling_prop.", 1:ints, "."),
               "rho", "mrca", "origin"),
  label    = c(paste("Diversification rate", 1:ints),
               paste("Turnover", 1:ints),
               paste("Sampling proportion", 1:ints),
               "Sampling at present prob", "MRCA", "Origin"),
  stringsAsFactors = FALSE
)

missing_cols <- setdiff(params$true_col, names(true.rates))
if (length(missing_cols))
  stop("Missing in true rates file: ", paste(missing_cols, collapse = ", "))

######################################################
# Helper functions
######################################################
test_vals <- function(prob, true, estimated) {
  hpd <- HPDinterval(as.mcmc(estimated), prob = prob)
  as.numeric(true >= hpd[1] & true <= hpd[2])
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# all summary statistics for one parameter in one simulation
summarise_draws <- function(true, draws) {
  hpd <- HPDinterval(as.mcmc(draws))
  coverage <- sapply(probs, test_vals, true = true, estimated = draws)
  names(coverage) <- paste0("test_", round(probs * 100))
  c(true           = true,
    estimated      = median(draws),
    upper          = unname(quantile(draws, 0.975)),
    lower          = unname(quantile(draws, 0.025)),
    rel.err.median = abs(median(draws) - true) / true,
    rel.err.mean   = abs(mean(draws) - true) / true,
    rel.err.mode   = abs(Mode(draws) - true) / true,
    cv             = sqrt(exp(sd(log(draws))^2) - 1),
    hpd.lower      = hpd[1],
    hpd.upper      = hpd[2],
    hpd.rel.width  = (hpd[2] - hpd[1]) / true,
    coverage)
}

######################################################
# Process each simulation
######################################################
results     <- list()
post_ess    <- rep(NA_real_, n)
remove_rows <- c()

for (i in 1:n) {
  cat("i =", i, " ")
  fl <- paste0(log_dir, i, ".log")
  if (!file.exists(fl)) {
    remove_rows <- c(remove_rows, i)
    cat("FILE DOES NOT EXIST", i, "\n")
    next
  }

  t <- read.table(fl, header = TRUE, sep = "\t")
  t <- t[-seq_len(ceiling(nrow(t) * burnin)), ]
  if (nrow(t) == 0) {
    remove_rows <- c(remove_rows, i)
    next
  }

  post_ess[i] <- as.numeric(effectiveSize(as.mcmc(t$posterior)))
  if (post_ess[i] < min_ess) {
    remove_rows <- c(remove_rows, i)
    cat("post ess for", i, "is too small (", post_ess[i], ")\n")
    next
  }

  truth <- as.numeric(unlist(true.rates[i, params$true_col, with = FALSE]))
  rows  <- lapply(seq_len(nrow(params)), function(k)
    summarise_draws(truth[k], t[[params$log_col[k]]]))

  results[[length(results) + 1]] <- data.frame(sim = i,
                                               parameter = params$name,
                                               do.call(rbind, rows))
}
cat("\n")

results <- do.call(rbind, results)
valid   <- length(unique(results$sim))
removed <- length(remove_rows)
cat("Valid runs:", valid, " Removed:", removed, "\n")

write.csv(results, file = paste0(fig_dir, "per_sim_summaries.csv"), row.names = FALSE)
save.image(file = paste0(fig_dir, "validationPlotEnvironment.RData"))

######################################################
# Plots
######################################################
test_cols <- paste0("test_", round(probs * 100))

plot_param <- function(name, label) {
  sub <- results[results$parameter == name, ]

  d  <- data.frame(p = probs, vals = colMeans(sub[, test_cols]))
  dd <- data.frame(p = rep(probs, each = nrow(sub)), vals = unlist(sub[, test_cols]))

  p_qq <- ggplot(d) +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    geom_point(aes(x = p, y = vals), size = 2) +
    stat_summary(data = dd, aes(x = p, y = vals), fun.data = "mean_cl_boot",
                 fun.args = list(conf.int = .95, B = 1000), colour = "blue") +
    xlab("Credible interval width") + ylab("Truth recovery proportion") +
    theme_minimal() + ggtitle(label)

  lims <- range(c(sub$true, sub$lower, sub$upper))
  p_est <- ggplot(sub) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    geom_errorbar(aes(x = true, ymin = lower, ymax = upper), colour = "grey", width = 0.01) +
    geom_point(aes(x = true, y = estimated), size = 2) +
    theme_minimal() + ggtitle(label)
  p_log <- p_est + scale_x_log10(limits = lims) + scale_y_log10(limits = lims)

  ggsave(paste0(fig_dir, "qq_", name, ".pdf"),       p_qq,  width = 6, height = 5)
  ggsave(paste0(fig_dir, name, ".pdf"),              p_est, width = 6, height = 5)
  ggsave(paste0(fig_dir, name, "_logScale.pdf"),     p_log, width = 6, height = 5)

  list(qq = p_qq, est = p_est, log = p_log)
}

plots <- setNames(lapply(seq_len(nrow(params)), function(k)
  plot_param(params$name[k], params$label[k])), params$name)

# view any plot interactively, e.g. plots$turnover2$qq

######################################################
# Summary tables
######################################################
param_mean <- function(col)
  sapply(params$name, function(p) mean(results[[col]][results$parameter == p]))

tables <- list(
  HPD_test             = c("HPD coverage",                     "test_95"),
  HPD_rel_width_score_ = c("Average relative HPD width",       "hpd.rel.width"),
  CV_score             = c("Average Coefficient of Variation", "cv"),
  rel_error_median     = c("Relative Error at Medians",        "rel.err.median"),
  rel_error_mean       = c("Relative Error at Means",          "rel.err.mean")
)

for (nm in names(tables)) {
  out <- data.frame(Parameter = params$name, value = param_mean(tables[[nm]][2]))
  names(out)[2] <- tables[[nm]][1]
  write.csv(out, file = paste0(fig_dir, nm, ".csv"), row.names = FALSE)
}