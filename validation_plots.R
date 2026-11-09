######################################################
# Coverage / accuracy validation for skyline SRFBD sims
# Based on code from Stolz et al. 2024: https://github.com/jugne/sRanges-material/blob/714cde4141bedb997a52ec1169b480b51823198e/validation/Figures2_A10_validationPlots.R#L4
# and https://github.com/jugne/sRanges-material/blob/714cde4141bedb997a52ec1169b480b51823198e/validation/FiguresA11_A12_validationPlotsRUV.R#L4

# Edits made by Kate Truman with Claude and ChatGPT
#
# Usage (from anywhere):
#   Rscript graphing/validation_plots.R 3      # analyse runs/run_3
#
# Reads   runs/run_NUM/true_rates.csv
#         runs/run_NUM/inf/<i>/*.log   (BEAST logs, read in place)
# Writes  runs/run_NUM/figures/
#           coverage plots, estimate plots, summary tables, and
#           rank uniformity validation (RUV / SBC) ECDF plots
# The number of skyline intervals is taken from true_rates.csv.
######################################################
library(ggplot2)
library(coda)
library(data.table)
library(Hmisc)   # needed for mean_cl_boot
library(SBC)     # ECDF plots for rank uniformity validation
options(digits = 7)

rm(list = ls())

base_dir <- "/home/ket581/skyline/SUMMER/skyline-SRFBD-validation/"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript validation_plots.R RUN_NUMBER")
run_num <- args[1]

run_dir <- paste0(base_dir, "runs/run_", run_num, "/")
inf_dir <- paste0(run_dir, "inf/")
fig_dir <- paste0(run_dir, "figures/")
if (!dir.exists(inf_dir)) stop("No inference folder at ", inf_dir)
dir.create(fig_dir, showWarnings = FALSE)

true.rates <- fread(paste0(run_dir, "true_rates.csv"), sep = ",", fill = TRUE)

n       <- nrow(true.rates)
ints    <- length(grep("^div_rate_[0-9]+$", names(true.rates)))
burnin  <- 0.1
min_ess <- 200
n_keep  <- 200   # posterior draws kept (evenly thinned) per simulation for rank calculation
probs   <- seq(0, 1, 0.05)
cat("run_", run_num, ": ", n, " simulations, ", ints, " interval(s)\n", sep = "")

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
  true_col = c(paste0("div_rate_", 1:ints),
               paste0("turnover_", 1:ints),
               paste0("sampling_prop_", 1:ints),
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

# BEAST may log a one-element parameter without the ".1" suffix
# (e.g. netDiversification instead of netDiversification.1); accept either
resolve_log_cols <- function(log_names) {
  sapply(params$log_col, function(col) {
    if (col %in% log_names) return(col)
    bare <- sub("\\.1$", "", col)
    if (ints == 1 && bare %in% log_names) return(bare)
    stop("Column ", col, " not found in log. Log columns: ", paste(log_names, collapse = ", "))
  }, USE.NAMES = FALSE)
}

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
ranks       <- list()
post_ess    <- rep(NA_real_, n)
remove_rows <- c()
log_cols    <- NULL

for (i in true.rates$sim) {
  cat("i =", i, " ")
  fl <- Sys.glob(paste0(inf_dir, i, "/*.log"))
  if (length(fl) == 0) {
    remove_rows <- c(remove_rows, i)
    cat("NO LOG IN", paste0(inf_dir, i), "\n")
    next
  }
  if (length(fl) > 1) {
    remove_rows <- c(remove_rows, i)
    cat("MULTIPLE LOGS IN", paste0(inf_dir, i), "- skipped:", basename(fl), "\n")
    next
  }

  t <- read.table(fl, header = TRUE, sep = "\t")
  if (is.null(log_cols)) log_cols <- resolve_log_cols(names(t))
  t <- t[-seq_len(ceiling(nrow(t) * burnin)), ]
  if (nrow(t) < n_keep) {
    remove_rows <- c(remove_rows, i)
    cat("only", nrow(t), "samples after burnin for", i, "( need", n_keep, ")\n")
    next
  }

  truth <- as.numeric(unlist(true.rates[sim == i, params$true_col, with = FALSE]))

  # rank of the true value among n_keep evenly thinned posterior draws (0 .. n_keep).
  # Computed for every run; low-ESS runs are removed later in the RUV section.
  keep <- round(seq(1, nrow(t), length.out = n_keep))
  ranks[[length(ranks) + 1]] <- data.frame(
    sim_id   = i,
    variable = params$name,
    rank     = sapply(seq_len(nrow(params)), function(k) sum(t[[log_cols[k]]][keep] < truth[k])))

  post_ess[i] <- as.numeric(effectiveSize(as.mcmc(t$posterior)))
  if (post_ess[i] < min_ess) {
    remove_rows <- c(remove_rows, i)
    cat("post ess for", i, "is too small (", post_ess[i], ")\n")
    next
  }
  rows  <- lapply(seq_len(nrow(params)), function(k)
    summarise_draws(truth[k], t[[log_cols[k]]]))

  results[[length(results) + 1]] <- data.frame(sim = i,
                                               parameter = params$name,
                                               do.call(rbind, rows))
}
cat("\n")

if (length(results) == 0) stop("No usable logs for run_", run_num)
results <- do.call(rbind, results)
valid   <- length(unique(results$sim))
removed <- length(remove_rows)
cat("Valid runs:", valid, " Removed:", removed, "\n")

ranks   <- do.call(rbind, ranks)

write.csv(results, file = paste0(fig_dir, "per_sim_summaries.csv"), row.names = FALSE)
ranks$post_ess <- post_ess[ranks$sim_id]
write.csv(ranks,   file = paste0(fig_dir, "ranks.csv"),             row.names = FALSE)
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
# Rank uniformity validation (RUV / SBC)
# If inference is calibrated, each parameter's ranks are uniform on 0..n_keep,
# so the ECDF stays inside the band. The origin is fixed at 4 in the simulation
# rather than drawn from its prior, so its ranks are not expected to be uniform.
######################################################
# runs whose posterior ESS is below min_ess are left out of the RUV plots
remove <- which(post_ess < min_ess)
cat("RUV: removing", length(remove), "run(s) with posterior ESS <", min_ess,
    if (length(remove)) paste0("(", paste(remove, collapse = ", "), ")") else "", "\n")
ranks_ruv <- ranks[!ranks$sim_id %in% remove, ]
writeLines(as.character(remove), paste0(fig_dir, "ruv_removed_low_ess.txt"))

ruv_theme <- ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none",
                 text = ggplot2::element_text(size = 30),
                 strip.text = ggplot2::element_blank(),
                 strip.background = ggplot2::element_blank())

for (k in seq_len(nrow(params))) {
  v   <- params$name[k]
  sub <- ranks_ruv[ranks_ruv$variable == v, c("sim_id", "variable", "rank")]

  p_ecdf <- SBC::plot_ecdf(sub, max_rank = n_keep) + ruv_theme +
    ggplot2::xlab("Normalised rank") + ggplot2::ylab("Cumulative Probability") +
    ggplot2::ggtitle(params$label[k])
  p_diff <- SBC::plot_ecdf_diff(sub, max_rank = n_keep) + ruv_theme +
    ggplot2::xlab("Normalised rank") + ggplot2::ylab("ECDF difference") +
    ggplot2::ggtitle(params$label[k])

  ggsave(paste0(fig_dir, "ecdf_", v, ".pdf"),      p_ecdf, width = 10, height = 10)
  ggsave(paste0(fig_dir, "ecdf_diff_", v, ".pdf"), p_diff, width = 10, height = 10)
}

# all parameters on one page, one panel each
p_all <- SBC::plot_ecdf_diff(ranks_ruv[, c("sim_id", "variable", "rank")], max_rank = n_keep) +
  ggplot2::theme_bw()
ggsave(paste0(fig_dir, "ecdf_diff_all.pdf"), p_all, width = 14, height = 10)

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

cat("Figures and tables written to", fig_dir, "\n")