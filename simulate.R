# Code to simulate 200 trees under a skyline FBD model to validate the sRanges package.
# Skyline interval start times and template paths are given as command-line arguments:
#
#   Rscript simulate.R                       # defaults: four intervals, times 0,1,2,3
#   Rscript simulate.R --times=0 \
#       --inference_template=templates/ssRanges_inference_template_one_int.xml
#
# Arguments (all optional, form --name=value):
#   --times               comma-separated interval start times, must start at 0 and be below the origin (4)
#   --inference_template  inference XML template
#   --dna_template        DNA simulation XML template
#   --morph_template      morphological simulation XML template
# Template paths are relative to the working directory (wd below) unless absolute.
# Based on https://github.com/jugne/sRanges-material/blob/main/validation/estimate_all_params_ext_dna/simulate.R

rm(list = ls())
options(digits = 16)

library(ape)
library(phytools)
library(ggplot2)
library(stringr)
library(beastio)
library(FossilSim)
library(R.utils)


## helper function
"%notin%" <- Negate("%in%")
wd = "/home/ket581/skyline/SUMMER/skyline-SRFBD-validation/"

setwd(wd)

## change paths if needed
templates_dir <- paste0(wd, "templates/")
beast_dir <- "~/beast/bin/beast"

######################################################
# Command-line arguments (defaults used when running interactively)
######################################################
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[1]) else default
}

times_arg          <- get_arg("times", "0,1,2,3")
inference_template <- normalizePath(get_arg("inference_template", paste0(templates_dir, "ssRanges_inference_template.xml")), mustWork = TRUE)
dna_template       <- normalizePath(get_arg("dna_template",       paste0(templates_dir, "sRanges_simDNA_template.xml")),       mustWork = TRUE)
morph_template     <- normalizePath(get_arg("morph_template",     paste0(templates_dir, "sRanges_simMorph_template.xml")),     mustWork = TRUE)

######################################################
# Output folder: runs/run_NUM
#   overwrite = FALSE -> create run_(highest + 1), or run_1 if none exist
#   overwrite = TRUE  -> delete and recreate the highest existing run_NUM
######################################################
overwrite <- FALSE
runs_root <- paste0(wd, "runs/")
dir.create(runs_root, showWarnings = FALSE, recursive = TRUE)

existing <- list.dirs(runs_root, full.names = FALSE, recursive = FALSE)
run_nums <- as.integer(sub("^run_", "", grep("^run_[0-9]+$", existing, value = TRUE)))

if (overwrite && length(run_nums) > 0) {
  run_num <- max(run_nums)
  unlink(paste0(runs_root, "run_", run_num), recursive = TRUE)
  cat("Overwriting run_", run_num, "\n", sep = "")
} else {
  run_num <- if (length(run_nums) > 0) max(run_nums) + 1 else 1
  cat("Creating run_", run_num, "\n", sep = "")
}

run_dir <- paste0(runs_root, "run_", run_num, "/")
dir.create(run_dir, recursive = TRUE)
out_dir <- run_dir
file.create(paste0(out_dir, "out"))


set.seed(5647829)


### Parameter transformation functions
f_lambda <- function(diversification, turnover) {
  return(diversification / (1 - turnover))
}

f_mu <- function(turnover, lambda) {
  return(turnover * lambda)
}

f_psi <- function(mu, sampling_prop) {
  return(mu * sampling_prop / (1 - sampling_prop))
}

## Set initial params to null, origin will be fixed (but still estimated)
div_rate = NULL
turnover = NULL
sampling_prop = NULL
sampl_extant_prob = NULL
lambda = NULL
mu = NULL
psi = NULL
redraws = 0
origin = 4

# Boundaries of skyline intervals, from the --times argument. This should always include zero, and not the origin.
# e.g. --times=0,1,2,3 gives four intervals (0 - 1, ..., 3 - 4 (the origin)); --times=0 gives one interval (0 - 4).
times = as.numeric(strsplit(times_arg, "[, ]+")[[1]])
if (any(is.na(times)) || times[1] != 0 || is.unsorted(times, strictly = TRUE) || max(times) >= origin)
  stop("--times must be increasing, start at 0 and be below the origin (", origin, "): got '", times_arg, "'")
ints= length(times)
cat("Intervals:", ints, " times:", paste(times, collapse = " "), "\n")
cat("Inference template:", inference_template, "\n")

## function to redraw parameters from their posteriors (uniform in all cases here)
redraw_params <- function() {
  div_rate <<- runif(ints, 0.7, 0.9)
  turnover <<- runif(ints, 0.2, 0.8)
  sampling_prop <<- runif(ints, 0.2, 0.8)
  sampl_extant_prob <<- runif(1, 0.7, 1.)

  lambda <<- f_lambda(div_rate, turnover)
  mu <<- f_mu(turnover, lambda)
  psi <<- f_psi(mu, sampling_prop)
}

# number of trees to simulate
ntrees <- 200 


## The following min and max values are set to ensure a small enough tree.
## If they are violated, the simulation will be rejected and new parameter values redrawn.
## This is a bit problematic, since we do not account in any way for rejected simulations and parameter combinations when doing inference.
## In inference we only condition on survival (so not on min or max number of samples).

# min for extant samples and max for total nodes
min_ext_samples = 5
max_nodes = 1000

# min and max value for fossils
min_fossils = 0
max_fossils = 10000

fix_singleton_last <- function(tree, newick) {
  labs <- tree$tip.label
  sp   <- sub("_.*", "", labs)
  for (s in unique(sp)) {
    idx <- which(sp == s)
    if (length(idx) == 1 && labs[idx] == paste0(s, "_last")) {
      old <- labs[idx]
      new <- paste0(s, "_first")
      tree$tip.label[idx] <- new
      # replace label in the newick string; labels sit between ( or , and :
      newick <- gsub(paste0("([(,])", old, ":"), paste0("\\1", new, ":"), newick)
    }
  }
  list(tree = tree, newick = newick)
}


true_rates_list <- list()

### Simulate Trees

# simulating trees and fossils with parameters
trees = list()
beast_trees = list()
fossils = list()
taxonomy = list()
samp_trees = list()
while(length(trees) < ntrees) {
  redraw_params()
  redraws = redraws + 1
  # Simulate tree using skyline parameters
  tree_tmp <- TreeSim::sim.rateshift.age(age = origin, numbsim = 1, lambda = lambda, mu = mu, times = times, mrca = FALSE, complete = TRUE)

  if (length(tree_tmp[[1]]) == 1 || tree_tmp[[1]]$Nnode > max_nodes) {
    next # reject and redraw parameters
  }
  mrca <- max(ape::node.depth.edgelength(tree_tmp[[1]]))
  if (length(which((mrca - ape::node.depth.edgelength(tree_tmp[[1]])) < 1e-7)) < min_ext_samples) {
    next # reject and redraw parameters
  }
    t <- tree_tmp[1][[1]]
    actual_origin <- tree.max(as.phylo(t))
  # Call to FossilSim needs to include the origin in the interval input  
  horizons <- c(times, actual_origin)
  taxonomy_tmp <- sim.taxonomy(tree_tmp[[1]], beta = 0, lambda.a = 0)
  # Simualte fossils using skyline rates.
  fossils_tmp = FossilSim::sim.fossils.intervals(rates=psi, taxonomy = taxonomy_tmp, interval.ages = horizons)
  beast_tree_tmp = beast.fbd.format(tree_tmp[[1]], fossils_tmp, rho = sampl_extant_prob, digits = 16)
  tree_after_rho = ape::read.tree(text = beast_tree_tmp)
  mrca = max(ape::node.depth.edgelength(tree_after_rho))
  n_ext <- length(which((mrca - ape::node.depth.edgelength(tree_after_rho)) < 1e-7))
  n_fossils <- length(tree_after_rho$tip.label) - n_ext
  if (n_ext < min_ext_samples || n_fossils < min_fossils || n_fossils > max_fossils){
    next # reject and redraw parameters
  }
  write(beast_tree_tmp, file = paste0(out_dir, "out"), append = TRUE)
  trees <- c(trees, tree_tmp)
  beast_trees[[length(trees)]] <- beast_tree_tmp
  fossils[[length(trees)]] <- fossils_tmp
  taxonomy <- c(taxonomy, taxonomy_tmp)
  vals <- c(
    setNames(div_rate,      paste0("div_rate_", 1:ints)),
    setNames(turnover,      paste0("turnover_", 1:ints)),
    setNames(sampling_prop, paste0("sampling_prop_", 1:ints)),
    rho = sampl_extant_prob,
    setNames(lambda,        paste0("birth_", 1:ints)),
    setNames(mu,            paste0("death_", 1:ints)),
    setNames(psi,           paste0("sampling_", 1:ints)),
    origin    = origin,
    mrca      = mrca,
    n_samples = n_fossils + n_ext,
    n_extant  = n_ext,
    n_ranges  = 0,
    draws     = redraws
  )
  true_rates_list[[length(trees)]] <- data.frame(as.list(vals),
                                                 times = paste(times, collapse = " "),
                                                 tree  = beast_tree_tmp,
                                                 stringsAsFactors = FALSE)
  redraws = 0
  i = length(trees)
  print(i)
}

true_rates <- do.call(rbind, true_rates_list)
true_rates <- data.frame(sim = seq_len(nrow(true_rates)), true_rates)
# Now create simulation xmls for each tree.

for (i in 1:ntrees) {
  beast_tree <- beast_trees[[i]]
  true_rates$tree[i] = beast_tree
  tmp_tree <- ape::read.tree(text = beast_tree)
  fixed      <- fix_singleton_last(tmp_tree, beast_tree)
  tmp_tree   <- fixed$tree
  beast_tree <- fixed$newick
  beast_trees[[i]]   <- beast_tree
  sample_times = ape::node.depth.edgelength(tmp_tree)
  sample_times <- max(sample_times) - sample_times
  mrca <- max(sample_times)  
  labs <- tmp_tree$tip.label
  sp   <- sub("_.*", "", labs)
  bad  <- Filter(function(s) {
  if (sum(sp == s) > 1) !(paste0(s, "_last") %in% labs) else !(paste0(s, "_first") %in% labs)
  }, unique(sp))
  if (length(bad)) { cat("Tree", i, "\n"); print(labs[sp %in% bad]) } 
  sample_times_round = round(sample_times, 15)
  sample_times_round[which(sample_times_round < 1e-10)] <- 0
  oldest_sample <- max(sample_times_round[seq_len(ape::Ntip(tmp_tree))])

  taxon =list()
  taxon_extant = list()
  strat_ranges = list()
  strat_ranges_refs = list()
  j = 0
  k = 1
  for (k in 1:length(unique(sub("_.*", "", tmp_tree$tip.label)))) {

    tip = unique(sub("_.*", "", tmp_tree$tip.label))[k]
    taxon <- append(taxon, paste0(tip, "_first"))
    if (length(which(sub("_.*", "", tmp_tree$tip.label) %in% tip)) > 1) {
      taxon <- append(taxon, paste0(tip, "_last"))
      strat_ranges <- append(strat_ranges,
        paste0('<stratigraphicRange id="r',
          j,'" spec="StratigraphicRange" firstOccurrence="@',
          paste0(tip, "_first"), '" lastOccurrence="@',
          paste0(tip, "_last"), '"/>'))
      if (sample_times_round[which(tmp_tree$tip.label == paste0(tip, "_last"))] == 0.0) {
        taxon_extant <- append(taxon_extant, paste0(tip, "_last"))
      }
    } else {
      strat_ranges <- append(
        strat_ranges,
        paste0(
          '<stratigraphicRange id="r',
          j, '" spec="StratigraphicRange" firstOccurrence="@',
          paste0(tip, "_first"), '" lastOccurrence="@',
          paste0(tip, "_first"), '"/>'
        )
      )
      if (sample_times_round[which(tmp_tree$tip.label == paste0(tip, "_first"))] == 0.0) {
        taxon_extant <- append(taxon_extant, paste0(tip, "_first"))
      }
    }
    strat_ranges_refs <- append(strat_ranges_refs, paste0('<stratigraphicRange idref="r', j, '"/>'))
    j = j + 1
  }
  true_rates$n_ranges[i] = j


  sim <- readLines(dna_template)
  sim  <- gsub(pattern = "insertNewick",
               replace = paste0("newick='", beast_tree, "'"), x = sim)
  taxon_extant_str = c()
  for (tx in taxon){
    taxon_extant_str <- c(taxon_extant_str, paste0("<sequence spec='Sequence' taxon='",tx,"' value='?'/>"))
  }
  taxon_extant_str <- paste0(taxon_extant_str, collapse = "\n\t\t\t")
  
  taxon_str = c()
  taxon_set_str = c()
  taxa_age_str = c()
  for (tx in taxon){
    taxon_str <- c(taxon_str, paste0("<sequence spec='Sequence' taxon='",tx,"' value='?'/>"))
    taxon_set_str <- c(taxon_set_str, paste0("<taxon spec='Taxon' id='",tx,"'/>"))
    taxa_age_str <- c(taxa_age_str, paste0(tx, "=", sample_times_round[which(tmp_tree$tip.label==tx)]))
  }
  taxon_str <- paste0(taxon_str, collapse = "\n\t\t\t")
  taxon_set_str <- paste0(taxon_set_str, collapse = "\n\t\t\t\t\t\t")
  taxa_age_str <- paste0(taxa_age_str, collapse = ", ")
  
  sim  <- gsub(pattern = "<insertSequence/>",
               replace = taxon_extant_str, x = sim)
  
  ##### run beast on the xml simulating DNA data
  sim_dir <- paste0(run_dir, "sim/", i)
  dir.create(sim_dir, recursive = TRUE)
  setwd(sim_dir)
  writeLines(sim, con='sRanges_simDNA.xml')
  cmd <- paste0(beast_dir," -seed 42 sRanges_simDNA.xml")
  system(cmd)
  
  ## remove DNA from extinct occurences
  extinct <- taxon[which(taxon %notin% taxon_extant)]
  dna_sim <- readLines("simulated_dna_alignment.xml")
  idx  <- grep("<data id=", dna_sim)
  dna_sim[idx] <- "<data id='dna_alignment' spec='beast.base.evolution.alignment.Alignment'>"
  for (tax in extinct){
    idx <- grep(paste0("taxon=['\"]", tax, "['\"]"), dna_sim)
    dna_sim[idx] <- paste0("    <sequence spec='beast.base.evolution.alignment.Sequence' taxon='",
                           tax,"' value='", paste0(rep("-", 1000),
                                                   collapse = ""), "'/>")
  }
  dna_sim  <- gsub(pattern = "id='Sequence",
                   replace = "id='dna", x = dna_sim)
  
  writeLines(dna_sim, "simulated_dna_alignment.xml")
  
  
  
  ####### simulate morph sequences (in the same sim folder as the DNA)
  sim <- readLines(morph_template)
  sim  <- gsub(pattern = "insertNewick",
               replace = paste0("newick='", beast_tree, "'"), x = sim)
  
  sim  <- gsub(pattern = "<insertMorphSequence/>",
               replace = taxon_str, x = sim)
  
  ##### run beast on the xml simulating morphological data
  writeLines(sim, con='sRanges_simMorph.xml')
  cmd <- paste0(beast_dir," -seed 42 sRanges_simMorph.xml")
  system(cmd)
  
  morph_sim <- readLines("simulated_morph_alignment.xml")
  idx  <- grep("<data id=", morph_sim)
  morph_sim[idx] <- "<data id='morph_alignment' spec='beast.base.evolution.alignment.Alignment'>"
  morph_sim  <- gsub(pattern = "id='Sequence",
                     replace = "id='morph", x = morph_sim)
  morph_sim  <- gsub(pattern = ",",
                     replace = "", x = morph_sim)
  setwd(wd)
  
  ####### now create inference xmls, with the previosuly simulated data
  sim <- readLines(inference_template)
  sim  <- gsub(pattern = "<insertStartMorphData/>",
               replace = paste0(morph_sim, collapse = "\n"), x = sim)
  sim  <- gsub(pattern = "<insertStartDNAData/>",
               replace = paste0(dna_sim, collapse = "\n"), x = sim)
  
  sim  <- gsub(pattern = "<insertMorphSequence/>",
               replace = taxon_str, x = sim)
  
  sim  <- gsub(pattern = "<inputTaxa/>",
               replace = taxon_set_str, x = sim)

  sim <- gsub(
    pattern = "</intervalTimes/>",
    replace = paste(times, collapse = " "), x = sim
  )
  
  sim  <- gsub(pattern = "<inputTaxaAge/>",
               replace = taxa_age_str, x = sim)
  
  sim  <- gsub(pattern = "<inputStratRanges/>",
               replace = paste0(strat_ranges, collapse = "\n\t\t\t\t"), x = sim)
  sim  <- gsub(pattern = "<inputStratRangesRef/>",
               replace = paste0(strat_ranges_refs, collapse = "\n\t\t\t"), x = sim)
  
    rnd_origin <- runif(1, max(2 * oldest_sample, 1), 500) # start above the oldest sample 
  rnd_div_rate <- runif(ints, 0.7, 0.9)
  rnd_turnover <- runif(ints, 0.2, 0.8)
  rnd_sampling_prop <- runif(ints, 0.2, 0.8)
  rnd_sampl_extant_prob <- runif(1, 0.7, 1.)

  sim <- gsub(
    pattern = "<initOrigin/>",
    replace = paste0("<parameter id='origin' lower='0.0'
                                name='stateNode'>", rnd_origin, "</parameter>"),
    x = sim
  )
  sim <- gsub(
    pattern = "<initDiversificationRate/>",
    replace = paste0("<parameter id='netDiversification' lower='0.0'
                                name='stateNode'>", paste(rnd_div_rate, collapse = " "), "</parameter>"),
    x = sim
  )
  sim <- gsub(
    pattern = "<initTurnover/>",
    replace = paste0("<parameter id='turnOver' lower='0.' upper = '1.'
                                name='stateNode'>", paste(rnd_turnover, collapse = " "), "</parameter>"),
    x = sim
  )
  sim <- gsub(
    pattern = "<initSamplingProportion/>",
    replace = paste0("<parameter id='samplingProportion' lower='0.0'
                                name='stateNode'>", paste(rnd_sampling_prop, collapse = " "), "</parameter>"),
    x = sim
  )

  sim <- gsub(
    pattern = "<initSamplingAtPresentProb/>",
    replace = paste0("<parameter id='rho' lower='0.0'
                                name='stateNode'>", rnd_sampl_extant_prob, "</parameter>"),
    x = sim
  )

  inf_dir <- paste0(run_dir, "inf/", i)
  dir.create(inf_dir, recursive = TRUE)

  sim <- gsub(
    pattern = "<logname/>",
    replace = paste0(
      "<logger logEvery=\"200\" fileName=\"", inf_dir,
      "/sRanges.$(seed).log\">"
    ),
    x = sim
  )

  sim <- gsub(
    pattern = "<logtreename/>",
    replace = paste0(
      "<logger logEvery=\"5000\" fileName=\"", inf_dir,
      "/sRanges.$(seed).trees\">"
    ),
    x = sim
  )

  writeLines(sim, con = paste0(inf_dir, "/ssRanges_inference.xml"))
}

save.image(file = paste0(run_dir, "simulation.RData"))
write.csv(true_rates, file = paste0(run_dir, "true_rates.csv"), row.names = FALSE)