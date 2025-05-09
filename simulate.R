rm(list = ls())
options(digits = 16)

library(ape)
library(phytools)
library(ggplot2)
library(stringr)
library(beastio)
library(FossilSim)
library(R.utils)

# Based on https://github.com/jugne/sRanges-material/blob/main/validation/estimate_all_params_ext_dna/simulate.R

## helper function
"%notin%" <- Negate("%in%")
setwd("/home/ket581/skyline/SUMMER/skyline-SRFBD-validation/")

## change paths if needed
templates_dir <- "templates/"
beast_dir <- "beast/bin/beast"
run_count = 0
while (dir.exists(paste("simulated_trees/", run_count, sep = ""))) {
    run_count <- run_count + 1
}
out_dir <- paste("simulated_trees/", run_count, "/", sep = "")
dir.create(paste("simulated_trees/", run_count, sep = ""))
file.create(paste(out_dir, "out", sep = ""))


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
div_rate <- NULL
turnover <- NULL
sampling_prop <- NULL
sampl_extant_prob <- NULL
lambda <- NULL
mu <- NULL
psi <- NULL
redraws <- 0
origin <- 4

# Boundaries of skyline intervals. This should always include zero, and not the origin.
# In this case there are four boundaries. The first boundary is from 0 - 1, and the last from 3 - 4 (the origin).
times = c(0,1,2,3)
ints= length(times)

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
min_ext_samples <- 5
max_nodes <- 1000

# min and max value for fossils
min_fossils <- 0
max_fossils <- 10000


# This structure is getting rather large as the number of intervals increase - it would be worth finding a dynamic way to define this if there are significantly more intervals.
true_rates <- data.frame(
  div_rate_1 = numeric(), div_rate_2 = numeric(), div_rate_3 = numeric(), div_rate_3 = numeric(), 
  turnover_1 = numeric(), turnover_2 = numeric(),   turnover_3 = numeric(), turnover_4 = numeric(),
  sampling_prop_1 = numeric(), sampling_prop_2 = numeric(), sampling_prop_3 = numeric(), sampling_prop_4 = numeric(),
  rho = numeric(),  birth_1=numeric(), birth_2=numeric(), birth_3=numeric(), birth_4=numeric(),
   death_1=numeric(), death_2=numeric(), death_3=numeric(), death_4=numeric(),  
   sampling_1=numeric(), sampling_2=numeric(), sampling_3=numeric(), sampling_4=numeric(),
   times=numeric(), origin = numeric(), mrca = numeric(), tree = character(),
  n_samples = numeric(), n_extant = numeric(), n_ranges = numeric(), draws = numeric()
)

### Simulate Trees

# simulating trees and fossils with parameters
trees <- list()
beast_trees <- list()
fossils <- list()
taxonomy <- list()
samp_trees <- list()
while (length(trees) < ntrees) {
  redraw_params()
  redraws <- redraws + 1
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
    origin <- tree.max(as.phylo(t))
  # Call to FossilSim needs to include the origin in the interval input  
  horizons <- c(times, origin)
  taxonomy_tmp <- sim.taxonomy(tree_tmp[[1]], beta = 0, lambda.a = 0)
  # Simualte fossils using skyline rates.
  fossils_tmp <- FossilSim::sim.fossils.intervals(rates=psi, taxonomy = taxonomy_tmp, interval.ages = horizons)
  beast_tree_tmp <- beast.fbd.format(tree_tmp[[1]], fossils_tmp, rho = sampl_extant_prob, digits = 16)
  write(beast_tree_tmp, file = paste(out_dir, "out", sep = ""), append = TRUE)
  tree_after_rho <- ape::read.tree(text = beast_tree_tmp)
  mrca <- max(ape::node.depth.edgelength(tree_after_rho))
  n_ext <- length(which((mrca - ape::node.depth.edgelength(tree_after_rho)) < 1e-7))
  n_fossils <- length(tree_after_rho$tip.label) - n_ext
  if (n_ext < min_ext_samples || n_fossils < min_fossils || n_fossils > max_fossils) {
    next # reject and redraw parameters
  }

  trees <- c(trees, tree_tmp)
  beast_trees[[length(trees)]] <- beast_tree_tmp
  fossils[[length(trees)]] <- fossils_tmp
  taxonomy <- c(taxonomy, taxonomy_tmp)
  true_rates_tmp <- data.frame(div_rate[1],  div_rate[2], div_rate[3],  div_rate[4], turnover[1], turnover[2], turnover[3], turnover[4], 
  sampling_prop[1], sampling_prop[2], sampling_prop[3], sampling_prop[4],
    rho = sampl_extant_prob, lambda[1], lambda[2], lambda[3], lambda[4], mu[1], mu[2], mu[3], mu[4], psi[1], psi[2], psi[3], psi[4], times=paste(times,collapse=" "),  origin, mrca = mrca, tree = beast_tree_tmp,
    n_samples = n_fossils + n_ext, n_extant = n_ext, n_ranges = 0, draws = redraws
  )
  true_rates <- rbind(true_rates, true_rates_tmp)
  redraws <- 0
  i <- length(trees)
    print(i)
}

# Now create simulation xmls for each tree. We don't currently include DNA or morphological characters.

for (i in 1:ntrees) {
  beast_tree <- beast_trees[[i]]
  mrca = true_rates$mrca[i]
  true_rates$tree[i] <- beast_tree
  tmp_tree <- ape::read.tree(text = beast_tree)
  sample_times <- ape::node.depth.edgelength(tmp_tree)
  sample_times <- max(sample_times) - sample_times
  sample_times_round <- round(sample_times, 15)
  sample_times_round[which(sample_times_round < 1e-10)] <- 0

  taxon <- list()
  taxon_extant <- list()
  strat_ranges <- list()
  strat_ranges_refs <- list()
  j <- 0
  k <- 1
  for (k in 1:length(unique(sub("_.*", "", tmp_tree$tip.label)))) {
    tip <- unique(sub("_.*", "", tmp_tree$tip.label))[k]
    taxon <- append(taxon, paste0(tip, "_first"))
    if (length(which(sub("_.*", "", tmp_tree$tip.label) %in% tip)) > 1) {
      taxon <- append(taxon, paste0(tip, "_last"))
      strat_ranges <- append(
        strat_ranges,
        paste0(
          '<stratigraphicRange id="r',
          j, '" spec="StratigraphicRange" firstOccurrence="@',
          paste0(tip, "_first"), '" lastOccurrence="@',
          paste0(tip, "_last"), '"/>'
        )
      )
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
    j <- j + 1
  }
  true_rates$n_ranges[i] <- j

  # sim <- readLines(paste0(templates_dir, "ssRanges_inference_template.xml"))
  # sim <- gsub(
  #   pattern = "insertNewick",
  #   replace = paste0("newick='", beast_tree, "'"), x = sim
  # )
  taxon_extant_str <- c()
  for (tx in taxon) {
    taxon_extant_str <- c(taxon_extant_str, paste0("<sequence spec='Sequence' taxon='", tx, "' value='?'/>"))
  }
  taxon_extant_str <- paste0(taxon_extant_str, collapse = "\n\t\t\t")

  taxon_str <- c()
  taxon_set_str <- c()
  taxa_age_str <- c()
  for (tx in taxon) {
    taxon_str <- c(taxon_str, paste0("<sequence spec='Sequence' taxon='", tx, "' value='?'/>"))
    taxon_set_str <- c(taxon_set_str, paste0("<taxon spec='Taxon' id='", tx, "'/>"))
    taxa_age_str <- c(taxa_age_str, paste0(tx, "=", sample_times_round[which(tmp_tree$tip.label == tx)]))
  }
  taxon_str <- paste0(taxon_str, collapse = "\n\t\t\t")
  taxon_set_str <- paste0(taxon_set_str, collapse = "\n\t\t\t\t\t\t")
  taxa_age_str <- paste0(taxa_age_str, collapse = ", ")

#   sim <- gsub(
#     pattern = "<insertSequence/>",
#     replace = taxon_extant_str, x = sim
#   )

  ##### run beast on the xml simulating DNA data
#   sim_dir <- paste0(wd, "run_", i, "/sim")
#   dir.create(sim_dir, recursive = T)
#   setwd(sim_dir)
#   writeLines(sim, con = "sRanges_simDNA.xml")
#   cmd <- paste0(beast_dir, " -seed 42 sRanges_simDNA.xml")
#   system(cmd)

#   ## remove DNA from extinct occurences
#   extinct <- taxon[which(taxon %notin% taxon_extant)]
#   dna_sim <- readLines("simulated_dna_alignment.xml")
#   idx <- grep("<data id=", dna_sim)
#   dna_sim[idx] <- "<data id='dna_alignment' spec='beast.base.evolution.alignment.Alignment'>"
#   for (tax in extinct) {
#     idx <- grep(tax, dna_sim)
#     dna_sim[idx] <- paste0(
#       "    <sequence spec='beast.base.evolution.alignment.Sequence' taxon='",
#       tax, "' value='", paste0(rep("-", 1000),
#         collapse = ""
#       ), "'/>"
#     )
#   }
#   dna_sim <- gsub(
#     pattern = "id='Sequence",
#     replace = "id='dna", x = dna_sim
#   )

#   writeLines(dna_sim, "simulated_dna_alignment.xml")



  ####### simulate morph sequences

#   sim <- readLines(paste0(templates_dir, "sRanges_simMorph_template.xml"))
#   sim <- gsub(
#     pattern = "insertNewick",
#     replace = paste0("newick='", beast_tree, "'"), x = sim
#   )
#   taxon_extant_str <- c()
#   for (tx in taxon) {
#     taxon_extant_str <- c(taxon_extant_str, paste0("<sequence spec='Sequence' taxon='", tx, "' value='?'/>"))
#   }
#   taxon_extant_str <- paste0(taxon_extant_str, collapse = "\n\t\t\t")

#   taxon_str <- c()
#   taxon_set_str <- c()
#   taxa_age_str <- c()
#   for (tx in taxon) {
#     taxon_str <- c(taxon_str, paste0("<sequence spec='Sequence' taxon='", tx, "' value='?'/>"))
#     taxon_set_str <- c(taxon_set_str, paste0("<taxon spec='Taxon' id='", tx, "'/>"))
#     taxa_age_str <- c(taxa_age_str, paste0(tx, "=", sample_times_round[which(tmp_tree$tip.label == tx)]))
#   }
#   taxon_str <- paste0(taxon_str, collapse = "\n\t\t\t")
#   taxon_set_str <- paste0(taxon_set_str, collapse = "\n\t\t\t\t\t\t")
#   taxa_age_str <- paste0(taxa_age_str, collapse = ", ")

#   sim <- gsub(
#     pattern = "<insertMorphSequence/>",
#     replace = taxon_str, x = sim
#   )

#   ##### run beast on the xml simulating morphological data
#   writeLines(sim, con = "sRanges_simMorph.xml")
#   cmd <- paste0(beast_dir, " -seed 42 sRanges_simMorph.xml")
#   system(cmd)

#   morph_sim <- readLines("simulated_morph_alignment.xml")
#   idx <- grep("<data id=", morph_sim)
#   morph_sim[idx] <- "<data id='morph_alignment' spec='beast.base.evolution.alignment.Alignment'>"
#   morph_sim <- gsub(
#     pattern = "id='Sequence",
#     replace = "id='morph", x = morph_sim
#   )
#   morph_sim <- gsub(
#     pattern = ",",
#     replace = "", x = morph_sim
#   )

  ####### now create inference xmls, with the previosuly simulated data
  sim <- readLines(paste0(templates_dir, "ssRanges_inference_template.xml"))
#   sim <- gsub(
#     pattern = "<insertStartMorphData/>",
#     replace = paste0(morph_sim, collapse = "\n"), x = sim
#   )
#   sim <- gsub(
#     pattern = "<insertStartDNAData/>",
#     replace = paste0(dna_sim, collapse = "\n"), x = sim
#   )

#   sim <- gsub(
#     pattern = "<insertMorphSequence/>",
#     replace = taxon_str, x = sim
#   )

  sim <- gsub(
    pattern = "<inputTaxa/>",
    replace = taxon_set_str, x = sim
  )

    sim <- gsub(
    pattern = "</intervalTimes/>",
    replace = paste(times, collapse = " "), x = sim
  )

  sim <- gsub(
    pattern = "<inputTaxaAge/>",
    replace = taxa_age_str, x = sim
  )

  sim <- gsub(
    pattern = "<inputStratRanges/>",
    replace = paste0(strat_ranges, collapse = "\n\t\t\t\t"), x = sim
  )
  sim <- gsub(
    pattern = "<inputStratRangesRef/>",
    replace = paste0(strat_ranges_refs, collapse = "\n\t\t\t"), x = sim
  )

  rnd_origin <- runif(1, mrca * 2, 500) # just so origin is not smaller than mrca
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
    cat("mrca", mrca, "\n")
    index = i - 1
    inf_dir <- paste0("inf/", index)
    inf_full <- paste0("/home/ket581/skyline/SUMMER/skyline-SRFBD-validation/inf/", (i-1))
    dir.create(inf_dir, recursive = T)

        sim <- gsub(
        pattern = "<logname/>",
        replace = paste0(
            "<logger logEvery=\"5000\" fileName=\"", inf_full,
            "/sRanges.$(seed).log\">"
        ),
        x = sim
    )

    sim <- gsub(
        pattern = "<logtreename/>",
        replace = paste0(
            "<logger logEvery=\"5000\" fileName=\"", inf_full,
            "/sRanges.$(seed).trees\">"
        ),
        x = sim
    )

  writeLines(sim, con = paste(inf_dir, "/sRanges_inference.xml", sep=""))
}
save.image(file = "simulation.RData")

write.csv(true_rates, "true_rates.csv")
