library(ape)
library(phytools)
library(ggplot2)
library(stringr)
library(beastio)
library(FossilSim)
library(R.utils)
"%notin%" <- Negate("%in%")


templates_dir <- "sRanges-material/validation/estimate_all_params_full_dna/templates/"
trees_dir <- "simulated_trees/0"
con <- file(paste(trees_dir, "/out_short", sep = ""), "r")
line <- readLines(con, n = 1)
tree_count <- 0
while (length(line) != 0) {
    tree <- ape::read.tree(text = line)
    tmp_tree <- ape::read.tree(text = line)
    sample_times <- ape::node.depth.edgelength(tmp_tree)
    mrca <- max(sample_times)
    sample_times <- max(sample_times) - sample_times
    sample_times_round <- round(sample_times, 15)
    sample_times_round[which(sample_times_round < 1e-10)] <- 0
    print(sample_times_round)

    taxon <- list()
    taxon_extant <- list()
    strat_ranges <- list()
    strat_ranges_refs <- list()

    j <- 0
    k <- 1
    for (k in 1:length(unique(sub("_.*", "", tmp_tree$tip.label)))) {
        tip <- unique(sub("_.*", "", tmp_tree$tip.label))[k]
        print(tip)
        taxon <- append(taxon, paste0(tip, "_first"))
        if (length(which(sub("_.*", "", tmp_tree$tip.label) %in% tip)) > 1) {
            taxon <- append(taxon, paste0(tip, "_last"))
            cat("\ntaxon\n")
            print(taxon)
            strat_ranges <- append(
                strat_ranges,
                paste0(
                    '<stratigraphicRange id="r',
                    j, '" spec="StratigraphicRange" firstOccurrence="@',
                    paste0(tip, "_first"), '" lastOccurrence="@',
                    paste0(tip, "_last"), '"/>'
                )
            )
            cat("sr\n")
            print(strat_ranges)
            if (paste0(tip, "_last") %in% tmp_tree$tip.label) {
                cat("\nMATCH\n")
                print(tmp_tree$tip.label)
                print(paste0(tip, "_last"))
                if (sample_times_round[which(tmp_tree$tip.label == paste0(tip, "_last"))] == 0.0) {
                    taxon_extant <- append(taxon_extant, paste0(tip, "_last"))
                }
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
            cat("sr\n")
            print(strat_ranges)
            if (paste0(tip, "_first") %in% tmp_tree$tip.label) {
                if (sample_times_round[which(tmp_tree$tip.label == paste0(tip, "_first"))] == 0.0) {
                    taxon_extant <- append(taxon_extant, paste0(tip, "_first"))
                }
            }
        }
        strat_ranges_refs <- append(strat_ranges_refs, paste0('<stratigraphicRange idref="r', j, '"/>'))
        j <- j + 1
    }
    n_ranges <- j
    sim <- readLines(paste0(templates_dir, "sRanges_simDNA_template.xml"))
    sim <- gsub(
        pattern = "insertNewick",
        replace = paste0("newick='", line, "'"), x = sim
    )
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

    sim <- gsub(
        pattern = "<insertSequence/>",
        replace = taxon_extant_str, x = sim
    )
    print(sim)
    writeLines(sim, con = "sRanges_simDNA.xml")
    cmd <- "'beast/bin/beast' -seed 42 sRanges_simDNA.xml"
    system(cmd)

    extinct <- taxon[which(taxon %notin% taxon_extant)]
    dna_sim <- readLines("simulated_dna_alignment.xml")
    idx <- grep("<data id=", dna_sim)
    dna_sim[idx] <- "<data id='dna_alignment' spec='beast.base.evolution.alignment.Alignment'>"
    dna_sim <- gsub(
        pattern = "id='Sequence",
        replace = "id='dna", x = dna_sim
    )

    writeLines(dna_sim, "simulated_dna_alignment.xml")


    ####### now create inference xmls, with the previosuly simulated data
    sim <- readLines(paste0(templates_dir, "sRanges_inference_template.xml"))
    sim <- gsub(
        pattern = "<insertStartDNAData/>",
        replace = paste0(dna_sim, collapse = "\n"), x = sim
    )

    sim <- gsub(
        pattern = "<insertMorphSequence/>",
        replace = taxon_str, x = sim
    )

    sim <- gsub(
        pattern = "<inputTaxa/>",
        replace = taxon_set_str, x = sim
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

    rnd_origin <- runif(1, mrca * 2, 500)
    rnd_div_rate <- runif(1, 0.7, 0.9)
    rnd_turnover <- runif(1, 0.2, 0.8)
    rnd_sampling_prop <- runif(1, 0.2, 0.8)
    rnd_sampl_extant_prob <- runif(1, 0.7, 1.)

    sim <- gsub(
        pattern = "<initOrigin/>",
        replace = paste0("<parameter id='origin' lower='0.0'
                                name='stateNode'>", rnd_origin, "</parameter>"),
        x = sim
    )
    sim <- gsub(
        pattern = "<initDiversificationRate/>",
        replace = paste0("<parameter id='diversificationRate' lower='0.0'
                                name='stateNode'>", rnd_div_rate, "</parameter>"),
        x = sim
    )
    sim <- gsub(
        pattern = "<initTurnover/>",
        replace = paste0("<parameter id='turnover' lower='0.' upper = '1.'
                                name='stateNode'>", rnd_turnover, "</parameter>"),
        x = sim
    )
    sim <- gsub(
        pattern = "<initSamplingProportion/>",
        replace = paste0("<parameter id='samplingProportion' lower='0.0'
                                name='stateNode'>", rnd_sampling_prop, "</parameter>"),
        x = sim
    )

    sim <- gsub(
        pattern = "<initSamplingAtPresentProb/>",
        replace = paste0("<parameter id='samplingAtPresentProb' lower='0.0'
                                name='stateNode'>", rnd_sampl_extant_prob, "</parameter>"),
        x = sim
    )
    cat("mrca", mrca, "\n")
    inf_dir <- paste0("inf/", tree_count)
    dir.create(inf_dir, recursive = T)
    writeLines(sim, con = paste(inf_dir, "/sRanges_inference.xml", sep = ""))
    line <- readLines(con, n = 1)
    tree_count <- tree_count + 1
}
close(con)
save.image("simulation.RData")
