library(ape)
library(data.table)
library(phytools)
library(ggplot2)
library(stringr)
library(beastio)
library(zeallot)
library(tibble)
library(rlist)
library(stringr)
# library(FossilSim)
library(R.utils)
"%notin%" <- Negate("%in%")


templates_dir <- "templates/"
trees_dir <- "simulated_trees/0"
pcon <- file(paste(trees_dir, "/params", sep = ""), "r")
l = 4
ntrees = 0
t <- fread(paste(trees_dir, "/params", sep = ""), sep=" ", fill = TRUE)
n<-dim(t)[1]
t<-t[1:(n-2),]
print(t)




con <- file(paste(trees_dir, "/out", sep = ""), "r")
line <- readLines(con, n = 1)
tree_count <- 0
while ((length(line) != 0) && (tree_count < 1)) {
    edited_line <- line
    tree <- ape::read.tree(text = line)
    tmp_tree <- ape::read.tree(text = line)
    renamed_tree = ape::read.tree(text = line)
    sample_times <- ape::node.depth.edgelength(tmp_tree)
    true_params = t[(tree_count*4+1):(tree_count*4+4),]
    cat("tree count ", tree_count, " ",tree_count*4+1, " ", tree_count*4+4,  "\n")
    print(true_params)
    mrca <- max(sample_times)
    sample_times <- max(sample_times) - sample_times
    sample_times_round <- round(sample_times, 15)
    sample_times_round[which(sample_times_round < 1e-5)] <- 0
    #print(sample_times_round)

    taxon <- list()
    taxon_extant <- list()
    taxon_death <- list()
    ages_death <- list()
    strat_ranges <- list()
    strat_ranges_refs <- list()
    ages <- list()
    extinct = list()

    j <- 0
    k <- 1
    done = list()
    for (k in 1:length(tmp_tree$tip.label)){
        tip = tmp_tree$tip.label[k]
        start = sub("_.*", "", tip)
        a = sample_times_round[k]
        # cat("a = ", a, "\n")
        # if (a == 0.0){
        #     extinct = append(extinct, tip)
        # }

        cat("k = ", k, "\n")
        if(str_detect(tip, ".*_\\d+") && start %notin% done){
            done = append(done, start)
            # cat("\n")
            # cat("Start is ", start, "\n")
            first = NA
            first_time = 0
            last = NA
            last_time = mrca*10
            # cat("last time ", last_time, "\n")
            # print(sample_times_round[tmp_tree$tip.label == tip])
            selected_index = str_detect(tree$tip.label, paste0(start,"_\\d+"))
            selected_tips = tree$tip.label[selected_index]
            all_ages = (sample_times_round[selected_index])
            selected_ages = all_ages[1:(length(all_ages)/2)]
            # cat("all ages", sample_times_round[selected_index],"\n")
            # cat("select ages", selected_ages,"\n")
            # cat("all tips", selected_tips,"\n")
            smallest_time =sample_times_round[selected_index]
            smallest_time = min(smallest_time[1:(length(smallest_time)/2)])
            # cat("smallest time ", smallest_time, "\n")
            smallest = selected_tips[selected_ages == smallest_time]
            if (smallest_time != 0.0){
                taxon_death = append(taxon_death, smallest)
                ages_death = append(ages_death, smallest_time)
            }
            cat("smallest ", smallest, "\n")
            # print(tree$tip.label[str_detect(tree$tip.label,paste0(start,"_\\d+"))])
            for (label in tree$tip.label[str_detect(tree$tip.label,paste0(start,"_\\d+"))]){
                age = sample_times_round[tree$tip.label == label][1]
                if (age == 0.0) {
                    taxon_extant <- append(taxon_extant, label)
                }
                else {
                    if (age > smallest_time){
                        # cat("l time is ", l_time, "\n")
                        # cat("other label ", label, "\n")
                        if(age < last_time){
                            last = label
                            last_time = age
                        }
                        if(age > first_time){
                            first = label
                            first_time = age
                        }

                    }
                    }
            }
            cat("first ", first, " last ", last,"\n")
            if((!(is.na(first))) && !(is.na(last))){
            tmp_tree$tip.label[tmp_tree$tip.label == first] = paste0(start, "_first")
            taxon <- append(taxon, paste0(start, "_first"))
            edited_line <- gsub(pattern = first, replace = paste0(start, "_first"), x = edited_line)
            ages <- append(ages, first_time)
            if (first_time != last_time){
                tmp_tree$tip.label[tmp_tree$tip.label == last] = paste0(start, "_last")
                strat_ranges <- append(
                strat_ranges,
                paste0(
                    '<stratigraphicRange id="r',
                    j, '" spec="StratigraphicRange" firstOccurrence="@',
                    paste0(start, "_first"), '" lastOccurrence="@',
                    paste0(start, "_last"), '"/>'
                )
            )
            taxon <- append(taxon, paste0(start, "_last"))
            ages <- append(ages, last_time)
            edited_line <- gsub(pattern = last, replace = paste0(start, "_last"), x = edited_line)
            }
            else {
                strat_ranges <- append(
                strat_ranges,
                paste0(
                    '<stratigraphicRange id="r',
                    j, '" spec="StratigraphicRange" firstOccurrence="@',
                    paste0(start, "_first"), '" lastOccurrence="@',
                    paste0(start, "_first"), '"/>'
                )
            )
            }
            strat_ranges_refs <- append(strat_ranges_refs, paste0('<stratigraphicRange idref="r', j, '"/>'))
            j <- j + 1
            }

        }
    }

    # print(tmp_tree$tip.label)
  
    # cat("AGES\n")
    # print(ages)
    #     cat("TAXON HERE\n")
    # print(taxon)
    cat("LINES\n")
    print(line)
    cat("lines2\n")
    print(edited_line)
    n_ranges <- j
    sim <- readLines(paste0(templates_dir, "ssRanges_simDNA_template.xml"))
    sim <- gsub(
        pattern = "insertNewick",
        replace = paste0("newick='", edited_line, "'"), x = sim
    )
    taxon_extant_str <- c()
    taxon_str <- c()
    taxon_set_str <- c()
    taxa_age_str <- c()
    txc = 1
    for (tx in taxon_extant) {
        taxon_extant_str <- c(taxon_extant_str, paste0("<sequence spec='Sequence' taxon='", tx, "' value='?'/>"))
        taxon_set_str <- c(taxon_set_str, paste0("<taxon spec='Taxon' id='", tx, "'/>"))
        taxa_age_str <- c(taxa_age_str, paste0(tx, "=", 0))
        txc = txc + 1
    }
    taxon_extant_str <- paste0(taxon_extant_str, collapse = "\n\t\t\t")


    txc = 1
    txcnow = 0
    for (tx in taxon) {
        taxon_str <- c(taxon_str, paste0("<sequence spec='Sequence' taxon='", tx, "' value='?'/>"))
        taxon_set_str <- c(taxon_set_str, paste0("<taxon spec='Taxon' id='", tx, "'/>"))
        taxa_age_str <- c(taxa_age_str, paste0(tx, "=", ages[txc]))
        txc = txc + 1
    }



    sim <- gsub(
        pattern = "<insertSequence/>",
        replace = taxon_extant_str, x = sim
    )
    # print(sim)
    writeLines(sim, con = "ssRanges_simDNA.xml")
    cmd <- "'/Applications/BEAST\ 2.7.7/bin/beast' -seed 42 ssRanges_simDNA.xml"
    system(cmd)
    # cat("taxon death\n")
    # print(taxon_death)

    # extinct <- taxon[which(taxon %notin% taxon_extant)]
    cat("extinct is\n")
    print(extinct)
 
    dna_sim <- readLines("simulated_dna_alignment.xml")
    txc = 1
      for (tax in taxon_death){
        cat("tax = ", tax, " age \n")
        print(ages_death[txc])
    idx  <- grep(paste0("'",tax), dna_sim)
    dna_sim[idx] <- paste0("    <sequence spec='beast.base.evolution.alignment.Sequence' taxon='",
                           tax,"' value='", paste0(rep("-", 1000),
                                                   collapse = ""), "'/>")
        taxon_set_str <- c(taxon_set_str, paste0("<taxon spec='Taxon' id='", tax, "'/>"))
        taxa_age_str <- c(taxa_age_str, paste0(tax, "=", ages_death[txc]))
        txc = txc + 1
  }
    idx <- grep("<data id=", dna_sim)
    dna_sim[idx] <- "<data id='dna_alignment' spec='beast.base.evolution.alignment.Alignment'>"
#     for (tax in extinct){
#     idx  <- grep(paste0("'",tax), dna_sim)
#     dna_sim[idx] <- paste0("    <sequence spec='beast.base.evolution.alignment.Sequence' taxon='",
#                            tax,"' value='", paste0(rep("-", 1000),
#                                                    collapse = ""), "'/>")
#   }
  dna_sim <- gsub(pattern = "id='Sequence",
                   replace = "id='dna", x = dna_sim)
  

    writeLines(dna_sim, "simulated_dna_alignment.xml")


    ####### now create inference xmls, with the previously simulated data
    sim <- readLines(paste0(templates_dir, "ssRanges_inference_template.xml"))
    sim <- gsub(
        pattern = "<insertStartDNAData/>",
        replace = paste0(dna_sim, collapse = "\n"), x = sim
    )

    # sim <- gsub(
    #     pattern = "<insertMorphSequence/>",
    #     replace = taxon_str, x = sim
    # )

    taxon_str <- paste0(taxon_str, collapse = "\n\t\t\t")
    taxon_set_str <- paste0(taxon_set_str, collapse = "\n\t\t\t\t\t\t")
    taxa_age_str <- paste0(taxa_age_str, collapse = ", ")

    sim <- gsub(
        pattern = "<inputTaxa/>",
        replace = taxon_set_str, x = sim
    )
    cat("SET STR\n")
        print(taxon_set_str)
    # print(taxa_age_str)
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
    rnd_div_rate <- runif(l, 0.7, 0.9)
    rnd_turnover <- runif(l, 0.2, 0.8)
    rnd_sampling_prop <- runif(l, 0.2, 0.8)
    rnd_sampl_extant_prob <- runif(1, 0.7, 1)

    # rnd_birth = rnd_div_rate/(1-rnd_turnover)
    # rnd_death = rnd_div_rate*rnd_turnover/(1-rnd_turnover)
    # rnd_sampling = rnd_sampl_extant_prob*rnd_death/(1-rnd_sampling_prop)

    
    cat("times\n")
    print(true_params$times)

    #         sim <- gsub(
    #     pattern = "<initTimes/>",
    #     replace = paste0("<parameter id='intervalTimes'
    #                             name='stateNode'>", true_params$times, "</parameter>"),
    #     x = sim
    # )

  sim  <- gsub(pattern = "<initOrigin/>",
               replace = paste0("<parameter id='origin' lower='0.0'
                                name='stateNode'>",rnd_origin,"</parameter>"),
               x = sim)
  sim  <- gsub(pattern = "<initDiversificationRate/>",
               replace = paste0("<parameter id='netDiversification' lower='0.0'
                                name='stateNode'>",rnd_div_rate,"</parameter>"),
               x = sim)
  sim  <- gsub(pattern = "<initTurnover/>",
               replace = paste0("<parameter id='turnOver' lower='0.' upper = '1.'
                                name='stateNode'>",rnd_turnover,"</parameter>"),
               x = sim)
  sim  <- gsub(pattern = "<initSamplingProportion/>",
               replace = paste0("<parameter id='samplingProportion' lower='0.0'
                                name='stateNode'>",rnd_sampling_prop,"</parameter>"),
               x = sim)
  
  sim  <- gsub(pattern = "<initSamplingAtPresentProb/>",
               replace = paste0("<parameter id='rho' lower='0.0'
                                name='stateNode'>",rnd_sampl_extant_prob,"</parameter>"),
               x = sim)

    #   sim  <- gsub(pattern = "<initOrigin/>",
    #            replace = paste0("<parameter id='origin' lower='0.0'
    #                             name='stateNode'>",rnd_origin,"</parameter>"),
    #            x = sim)



    # sim <- gsub(
    #     pattern = "<initBirth/>",
    #     replace = paste0("<parameter id='birthRate' lower='0.0'
    #                             name='stateNode'>", rnd_birth, "</parameter>"),
    #     x = sim
    # )
    # sim <- gsub(
    #     pattern = "<initDeath/>",
    #     replace = paste0("<parameter id='deathRate' lower='0.0'
    #                             name='stateNode'>", rnd_death, "</parameter>"),
    #     x = sim
    # )
    # sim <- gsub(
    #     pattern = "<initSampling/>",
    #     replace = paste0("<parameter id='samplingRate' lower='0.' upper = '1.'
    #                             name='stateNode'>", rnd_sampling, "</parameter>"),
    #     x = sim
    # )
    # sim <- gsub(
    #     pattern = "<initRho/>",
    #     replace = paste0("<parameter id='m_rho' lower='0.0'
    #                             name='stateNode'>", rnd_sampl_extant_prob, "</parameter>"),
    #     x = sim
    # )
    cat("mrca", mrca, "\n")
    inf_dir <- paste0("inf/", tree_count)
    dir.create(inf_dir, recursive = T)
    writeLines(sim, con = paste(inf_dir, "/ssRanges_inference.xml", sep = ""))
    line <- readLines(con, n = 1)
    tree_count <- tree_count + 1
}
close(con)
save.image("simulation.RData")
