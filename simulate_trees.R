# Generates 200 trees of sufficient size with an origin age of 4 and piecewise-constant
# diversification rates.
# Following the work of Stolz et al., trees with more than 1000 total samples or less than
# 5 extant samples are excluded, along with extinct trees.

library("TreeSim")
library("FossilSim")
library("stringr")
set.seed(321)
"%notin%" <- Negate("%in%")

good <- 0
# Use 4 intervals
times <- c(0, 2) # c(0, 1, 2, 3)
l <- length(times)
ntrees <- 200
run_count <- 0
if (!(dir.exists("simulated_trees"))) {
    dir.create("simulated_trees")
}

while (dir.exists(paste("simulated_trees/", run_count, sep = ""))) {
    run_count <- run_count + 1
}
out_dir <- paste("simulated_trees/", run_count, "/", sep = "")
dir.create(paste("simulated_trees/", run_count, sep = ""))

file.create(paste(out_dir, "params", sep = ""))
params_f <- paste(out_dir, "params", sep = "")
pf <- file(params_f, "wb")
sample_reject <- 0
extinction_reject <- 0
while (good < ntrees) {
    origin <- 4
    cat("Valid trees generated = ", good, "\n")
    # Distributions for parameters (diversification rate, turnover, rho, sampling proportion)
    d <- runif(l, min = 0.7, max = 0.9)
    v <- runif(l, min = 0.2, max = 0.8)
    rho <- runif(1, min = 0.7, max = 1)
    s <- runif(l, min = 0.2, max = 0.8)

    birth <- d / (1 - v)
    death <- d * v / (1 - v)
    sampling <- s * death / (1 - s)
    cat("birth ", birth, "  death ", death, " sampling ", sampling, "\n")
    tree <- TreeSim::sim.rateshift.age(age = origin, numbsim = 1, lambda = birth, mu = death, times = times, mrca = FALSE, complete = TRUE)
    if ((is.numeric(tree[[1]]))) {
        tnum <- as.numeric(tree[[1]])
        if (tnum == 0) {
            extinction_reject <- extinction_reject + 1
            cat("Tree went extinct\n")
        } else {
            sample_reject <- sample_reject + 1
            cat("Tree had 1 sample only\n")
        }
    } else {
        t <- tree[1][[1]]
        origin <- tree.max(as.phylo(t))
        horizons <- c(times, origin)
        tax <- sim.taxonomy(tree = t)
        f <- sim.fossils.intervals(taxonomy = tax, interval.ages = horizons, rates = sampling)

        f_ss <- FossilSim::subsample.fossils.uniform.intervals(f, s, times)
        f_ss <- sim.extant.samples(f_ss, t, rho = rho)
        pdf(paste(out_dir, good, ".fss.pdf", sep = ""))
        plot(f_ss, tree = t, taxonomy = tax, show.taxonomy = TRUE)
        dev.off()

        oy <- subsample.fossils.oldest.and.youngest(f_ss, t)
        tree <- FossilSim::SAtree.from.fossils(tree = t, fossils = oy)


        fossils <- tree$fossils
        extant_samples <- fossils$tip.label[fossils$hmin < 0.00005]
        tree <- tree$tree
        tree <- sampled.tree.from.combined(tree, sampled_tips = extant_samples)

        nfossils <- length(tree$tip.label)

        if ((nfossils > 1000) || (length(extant_samples) < 5)) {
            sample_reject <- sample_reject + 1
            next
        }

        write.table(cbind(d, v, rho, s, birth, death, sampling, times, origin), pf,
            append = file.exists(params_f),
            row.names = F, quote = F, col.names = (good == 0)
        )

        print(write.tree(tree, file = paste(out_dir, "out", sep = ""), append = TRUE))
        good <- good + 1
    }
}
write.table(cbind(extinction_reject, sample_reject), pf, append = file.exists(params_f), col.names = T, row.names = F, quote = F)
