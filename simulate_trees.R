# Generates 200 trees of sufficient size with an origin age of 4 and piecewise-constant
# diversification rates.
# Following the work of Stolz et al., trees with more than 1000 total samples or less than
# 5 extant samples are excluded, along with extinct trees.

library("TreeSim")
library("FossilSim")
set.seed(321)
good <- 0
# Use 4 intervals
times <- c(0, 1, 2, 3)
l <- length(times)
ntrees <- 200
run_count <- 0
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
    # Distribitions for parameters (diversification rate, turnover, rho, sampling proportion)
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
        cat("origin as tree max is ", origin)
        horizons <- c(times, origin)
        tax <- sim.taxonomy(tree = t)
        print(t)
        print(tax)
        # pdf(paste(out_dir, good, ".pdf", sep = ""))
        # plot(tax, tree = t)
        # dev.off()

        f <- sim.fossils.intervals(taxonomy = tax, interval.ages = horizons, rate = sampling)
        # pdf(paste(out_dir, good, ".f.pdf", sep = ""))
        # plot(f, tree = t, taxonomy = tax, show.taxonomy = TRUE)
        # dev.off()

        # print(write.tree(t))

        f_ss <- FossilSim::subsample.fossils.uniform.intervals(f, s, times)
        f_ss <- sim.extant.samples(f_ss, t, rho = rho)
        print(f_ss)
        pdf(paste(out_dir, good, ".fss.pdf", sep = ""))
        plot(f_ss, tree = t, taxonomy = tax, show.taxonomy = TRUE)
        dev.off()

        # fn <- paste(out_dir, good, ".fss", sep = "")
        # write.table(f_ss, fn, row.names = F, append = file.exists(fn), quote = F, sep = "\t")

        oy <- subsample.fossils.oldest.and.youngest(f_ss, t)
        extant_samples <- oy$hmin[oy$hmin < 0.00005]
        if (length(extant_samples) < 5) {
            unlink(paste(out_dir, good, sep = ""))
            cat("\nDELETING DIRECTORY\n")
            sample_reject <- sample_reject + 1
            next
        }

        pdf(paste(out_dir, good, ".old_young.pdf", sep = ""))
        plot(oy, tree = t, taxonomy = tax, show.taxonomy = TRUE)
        dev.off()
        # pdf(paste(out_dir, good, ".SR.pdf", sep = ""))
        # dev.off()

        # fn <- paste(out_dir, good, ".oy", sep = "")
        # write.table(oy, fn, append = file.exists(fn), row.names = F, quote = F, sep = "\t")

        write.table(cbind(d, v, rho, s, birth, death, sampling, times, origin), pf,
            append = file.exists(params_f),
            row.names = F, quote = F, col.names = (good == 0)
        )

        tree <- FossilSim::SAtree.from.fossils(tree = t, fossils = oy)
        fossils <- tree$fossils
        tree <- tree$tree
        print(fossils)
        nfossils <- nrow(fossils)
        if (nfossils > 1000) {
            unlink(paste(out_dir, good, sep = ""))
            cat("\nDELETING DIRECTORY\n")
            sample_reject <- sample_reject + 1
            next
        }
        origin <- FossilSim::tree.max(tree) + tree$root.edge
        tree <- SAtree(tree, TRUE)

        cat("\nFINAL TREE\n")
        print(write.tree(tree))
        print(write.tree(tree, file = paste(out_dir, "out", sep = ""), append = TRUE))


        good <- good + 1
    }
}
write.table(cbind(extinction_reject, sample_reject), pf, append = file.exists(params_f), col.names = T, row.names = F, quote = F)
