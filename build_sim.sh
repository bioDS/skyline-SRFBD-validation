R CMD build ../../TreeSim
R CMD INSTALL TreeSim_2.4.tar.gz
R CMD build ../../fossilsim
R CMD INSTALL FossilSim_2.4.1.tar.gz
Rscript simulate_trees.R
