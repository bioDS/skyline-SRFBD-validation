This repository contains validation code for the [ssRanges BEAST 2 pacage](https://github.com/bioDS/skyline-stratigraphic-ranges), and is based on similar validations for the sRanges package, see the [sRanges material](https://github.com/jugne/sRanges-material) repository.

The code uses a [fork of TreeSim](https://github.com/bioDS/TreeSim) suitable for skyline diversification rates.

`simulate.R` contains code to simulate 200 trees using four equal length skyline intervals, and to produce XML files to infer the diversification rates for each tree, while `simulateOne.R` is the equivalent script for 200 trees with only one skyline inteval (equivalent to constant rates). The scripts use the template XML `templates/ssRanges_inference_template` to produce the XML for each simulation.

Use `slurm-inferences.sh` to run the BEAST 2 analyses for each simulated tree. The script uses the `Intellij_module_versions` folder to control the versions of BEAST 2 packages used.

Code to graph the log files from analyses is included in the `graphing` folder.
