This repository contains validation code for the [ssRanges BEAST 2 pacage](https://github.com/bioDS/skyline-stratigraphic-ranges), and is based on similar validations for the sRanges package, see the [sRanges material](https://github.com/jugne/sRanges-material) repository.

The code uses a [fork of TreeSim](https://github.com/bioDS/TreeSim) suitable for skyline diversification rates.

Use `bash pipeline.sh two` to run the complete validation process for two analyes, one for 200 simulated trees on four equal length skyline intervals, and one for 200 simulated trees with one skyline interval (equiavalent to using constant rates). Calling `bash pipeline.sh` without arguments runs these analyses, as well as a third analysis which uses the *sRanges* package on one interval in order to compare results. 

The four interval script use the template XML `templates/ssRanges_inference_template.xml` to produce the XML for each simulation, while the one interval script using the *ssRanges* package relies on `ssRanges_inference_template_one_int.xml`. All analyses utilise the DNA and morphological data xml files in the [templates](templates/) folder.

The pipeline script also calls `validation_plots.R` to produce qq plots and rank uniformity plots from the finished analyses.
