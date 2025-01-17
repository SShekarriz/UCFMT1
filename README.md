Investigating Engraftment in FMT
================================

This is the github repository accompanying the paper [Detecting microbial
engraftment after FMT: sorting signals from
noise](https://www.biorxiv.org/content/10.1101/2024.09.11.612315v1), currently
in review with Nature Communications.

This repository is limited to code for identifying and counting engraftment
events, and analysis downstreatm of that. Read assembly, bin and MAG calling,
taxonomic assignment, mapping reads to assemblies or marker databases is
upstream of this repository.

## Workflow

The scripts in the [scripts](./scripts/) directory are run in the following
order:

0. `permutation_functions.R` & `paperfig_functions.R`
	* contain the functions used to process the data, run the analysis, and
	generate the figures. Sourced by the other scripts.
1. `data_processing.Rmd`
	* Determines cutoffs for presence and absence by profiling donor samples.
	* Identifies and counts engraftment events for all feature types and
	profiles patient samples for all event types.
	* writes files to `processed_data/`, which is not tracked by git, but should
	be completely generateable from the scripts and data files in this
	repository.
2. `test_permutation.R`
	* uses data from [test_data](./test_data/) to test permutation functions and
	make sure they are acting as expected.


## Directories

* `scripts/`
	* `data_processing.Rmd`
	* `obj_and_desc.Rmd` & `obj_and_desc.html` 
		* an attempt to track the data files not in git and make sure they are
		well-behaved.
* `data/`
* `logs/`
* `plots/`
* `figs/`
* `permut_data/`
* `plots_paper/`
* `results/`
* `test_data/`
