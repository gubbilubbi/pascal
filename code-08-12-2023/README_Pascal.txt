
To do before real run

0. run addPaths.m (fix name to pascal instead of scilife and add spm to the path.)

Settings:
1. Set correct idsPD and idsHC

	- Normalisation (zscore/eucnorm)
	- FC -> Adj (abs, geq0)
	- *groupFC (One FC per group or one per subject)
	- *NNcutoff based on signals
	- capOutliers of the zscored data (didn't seem to change much at first try)
	- timepoints
	- boolFC (FC or SC)

2. run create_codebook by uncommenting in plot on brain, comment back afterwards because
   it is annoying with the warning


Non-mandatory changes for running

	- Plot the PSD for each subject/group

	- Plot the variance of the low/high U/X (Barplot?) for certain regions of 		  interest?

	- median_s vs median_g is there a difference between these?


