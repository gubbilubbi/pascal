# GSP_StructuralDecouplingIndex
Code to compute and test structural decoupling index (Preti and Van De Ville 2019)


The code includes two parts: the first part runs in Matlab (folder Matlab), while the second part (Neurosynth analysis) is implemented in Python (folder Python, code adapted from Margulies PNAS 2016).


*******SYSTEM REQUIREMENTS:

The following softwares need to be installed:

- For the first part of the code: MATLAB (originally created on version R2016B), spm12 (please download spm12 from here: https://www.fil.ion.ucl.ac.uk/spm/software/download/spmreg.php and set it under the known Matlab paths with Matlab --> Set path --> Add folder with subfolders).

- For the second part of the code: PYTHON version 2.7.15, Anaconda distribution (including Jupyter Notebook App), Neurosynth (you can follow instructions reported here to easily install Neurosynth: https://github.com/neurosynth/neurosynth).

Code tested on macOS High Sierra 10.13.6.


*******INSTALLATION GUIDE:

- It will be enough to download the code folder and set the Matlab path to the folder (Matlab --> Set path --> Add folder with subfolders). To avoid conflicts with existing toolboxes / functions, it is suggested to include ONLY the paths of this folder and of the spm folder in the Matlab paths list, before running the code (by setting the path to 'default' before adding the new paths).

- In the folder Python/database_feb_2015, please unzip the two files databasetxt and features.txt. so that they are available for the meta-analysis.

Typical install time: within minutes.


*******DEMO:

///PART 1 (Matlab):

Instructions: After having set the Matlab path (see Installation), run the main script by typing GSP_FullPipeline in the Matlab command window.

Expected output: all computations described in the paper are performed and all paper figures (except for Figure 3, produced in part 2) are created for the demo dataset (10 subjects). The structural decoupling index percentile masks (which will serve as input for the following Python meta-analysis) are created and saved in the folder results/my_masks. 

Expected run time for demo on a "normal" desktop computer: around 10 minutes in total.


///PART 2 (Python):

Instructions: open and run the script 05_metaanalysis_neurosynth_myanalysis.ipynb (included in the folder Python) in the Jupyter Notebook. 

Expected output: The script runs the Neurosynth meta-analysis by taking as inputs the maps in Python/my_masks and produces the plot of figure 3 (left). 

Expected run time for demo on a "normal" desktop computer: around 4 minutes.


*******INSTRUCTIONS FOR USE: How to run the software on your data

///PART 1 (Matlab):

replace the 2 input mat files inside the folder Matlab/data, and their names in the code (script GSP_Laplacian). 

///PART 2 (Python):

replace the Neurosynth input masks inside the folder Python/my_masks and their names in the code (script 05_metaanalysis_neurosynth_myanalysis.ipynb).

Reproduction instructions: all computations described in the paper are performed in the code.
