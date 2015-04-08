## ica+ PACKAGE FOR DUAL REGRESSION IN PYTHON
Vanessa Sochat

Currently being updated for Stanford, Sherlock Cluster
Updated: September 21st, 2011

### Overview

This package derives functional networks from resting bold data, filters, performs quality analysis, matches to a user specified template image, and produces AIM xml files for use with Annotated Image Markup. MRtools.py is a python module that can be used for basic image manipulation, matching, and filtering, and the matching (pyMatch), AIM Template generation (AIMTemp), Quality analysis (melodic_qa), and filtering (melodic_hp) can be run by the package, or standalone on command line.


### ica+ Package Contents:

ica+.py: The submission script, which takes user input and either batch runs single subject ICAS, does a group ICA, checks QA of single subject ICA reports, or runs a dual regression with a group ICA. 
MRtools.py A python module that allows for reading in image files and simple coordinate lookup, matching, and filtering. The scripts pyMatch.py, and AIMTemp.py use this module for their operations. More details to come soon! 
melodic_ss.sh: a bash script that takes raw anatomical and functional data for one subject through preprocessing for group ICA (see script for details). The preprocessing is the equivalent of what would be done by selecting “multi-session temporal concatenation” in the GUI. The functional data is preprocessed, and registered to the standard MNI152_2mm_brain template with the extracted anatomical data as an intermediate registration. This script is submitted in batch via the submission python ica+.py. The current method requires using matlab for bandpass filtering with bandpass.m, you will need to customize this to work with your particular setup. 
melodic_gp.sh: a bash script that combines a user specified list of ica directories to do a group ICA (multi-session temporal concatenation) analysis. This script is also run via the submission python ica+.py. 
melodic_qa.py: A python script that can be run on its own, or submit through the submission python. It takes a list of single subject ICA directories, as well as a user specified translational motion benchmark (in mm), and a rotational motion benchmark (in degrees) and checks mcflirt output in the single subject ICA folders. It produces an HTML output page that displays flagged subjects with motion charts, as well as a _qa.txt file under “list” with subjects that pass QA, that can be used as an input file to a group ICA analysis. 
melodic_dr.sh: Performs dual regression. 
melodic_hp.sh: Runs a highpass filter over all group network gica results to produce a list of “good” ones. Can be run with ica+ (not yet tested) or standalone on command line. 
pyMatch.py: Takes list of “good” images (can be output from melodic_hp) and produces list of best matches to user specified template. 
AIMTemp.py Create an AIM xml template for the creation of xml files for AIM (Annotated Image Markup)

The entire project is hosted at github-ica+


### ica+ Add Ons:

The following scripts are considered add-ons because they are not controlled by the submission script, but still created for use with the package. 
MRlog: uses fsl or MRtools (with nibabel) to print a log of all data below a certain top directory, of type nifti.  This is is useful for distinguishing images based on timepoints or dimensions.
melodic_ts takes dual regression output, extracts the top six voxels, and then goes back to the raw data to extract the complete timeseries per each max voxel per each subject. Output is in a text file that can be used however you like. A sample use is with… 
ica_corr.m: Takes melodic_ts output, calculates the correlation matrix for the max voxels, and outputs a list of these correlations for usage with whatever your heart desires (machine learning algorithms, covariates, etc). 
pyCorr: A separate package for matching group networks to individual networks, or any sort of template to input images using nibabel and numpy in python. Includes a main script for doing the analysis (pyCorr.py) and a separate script for printing HTML report results (pyCorrPrint.py). Please note that this functionality to match is now incorporated into the MRtools module, and run with pyMatch.py!


### Requirements

python - with scipy, numpy, and nibabel modules to do matching, filtering, or AIM template generation
fsl

### Details

Using the Melodic GUI in FSL is not a transparent process, given that you input some parameters, press Go, and don't really have a good idea of what is going on. Additionally, when you perform group melodic analysis, it automatically does single subject runs and moves right into the group in a serial fashion, which is incredibly slow. It does not give you the option to use intermediate data (single subject ICA directories), so if your group ICA run has an error, you are forced to start from scratch, and wait for all the intermediate data to be re-created. Additionally, if you were to use the melodic command line utility, this assumes already preprocessed data. At the end of the day, it's very piecewise, slow, and messy to do a multi-session temporal concatenation with Melodic. 

I have created a package for running melodic and then dual regression, starting with raw data, and ending with the running of the dual regression. The user then can match this output to an image of his/her choice and produce AIM files for web query. To do this, I simply looked through FSL's source code, and wrote out the commands in bash scripts, so the entire process is completely transparent, and most importantly, easy for the user to customize. The package includes 1) a master python submission script, 2) a single subject bash ICA script, 3) a group ICA bash script, 4) a motion parameters (QA) checking python script, and 5) a dual regression script (not written by me, just modified to work with my submission script!). This package is by no means anything official, however I believe in the sharing of knowledge and tools for research, so I have chosen to include the scripts here for learning and usage.


### ica+ Usage

# Produced by typing:  python ica+.py --help
 
     Melodic+ Dual Regression Package for Python
 
     This python submission script allows for either single subject ICA analyses (ica), 
     checking quality analysis after ica, group ICA analysis (gica), and finally,
     dual regression based on a group ICA.  The script sets up sub directories based
     on a user specified experiment folder, and multiple gica and dual regression runs 
     can be done from one common set of single subject ica runs, all within the same
     experiment folder.
 
     Usage: python ica+.py <-o /fullpath/here> [options]
          see below for ica, gica, qa, and dr use cases
 
     Main Options:
       -h, --help             show this help  
       --ica                  do single subject ICA
       --gica                 group ICA
       --qa                   run QA
       --dr                   dual regression
       --match                find top DR result matches to a template image
       --aim                  create AIM xmls for top matches
       --queue                name of slurm queue to run jobs (uses sbatch)
 
     ALWAYS REQUIRED:
       -o ..., --output=...   experiment folder where the following subdirectories will be created:
			   ica: (sub1.ica, sub2.ica, sub3.ica...)
			   qa: 
			   gica: (group1.gica, group2.gica, group3.gica...)
			   dr: (dr1, dr2, dr3...)
 
      TO RUN SINGLE SUBJECT ICA:
        python ica+.py -o /my/experiment --ica=input.txt
          input.txt is a three column csv file with one row per subject
          each row contains: subjectID,full/path/highres.nii.gz,full/path/functional.nii.gz
          output goes to /my/experiment/ica/sub1.ica, /my/experiment/ica/sub2.ica...
 
     TO RUN QUALITY ANALYSIS:
        python ica+.py -o /my/experiment --qa=ica_dirs.txt --name=qarun
          (must be done after single subject ica is complete)
          ica_dirs.txt is a csv file with list of full paths to ica directories
          you can use the *_ica.txt output file (from ica) under /my/experiment/list
          qarun is a name for the unique run
          output goes to /my/experiment/qa/
 
     TO RUN GROUP ICA:
        python ica+.py -o /my/experiment --gica=ica_dirs.txt --name=run_name
          ica_dirs.txt is a single column text file listing full paths to ica directories
          you can use the *_ica.txt output file (from ica) under /my/experiment/list
          output goes to /my/experiment/gica/run_name.gica
 
     TO RUN DUAL REGRESSION:
        python ica+.py -o /my/experiment --dr=group.gica --name=dr_name --con=design.con --mat=design.mat --iter=500
          (must be done after group ICA is complete)
          group.gica is the group gica directory under /my/experiment/gica/
          design.con and design.mat are created in the FSL GUI with GLM (general linear model)
          500 is the number of iterations to run
          output goes to /my/experiment/dr/dr_name
 
     TO FIND TOP MATCHES OF DUAL REGRESSION RESULTS TO A TEMPLATE IMAGE:
        python ica+.py -o /my/experiment --match=dr_name --template=template.nii.gz
          (must be done after dual regression is complete)
          match should specify the dual regression run to use in the experiment folder. All dr_stage3_corrp* will 
          be filtered and used, and output will go to /my/experiment/match as a text file of the 
          name "dr_name-template_bestcomps.txt." This file can be specified as input to create AIM xml for AIM files
 
     TO CREATE AIM TEMPLATES FROM TOP MATCHES
        python ica+.py -o /my/exp --aim=/path/template_bestcomps.txt --ics=/path/match/template_original-ic.txt --name=dr_name
          (must be done after match is complete)
          --aim specifies a single column text file produced by --match that lists full image paths to a template and dual 
          regression images.  Should be in the format: /path/to/image.nii.gz:score.  
          --name specifies the name of the dr_run, for naming the output files
          --ics [OPTIONAL] specifies a single column text file produced by --match that lists corresponding "good" networks
          Will output AIM xml files under /my/exp/aim/dr_name-template/ for creation of AIM files for web query of results.  
 
     To use the AIMTemp.py script standalone to create an AIM template
     for one image:  python AIMTemp.py -o /fullpath/here --input=myimage.nii.gz --name=outname


### Instructions for using ica+ for multi-session temporal concatenation


#### Setup

- Put the files on your server, where you can submit jobs to a cluster, and you have access to the data 
- Put all script files in the same folder, the submission python will not run if it can't find it's buddy scripts. 
- You will want to look at the submission command, as well as the method and variables. 
- Create a text input file for your raw data, with one row per subject. Each row should have: 
the subject ID (what the ica folder will be called) followed by a comma, the FULL path to the raw anatomical data, (comma), and FULL path to raw functional data (comma). For example:

     111111,/path/to/111111/anat/mprage.nii.gz,path/to/111111/func/rest.nii.gz
     222222,/path/to/222222/anat/mprage.nii.gz,path/to/222222/func/rest.nii.gz
     333333,/path/to/333333/anat/mprage.nii.gz,path/to/333333/func/rest.nii.gz
     FOR ALL INPUT FILES: USE FULL PATHS!


### Run Single Subject ICA

Note that you will need at least python 2.4 for this to work, which has the subprocess module. You will also need an installation of matlab to either run bandpass.m with matlab directly, or compile it to create an executable for your cluster. If you type “deploytool” in matlab there is a GUI for doing this, and the runtime that you use to run it on your cluster must be from the same version that it was compiled with. (For example, deplytool via a windows matlab installation will produce an .exe file, matlab on a mac will produce a .app, and matlab on linux will produce a .sh and binary executable. To do single subject ICAs, this is the general command you will use:

     python ica+.py -o /my/experiment --ica=input.txt
 
     # -o specifies the FULL PATH to what you want for your output directory.  It can either exist, or if it doesn't, it will be created..
     # --ica=input.txt is the text file prepared above.  This tells the script where the raw, unprocessed data is!

#### Output
- ICA Directories: will be be under /my/experiment/ica/sub1.ica … sub2.ica … sub3.ica, etc. I decided to put all the single subject ica directories in the same place, because I found the FSL standard of creating them with the raw data made it very hard to delete errored runs, and find them generally.
- List of ICA Directories: under /my/experiment/list you will find a *_ica.txt file, which is basically a list of the ica directories that you just created. This file is automatically created so you can use it as the INPUT FILE for checking QA, the next step. It is INCREDIBLY IMPORTANT to note that this file is ordered by the sorted subject ID. This happens because I use a python dictionary to store paths based on the ID as a lookup, and then print it sorted. If you have some reason to use a differently ordered list, then you will need to modify the script! This only gets important if you set different covariates for the dual regression, when the subject order for subjects that go into the analysis matters. Since this is a potential issue, as a good measure, it's ALWAYS pertinent that you double check the order in the _ica list, the _qa list, and then finally, the dual regression GLM design. Minimally, you can just modify the input text file to meet your needs.
- Single Subject Log Folder: /my/experiment/ica/sub1.ica/log: Contains the (cluster specific) output file (.out) and error file (.err) in case you have an error with a specific run. Where these files are sent, and that they are created period, is one of the options that comes with the bsub command. Again, this is something you may need to modify for your specific cluster environment.


### Quality Analysis
After a single subject ICA is run, since melodic is used for motion correction, we should have a good way to cycle through all motion parameters (rotational and translational for x,y,z) for all subjects, and flag a subject with motion above this benchmark. The Quality Analysis portion of this package does exactly that, and outputs an HTML page with a summary of passing and flagged subjects. Right now, it seems pretty standard that we would want to check for 2mm translational motion, or 2 degrees rotational motion, so these parameters are hard coded into the script. I should note that melodic measures rotation in radians, however we are more comfortable with degrees, so the script takes degrees as input, and converts to radians. If you want to change the standard from 2mm or 2 degrees, feel free to edit the python script at the quality analysis submission command.

     python ica+.py -o /my/experiment --qa=ica_dirs_ica.txt --name=qarun
 
     # -o is the FULL PATH to the output experiment directory (same as for ica) - since we are checking QA this directory MUST exist!
     # (must be done after single subject ica is complete)
     # ica_dirs_ica.txt is the single column file with list of full paths to ica directories
     # again, you can use the *_ica.txt output file (from ica) under /my/experiment/list
     # qarun is a name for the unique run, so you could do multiple quality analysis checks for additional subjects
     # output goes to /my/experiment/qa/

#### Output
- HTML Report: will be be under /my/experiment/qa/qa_run.html. It will display a list of passing and flagged subjects, and for each flagged subject, a list of what they were flagged for, in what direction. It will also show, for each flagged subject, rotational, translational, and displacement charts, to help you decide to use the data.
- Passing Subject Input: The script assumes that we would want to use all single subject ICA directories that pass QA for a higher level (group) ICA analysis. To help facilitate this, it creates another input text file, under /my/experiment/list/qaname_qa.txt. This file will only contain the list of ICA directories for subjects that pass QA. If you want to re-add someone, you will need to add their directory to this file. For your own organizational sanity, it makes sense to keep this list ordered with whatever format you are using for your data - since when you create a design matrix with covariates, the order will need to coincide with this list. More details under Dual Regression.
- Flagged Subject List: Simply a text file of flagged subjects, /my/experiment/qa/qaname.flag used for generation of the HTML report. I created this with the mindset that I might want to implement other functionality for the flagged subjects at some future point.
- QA Submission Log Folder: /my/experiment/qa/log/qarun.out and qarun.err are the output and error logs from the cluster, in case you need to troubleshoot a job submission. Again, these are set with the bsub command (the job submission command) and may need to be modified for your cluster environment.


### Group ICA
After single subject ICAs have been run, and motion parameters have been checked, we are ready to run the multi-session temporal concatenation, or group ICA. You need a single column file of the single subject ICA directories to use, which means you can use /my/experiment/list/qarun_qa.txt for your input.

     python ica+.py -o /my/experiment -gica=ica_dirs_qa.txt -n run_name
 
     # ica_dirs_qa.txt is the FULL PATH to the single column text file listing full paths to ica directories
     # output goes to /my/experiment/gica/run_name.gica

#### Output
Group .gica Folder: will be under /my/experiment/gica/ with name groupname.gica. Within this folder, you will find a “log” folder with the cluster output and error logs, an “inputreg” folder with the compiled group masks and background image, and the meat of the analysis, the multi-session temporal concatenation results, are under groupmelodic.ica. See FSL's page on this type of analysis for details to the output.


### Dual Regression
Lastly, we will want to run dual regression. I have not modified this script beyond changing the name and including it with the package, the original documentation is here: http://www.fmrib.ox.ac.uk/analysis/dualreg/.

- You will first want to set up your general linear model to produce the design.con and design.mat files for your study. Open the FSL GUI, and click on “Misc” in the bottom left, and then “GLM Setup.” If you've used FSL before, this will be familiar. What the dual regression script does is basically run the randomise command on the data with your specified .con and .mat files. More details on GLM can be found at: http://www.fmrib.ox.ac.uk/fsl/feat5/glm.html

     python ica+.py -o /my/experiment --dr=group.gica --name=dr_name --con=design.con --mat=design.mat --iter=500
 
     # (must be done after group ICA is complete)
     # -o is the FULL PATH to the ica+ experiment directory
     # group.gica is the group gica directory under /my/experiment/gica/ - you do NOT need a full path here!
     # design.con and design.mat are created in the FSL GUI with GLM (general linear model)
     # dr_name is what you want to call it, which will be the name of the folder under /my/experiment/qa/
     # 500 is the number of iterations to run
     # output goes to /my/experiment/dr/dr_name

#### Output

Dual Regression Output Folder: will be under /my/experiment/dr/ with name dr_name (that you have specified). Within this folder, there is also a “log” folder with output and error files from the cluster submission, and the “meat” of your dual regression analysis is within the “scripts+logs” folder.

Keep in mind that running 500 iterations can take part of a day, but adding zeroes to that number can add days to the run! It's best to start conservatively (500) and increase as you see fit.
