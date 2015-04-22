#!/usr/bin/env python2
"""Melodic+ Dual Regression Package for Python

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


ALWAYS REQUIRED:
  -o ..., --output=...   experiment folder where the following subdirectories will be created:
			   ica: (sub1.ica, sub2.ica, sub3.ica...)
			   qa: 
			   gica: (group1.gica, group2.gica, group3.gica...)
			   dr: (dr1, dr2, dr3...)
  --queue                name of slurm queue to run jobs (uses sbatch)
  --tr=                  TR (in seconds)

TO RUN SINGLE SUBJECT ICA:
   python ica+.py -o /my/experiment --ica=input.txt --queue=russpold --tr=0.72
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
   python ica+.py -o /my/experiment --gica=ica_dirs.txt --name=run_name --tr=0.72
     ica_dirs.txt is a single column text file listing full paths to ica directories
     you can use the *_ica.txt output file (from ica) under /my/experiment/list
     output goes to /my/experiment/gica/run_name.gica

TO RUN DUAL REGRESSION:
   python ica+.py -o /my/experiment --dr=group.gica --name=dr_name --con=design.con --mat=design.mat --iter=500 --queue=russpold
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

"""

__author__ = "Vanessa Sochat (vsochat@stanford.edu)"
__version__ = "$Revision: 3.1 $"
__date__ = "$Date: 2015/04/08 $"
__license__ = "Python"

from xml.dom import minidom
import random
import sys
import os
import datetime
import time
import getopt
import subprocess

class NoSourceError(Exception): pass

#----MELODIC-------------------------------------------------------------------------------
# created as an object to allow for future functionality to perform single and group in same run
class Melodic:
    def __init__(self):
        self.subs = []
        self.icas = []                   # Single subject ica directories created
	self.gica = []			 # Single subject ica directory input
        self.anat = []
	self.func = []
        self.timepoints = 0
	self.outdir = None

    def group(self,icadirs,gicaname,outdir,scriptinput,pyexec,filterscript,tr,queue="bigmem",mask=None):
        self.outdir = outdir
	
	# Create an output name based on date / time, if one not specified
	if not gicaname:		
	    now = datetime.datetime.now()
            gicaname = now.strftime("%Y-%m-%d_%H_%M")

	# Read in all single subject ica directories, make sure each exists
        print "Reading ICA directory input..."
        try:
	    readdata = open(icadirs,'r')
	    for line in readdata:
	        line = line.rstrip("\n").rstrip()
                if os.path.exists(line):
                    self.gica.append(line)
                else:
                    print "Cannot find ica directory %s. Exiting." %(line)
                    sys.exit()
            readdata.close
	except:
            print "Cannot open file %s. Exiting" %(icadirs)
            sys.exit()
	
	# Create output directories, exit if the run name already exists
	self.createGPout("%s/gica" %(outdir))
        gpout = "%s/gica/%s.gica" %(outdir,gicaname)
	if not os.path.exists(gpout):
            self.createGPout(gpout)
            self.createGPout(gpout + "/log")	
	else:
            print "Run %s already exists, and will not be overwritten." %(gpout)
            print "Specify a new name and re-run, or delete old run.  Exiting!"
            sys.exit()

	# Prepare list of input directories into one string
	# Also print list of data paths to file, for use with dual regression
        subinput = ''
        gica_file = open("%s/.filelist" %(gpout),"w")
        for gicapath in self.gica:
            subinput = subinput + " " + gicapath
            gica_file.write("%s/reg_standard/filtered_func_data\n" %(gicapath))
        gica_file.close()

        # Submit GICA script
        filey = "%s/log/gica.job" %(gpout)
        filey = open(filey,"w")
        filey.writelines("#!/bin/bash\n")
        filey.writelines("#SBATCH --job-name=%s_gica\n" %(gicaname))
        filey.writelines("#SBATCH --output=%s/log/gica.out\n" %(gpout))
        filey.writelines("#SBATCH --error=%s/log/gica.err\n" %(gpout))
        filey.writelines("#SBATCH --time=24:00:00\n")
        filey.writelines("#SBATCH --mem=48000\n")
        filey.writelines('%s %s %s %s %s %s "%s"'
%(scriptinput,gpout,pyexec,filterscript,tr,mask,subinput))
        filey.close()
        os.system("sbatch -p %s --qos=%s %s/log/ica.job" %(queue,queue,gpout))


    def single(self,anat,func,timepoints,outdir,script,tr,queue="normal",mask=None):
        self.outdir = outdir
        self.anat = anat
        self.func = func
        self.createGPout(outdir + "/ica")

        for sub,funcdata in sorted(self.func.items()):
            if self.createSSOut(sub):
              self.printSS(sub,self.anat[sub],funcdata)
              self.runICA(sub,funcdata,self.anat[sub],script,
                          "%s/ica/%s.ica" %(outdir,str(sub)),queue,tr,mask)
              # Save the full single subject directory path to the icas list
	      self.icas.append("%s/ica/%s.ica" %(outdir,str(sub)))
        
        # Format the date for printing
        now = datetime.datetime.now()
        now = now.strftime("%Y-%m-%d_%H_%M")

	# Print the icas list to the list folder, for group melodic input
	GPfile = open(self.outdir + "/list/" + now + "_ica.txt","w")
        for eachica in self.icas[:-1]:
            GPfile.writelines(eachica)
            GPfile.writelines("\n")
        GPfile.writelines(self.icas[-1])
        GPfile.close()
        
    def runICA(self,subject,funcinput,anatinput,scriptinput,ssoutdir,queue,tr,mask):
        # Write job file
        filey = "%s/log/ica.job" %(ssoutdir)
        filey = open(filey,"w")
        filey.writelines("#!/bin/bash\n")
        filey.writelines("#SBATCH --job-name=%s_ica\n" %(str(subject)))
        filey.writelines("#SBATCH --output=%s/log/ica.out\n" %(ssoutdir))
        filey.writelines("#SBATCH --error=%s/log/ica.err\n" %(ssoutdir))
        filey.writelines("#SBATCH --time=2-00:00\n")
        filey.writelines("#SBATCH --mem=64000\n")
        filey.writelines("%s %s %s %s %s %s" %(scriptinput,ssoutdir,funcinput,anatinput,tr,mask))
        filey.close()
        os.system("sbatch -p %s %s/log/ica.job" %(queue,ssoutdir))
        time.sleep(2) 

    def createGPout(self,outdir):
	if not os.path.exists(outdir):
    	    print "Creating output directory %s..." %(outdir)
            os.makedirs(outdir)	    	
	else:
	    print "Output directory %s already created." %(outdir)

    def createSSOut(self,subject):
        if not os.path.exists("%s/ica/%s.ica" %(self.outdir,str(subject))):
            os.makedirs("%s/ica/%s.ica" %(self.outdir,str(subject)))
            os.makedirs("%s/ica/%s.ica/log" %(self.outdir,str(subject)))	    	
	    return True
	else:
	    print "Output ica directory already exists for %s. Will not run!" %(subject)
	    return False

    def printSS(self,subject,anatprint,funcprint):
        # Print raw subject data to log txt file, for later use 
        SSfile = open("%s/ica/%s.ica/log/input.txt" %(self.outdir,str(subject)),"w")
        SSfile.write("%s,%s,%s" %(str(subject),anatprint,funcprint))
        SSfile.close()
    
	 
#-----SETUP--------------------------------------------------------------------------------
class Setup:
    def __init__(self):
	self.anat = {}
	self.func = {}
        self.subs = []        
	self.outdir = None
        self.timepoints = 0
	    
    def ica(self,scans,outdir):
	self.readData(scans,outdir)
        self.checkData()
        melodicinput = self.returnData()
        return melodicinput

    def readData(self,scans,outdir): 
        # Read text file of sub IDs, anatomical data, functional data
	print "Reading input data file..."
        try:
	    readdata = open(scans,'r')
	    for line in readdata:
	        subss,anatss,funcss = line.split(",")
                self.anat[subss.rstrip()] = anatss.rstrip("\n").rstrip()
		self.func[subss.rstrip()] = funcss.rstrip("\n").rstrip()
            readdata.close()

	    # Copy scans file to list folder, in case we want it again
            now = datetime.datetime.now()
            now = now.strftime("%Y-%m-%d_%H_%M")
            os.system("cp %s %s/list/%s_raw.txt" %(scans,outdir,now))
       
        except:
            print "Cannot open file %s. Exiting" %(scans)
            sys.exit()

    def checkData(self):
        # Check that all files exist for anat, exit if one missing
	print "Checking for all anatomical raw data..."
	for k,anatfile in self.anat.items():
            if not os.path.isfile(anatfile):
	        print "Cannot find %s. Exiting" %(anatfile)
                sys.exit()

        # Check that all files exist for func, exit if one missing
	print "Checking for all functional raw data..."
        for k,funcfile in self.func.items():
            if not os.path.isfile(funcfile):
                print "Cannot find %s. Exiting" %(funcfile)
                sys.exit()

	# Check that the orientation is LAS for the anatomical data
	print "Checking for LAS orientation of anatomical data..."
	for k,i in self.anat.items():
	    anat_orient = subprocess.Popen(['fslorient','-getorient',i],stdout=subprocess.PIPE)
            if anat_orient.stdout.read() != "RADIOLOGICAL\n":
		print "Anatomical %s is not in radiological (LAS) orientation. Exiting." %(i)
		sys.exit(2)
	
	# Check that the number of timepoints is equal for all functional runs, and check orientation
	print "Checking for equal timepoints between functional input..."
        func_standard = subprocess.Popen(['fslval',self.func.items()[0][1],'dim4'],stdout=subprocess.PIPE)
	func_standard = func_standard.stdout.read()

	for k,i in self.func.items():
	    func_check = subprocess.Popen(['fslval',i,'dim4'],stdout=subprocess.PIPE)
            func_check = func_check.stdout.read()
	    if func_check != func_standard:
		print "First subject resting data has %s timepoints" %func_standard
                print "Functional %s has %s timepoints." %(i,func_check)
                print "Timepoints must be equal! Exiting"
                sys.exit(2)
        print "The number of timepoints for all runs is %s" %func_standard
	self.timepoints = func_standard      # save this value to print to each ica folder

	print "Checking for LAS orientation of functional data..."		
	for k,i in self.func.items():
            func_orient = subprocess.Popen(['fslorient','-getorient',i],stdout=subprocess.PIPE)
            if func_orient.stdout.read() != "RADIOLOGICAL\n":
		print "Functional %s is not in radiological (LAS) orientation. Exiting." %i
		sys.exit(2)

	# return the number of timepoints
	self.timepoints = func_standard	
	
    def returnData(self):
	# Make sure the lists are the same length, as a first check
        if len(self.anat) == len(self.func):	
	    # Return list of anatomical files, list of functional files, and timepoints
	    return (self.anat,self.func,self.timepoints)
        else: 
            print "Error: There are %s anatomical input and %s functional paths.  Check input file and rerun." %(str(len(anat)),str(len(func)))
            sys.exit()

#----QUALITY ANALYSIS---------------------------------------------------------------
class QualityAnalysis:
    def __init__(self):
        self.outdir = None
        self.qaname = None
        self.data = None

    def createGPout(self,outdir):
        if not os.path.exists(outdir):
    	    print "Creating output directory " + outdir + "..."
            os.makedirs(outdir)	    	
	else:
	    print "Output directory " + outdir + " already created."

    def setup(self,icas,qaname,outdir,scriptinput):
        self.qaname = qaname

        # Create output QA directory
        print "Checking for output directory..."
        self.createGPout(outdir)
 	self.createGPout(outdir + "/qa")
        self.createGPout(outdir + "/qa/log")
        self.outdir = outdir

        # Check for all input files in ica directories...
        try:
            fopen = open(icas,'r')
            for line in fopen:
                line = line.rstrip("\n").rstrip()
                if not os.path.isfile(line + "/mc/prefiltered_func_data_mcf.par"):
                    print "Cannot find mc/prefiltered_func_data_mcf.par in " + line + ". Exiting."
                    sys.exit()

        except:
            print "Cannot open " + icas + " for reading.  Exiting"
            sys.exit()

        print "Found all mc/prefiltered_func_data_mcf.par files..."
        self.data = icas
        fopen.close

        # Submit QA python script to run QA, send same input text file
        subprocess.Popen(['bsub','-J',self.qaname + "_qa",'-o',self.outdir + "/qa/log/" + self.qaname + ".out",'-e',self.outdir + "/qa/log/" + self.qaname + ".err","python",scriptinput,"-o",self.outdir,"--icas=" + icas,"--rot=2.0","--tran=2.0","--name=" + self.qaname])
        time.sleep(2) 

#----DUAL REGRESSION---------------------------------------------------------------
class DualRegression:
    def __init__(self):
      self.gicadir = None
      self.drname = None
      self.outdir = None
      self.fullout = None
      self.iters = 0
      self.mat = None
      self.con = None
    
    def createGPout(self,outdir):
        if not os.path.exists(outdir):
    	    print "Creating output directory " + outdir + "..."
            os.makedirs(outdir)	    	
	else:
	    print "Output directory " + outdir + " already created."


    def runDR(self,gicadir,drname,outdir,con,mat,iters,scriptinput,queue):
         if os.path.exists(outdir):
             self.outdir = outdir
             self.createGPout(outdir + "/dr")
             if drname:
                 if not os.path.exists(outdir + "/dr/" + drname):
                     self.drname = drname
                     self.fullout = outdir + "/dr/" + drname
                     self.createGPout(self.fullout)
                     self.createGPout(self.fullout + "/log")
                 else: 
                     print outdir + "/dr/" + drname + " already exists"
                     self.dtOut()
             else:
                     self.dtOut()
         else: 
             print outdir + "does not exist.  Check path and rerun.  Exiting." 
             sys.exit()
         self.checkGICA(gicadir)
         self.checkDesign(con,mat,iters)
         self.submitDR(scriptinput,queue)

    def checkGICA(self,gicadir):
        if os.path.isfile(self.outdir + "/gica/" + gicadir + "/groupmelodic.ica/melodic_IC.nii.gz"):
            self.gicadir = self.outdir + "/gica/" + gicadir
        else:
            print "Cannot find" + self.outdir + "/gica/" + gicadir + "/groupmelodic.ica/melodic_IC.nii.gz.  Exiting."
            sys.exit()

    def checkDesign(self,con,mat,iters):
        # Check for contrast file
        if os.path.isfile(con):
            self.con = os.path.abspath(con)
        else:
            print "Cannot find " + con + ". Make sure it is in the same directory, OR specify full path.  Exiting."
            sys.exit()  

        # Check for design matrix
        if os.path.isfile(mat):
            self.mat = os.path.abspath(mat)
        else:
            print "Cannot find " + mat + ". Make sure it is in the same directory, OR specify full path.  Exiting."
            sys.exit()

        # Input iterations
        try:
            self.iters = int(iters)
            print "Dual regression will be run with " + str(self.iters) + " iterations."
        except:
            "Error setting iterations to " + iters + ". Exiting"
        

    def dtOut(self): 
         # If no name specified, or name in use, name directory by date/time
	 now = datetime.datetime.now()
         self.drname = now.strftime("%Y-%m-%d_%H_%M")
         self.fullout = self.outdir + "/dr/" + self.drname
         self.createGPout(self.fullout)
         self.createGPout(self.fullout + "/log")

    def submitDR(self,scriptinput,queue):
        # Get list of files from group run
        gica_filelist = open(self.gicadir + "/.filelist","r")
        sublist = ''
        for funcdata in gica_filelist:
            sublist = sublist + " " + funcdata.rstrip("\n").rstrip()
        sublist = '"' + sublist + '"'
        gica_filelist.close()

        # Print group gica directory to log, for later use with match
        gica_log = open(self.fullout + "/gica_name.txt","w")
        gica_log.write(self.gicadir)
        gica_log.close()

        melodic_mix = "%s/groupmelodic.ica/melodic_IC 1" %(self.gicadir)

        # Submit script for dual regression
        filey = "%s/log/dr.job" %(self.fullout)
        filey = open(filey,"w")
        filey.writelines("#!/bin/bash\n")
        filey.writelines("#SBATCH --job-name=%s_dr\n" %(self.drname))
        filey.writelines("#SBATCH --output=%s/log/dr.out\n" %(self.fullout))
        filey.writelines("#SBATCH --error=%s/log/dr.err\n" %(self.fullout))
        filey.writelines("#SBATCH --time=24:00:00\n")
        filey.writelines("#SBATCH --mem=48000\n")
        filey.writelines('%s %s %s %s %s %s %s %s'
%(scriptinput,self.gicadir,melodic_mix,self.mat,self.con,str(self.iter),"%s/result" %(self.fullout),self.gicadir))
        filey.close()
        os.system("sbatch -p %s %s/log/dr.job" %(queue,self.fullout))
        
#----MATCH-------------------------------------------------------------------------------
# do prep to submit pyMatch.py with user specified dual regression results and template image
class Match:
    def __init__(self,outdir,drname,template):
      self.drname = None              # Dual Regression run name under /experiment/dr/
      self.outdir = None              # Experiment output directory
      self.fullout = None             # Full output path /experiment/match/drname
      self.template = None            # template image file to match to
      self.images = None              # text file created for pyMatch with list of dr_stage3_corrp image names (under gica/run1/filter)
      self.drimfolder = None          # Dual regression image folder with stage3_corrp* images
      self.tempname = None            # Name of template image
      self.subs = None                # text file for pyMatch with directory name containing contender images, should be dual regression folder
      self.checkMatch(drname)
      self.setupDir(outdir,template)      
      self.setupMatch()

    def checkMatch(self,drrun):
    # Make sure dr_run exists
        if not os.path.exists(self.outdir + "/dr/" + drrun):
            print "Cannot find dual regression run " + drrun + ". Exiting!"
            sys.exit()         
        else: 
            self.drname = drrun
            self.drimfolder = self.outdir + "/dr/" + self.drname + "/result"
   
    # Create output folder
    def createDir(self,dirname):
        if not os.path.exists(dirname):
    	    print "Creating output directory " + dirname + "..."
            os.makedirs(dirname)	    	
	else:
	    print "Output directory " + outdir + " already created."

    
    def setupDir(self,outdir,template):
        # Check for output directories, create subdirectories 
        if os.path.exists(outdir):
            self.outdir = outdir
            self.createDir(outdir + "/match")
            if not os.path.exists(self.outdir + "/match/" + self.drname):
                self.fullout = self.outdir + "/match/" + self.drname
                self.createDir(self.fullout)
                self.createDir(self.fullout + "/log")
            else: 
                print "Output directory " + self.outdir + "/match/" + self.drname + " already exists."
                self.fullout = self.outdir + "/match/" + self.drname

            # Check for template, copy to output folder if it exists
            if os.path.isfile(template):
                if not os.path.isfile(self.fullout + "/" + drname + "-" + os.path.basename(template)):
                    self.tempname = os.path.basename(template)
                    shutil.copy(template,self.fullout + "/" + drname + "-" + self.tempname )
                    self.template = self.fullout + "/" + drname + "-" + self.tempname
                else:
                    print "Template already been used for this dr_run, please delete old results, or use a different template."
                    sys.exit()
            else:
                print "Cannot find " + template + ". Make sure path is correct, and re-run."
                sys.exit()
        else:
            print "Output directory " + self.outdir + " not found! Exiting."
            sys.exit()

    def setupMatch(self):
        # Find file of "good" dual regression images from original gica directory
        # We know the path of the dual regression run, and from these can load the original
        # gica run that the dr was done with from file...
        print "Finding group melodic directory that dual regression was run from..."
        try:
            gicanamefile = open(self.outdir + "/dr/" + self.drname + "/log/gica_name.txt","r")
            gicadir = gicanamefile.readline().rstrip()
            print "Group melodic directory is " + gicadir
        except:
            print "Cannot find gica_name.txt in " + self.outdir + "/dr/" + self.drname + "/log/.  Make sure this"
            print "text file exists with the full path to the gica directory used for the dual regression, and re-run."
            sys.exit()
        
        # DUAL REGRESSION IMAGES 
        # Read in DR image names from gica_DR-hpfilter-good.txt file in gicadir...
        print "Reading in dual regression images for networks that passed high pass filtering..."
        drgood = open(gicadir + "/filter/gica_DR-hpfilter-good.txt","r")    
        drgoodlist = []    
        for line in drgood:
            print "Adding " + line.rstrip() + " as a contender to be matched..."
            drgoodlist.append(line.rstrip())

        # Write these images and full path to file, for input to pyMatch.py
        print "Writing input file " + self.fullout + "/" + self.tempname + "-images.txt for pyMatch.py..."
        drfile = open(self.fullout + "/" + self.tempname + "-images.txt","w")      # Image input file for pyMatch.py
        for drimg in drgoodlist:
            drfile.write(drimg + "\n")
        drfile.close()

        # ICA IMAGES
        # Read in ICA image names from gica_IC-hpfilter-good.txt file in gicadir...
        print "Reading in gica IC images for later use to create AIM template..."
        icgood = open(gicadir + "/filter/gica_IC-hpfilter-good.txt","r")    
        icgoodlist = []    
        for line in icgood:
            print "Adding " + line.rstrip() + " to list of good networks..."
            icgoodlist.append(gicadir + "/groupmelodic.ica/stats/" + line.rstrip())

        # Write these images and full path to file, for input to AIMTemp.py
        print "Writing input file " + self.fullout + "/" + self.tempname + "-original-ic.txt for AIMTemp.py..."
        icfile = open(self.fullout + "/" + self.tempname + "-original-ic.txt","w")   # Input file for AIMTemp.py
        for icimg in icgoodlist:
            icfile.write(icimg + "\n")
        icfile.close()

        # set self.images to this full path
        self.images = self.fullout + "/" + self.tempname + "-images.txt"

        # Write dr_folder path to file for subs list input to pyMatch.py
        subfile = open(self.fullout + "/" + self.tempname + "-subs.txt","w")      # Subs/Output folder(s) list for pyMatch.py
        subfile.write(self.drimfolder + "\n")
        subfile.close()
        self.subs = self.fullout + "/" + self.tempname + "-subs.txt"

    def runMatch(self,scriptinput,pyexec):
        
        # Submit script to run Match
        subprocess.Popen(['bsub','-J',self.drname + "_match",'-o',self.fullout + "/log/" + self.tempname + ".out",'-e',self.fullout + "/log/" + self.tempname + ".err",'-W 99:30',pyexec,scriptinput,"--output=" + self.fullpath,"--subs=" + self.subs,"--template=" + self.template,"--images=" + self.images])
        time.sleep(2)


#----AIM-------------------------------------------------------------------------------
# create AIMTemp files from a Match run, including original template, and top three IC
# component network images (under groupmelodic.ica/stats/IC_* images) and dual regression
# results (dr_stage3_corrp* images under dr/results)

class AIM:
    def __init__(self,outdir,drname):
      self.drname = None              # Dual Regression run name under /experiment/dr/
      self.expdir = None              # Experiment output directory
      self.outdir = None              # aim/template output directory
      self.fullout = None             # Full output path /experiment/aim/drname
      self.template = None            # template image file to match to
      self.tempname = None            # Name of template image
      self.inputs = {}                # Dictionary of input image paths indexed by image name
      self.scores = {}                # Dictionary of match scores, if we ever need them
      self.setupDir(outdir,drname)      
     
    def setupDir(self,outdir,drname):
        # Make sure experiment exists
        if os.path.exists(outdir):
            self.expdir = outdir 
            self.createDir(outdir + "/aim")
        else:
            print "Cannot find " + outdir + ". Exiting!"
            sys.exit()

        # Make sure dr-run exists
        if os.path.exists(outdir + "/dr/" + drname):
            self.drname = drname
            print "Found dual regression run " + outdir + "/dr/" + drname
        else:
            print "Cannot find dual regression run " + outdir + "/dr/" + drname + ". Exiting!"
            sys.exit()

    # Create output folders
    def createDir(self,dirname):
        if not os.path.exists(dirname):
    	    print "Creating output directory " + dirname + "..."
            os.makedirs(dirname)	    	
	else:
	    print "Output directory " + outdir + " already created."

    # Check if image exists
    def exists(self,imagetocheck,infile):
        if os.path.exists(imagetocheck):
            print "Found image " + imagetocheck
            return True
        else:
            print "Cannot find image " + imagetocheck + " in file " + infile + ". Check path.  Exiting!"
            sys.exit()

    def DR(self,drinputfile):
        # Read dr_input file - should be single column with full path to template, followed by full paths to each dr image
        try:
            drinfile = open(drinputfile,"r")

            # First find the template image
            firstline = drinfile.readline().rstrip()
            template = firstline.split(':')[0]
            if self.exists(template,drinputfile):
                self.template = os.path.abspath(template)
                self.tempname = os.path.basename(template).split('.')[0]
                self.inputs[self.tempname + self.drname + '_TEMPLATE'] = self.template
                # The template is a perfect match to itself, so give it a ridiculously high score
                self.scores[self.tempname + self.drname + '_TEMPLATE'] = 99999
                # Make output folder to correspond to dual regression name
                createDir(self.expdir + "/" + self.drname)
                createDir(self.expdir + "/" + self.drname + "/log")
                self.fullout = self.expdir + "/" + self.drname
                
            # Now find the remaining dual regression images
            for line in drinfile:
                drimres = line.rstrip().split(':')[0]
                if self.exists(drimres):
                    # Add to list of input, indexed by image and dr_run name, if we find them
                    self.inputs[self.drname + "-" + os.path.basename(drimres).split('.')[0]] = os.path.abspath(drimres)
                    self.scores[self.drname + "-" + os.path.basename(drimres).split('.')[0]] = line.rstrip().split(':')[1]
            drinfile.close()
    
        except:
            print "Cannot read input file " + drinputfile + ". Exiting!"


    def DRandIC(self,drinputfile,icinputfile):
        # First read in dual regression results
        self.DR(drinputfile)

        # Now read in corresponding IC networks 
        try:
            icinfile = open(icinputfile,"r")
            for line in icinfile:
              icimres = line.rstrip()
              if self.exists(icimres,icinputfile):
              # Add to list of inputs, indexed by dr_run name and IC image name, if we find them
              # An entire network image does not have a match score
                  self.inputs[self.drname + "-" + os.path.basename(icimres).split('.')[0]] = os.path.abspath(icimres)
        except:
            print "Cannot read input file " + icinputfile + ". Exiting!"

    def runAIM(self,scriptinput,pyexec):
        print "Submitting individual AIM template jobs..."
        for outname,fullpath in self.images.iteritems():

            # Submit script to create AIM templates...
            print "Submitting job for " + outname + "..."
            subprocess.Popen(['bsub','-J',outname + "_aim",'-o',self.fullout + "/log/" + outname + ".out",'-e',self.fullout + "/log/" + outname + ".err",'-W 99:30',pyexec,scriptinput,"-o",self.fullpath,"--input=" + fullpath,"--name=" + outname])
            time.sleep(2)
            # Usage: python AIMTemp.py -o /fullpath/here --input=myimage.nii.gz --name=outname

#-----------------------------------------------------------------------------------
def usage():
    print __doc__

# Checks for fsl installation, exits if not found
def fslcheck():
    print "Checking for FSL installation..."
    dirfsl = subprocess.Popen(['which','fsl'],stdout=subprocess.PIPE)
    if len(dirfsl.stdout.read()) is 0:
        print "Cannot find FSL installation.  Exiting."
        sys.exit() 

# Finds .sh scripts to run jobs and remembers full path
def scriptcheck():
    scriptdict = {}
    for scriptname in ("melodic_ss.sh","melodic_gp.sh","melodic_dr.sh","melodic_qa.py","run_Bandpass.sh","Bandpass","AIMTemp.py","MRtools.py","pyMatch.py","melodic_hp.py"):
        if os.path.isfile(scriptname):
	    scriptdict[scriptname] = os.path.abspath(scriptname)
        else:
            print "Cannot find " + scriptname + ". Make sure it is in the same directory, and re-run."
            sys.exit()
    return scriptdict

# Setup experiment directories
def setupout(useroutdir):
    # Make sure we don't end in a slash
    if useroutdir[-1] == "/":
        useroutdir = useroutdir[:-1]

    # If the directory doesn't exist, make it.
    if not os.path.exists(useroutdir):
    	print "Creating output directory..."
        os.makedirs(useroutdir)
        os.makedirs(useroutdir + "/list")
        return useroutdir
    else:
	print "Output directory already exists."
        if not os.path.exists(useroutdir + "/list"):
            os.makedirs(useroutdir + "/list")
        return useroutdir	

# Checks that variable has been defined, exits if has not
def varcheck(vartocheck):
    for varname, desname in vartocheck.items():
        if not varname:
            print "Missing variable " + desname + ".  Please specify and re-run!" 
            sys.exit()

#----------------------------------------------
# MAIN SCRIPT STARTS RUNNING HERE  
#----------------------------------------------
def main(argv):
    try:
        opts, args = getopt.getopt(argv, "ho:", ["help","output=","qa=","ica=","gica=",
                           "dr=",'name=',"con=","mat=","iter=","match=","tr=","queue=","mask="])

    except getopt.GetoptError:
        usage()
        sys.exit(2)
    
    # First cycle through the arguments to collect user variables
    runtype = None
    outdir = None
    scans = None
    runname = None
    gicadir = None
    matchdrname = None
    mask = None

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()   
 
        if opt in ("--ica") and not runtype:
            scans = arg
            runtype = "ica"
	if opt in ("--gica") and not runtype:
            icadirs = arg
            runtype = "gica"
	if opt in ("--qa") and not runtype:
	    icadirs = arg
            runtype = "qa"
        if opt in ("--dr") and not runtype:
	    gicadir = arg
            runtype = "dr"
        if opt in ("--match") and not runtype:
	    matchdrname = arg
            runtype = "match"
        if opt in ("--aim") and not runtype:
	    aim = arg
            runtype = "aim" 

        if opt in ("-o", "--output"):
            outdir = arg        
        if opt in ("--tr"):
            tr = arg        
        if opt in ("--name"):
            runname = arg
        if opt in ("--template"):
            aimtemp = arg
        if opt in ("--con"):
            con = arg
        if opt in ("--mat"):
            mat = arg
        if opt in ("--iter"):
            iters = arg
        if opt in ("--ics"):
            ics = arg
        if opt in ("--queue"):
            queue = arg
        if opt in ("--mask"):
            mask = arg


    fslcheck()            	  # Check to make sure fsl is installed!
    scriptdict = scriptcheck()	  # look for required scripts
    outdir = setupout(outdir)     # setup output directory
    pyexec = sys.executable       # save path to current python executable, to pass along

    # SINGLE SUBJECT ICA
    if runtype is "ica":
        # Make sure we have an output directory and scans list
        try:
            varcheck({scans:"scan list file (--ica=file.txt)",outdir:"experiment output directory (-o)"})
        except:
            usage()
            sys.exit()
        icaRun = Setup()
        anat,func,timepoints = icaRun.ica(scans,outdir)  
        melRun = Melodic()
        melRun.single(anat,func,timepoints,outdir,scriptdict["melodic_ss.sh"],tr,queue,mask)
        print "Done submitting ICA jobs."
        print "Follow output at " + outdir + "/ica/"
        print "When complete, use " + outdir + "/list/*_ica.txt for qa or gica input file."
        sys.exit()

    # GROUP ICA
    elif runtype is "gica":
        print "Preparing to run GICA"
        # Make sure we have an output directory and single subject ica scan list
        try: varcheck({icadirs:"ica directory file (--gica=icadirs.txt)",outdir:"experiment output directory (-o)",runname:"name for run (--name=run_name)"})
        except:
             usage()
             sys.exit()
        melGP = Melodic()
        try:
            melGP.group(icadirs,runname,outdir,scriptdict["melodic_gp.sh"],pyexec,scriptdict["melodic_hp.py"],tr,queue,mask)
        except:
            melGP.group(icadirs,runname,outdir,scriptdict["melodic_gp.sh"],pyexec,scriptdict["melodic_hp.py"],tr,queue)
        print "Done submitting GICA job."
        print "Follow output at %s/gica/" %(outdir)
        print "HP Filtering will happen at end, good IC and DR image lists will be under gica/filter"
        
    # QUALITY ANALYSIS
    elif runtype is "qa":
        try: varcheck({icadirs:"ica directory file (--qa=icadirs.txt)",outdir:"experiment output directory (-o)",runname:"name for qa run (--name=qa_run)"})
        except:
            usage()
            sys.exit()
        qaRun = QualityAnalysis()
        qaRun.setup(icadirs,runname,outdir,scriptdict["melodic_qa.py"])
        
    # DUAL REGRESSION
    elif runtype is "dr":
        try: varcheck({gicadir:"gica directory (--dr=group.gica)",outdir:"experiment output directory (-o)",con:"design contrasts (--con=design.con",mat:"design matrix (--mat=design.mat)",iters:"iterations (--iter=500)"})
        except:
            usage()
            sys.exit()
        drRun = DualRegression()
        drRun.runDR(gicadir,runname,outdir,con,mat,iters,scriptdict["melodic_dr.sh"],queue)
        print "Dual Regression job submit."
        print "Output will be in /dr/" + runname
        
    # MATCH
    elif runtype is "match":
        import shutil
        import dircache
        try: varcheck({matchdrname:"dual regression results folder to match (--match=dr_name)",outdir:"experiment output directory (-o)",aimtemp:"template image (--template=image.nii.gz"})
        except:
            usage()
            sys.exit()
        print "Preparing Match Object to perform matching..."
        MatchRun = Match(outdir,matchdrname,aimtemp)
        MatchRun.runMatch(scriptdict["pyMatch.py"],pyexec)
        print "Match job submit."
        print "Output will be in /match/" + runname
        print "Top three matches (for excel import) in file " + aimtemp + "-beststats.txt"
        print "File (for ica+ --aim= input) is " + aimtemp + "-bestcomps.txt"  

    # AIM
    elif runtype is "aim":
        print "Preparing to run ica+AIM..."
        try: varcheck({aim:"*-bestcomps.txt file input under match/drname/ (--aim=template_bestcomps.txt)",outdir:"experiment output directory (-o)",ics:"original filtered network components (--ics=/path/to/match/template-original-ics.txt)",runname:"dual regression run name (--name=dr_run"})
        except:
            usage()
            sys.exit()
        AIMRun = AIM(outdir,runname)             # Setup AIM instance
        if not ics: AIMRun.DR(aim)       # If user didn't specify original IC networks to create AIM instances for 
        else: AIMRun.DRandIC(aim,ics)    # If user specified IC networks
        AIMRun.runAIM(scriptdict["AIMTemp.py"],pyexec)
        print "Finished submitting AIM xml template jobs."
        print "Look for output in " + outdir + "/aim/" + runname

    # USER FAIL
    else: 
        print "Error, you must specify a run type!"
        usage()
        sys.exit()
    
if __name__ == "__main__":
    main(sys.argv[1:])
