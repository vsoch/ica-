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
   python ica+.py -o /my/experiment -gica=ica_dirs.txt -n run_name
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

"""

__author__ = "Vanessa Sochat (vsochat@stanford.edu)"
__version__ = "$Revision: 1.0 $"
__date__ = "$Date: 2011/07/18 $"
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

    def group(self,icadirs,gicaname,outdir,scriptinput):
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
                    print "Cannot find ica directory " + line + ". Exiting."
                    sys.exit()
            readdata.close
	except:
            print "Cannot open file " + icadirs + " . Exiting"
            sys.exit()
	
	# Create output directories, exit if the run name already exists
	self.createGPout(outdir + "/gica")
        gpout = outdir + "/gica/" + gicaname + ".gica"
	if not os.path.exists(gpout):
            self.createGPout(gpout)
            self.createGPout(gpout + "/log")	
	else:
            print "Run " + gpout + " already exists, and will not be overwritten."
            print "Specify a new name and re-run, or delete old run.  Exiting!"
            sys.exit()

	# Prepare list of input directories into one string
	# Also print list of data paths to file, for use with dual regression
        subinput = ''
        gica_file = open(gpout + "/.filelist","w")
        for gicapath in self.gica:
            subinput = subinput + " " + gicapath
            gica_file.write(gicapath + "/reg_standard/filtered_func_data\n")
        gica_file.close()

        # Submit GICA script
        scommand = "bsub -J " + gicaname + "_gica -M 8000000 -o " + gpout + "/log/gica.out -e " + gpout + "/log/gica.err -R \"rusage[mem=8192]\" -W 50:30 " + scriptinput + " " + gpout + " \"" + subinput + "\""
        # R \"rusage[mem=6144]\"
        # print scommand
	subprocess.Popen(['%s' % scommand], shell=True, executable = "/bin/bash")
        time.sleep(2)	

    def single(self,anat,func,timepoints,outdir,script):
        self.outdir = outdir
        self.anat = anat
        self.func = func
        self.createGPout(outdir + "/ica")

        for sub,funcdata in sorted(self.func.items()):
            if self.createSSOut(sub):
              self.printSS(sub,self.anat[sub],funcdata)
              self.runICA(sub,funcdata,self.anat[sub],script,outdir + "/ica/" + str(sub) + ".ica")
              # Save the full single subject directory path to the icas list
	      self.icas.append(outdir + "/ica/" + str(sub) + ".ica")
        
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
        
    def runICA(self,subject,funcinput,anatinput,scriptinput,ssoutdir):
        # Submit subject ICA script	
	subprocess.Popen(['bsub','-J',str(subject) + "_ica",'-o',ssoutdir + "/log/ica.out",'-e',ssoutdir + '/log/ica.err',scriptinput,ssoutdir,funcinput, anatinput])
        time.sleep(2) 

    def createGPout(self,outdir):
	if not os.path.exists(outdir):
    	    print "Creating output directory " + outdir + "..."
            os.makedirs(outdir)	    	
	else:
	    print "Output directory " + outdir + " already created."

    def createSSOut(self,subject):
        if not os.path.exists(self.outdir + "/ica/" + str(subject) + ".ica"):
            os.makedirs(self.outdir + "/ica/" + str(subject) + ".ica")
            os.makedirs(self.outdir + "/ica/" + str(subject) + ".ica/log")	    	
	    return True
	else:
	    print "Output ica directory already exists for " + subject + ". Will not run!"
	    return False

    def printSS(self,subject,anatprint,funcprint):
        # Print raw subject data to log txt file, for later use 
        SSfile = open(self.outdir + "/ica/" + str(subject) + ".ica/log/input.txt","w")
        SSfile.write(str(subject) + "," + anatprint + "," + funcprint)
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
            os.system("cp " + scans + " " + outdir + "/list/" + now + "_raw.txt")
       
	except:
            print "Cannot open file " + scans + ". Exiting"
            sys.exit()

    def checkData(self):
        # Check that all files exist for anat, exit if one missing
	print "Checking for all anatomical raw data..."
	for k,anatfile in self.anat.items():
            if not os.path.isfile(anatfile):
	        print "Cannot find " + anatfile + ". Exiting"
                sys.exit()

        # Check that all files exist for func, exit if one missing
	print "Checking for all functional raw data..."
        for k,funcfile in self.func.items():
            if not os.path.isfile(funcfile):
                print "Cannot find " + funcfile + ". Exiting"
                sys.exit()

	# Check that the orientation is LAS for the anatomical data
	print "Checking for LAS orientation of anatomical data..."
	for k,i in self.anat.items():
	    anat_orient = subprocess.Popen(['fslorient','-getorient',i],stdout=subprocess.PIPE)
            if anat_orient.stdout.read() != "RADIOLOGICAL\n":
		print "Anatomical " + i + " is not in radiological (LAS) orientation. Exiting."
		sys.exit(2)
	
	# Check that the number of timepoints is equal for all functional runs, and check orientation
	print "Checking for equal timepoints between functional input..."
        func_standard = subprocess.Popen(['fslval',self.func.items()[0][1],'dim4'],stdout=subprocess.PIPE)
	func_standard = func_standard.stdout.read()

	for k,i in self.func.items():
	    func_check = subprocess.Popen(['fslval',i,'dim4'],stdout=subprocess.PIPE)
            func_check = func_check.stdout.read()
	    if func_check != func_standard:
		print "First subject resting data has " + func_standard + " timepoints"
                print "Functional " + i + " has " + func_check + " timepoints."
                print "Timepoints must be equal! Exiting"
                sys.exit(2)
        print "The number of timepoints for all runs is " + func_standard
	self.timepoints = func_standard      # save this value to print to each ica folder

	print "Checking for LAS orientation of functional data..."		
	for k,i in self.func.items():
            func_orient = subprocess.Popen(['fslorient','-getorient',i],stdout=subprocess.PIPE)
            if func_orient.stdout.read() != "RADIOLOGICAL\n":
		print "Functional " + i + " is not in radiological (LAS) orientation. Exiting."
		sys.exit(2)

	# return the number of timepoints
	self.timepoints = func_standard	
	
    def returnData(self):
	# Make sure the lists are the same length, as a first check
        if len(self.anat) == len(self.func):	
	    # Return list of anatomical files, list of functional files, and timepoints
	    return (self.anat,self.func,self.timepoints)
        else: 
            print "Error: There are " + str(len(anat)) + " anatomical input and " + str(len(func)) + " functional paths.  Check input file and rerun."
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


    def runDR(self,gicadir,drname,outdir,con,mat,iter,scriptinput):
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
         self.checkDesign(con,mat,iter)
         self.submitDR(scriptinput)

    def checkGICA(self,gicadir):
        if os.path.isfile(self.outdir + "/gica/" + gicadir + "/groupmelodic.ica/melodic_IC.nii.gz"):
            self.gicadir = self.outdir + "/gica/" + gicadir
        else:
            print "Cannot find" + self.outdir + "/gica/" + gicadir + "/groupmelodic.ica/melodic_IC.nii.gz.  Exiting."
            sys.exit()

    def checkDesign(self,con,mat,iter):
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
            self.iter = int(iter)
            print "Dual regression will be run with " + str(self.iter) + " iterations."
        except:
            "Error setting iterations to " + iter + ". Exiting"
        

    def dtOut(self): 
         # If no name specified, or name in use, name directory by date/time
	 now = datetime.datetime.now()
         self.drname = now.strftime("%Y-%m-%d_%H_%M")
         self.fullout = self.outdir + "/dr/" + self.drname
         self.createGPout(self.fullout)
         self.createGPout(self.fullout + "/log")

    def submitDR(self,scriptinput):
        # Get list of files from group run
        gica_filelist = open(self.gicadir + "/.filelist","r")
        sublist = ''
        for funcdata in gica_filelist:
            sublist = sublist + " " + funcdata.rstrip("\n").rstrip()
        sublist = '"' + sublist + '"'
        gica_filelist.close()

        # Submit script for dual regression
        subprocess.Popen(['bsub','-J',self.drname + "_dr",'-o',self.fullout + "/log/dr.out",'-e',self.fullout + '/log/dr.err','-W 99:30',scriptinput,self.gicadir + "/groupmelodic.ica/melodic_IC","1",self.mat,self.con,str(self.iter),self.fullout + "/result",self.gicadir])
        time.sleep(2)

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
    for scriptname in ("melodic_ss.sh","melodic_gp.sh","melodic_dr.sh","melodic_qa.py","run_Bandpass.sh","Bandpass"):
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
        opts, args = getopt.getopt(argv, "ho:", ["help","output=","qa=","ica=","gica=","dr=",'name=',"con=","mat=","iter="])

    except getopt.GetoptError:
        usage()
        sys.exit(2)
    
    # First cycle through the arguments to collect user variables
    runtype = None
    outdir = None
    scans = None
    runname = None

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()    
        if opt in ("--ica") and not runtype:
            scans = arg
            runtype = "ica"
        if opt in ("-o", "--output"):
            outdir = arg
	if opt in ("--gica") and not runtype:
            icadirs = arg
            runtype = "gica"
	if opt in ("--qa") and not runtype:
	    icadirs = arg
            runtype = "qa"
        if opt in ("--dr") and not runtype:
	    gicadir = arg
            runtype = "dr"
        if opt in ("--name"):
            runname = arg
        if opt in ("--con"):
            con = arg
        if opt in ("--mat"):
            mat = arg
        if opt in ("--iter"):
            iter = arg

    fslcheck()            	  # Check to make sure fsl is installed!
    scriptdict = scriptcheck()	  # look for required scripts
    outdir = setupout(outdir)     # setup output directory

    # SINGLE SUBJECT ICA
    if runtype is "ica":
        # Make sure we have an output directory and scans list
        try:
            varcheck({scans:"scan list file (--ica=file.txt)",outdir:"experiment output directory (-o)"})
        except:
            usage()
            sys.exit()
        icaRun = Setup()
        melodicinput = icaRun.ica(scans,outdir)  
        anat = melodicinput[0]
        func = melodicinput[1]
        timepoints = melodicinput[2]
        melRun = Melodic()
        melRun.single(anat,func,timepoints,outdir,scriptdict["melodic_ss.sh"])
        print "Done submitting ICA jobs."
        print "Follow output at " + outdir + "/ica/"
        print "When complete, use " + outdir + "/list/*_ica.txt for qa or gica input file."
        sys.exit()

    # GROUP ICA
    elif runtype is "gica":
        # Make sure we have an output directory and single subject ica scan list
        try: varcheck({icadirs:"ica directory file (--gica=icadirs.txt)",outdir:"experiment output directory (-o)",runname:"name for run (--name=run_name)"})
        except:
             usage()
             sys.exit()
        melGP = Melodic()
        melGP.group(icadirs,runname,outdir,scriptdict["melodic_gp.sh"])
        print "Done submitting GICA job."
        print "Follow output at " + outdir + "/gica/"
        
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
        try: varcheck({gicadir:"gica directory (--dr=group.gica)",outdir:"experiment output directory (-o)",con:"design contrasts (--con=design.con",mat:"design matrix (--mat=design.mat)",iter:"iterations (--iter=500)"})
        except:
            usage()
            sys.exit()
        drRun = DualRegression()
        drRun.runDR(gicadir,runname,outdir,con,mat,iter,scriptdict["melodic_dr.sh"])
        print "Dual Regression job submit."
        print "Output will be in /dr/" + runname

    # USER FAIL
    else: 
        print "Error, you must specify a run type!"
        usage()
        sys.exit()
    
if __name__ == "__main__":
    main(sys.argv[1:])
