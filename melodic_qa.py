#!/usr/bin/env python2
"""
 Melodic_qa: QA checking portion of ica+ package
 
 This python script does a quick quality analysis - it simply reads the mcflirt output and flags
 subjects with rotational or translational motion that exceeds a user set threshold, set below.
 It outputs a *_qa.txt file with a list of subjects ica directories that pass, for gica use.

McFlirt output utilized (under sub.ica/mc folder)
  prefiltered_func_data_mcf.par: contains rotation and translation motion parameters estimated by MCFLIRT, one row per volume
  This is the file that the script will utilize
  rot_x rot_y rotz tran_x tran_y tran_z

Usage: python melodic_qa.py -o /my/experiment --icas=input.txt --name=subgroup --rot=2.0 --tran=2.0
 
Options:
  -h, --help             show this help  
  -o, --output           experiment folder
  --icas                 single column file with list of ica directories to check qa for
  --name                 name for qa run (for output file under /my/experiment/qa
  --rot                  rotational motion benchmark (degrees)
  --tran                 translational motion benchmark (mm)

"""
import os
import sys
import getopt
import csv
import math

# ------------------------------------------------------------------------------------
class melodic_qa(Exception): pass
   
def usage():
    print __doc__

# Checks that variable has been defined, exits if has not
def varcheck(vartocheck):
    for varname, desname in vartocheck.items():
        if not varname:
            print "Missing variable " + desname + ".  Please specify and re-run!" 
            sys.exit()

# Checks again for each data file, in case script run separately from submission python
def checkData(inputfile):
    data = []
    fopen = open(inputfile,'r')
    for line in fopen:
        line = line.rstrip("\n").rstrip()
        if os.path.isfile(line + "/mc/prefiltered_func_data_mcf.par"):
            data.append(line)
        else: 
            print "Cannot find mc/prefiltered_func_data_mcf.par in " + line + ". Exiting."
            sys.exit()

    print "Found all mc/prefiltered_func_data_mcf.par files..."
    fopen.close
    return data
	    

# Read each file, compare against benchmark
def readData(inputdata,rotation,translation,outdir,qa_name):
    passing = []
    flagged = []

    # Convert degrees to radians
    pi = math.pi
    radians = pi * float(rotation) / 180
    print rotation + " degrees has been converted to " + str(radians) + " radians."	

    # Read from input file and compare motion parameters to benchmarks    
    for icadir in inputdata:
	rot = {"x":[],"y":[],"z":[]}
	tran = {"x":[],"y":[],"z":[]}
	f = open(icadir + "/mc/prefiltered_func_data_mcf.par", 'rb')
	reader = csv.reader(f, delimiter='\t')
        
        # Read all data from file into arrays, converting to float
        for row in reader:
            rot_x,rot_y,rot_z,tran_x,tran_y,tran_z = row[0].split()
            rot_x = float(rot_x)
            rot_y = float(rot_y)
            rot_z = float(rot_z)
            tran_x = float(tran_x)
            tran_y = float(tran_y)
            tran_z = float(tran_z)
            rot["x"].append(rot_x)
            rot["y"].append(rot_y)
            rot["z"].append(rot_z)
            tran["x"].append(tran_x)
            tran["y"].append(tran_y)
            tran["z"].append(tran_z)
	
        # Close the file     
        f.close()
    
        # Rotation: go through x,y,z, flag subject if greater than bench, add to passing list otherwise
	if flagCheck(rot,icadir,"rotation",radians,outdir,qa_name) and flagCheck(tran,icadir,"translation",translation,outdir,qa_name):
            passing.append(icadir)
        else:
            flagged.append(icadir)               

    # Combine flagged and passing list into a dictionary to return
    qa_output = {"passing":passing,"flagged":flagged}
	
    # Return list of passing ica directories, to be printed to file    
    return qa_output

# Flagcheck takes a dictionary of lists, returns False if a value exceeds benchmark
def flagCheck(timeseries,icadir,descriptor,benchmark,outdir,qa_name):
    for direction,valuelist in timeseries.items():
        for value in valuelist: 
            if abs(value) > float(benchmark): 
                flagSubject(icadir,descriptor,direction,outdir,qa_name)
                return False
    return True

# Prepare the output file
def setupOutfile(output,qa_name):
    if os.path.isfile(output + "/qa/" + qa_name + ".flag") or os.path.isfile(output + "/qa/" + qa_name + ".html") or os.path.isfile(output + "/list/" + qa_name + "_qa.txt"):
         print "QA output files with name " + qa_name + " already exist under " + output
         print "Delete old files or choose a different name. Exiting."
         sys.exit()
    else:    
        print "Creating new flagged subjects file..."
        flagfile = open(output + "/qa/" + qa_name + ".flag",'w').close()
        print "Creating new passing subjects file..."
        passfile = open(output + "/list/" + qa_name + "_qa.txt",'w').close()

# Check output directory
def setupout(outdir):
    # Make sure we don't end in slash
    if outdir[-1] == "/":
        outdir = outdir[:-1]
    # Check that experiment exists
    if not os.path.exists(outdir):
        print "Cannot find experiment directory " + outdir
        print "Check name and rerun! Exiting."
        sys.exit()
    else:
        if os.path.exists(outdir + "/qa"):
            print "QA output directory already exists."
        else:
            print "Creating output qa directory..."
            os.makedirs(outdir + "/qa")
    return outdir

# Print flagged subject
def flagSubject(icadir,flagtype,flagdirection,output,qa_name):
    flagfile = open(output + "/qa/" + qa_name + ".flag",'a')
    flagfile.writelines(flagtype + " " + flagdirection + ": " + icadir + "\n")
    flagfile.close()

# Print HTML report with motion charts for flagged subjects
def printReport(qa_output,output,qa_name,RM,TM):
    if not os.path.isfile(output + "/qa/" + qa_name + ".html"):
        print "Creating flagged subjects HTML report..."
        flagreport = open(output + "/qa/" + qa_name + ".html",'w')
        flagreport.write("<html>\n<body>\n<h1>ica+ Motion Report</h1>\n")
        flagreport.write("<p><strong>Experiment</strong>: " + output + "</p>\n")
	flagreport.write("<p><strong>QA Report Name</strong>: " + qa_name + "</p>\n")
        flagreport.write("<p><strong>Total Subjects</strong>: " + str(len(qa_output["passing"])+len(qa_output["flagged"])) + "</p>\n")
	flagreport.write("<p><strong>Flagged Subjects</strong>: " + str(len(qa_output["flagged"])) + "</p>\n")
	flagreport.write("<p><strong>Passing Subjects</strong>: " + str(len(qa_output["passing"])) + "</p>\n")
	flagreport.write("<p><strong>Rotational Motion Benchmark (deg)</strong>: " + str(RM) + "</p>\n")
	flagreport.write("<p><strong>Translational Motion Benchmark (mm)</strong>: " + str(TM) + "</p>\n")
	
	# Open flagged subject summary file, print to page
	if len(qa_output["flagged"]) > 0:
            flagreport.write("<h1>Flagged Subject Summary</h1>\n<p>")
	    flagfile = open(output + "/qa/" + qa_name + ".flag",'r')
	    for line in flagfile:
	        flagreport.write(line + "<br /><br />")
	    flagfile.close		
	    flagreport.write("</p>\n")


	    # Cycle through list of flagged ica directories, print name and links to motion images:
	    flagreport.write("<h1>Flagged Subjects</h1>\n<p>")
	    for flagged in qa_output["flagged"]:
   	        icadir_name = os.path.basename(flagged)
                flagreport.write("<p>" + icadir_name + "</p>\n")
                flagreport.write("<img src=\"../ica/" + icadir_name + "/mc/rot.png\" />\"\n")
                flagreport.write("<img src=\"../ica/" + icadir_name + "/mc/trans.png\" />\"\n")
                flagreport.write("<img src=\"../ica/" + icadir_name + "/mc/disp.png\" />\"\n")
  
	# Print list of passing subjects
        if len(qa_output["passing"]) > 0: 
	    flagreport.write("<h1>Passing Subjects</h1>\n<p>")
	    for passing in qa_output["passing"]:
   	        flagreport.write(passing + "<br /><br />\n")
	

        flagreport.write("</body>\n</html>")
        flagreport.close()

# Print passing subjects to file
def printPassing(output,qa_name,passingData):
    PASSfile = open(output + "/list/" + qa_name + "_qa.txt","w")
    for eachica in passingData[:-1]:
        PASSfile.writelines(eachica)
        PASSfile.writelines("\n")
    PASSfile.writelines(passingData[-1])
    PASSfile.close()

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "o:", ["output=","icas=","name=","rot=","tran="])

    except getopt.GetoptError:
        usage()
        sys.exit(2)
    
    # First cycle through the arguments to collect user variables
    subs = []
    outdir = None
    name = None
    TM = None	# translational motion benchmark (mm)
    RM = None   # rotational motion benchmark (degrees)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()    
        if opt in ("--icas"):
            icas = arg            
        if opt in ("-o", "--output"):
            outdir = arg
	if opt in ("--rot"):
            RM = arg
      	if opt in ("--tran"):
	    TM = arg
        if opt in ("--name"):
            runname = arg

    varcheck({icas:"input icas (--icas=dirs.txt)",outdir:"experiment output directory (-o)",RM:"rotation benchmark, degrees (--rot=2.0)",TM:"translation benchmark, mm (--tran=2.0)",runname:"name for qa run (--name=run_name)"})
    outdir = setupout(outdir)        # setup output directory
    setupOutfile(outdir,runname)     # ready output files
    data = checkData(icas)           # make sure files exist in all ica directories
    qa_output = readData(data,RM,TM,outdir,runname) # Read each file, compare against benchmarks
    printPassing(outdir,runname,qa_output["passing"])   # print passing subject IDs to ica list output file
    printReport(qa_output,outdir,runname,RM,TM) # print HTML page with motion charts for flags  
    print "Done checking QA."
    print "See /qa/" + runname + ".html for flagged subject overview."
    print "See /qa/" + runname + ".flag for flagged subject list."
    print "Use /list/" + runname + "_qa.txt for gica input file."
    sys.exit()

if __name__ == "__main__":
    main(sys.argv[1:])
