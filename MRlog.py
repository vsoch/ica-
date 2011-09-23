#!/usr/bin/env python2
"""
MRLog: prints log of imaging data information found in a particular directory.
FSL is used by default, expected to execute on command line as "fsl." To use
nibabel, make sure MRtools is in same folder, and specify --sof=nibabel

 
OPTIONS:
  -h, --help             show this help  
  -d, --dir              path to top data directory
  -s, --sof              software to read (fsl or nibabel)
  -o, --out              output folder for log
                         if not specified, uses pwd
USAGE:
python MRLog.py --out=/path/to/out --dir=/path/to/Data

OUTPUT: 
image,path,Ydim,file_type,timepoints,dims,Xpixdim,Zdim,Xdim,Ypixdim,descrip,units,Zpixdim, (fsl)
image,path,dims,units,Xpixdim,ydim,Zpixdim,file_type,zdim,timepoints,xdim,Ypixdim, (nibabel)

"""

__author__ = "Vanessa Sochat (vsochat@stanford.edu)"
__version__ = "$Revision: 2.1 $"
__date__ = "$Date: 2011/09/23 $"
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

#-----------------------------------------------------------------------------------
def usage():
    print __doc__

# SETUP FUNCTIONS
# Makes sure directory path doesn't end in slash
def checkSlash(indir):
    if indir[-1] == "/":
        indir = indir[:-1]
    return indir

# Check for directory (either data or output)
def checkDir(userdir):
    userdir = checkSlash(userdir)
    # If the directory doesn't exist, exit
    if not os.path.exists(userdir):
    	print "Directory " + userdir + " not found.  Exiting."
        sys.exit()
    return userdir

# Sets the software to use, fsl or nibabel
def setSoft(software):
    if software not in ("fsl","nibabel"):
        print "Error: software " + software + " is not supported."
        print "Supported types include fsl and nibabel.  Exiting!"
        sys.exit()
    else:
        print "Software specified is " + software
        return software    

# Check to see if the user has software to extract header data (FSL or nibabel)
def sofCheck(softocheck):
    found = False
    for soft in softocheck:
        print "Checking for " + soft + "..."
        if soft == "fsl":
            dirfsl = subprocess.Popen(['which','fsl'],stdout=subprocess.PIPE)
            if len(dirfsl.stdout.read()) is 0:
                print "Cannot find fsl installation."
            else:
                print "Found fsl installation."
                return "fsl"    
        elif soft == "nibabel":
            try:
                import MRtools
                print "Found MRtools module with nibabel."
                return "nibabel"
            except:
                print "Cannot find nibabel module."    
    print "Not found!  FSL or nibabel is required to create MRlog.  Exiting."
    sys.exit()

# Checks that variable has been defined, exits if has not
def varcheck(vartocheck):
    for varname, desname in vartocheck.items():
        if not varname:
            print "Missing variable " + desname + ".  Please specify and re-run!" 
            sys.exit()

# Prints headers based on software selection and fields to be extracted
def setupOutFile(outdir,soft,vals):
    # Setup output file name based on date and time
    mrdate = datetime.datetime.now()
    now = mrdate.strftime("%Y-%m-%d_%H_%M")
    logname = "MRlog-" + now + ".txt"
    
    try:
        outfile = open(outdir + "/" + logname,"w")
        outfile.write("image,path,")

        for descrip,val in vals[soft].iteritems():
            outfile.write(descrip + ",")
        outfile.write("\n")
        outfile.close()
        
        # Return the full path of the log
        return os.path.abspath(outdir + "/" + logname)
    except:
        print "Cannot open " + outdir + "/" + logname + " for writing.  Exiting!"
        sys.exit()

# Prints a completed entry to the output file
def addOutFile(line,outfile):
    try:
        fout = open(outfile,'a')
        fout.write(line + "\n")
        fout.close()
    except:
        print "Cannot print output to " + outfile + ".  Exiting!"
        sys.exit()

# EXTRACTION FUNCTIONS

def extractFSL(imagefiles,vals,outdir):
    for img in imagefiles:
        errorreading = False
        valtoget = None

        # Array to hold all values, first two entries are image name and path
        imgvals = os.path.basename(img) + "," + img

        for desc,fslval in vals.iteritems():
            try:
                valtoget = subprocess.Popen(['fslval',img,fslval],stdout=subprocess.PIPE)
                valtoget = valtoget.stdout.read().rstrip()
                # If it's empty, we should still print some empty space to the file
                if not valtoget:
                    valtoget = ' '
                imgvals = imgvals + "," + str(valtoget)
            except:
                errorreading = True
                continue
        if errorreading: print "Error reading image " + img + ". Skipping."
        # Print the line of data to the output file
        addOutFile(imgvals,outdir)

# Extract header data using nibabel
def extractNIB(imagefiles,vals,outdir):
    import MRtools
    for img in imagefiles:
        # Read in image to MRtools Data object
        image = MRtools.Data(img)
        errorreading = False
        valtoget = None

        # Array to hold all values, first two entries are image name and path
        imgvals = os.path.basename(img) + "," + img

        for desc,nibval in vals.iteritems():
            nibval = nibval.split(':')
            
            try:
                valtoget = image.getMeta(nibval[0])
                if isinstance(valtoget,(list)):
                    valtoget = valtoget[int(nibval[1])]

                # If it's empty, we should still print some empty space to the file
                if not valtoget:
                    valtoget = ' '
                imgvals = imgvals + "," + str(valtoget).rstrip()
            except:
                errorreading = True
                continue
        if errorreading: print "Error reading image " + img + ". Skipping."
        # Print the line of data to the output file
        addOutFile(imgvals,outdir)

# Get list of all image files below specified data directory
def getFiles(startdir,extensions):
    # Make sure we have read access
    if os.access(startdir,os.R_OK):
        imagefiles = []

        # Get list of files:
        for r, d, f in os.walk(startdir, topdown=True):
            for filey in f:
                # Look at the file extension to check for .nii or .img.  Do it twice in the case of .nii.gz
                if os.path.splitext(filey)[1] in extensions or os.path.splitext(os.path.splitext(filey)[0])[1] in extensions:
                    imagefiles.append(os.path.abspath(r  + "/" + filey))
                    print "Adding image file " + os.path.abspath(r + "/" + filey)
        return imagefiles

    else:
        print "Do not have permissions to access " + startdir + ".  Exiting!"
        sys.exit()

#----------------------------------------------------------------------------------------
def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hd:spo", ["help","dir=","sof=","pre=","out="])

    except getopt.GetoptError:
        usage()
        sys.exit(2)
    
    # PATH VARIABLES
    outdir = os.getcwd()               # if not specified, output will go to pwd
    soft = None                        # if not specified, will try fsl first
    extensions = ('.nii','.img')          # file extensions to include, currently hardcoded
    vals = dict()                      # Do not need to specify ".nii.gz"

    # OUTPUT VARIABLES
    # Variables to read from header using FSL.  Dictionary should be in form ('descriptor':'fslval').  Add more as you need!
    # Descriptors will be headers, separated by commas, in output file, fslvals should be fields in "fslinfo image.nii.gz"
    vals['fsl'] = {'timepoints':'dim4','dims':'dim0','xdim':'dim1','ydim':'dim2','zdim':'dim3','file_type':'file_type','descrip':'descrip','Xpixdim':'pixdim1','Ypixdim':'pixdim2','Zpixdim':'pixdim3','units':'vox_units'}   
    
    # Since nibabel will pull a list of raw data for one field (dim, for example) the format here is {'descriptor':'meta:list_location'}
    vals['nibabel'] = {'timepoints':'dim:4','dims':'dim:0','xdim':'dim:1','ydim':'dim:2','zdim':'dim:3','file_type':'magic:0','Xpixdim':'dim:5','Ypixdim':'dim:6','Zpixdim':'dim:7','units':'xyzt_units:0'}   

    # First cycle through the arguments to collect user variables
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()   
        if opt in ("--dir","-d"):
            dirtop = checkDir(arg)    
        if opt in ("--sof","-s"):
            soft = setSoft(arg)    
        if opt in ("--pre"):
            # Pre will allow user to specify prefix of file to include, not yet implemented
            # ...to be added when I need the functionality!
	    pre = arg
        if opt in ("--out","-o"):
	    outdir = checkDir(arg)
       
    # if user specified software to use, check for that one first, otherwise check both
    if not soft:
        soft = sofCheck(["fsl","nibabel"])
    else:
        soft = sofCheck([soft.rstrip()])    

    # Setup output file for printing results
    fullout = setupOutFile(outdir,soft,vals)
    
    # Get list of imaging files with full paths
    imagefiles = getFiles(dirtop,extensions)

    # Perform extraction...
    if soft is "fsl": extractFSL(imagefiles,vals['fsl'],fullout)
    if soft is "nibabel": extractNIB(imagefiles,vals['nibabel'],fullout)
 
    print "Done creating MRlog, located at " + fullout
    sys.exit()

if __name__ == "__main__":
    main(sys.argv[1:])
