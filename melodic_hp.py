#!/usr/bin/env python2

"""

melodic_hp - Highpass filtering of group (gica) network results for ica+

This python script uses MRtools to read in a list of network images (gica results) and
print out a report of images that pass high pass filtering based on user specified criteria.
The script is submit automatically by melodic_gp.sh to prepare a list of both acceptable
gica networks and subsequent dual regression component images (which represent significant
differences between two groups for the functional networks).  The images are expected to be
named thresh_zstat*.nii.gz and the timeseries files as ts*.txt. The script can also be run
on command line with the following functionality:

INPUT:
-h, --help      Print this usage

REQUIRED
-o  --output    Name of output folder.  For melodic ica will be /experiment/list/
    --ts        Full path to folder with corresponding timeseries for components
    --gica      Full path to gica directory with zstat images
    --name      Gica run name for output prefix (*_IC-hpfilter-good.txt) & (*_DR-hpfilter-good.txt)

OPTIONAL
    --nframes   Signal length time (nframes), do not specify to use all timepoints
    --fthresh   High frequency noise thresh, > this % of total energy signal is noise
                The default value is 50 - leave blank to use default

USAGE: 

python melodic_hp.py -o /exp/list --name=run1 --ts=/exp/gica/run1/groupmelodic.ica/report --gica=/exp/run1/gica/groupmelodic.ica/stats

OUTPUT: (name_IC-hpfilter-good.txt) and (name_DR-hpfilter-good.txt)

"""

__author__ = "Vanessa Sochat (vsochat@stanford.edu)"
__version__ = "$Revision: 1.0 $"
__date__ = "$Date: 2011/09/12 $"
__license__ = "Python"

import os
import sys
import MRtools # includes classes Data, Filter, and Match 
import operator
import getopt
import re
import dircache
import time


# Return list of thresh_zstat*.nii.gz image files from gica image directory
def getImages(gicadir):
        
        # Get all names of thresh_zstat*.nii.gz result images, the networks to be filtered, from gica directory
        gicaexpr = re.compile('thresh_zstat\d+[.]nii[.]gz')
        zstats = []
        
        # set self.images to this full path
        if os.path.exists(gicadir):
            for img in dircache.listdir(gicadir):
                if gicaexpr.match(img):
                    zstats.append(img)
                    print "Found image " + img + ". Added to list."
                    time.sleep(.2)
        else:
            print "Cannot find gica image directory " + gicadir + ". Exiting!"
            sys.exit()

        # If no thresh_zstat images found, exit
        if not zstats: 
            print "Error: no thresh_zstat*.nii.gz images found in " + gicadir + ". Exiting!"
            sys.exit()

        return zstats

# Check that directory exists, and doesn't end in "/"
def nixEndSlash(dirname):
    # Check that directory exists, period
    if not os.path.exists(dirname):
        print "Cannot find " + dirname + ". Exiting!"

    # Make sure we don't end in a slash
    if dirname[-1] == "/": return dirname[:-1]
    else: return dirname


# Print name_bestcomps.txt and name_drbestcomps.txt in output folder
def printRes(comps,output,outname):
    goodcompfile = open(output + "/" + outname + "_IC-hpfilter-good.txt","w") 
    gooddrfile = open(output + "/" + outname + "_DR-hpfilter-good.txt","w")       
    for threshnum,img in sorted(comps.iteritems()):
        goodcompfile.write(img + "\n")
        # The dual regression result count starts at 0, and IC images start at 1, so we subtract one from thresh image to get dr result
        if len(str(int(threshnum)-1)) is 1:
            gooddrfile.write('dr_stage3_ic000' + str(int(threshnum)-1) + '_tfce_corrp_tstat1.nii.gz\n')     
        elif len(str(int(threshnum)-1)) is 2:
            gooddrfile.write('dr_stage3_ic00' + str(int(threshnum)-1) + '_tfce_corrp_tstat1.nii.gz\n')
    goodcompfile.close()
    gooddrfile.close()


# MAIN ----------------------------------------------------------------------------------
def usage():
    print __doc__

def main(argv):

    # Prepare MRtools Filter object to do filtering...
    Filter = MRtools.Filter()

    try:
        opts, args = getopt.getopt(argv, "ho:", ["help","output=","ts=","gica=","name=","nframes","fthresh"])

    except getopt.GetoptError:
        usage()
        sys.exit(2)
    
    # First cycle through the arguments to collect user variables
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()    
        if opt in ("--ts"):
            timepath = nixEndSlash(arg)
        if opt in ("--name"):
            outname = arg
        if opt in ("--gica"):
            gicapath = nixEndSlash(arg)
	if opt in ("--nframes"):
            Filter.setSignalLength(arg)
        if opt in ("--fthresh"):
            Filter.setHFNthresh(arg)
        if opt in ("-o","--output"):
            output = nixEndSlash(arg)


    # FILTER COMPONENTS TO USE IN MATCHING
    # Get list of thresh_zstat*.nii.gz images in user specified gica image directory
    print "\nGetting list of images to filter in " + gicapath + "..."
    images = getImages(gicapath)
    
    # Dictionary "good" will index by zstat number
    good = {}

    # Create a MRtools Data object for each component to filter, add to good dictionary if passes
    # Also add to goodlist or badlist to summarize for user at the end!
    goodlist = []
    badlist = []

    for img in images:
        img_current = gicapath + "/" + img
        zstatnum = img.split('zstat')[1].split('.nii.gz')[0]
        ts_current = timepath + "/t" + zstatnum + ".txt"
        freq_current = timepath + "/f" + zstatnum + ".txt"
        try:
            # Use MRtools Filter class to determine if this component is "good"
            Contender = MRtools.Data(img_current)
            
            # If it's good, add to dictionary to print
            if Filter.isGood(Contender,ts_current,freq_current):
                print "GOOD: " + img + "\n"
                good[zstatnum] = img_current
                goodlist.append(zstatnum)
            else:
                print "BAD: " + img + "\n"
                badlist.append(zstatnum)
        except: 
	    print "Problem with reading " + img + " with MRtools for Filtering.  Exiting!"
            sys.exit()    
               
    # PRINT RESULTS
    print "Printing results to " + output
    if badlist: print "Bad components include: " + str(badlist)
    if goodlist: print "Good components include: " + str(goodlist)
    printRes(good,output,outname)

if __name__ == "__main__":
    main(sys.argv[1:])
