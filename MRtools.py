#!/usr/bin/env python2

"""

MRtools.py - Python tools for image manipulation and filtration

MRtools.Data:   Translate between images of different dimensions and formats
MRtools.Filter: Determine goodness of an input image and a frequency timeseries
MRtools.Match:  Return match score for two MRtools Data objects

Class to create a nifti image object that can be queued for values in raw coordinate 
space as well as MNI space.  Intended use is for a translation between
a query for a value in MNI space to the stored image data.  The class uses nibabel
and numpy to read and manipulate the data.  Currently only "result" images with 1
timepoint are supported, as the script is intended for comparison between 3D images.

SAMPLE USAGE: Make sure script is somewhere on your path

To use Data class:
>> import MRtools
>> Image = Mrtools.Data('myimage.nii.gz')
>> Image.mritoRCP([x,y,z])
>> Image.getValMNI([x,y,z])

To use Filter class with Group Network Image:
>> import MRtools
>> Image = MRtools.Data('myimage.nii.gz')
>> Filter = MRtools.Filter()
>> Filter.isGood(Image,'timeseries.txt','frequency.txt')

To use Match class with Image:
>> import MRtools
>> Template = Mrtools.Data('myimage.nii.gz')
>> Match = MRtools.Match(Template)
>> Match.setIndexCrit(">",0)
>> Match.genIndexMNI()
>> Contender = MRtools.Data('contender.nii.gz')
>> Match.addComp(Contender)
>> Match.doTemplateMatch()

"""

__author__ = "Vanessa Sochat (vsochat@stanford.edu)"
__version__ = "$Revision: 1.0 $"
__date__ = "$Date: 2011/09/09 $"
__license__ = "Python"

import os
import sys
import nibabel as nib
import scipy
import scitools.numpytools as scinu
import numpy as np
import operator
import getopt

# Data------------------------------------------------------------------------------
class Data:
    def __init__(self,imname):
        self.name = imname  # name of the image, as user has input
        self.path = None    # Full path to the image        
        self.img = None     # a nibabel Nifti object to hold image
        self.go = self.checkFile()

        self.xdim = 0
        self.ydim = 0
        self.zdim = 0       
        if self.go:
            self.readDim()
            self.data = []      # the raw Y data 
            self.aff = []       # Affine transformation matrix
            self.readData()
            self.readAff()
                        
            self.XYZ = []       # XYZ coordinates to match raw data
            self.RCP = []       # "raw coordinate points"
            self.readXYZ()

    def __repr__(self):
        return self.name

# READING DATA FUNCTIONS -------------------------------------------------------------
# CHECK FILE 
    def checkFile(self):
        if os.path.isfile(self.name):
            print "Reading in data " + self.name + "..."
            self.path = os.path.abspath(self.name)
            try: # Read in template image:
                self.img = nib.load(self.path)
            except:
                 print "Cannot read image data " + self.name
                 return False
        else:
            print "Cannot find " + self.name + ". Check the path."
            return False
        return True

# READ DATA
    def readDim(self):
        self.xdim = self.img.get_shape()[0]
        self.ydim = self.img.get_shape()[1]
        self.zdim = self.img.get_shape()[2]

# READ AFFINE TRANSFORMATION MATRIX
    def readAff(self):
        self.aff = scinu.mat(self.img.get_affine())
    

# Read raw image data
    def readData(self):
        # Read in all data to a temporary variable
        dataTEMP = self.img.get_data()
        # Check to see if we have a 4D image, represented by a 4th dimension
        if len(self.img.get_shape()) > 3:
            if self.img.get_shape()[3] > 1:
                print self.name + " has more than one timepoint, will try using first by default."        		
                if self.notEmpty(dataTEMP[:,:,:,0:1]):	
                    self.data = dataTEMP[:,:,:,0:1]
                else:
                    # Once or twice I've seen melodic spit out a component with 2 TPs, the first empty, and the second correct.
                    print "Data at timepoint 1 is an empty or errored image.  Trying next timepoints..."
                    if self.img.get_shape()[3] < 2:
                        print "The template is a 4D file but only has one timepoint that cannot be used.  Exiting!"
                        sys.exit()
                    else:
                        # Here we are checking timepoint 2, which likely has the map.  We could continue checking timepoints
		        # if two is empty, but this probably means something hugely wrong with the image, and we should stop
		        # and alert the user
                        if self.notEmpty(dataTEMP[:,:,:,1:2]):	
                            print self.name + " has empty first timepoint, using second."
                            self.data = dataTEMP[:,:,:,1:2]
                        else:
                            print self.name + " is empty at both timepoints 1 and 2, and we cannot use it.  Exiting!"
	    else:
                # Make sure image isn't empty, and set to dataTEMP
                if self.notEmpty(dataTEMP):
                    self.data = dataTEMP
                else:
                    print self.name + " is empty and cannot be used as a template.  Exiting."
                    sys.exit()

        # Otherwise, we have a 3D image and only one set of datapoints to try	
        else:
	    # Make sure that we don't have an empty image
            if self.notEmpty(dataTEMP):	
                self.data = dataTEMP
            else:
                print self.name + " is empty and cannot be used as a template!  Exiting."
                sys.exit()

# Check if data is empty
    def notEmpty(self,data):
        for p in range(0,np.shape(data)[0]-1):
            for o in range(0,np.shape(data)[1]-1):
                for d in range(0,np.shape(data)[2]-1):
                    if data[p,o,d] != 0:
                        return True
        return False


# Read XYZ Coordinates 
    def readXYZ(self):
        
        # Create coordinate space based on x,y,z dimensions, and multiply by affine matrix
        # Examples shown if xdim = 3, ydim=4, zdim=5
        # Create R row variable [1 2 3 1 2 3...] ydim * zdim times
        Rrow = list(scinu.seq(1,self.xdim)) * (self.ydim*self.zdim)

        # Create C row variable [1 1 1 2 2 2 3 3 3 4 4 4 1 1 1 2 2 2...] zdim X
        Crow = []
        for y in range(1,self.ydim+1):
            for x in range(0,self.xdim):
                Crow.append(y)
        Crow = Crow * self.zdim
	
        # Create P row variable [ each of 1:zdim xdim*ydim times ]
        Prow = []
        for z in range(1,self.zdim+1):
            holder = ([z] * self.xdim*self.ydim)
            for i in holder:
                Prow.append(i)

        # Create row of 1s of length zdim*xdim*ydim so we can multiply matrices
        onedim = [1] * self.xdim * self.ydim * self.zdim

        # Stack each row on top of one another
        self.RCP = np.vstack((Rrow,Crow,Prow,onedim))

        # Make it into a matrix
        self.RCP = scinu.mat(self.RCP)

        # Grab the first three rows of the affine transformation matrix (4th is for time)
        affXYZ = self.aff[0:3]

        # Multiply affine transformation matrix by coordinate data to go from coordinate --> MNI space
        self.XYZ = affXYZ * self.RCP

        # self.XYZ[:,0] contains first set of xyz coordinates (in column)
        # It is a 3xn matrix of XYZ locations returned (in mm) 

# HEADER DATA RETURN
    def getMeta(self,field):
        # Retrieve all raw header data from the image
        rawheader = self.img.get_header().structarr

        # See if the user specified value exists
        try:
            headerval = rawheader[field].tolist()
            return headerval
        except:
            print "Cannot find field " + field + " in header data for " + self.name

# IMAGE DATA RETURN - One coordinate

    def getValMNI(self,MNIcoord):
        '''Image.getValMNI([MNIx,MNIy,MNIz]) returns data value from MNI coordinate input'''
        # First convert from MNI to the images raw coordinate space
        MNIxyz = self.mnitoRCP(MNIcoord)
        try:
            return str(self.data[MNIxyz[0],MNIxyz[1],MNIxyz[2]])
        except:
            return 0
        #else:
        #    return self.data[MNIxyz[0],MNIxyz[1],MNIxyz[2]].tolist()[0]
        
    def getValRCP(self,RCPxyz):
        '''Image.getValRCP([RCPx,RCPy,RCPz]) returns data value from RCP coordinate input'''
        try: 
            return self.data[RCPxyz[0],RCPxyz[1],RCPxyz[2]]
        except:
            return 0 
            # print "RCP coordinate " + str(RCPxyz) + " undefined in " + self.name


# IMAGE DATA RETURN - All data
    def getData(self):
        '''Image.getData() returns entire raw data from nibabel object (in RCP space)'''
        return self.data

    def getXYZArray(self):
        '''Image.getXYZArray() returns entire coordinate matrix (as an array) (in MNI space)'''
        return np.array(self.XYZ)

    def getXYZMatrix(self):
        '''Image.getXYZMatrix() returns entire coordinate matrix (as a matrix) (in MNI space)'''
        return self.XYZ


# TRANSLATION MATRIX RETURN

    def getAffArray(self):
        '''Image.getAffArray() returns the affine transformation matrix for the image as a numpy array'''
        return np.array(self.aff)

    def getAffMatrix(self):
        '''Image.getAffMatrix() returns the affine transformation matrix for the image as a matrix'''
        return self.aff


# ALL COORDINATES RETURN

    def getRCPArray(self):
        '''Image.getRCPArray() returns all RCP coordinates xyz sets (as an array)'''
        return np.array(self.RCP)

    def getRCPMatrix(self):
        '''Image.getRCPMatrix() returns all RCP coordinates xyz sets (as a matrix)'''
        return np.array(self.RCP)


# COORDINATE TRANSLATIONS

    def mnitoRCP(self,coord):
        '''Image.mnitoRCP([x,y,z]) returns an RCP from an MNI coordinate input'''
        xcor = (coord[0] - self.aff[0,3]) / self.aff[0,0]
	ycor = (coord[1] - self.aff[1,3]) / self.aff[1,1]
	zcor = (coord[2] - self.aff[2,3]) / self.aff[2,2]
	return [xcor,ycor,zcor]

    def rcptoMNI(self,coord):
        '''Image.rcptoMni([x,y,z]) returns the an MNI coordinate from an RCP input'''
        xcor = (coord[0] * self.aff[0,0]) + self.aff[0,3]
	ycor = (coord[1] * self.aff[1,1]) + self.aff[1,3]
	zcor = (coord[2] * self.aff[2,2]) + self.aff[2,3]
	return [xcor,ycor,zcor]

# THRESHOLDING AND FILTERING

    def threshmin(self,threshmin):
        '''Image.thresh(threshval) returns an array of MNI coordinate above a particular threshold'''
        '''For MRtools Match object filter, use Match.setIndexCrit('>",0) and then Match.genIndexMNI()'''
        coords = []
        indexes = np.nonzero( self.data > threshmin )    # indexes[]
        # indexes[n][0] is x coordinate in mm, indexes[n][1] is y coordinate in mm, indexes[n][2] is z coordinate in mm
        
        # Cycle through all indexes and save list of MNI coordinates to return
        for i in range(0,len(indexes[0])):
	    # The index corresponds with raw coordinate space, so we need to convert each to MNI
	    coords.append(self.rcptoMNI([indexes[0][i],indexes[1][i],indexes[2][i]]))
        return coords

    def getMax(self):
        '''Image.getMax() returns the maximum activation value in an image'''
        # NOTE: this function has not yet been tested thoroughly
        maxval = 0
	for p in range(0,self.data.shape[0]-1):
            for o in range(0,self.data.shape[1]-1):
                for d in range(0,self.data.shape[2]-1):
                    if self.data[p,o,d] >= maxval:
		        maxval = self.data[p,o,d]
        return maxval

    def getUniqueIDs(self):
        '''Image.getUniqueIDs() returns all unique activation values in an image, likely corresponding to an atlas ID, NOT including 0'''
        uniques = []
	for p in range(0,self.data.shape[0]-1):
            for o in range(0,self.data.shape[1]-1):
                for d in range(0,self.data.shape[2]-1):
                    if self.data[p,o,d] not in uniques:
		        uniques.append(self.data[p,o,d])
        if 0 in uniques: uniques.remove(0) 
        return uniques


# Filter------------------------------------------------------------------------------
class Filter:
    '''High frequency filter'''
    def __init__(self):
        self.nframes = 0    # Signal length time (the number of frames)
        self.freq = 0       # Frequency is half the signal length time
        self.th_hfn = 50    # Threshold for high frequency noise, if frequency noise energy is 
                            # more than th_hfn % of the total energy signal, we have noise component
        self.hf_init_index = 25  # Chosen to remove components having most energy in range f > 0.1 Hz

    def __repr__(self):
        return "<Filter> " + self.Data

    # SET PARAMETERS
    def setSignalLength(self,signal):
        '''Filter.setSignalLength(240) sets the number of timepoints, starting with first, to use'''
        
        print "Signal length time (nframes) currently set at " + str(self.nframes)
        self.nframes = signal
        self.freq = self.nframes / 2
        print "<Signal: set at " + str(self.nframes) + "><Frequency: set at " + str(self.freq) + ">"

    def setHFNthresh(self,HFNthresh):
        '''Filter.setHFNthresh(60) sets the high frequency noise threshold, more than this % of total energy signal is a noise component'''
        ''' the default is set at 50 upon initialization of Filter'''
        
        print "High frequency noise threshold % currently set at " + str(self.th_hfn)
        self.th_hfn = HFNthresh
        print "<HFN Thresh: set at " + str(self.th_hfn) + ">"


    def setHFinitIndex(self,HFinitIndex):
        '''Filter.setHFinitIndex(25) sets the high frequency initial index - the range f > some Hz that you want to remove'''
        ''' the default is set at 25 = 0.1 Hz, upon initialization of Filter'''
        
        print "High frequency initial index currently set at " + str(self.hf_init_index)
        self.hf_init_index = HFinitIndex
        print "<HF Init Index: set at " + str(self.th_hfn) + ">"


    # RUN FILTER AND RETURN TRUE (GOOD) OR FALSE (BAD)        
    def isGood(self,MRData,tsfile,fqfile):
        '''Filter.isGood(MRData,tsfile,fqfile) returns true if a component passes high frequency filter'''

        # Check for ts and fq files... we don't currently use fqfile, but added in case algorithm changes
        for fcheck in (tsfile,fqfile):
            if not os.path.exists(fcheck):
                print "Cannot find timeseries file " + tsfile + ". Component will be skipped."
                return false;
            
        # TIMESERIES DATA
        # Read in the TS file, get the number of timepoints, save values
        tfile = open(tsfile,'r')
        self.nframes = sum(1 for line in tfile)
        tfile.close

        # Read in lines of data.  A group gica ts*.txt file will have multiple columns, one/subject
        # An individual ica ts*.txt file will just have one column.  Regardless, we use the mean
        time_data = np.zeros((1,self.nframes))
        counter = 0
        tfile = open(tsfile,'r')
        for line in tfile.readlines():
            tpstring = line.rstrip().split()
            floats = [float(x) for x in tpstring]
            tpmean = sum(floats) / len(floats)
            if counter < self.nframes:
                time_data[0,counter] = tpmean
                counter = counter + 1
        tfile.close()
                
        # FREQUENCY DATA
        # The frequency is half the number of timepoints 
        self.freq = (self.nframes / 2)
        freq_data = np.zeros((1,self.freq))
                
        # Calculate FFT (discrete fourier transform) of time data
        temp = abs(scipy.fft(time_data[0,:],self.nframes))
        #print "Temp FFT is " + str(temp)
        counter = 0
        for i in range(self.freq):
            freq_data[0,counter] = temp[counter]
            counter = counter + 1

        # Keep track of the following:
        total_energy = 0;
        high_freq_noise = 0;

        # Find the max value of the frequency data
        max_value = max(freq_data[0,:])
        print "Max value of frequency data is " + str(max_value)

        # Determine high frequency noise energy
        for i in range(self.hf_init_index-1,freq_data.shape[1]):
            high_freq_noise = high_freq_noise + (freq_data[0,i] * freq_data[0,i])
        
        print "High frequency noise energy is " + str(high_freq_noise)
       
        # Determine total energy and energy percent                
        for i in range(self.hf_init_index-1):
            total_energy = total_energy + (freq_data[0,i] * freq_data[0,i])

        total_energy = total_energy + high_freq_noise
        energy_percent = 100 * (high_freq_noise / total_energy)
     
        print "Total energy is " + str(total_energy)
        print "Energy percent is " + str(energy_percent)
        print "Acceptable is under " + str(self.th_hfn) + "%"

        # Compare energy percentage to user specified threshold
        if energy_percent > self.th_hfn:
            print "Noise component found! " + MRData.name + " will not be used."
            return False
        else:
            return True


# Match------------------------------------------------------------------------------
class Match:
    # Methods for matching entire images or overlays.
    # Components should be filtered for high frequency

    def __init__(self,Template):
        self.Data = Template    # MRtools Data object
        self.path = None        # Full path to the timeseries file        
        self.thresh = 0         # Set to zero by default
        self.filter = '>'   
        self.indexes = []
        self.coordsMNI = []
        self.coordsRCP = []
        self.components = []                  # List of components (MRtools Data objects) to check

        # Dictionaries to hold all results for one template across components
        self.activation_difference = {}       # Holds score with direction (+/-)
        self.activation_differenceabs = {}    # Holds absolute value of score, for ranking

    def __repr__(self):
        return "<Match>" + self.Data

    def addComp(self,MRData):
        '''Match.addComp(MRDataObj) adds a component to the list to be matched'''
        self.components.append(MRData)        

    def clearComp(self):
        '''Match.clearComp() clears component list'''
        self.components = []        

    def reset(self):
        '''Clears all components, results, and activation scores to prepare for next subject or set of component images'''
        self.activation_difference = {}
        self.activation_differenceabs = {}
        self.components = []
     
    def setIndexCrit(self,filt,thresh):
        '''setIndexCrit(filter,thresh) Set filter threshold (ie, 0) and filter (ie, <,>,=)'''
        self.thresh = thresh     
        # VANESSA add check here for filter input
        self.filter = filt        

    def genIndexMNI(self): 
        '''getIndex() sets the filter type and threshold, and calculates indices, converting to MNI coordinates'''
        '''For simple Data indexing outside of matching, use MRtools Data.threshmin(0)'''
        criteria = "self.Data.getData() " + str(self.filter) + " " + str(self.thresh)
        print "Filtering with criteria " + str(self.filter) + " " + str(self.thresh) + "..."
        self.indexes = np.nonzero( eval(criteria) )    
        if not self.indexes:
            print "Warning: No indexes found to match filter criteria!"
        else:
            # Save coordinates in both MNI and RCP space
            # NOTE - coordinates lookup in tempXYZ also tested, results were equivalent to 11th decimal point! 		
            for i in range(0,len(self.indexes[0])):
                self.coordsMNI.append(self.Data.rcptoMNI([self.indexes[0][i],self.indexes[1][i],self.indexes[2][i]]))
                self.coordsRCP.append([self.indexes[0][i],self.indexes[1][i],self.indexes[2][i]])

    def doTemplateMatch(self):
        '''doTemplateMatch() performs matching with Match.components, and coordinates Match.coordsMNI, for a specified subject ica directory'''
        '''Two dictionaries are returned containing activation difference scores (and absolute value of the scores) with component names as keys'''
        '''Algorithm calcs absolute value of (mean activation per voxel inside mask) - (mean activation / voxel outside mask)'''
        '''Good for finding voxel to voxel matching, not good for finding overlap'''
        # (Adapted from Kaustubh Supekar doTemplateMatching.m, 2008)
        # Set activation_difference dictionaries to keep track of absolute and absolute value scores, indexed by component name
        activation_difference = {}
        activation_differenceabs = {}

        print "\nCalculating shared and unshared activation per voxel for each contender image..."
        # Cycle through components and...
        for com in self.components:		      
            # Set activation counter variables to zero
            activation_in_roi = 0
            voxel_in_roi = 0
            activation_out_roi = 0
            voxel_out_roi = 0
            # Get the data in coordinate space
            data = com.getData()  
    
            # For each, take the coordinate list (in MNI) and convert to the raw coordinate space of the image
            coordsRCP = []
            for point in self.coordsMNI:
                coordsRCP.append(com.mnitoRCP(point))	     

            # SHARED ACTIVATION
            # For each point, try to look it up.  If we query an index that doesn't exist, this means
            # we don't have data for that point, and we don't use it in our similarity calculation.
            for point in coordsRCP:
                try:
                # Get the activation value at the point
                    if data[point[0],point[1],point[2]] != 0:
                        # If it isn't zero, then add to shared scoring
                        activation_in_roi = activation_in_roi + data[point[0],point[1],point[2]]
                        voxel_in_roi = voxel_in_roi + 1
                except:
                    print "Coordinate " + str(point) + " is not in " + com.name
                    print "...will not be included in similarity calculation!"
            
            # ACTIVATION IN IMAGE NOT IN TEMPLATE
            # Set the activation of each coordinate that we found to overlap as 0
            for point in coordsRCP:
                try:
                    data[point[0],point[1],point[2]] = 0
                except:
                    print "Coordinate " + str(point) + " is not in " + com.name
                    print "This coordinate will not be included in similarity calculation!"
              
            # Cycle through the data, and find voxels that still have activation
            for p in range(0,np.shape(data)[0]-1):
                for o in range(0,np.shape(data)[1]-1):
                    for d in range(0,np.shape(data)[2]-1):
                        if data[p,o,d] != 0:
                            # Make sure the template includes the point before including it
                            try:
                                checkpoint = self.Data.getData()[p,o,d]
                                activation_out_roi = activation_out_roi + data[p,o,d]
                                voxel_out_roi = voxel_out_roi + 1
                            except:
                                print "Point " + str(p) + " " + str(o) + " " + str(d) + " has activation, but not present in template ROI.  Will not count!"
		                     
            # Each subject will have an activation difference score for each component to the template.
            comname = os.path.basename(com.name.split('.')[0]) 
            if (voxel_in_roi == 0) and (voxel_out_roi == 0):
                activation_difference[com.name] = 0
                activation_differenceabs[com.name] = 0
                print comname + " does not have voxels with activation within template."
                print comname + " does not have voxels with activation outside of template."
            elif voxel_out_roi == 0:
                activation_difference[com.name] = (activation_in_roi/voxel_in_roi)
                activation_differenceabs[com.name] = abs(activation_in_roi/voxel_in_roi)
                print comname + " mean activation/voxel within template is " + str(activation_in_roi/voxel_in_roi)
                print comname + " does not have voxels with activation outside of template."
            elif voxel_in_roi == 0:
                activation_difference[com.name] = (0 - (activation_out_roi/voxel_out_roi))          
                activation_differenceabs[com.name] = abs(0 - (activation_out_roi/voxel_out_roi))
                print comname + " does not have voxels with activation within template."
                print comname + " mean activation/voxel outside of template is " + str(activation_out_roi/voxel_out_roi)
            else:
                activation_difference[com.name] = ((activation_in_roi/voxel_in_roi) - (activation_out_roi/voxel_out_roi))
                activation_differenceabs[com.name] = abs((activation_in_roi/voxel_in_roi) - (activation_out_roi/voxel_out_roi))
                print comname + " mean activation/voxel within template is " + str(activation_in_roi/voxel_in_roi)
                print comname + " mean activation/voxel outside of template is " + str(activation_out_roi/voxel_out_roi)
	
            print comname + " activation difference score: " + str(activation_difference[com.name])
            print comname + " absolute activation difference score: " + str(activation_differenceabs[com.name]) + "\n"
        return activation_difference,activation_differenceabs
           
    def doTemplateMatchV(self):
        '''doTemplateMatchV() performs matching with Match.components, and coordinates Match.coordsMNI, for a specified subject ica directory'''
        '''Two dictionaries are returned containing activation difference scores (and absolute value of the scores) with component names as keys'''
        '''Algorithm calcs (absolute mean activation per voxel inside mask) - (absolute mean activation / voxel outside mask)'''
        # Set activation_difference dictionaries to keep track of absolute and absolute value scores, indexed by component name
        activation_difference = {}
        activation_differenceabs = {}

        print "\nCalculating shared and unshared activation per voxel for each contender image..."
        # Cycle through components and...
        for com in self.components:		      
            # Set activation counter variables to zero
            activation_in_roi = 0
            activation_in_roiabs = 0
            voxel_in_roi = 0
            activation_out_roi = 0
            activation_out_roiabs = 0
            voxel_out_roi = 0
            # Get the data in coordinate space
            data = com.getData()  
    
            # For each, take the coordinate list (in MNI) and convert to the raw coordinate space of the image
            coordsRCP = []
            for point in self.coordsMNI:
                coordsRCP.append(com.mnitoRCP(point))	     

            # SHARED ACTIVATION
            # For each point, try to look it up.  If we query an index that doesn't exist, this means
            # we don't have data for that point, and we don't use it in our similarity calculation.
            for point in coordsRCP:
                try:
                # Get the activation value at the point
                    if data[point[0],point[1],point[2]] != 0:
                        # If it isn't zero, then add to shared scoring
                        activation_in_roi = activation_in_roi + data[point[0],point[1],point[2]]
                        activation_in_roiabs = activation_in_roiabs + abs(data[point[0],point[1],point[2]])
                        voxel_in_roi = voxel_in_roi + 1
                except:
                    print "Coordinate " + str(point) + " is not in " + com.name
                    print "...will not be included in similarity calculation!"
            
            # ACTIVATION IN IMAGE NOT IN TEMPLATE
            # Set the activation of each coordinate that we found to overlap as 0
            for point in coordsRCP:
                try:
                    data[point[0],point[1],point[2]] = 0
                except:
                    print "Coordinate " + str(point) + " is not in " + com.name
                    print "This coordinate will not be included in similarity calculation!"
              
            # Cycle through the data, and find voxels that still have activation
            for p in range(0,np.shape(data)[0]-1):
                for o in range(0,np.shape(data)[1]-1):
                    for d in range(0,np.shape(data)[2]-1):
                        if data[p,o,d] != 0:
                            # Make sure the template includes the point before including it
                            try:
                                checkpoint = self.Data.getData()[p,o,d]
                                activation_out_roi = activation_out_roi + data[p,o,d]
                                activation_out_roiabs = activation_out_roiabs + abs(data[p,o,d])
                                voxel_out_roi = voxel_out_roi + 1
                            except:
                                print "Point " + str(p) + " " + str(o) + " " + str(d) + " has activation, but not present in template ROI.  Will not count!"
		                     
            # Each subject will have an activation difference score for each component to the template.
            comname = os.path.basename(com.name.split('.')[0]) 
            if (voxel_in_roi == 0) and (voxel_out_roi == 0):
                activation_difference[com.name] = 0
                activation_differenceabs[com.name] = 0
                print comname + " does not have voxels with activation within template."
                print comname + " does not have voxels with activation outside of template."
            elif voxel_out_roi == 0:
                activation_difference[com.name] = (activation_in_roi/voxel_in_roi)
                activation_differenceabs[com.name] = abs(activation_in_roi/voxel_in_roi)
                print comname + " mean activation/voxel within template is " + str(activation_in_roi/voxel_in_roi)
                print comname + " absolute activation/voxel within template used for scoring is " + str(activation_in_roiabs/voxel_in_roi)
                print comname + " does not have voxels with activation outside of template."
            elif voxel_in_roi == 0:
                activation_difference[com.name] = (0 - (activation_out_roi/voxel_out_roi))          
                activation_differenceabs[com.name] = (0 - (activation_out_roiabs/voxel_out_roi))
                print comname + " does not have voxels with activation within template."
                print comname + " mean activation/voxel outside of template is " + str(activation_out_roi/voxel_out_roi)
                print comname + " absolute activation/voxel outside of template used for scoring is " + str(activation_out_roiabs/voxel_out_roi)
            else:
                activation_difference[com.name] = ((activation_in_roi/voxel_in_roi) - (activation_out_roi/voxel_out_roi))
                activation_differenceabs[com.name] = ((activation_in_roiabs/voxel_in_roi) - (activation_out_roiabs/voxel_out_roi))
                print comname + " mean activation/voxel within template is " + str(activation_in_roi/voxel_in_roi)
                print comname + " absolute activation/voxel within template used for scoring is " + str(activation_in_roiabs/voxel_in_roi)
                print comname + " mean activation/voxel outside of template is " + str(activation_out_roi/voxel_out_roi)
	        print comname + " absolute activation/voxel outside of template used for scoring is " + str(activation_out_roiabs/voxel_out_roi)
            print comname + " activation difference score: " + str(activation_difference[com.name])
            print comname + " absolute activation difference score: " + str(activation_differenceabs[com.name]) + "\n"
        return activation_difference,activation_differenceabs

    def matchOverlap(self):
        '''matchOverlap() performs matching with Match.components, and coordinates Match.coordsMNI, for a specified subject ica directory'''
        '''Two dictionaries are returned containing activation difference scores (and absolute value of the scores) with component names as keys'''
        '''Algorithm chooses "best" as having highest mean activation per voxel shared, not caring about activation in voxels not shared'''
        # Set activation_difference dictionaries to keep track of absolute and absolute value scores, indexed by component name
        activation_overlap = {}
        activation_overlapabs = {}

        print "\nCalculating shared activation per voxel for each contender image..."
        # Cycle through components and...
        for com in self.components:		      
            # Set activation counter variables to zero
            activation_in_roi = 0
            activation_in_roiabs = 0
            voxel_in_roi = 0
            # Get the data in coordinate space
            data = com.getData()  
    
            # For each, take the coordinate list (in MNI) and convert to the raw coordinate space of the image
            coordsRCP = []
            for point in self.coordsMNI:
                coordsRCP.append(com.mnitoRCP(point))	     

            # SHARED ACTIVATION
            # For each point, try to look it up.  If we query an index that doesn't exist, this means
            # we don't have data for that point, and we don't use it in our similarity calculation.
            for point in coordsRCP:
                try:
                # Get the activation value at the point
                    if data[point[0],point[1],point[2]] != 0:
                        # If it isn't zero, then add to shared scoring
                        activation_in_roi = activation_in_roi + data[point[0],point[1],point[2]]
                        activation_in_roiabs = activation_in_roiabs + abs(data[point[0],point[1],point[2]])
                        voxel_in_roi = voxel_in_roi + 1
                except:
                    print "Coordinate " + str(point) + " is not in " + com.name
                    print "...will not be included in similarity calculation!"
            
            # ACTIVATION IN IMAGE NOT IN TEMPLATE
            # We don't care for this algorithm
                     
            # Each subject will have an activation overlap score for each component to the template.
            comname = os.path.basename(com.name.split('.')[0]) 
            if (voxel_in_roi == 0):
                activation_overlap[com.name] = 0
                activation_overlapabs[com.name] = 0
                print comname + " does not have voxels with activation within template."
            else:
                activation_overlap[com.name] = (activation_in_roi/voxel_in_roi)
                activation_overlapabs[com.name] = abs(activation_in_roiabs/voxel_in_roi)
                print comname + " mean activation/voxel within template is " + str(activation_in_roi/voxel_in_roi)
                print comname + " absolute activation/voxel within template used for scoring is " + str(activation_in_roiabs/voxel_in_roi)
                print comname + " activation overlap score: " + str(activation_overlap[com.name])
            print comname + " absolute activation overlap score: " + str(activation_overlapabs[com.name]) + "\n"
        return activation_overlap,activation_overlapabs


# MAIN ----------------------------------------------------------------------------------
def main():
    print __doc__

if __name__ == "__main__":
    main()
