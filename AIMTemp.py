#!/usr/local/bin/python

"""AIMTemplate for Python

This script takes a nifti image as input and using the MRtrans class, creates
a template XML file for use with AIM (Annotation and Image Markup).  The 
output image name, if not specified, will be the same as the input image.xml
To produce standard output, the output will be of the same dimensions as the
aalMNI152 template.  If the input image is missing data for a coordinate, or
if the value is zero, it will not be written to the xml file.  To have consistent
image dimensions for results, the input image will be labeled based on the
dimensions of the aal template.  If the input image goes outside of the aal 
template, this area will not be included.

Usage: python AIMTemp.py -o /fullpath/here --input=myimage.nii.gz --name=outname
 
Main Options:
  -h, --help             show this help  
  -i, --input            input nifti file
  -o, --out              output folder (created if doesn't exist)
  -n, --name             output name

Optional:
  --rdf                  full path to rdf file to map AALID to FMAID
  --url                  url to map this path
"""

__author__ = 'Nolan Nichols'
__credits__ = ["Vanessa Sochat"]
__version__ = "$Revision: 1.0 $"
__date__ = "$Date: 2011/09/04 $"
__license__ = "Python"

from lxml import etree                          # module to parse xml
from urllib2 import urlopen                     # module to open rdf at QueryService url
from rdflib import ConjunctiveGraph as Graph    # module to work with rdf mapping between AAL and FMAIDs
import sys
import os
import getopt
import scitools.numpytools as scinu
import MRtools

#----AIM-TEMPLATE-------------------------------------------------------------------------------
class AIMTemplate:
    def __init__(self,infile):
        self.infile = None                                          
        self.FMRI = MRtools.Data(infile)              # Read the input file into a MRTools Data object, for easy query 
        self.AAL = MRtools.Data('MR/aal2mni152.nii.gz')  # Read the MNI152 template with labels into MRTrans object
                                                # single analyze volume in MNI152 space with integers 1-116 for anat labels
                                                # In future this can come from online location, now is hard coded file
        self.aalID = self.getAALs()             # An array of all IDs found in the AAL input image
        self.xyzlabels = self.voxelsByLabel()   # Return an array of lists of xyzlabels, one for each label found in the atlas
        self.aimTree = None
        
    # Read aal image to get array of aalIDs
    def getAALs(self):
        '''returns a list of unique values in the atlas image - NOT including 0'''
        print "Atlas image has unique values: " + str(sorted(self.AAL.getUniqueIDs()))
        return self.AAL.getUniqueIDs()        

    # Return a list of inputs to configure aimInstance based on anatomical entity (i.e., aalID)
    def voxelsByLabel(self):
        '''atlas is a NiBabel 3D array in MNI space. aalID is an integer 1-116 corresponding to a brain structure'''
        print "Labeling voxels in fmri image..."
        xyzlabels = [] # list that stores 5 items... (x,y,z coordinates, aalID, and zScore)
        
        # Go through the entire dimensions of the aal template, and query the data for a value  
        for x in range(self.AAL.xdim):
            for y in range(self.AAL.ydim):
                for z in range(self.AAL.zdim):
                    aalID = self.AAL.getValRCP([x,y,z])                 # Get the aalID value from the RCP coordinate

                    if aalID in self.aalID:    # Ensures we don't look at values of zero:
                        try:                                                # Add fMRI coordinate value only if it exists
                            MNIX,MNIY,MNIZ = self.AAL.rcptoMNI([x,y,z])     # Get the MNI coordinate for this point
                            fmriVal = self.FMRI.getValMNI([MNIX,MNIY,MNIZ]) # Get the value for this MNI coordinate from input data
                            if fmriVal not in ('0.0'):
                                print "MNI Coordinates " + str(MNIX) + ", " + str(MNIY) + ", " + str(MNIZ) + " have value " + str(fmriVal)
                                xyzlabels.append((MNIX, MNIY, MNIZ, aalID, fmriVal))                                  
                        except:
                            # print "fmri data is undefined at " + str(MNIX) + ", " + str(MNIY) + ", " + str(MNIZ) + " and will not be added"
                            continue
        print "Found " + str(len(xyzlabels)) + " coordinates in image input with activation." 
        return xyzlabels

    # Iterate through voxelsByLabel using the lookup indices as input to configure AIM instance
    def aimGen(self,aalDict):
        print "Configuring AIM Instance..."
        AIM_ROOT = etree.Element('AIM-ROOT') # create a root AIM element for ImageAnnotations returned by aimInstance
        undefinedAAL = []                    # keep a list of aalIDs undefined and defined in rdf...
        definedAAL = []
        recordcount = 0                      # Keep a count of the number of voxels in the image to label each one
        recordsize = len(self.xyzlabels)
        while self.xyzlabels:

            # Each entry in xyzlabels looks like: ([MNIX, MNIY, MNIZ, aalID, fmriVal])
            recordcount = recordcount + 1
            recordindex = str(recordcount) + "/" + str(recordsize)
            record = self.xyzlabels.pop()
            x = record[0]          # MNIX
            y = record[1]          # MNIY
            z = record[2]          # MNIZ
            aalID = record[3] # aalDict key is a string aalID
            Zscore = record[4]

            # Each entry in aalDict looks like: aalDict['aalID'] = ('fmaName', 'aalName', 'FMAID'))
            #                                   aalDict['77'] = ('Left thalamus', 'Thalamus_LEFT', '258716')


	    # Prepare FMA / ontology information
            cagridId = '%s-%s-%s' % (x,y,z)

            try:
                # We go to "except" if the aalID isn't a valid key, meaning it's not in the dictionary
                fmaid = aalDict[str(aalID)][2] # mapping between aalID and fmaid
                fmaLabel = aalDict[str(aalID)][0]
                ImageAnnotation = self.aimInstance(x, y, z, fmaid, fmaLabel, Zscore, cagridId, recordindex)
                AIM_ROOT.append(ImageAnnotation)
                if str(aalID) not in definedAAL: definedAAL.append(str(aalID))
            except:
                if str(aalID) not in undefinedAAL: undefinedAAL.append(str(aalID))
                continue

        self.aimTree = etree.ElementTree(AIM_ROOT) # create an element tree from AIM_ROOT returned from aimInstance
        print "AIM Instance configuration complete."
        
        # Report to the user the aalIDs that were defined, and not defined
        if definedAAL:
            print "AALIDs with activation, in atlas, found in AAL dictionary:"
            print definedAAL    

        if undefinedAAL:
            print "AALIDs with activation, in atlas image, but not in the dictionary, NOT added to AIM:"
            print undefinedAAL    
        return self.aimTree

    # return an AIM instance based on the required coordinates, labels, and statistics
    def aimInstance(self,x, y, z, fmaid, fmaLabel, zScore, cagridId, record):
        ''' '''
        # Global variables for AIM template
        XSI_NS = 'http://www.w3.org/2001/XMLSchema-instance'
        XSI_TYPE = '{%s}type' % XSI_NS
        # create an AIM ImageAnnotation template, ImageAnnotation sub-tree w/attributes
        ImageAnnotation = etree.Element('ImageAnnotation')
        ImageAnnotation.attrib['aimVersion'] = '3.0'
        ImageAnnotation.attrib['cagridId'] = cagridId
        ImageAnnotation.attrib['codeMeaning'] = 'MRI 3D Image'
        ImageAnnotation.attrib['codeValue'] = 'birnlex_2033'
        ImageAnnotation.attrib['codingSchemeDesignator'] = 'NeuroLEX'
        ImageAnnotation.attrib['dateTime'] = ''
        ImageAnnotation.attrib['name'] = record
        ImageAnnotation.attrib['uniqueIdentifier'] = ''
        # configure namespaces
        schemaLocation = '{%s}schemaLocation' % XSI_NS
        ImageAnnotation.attrib[schemaLocation] = 'gme://caCORE.caCORE/3.2/edu.northwestern.radiology.AIM AIM_v3_rv11_XML.xsd'
        # create an AIM ImageAnnotation template - calculationCollection sub-tree w/attributes
        calculationCollection = etree.SubElement(ImageAnnotation, 'calculationCollection')
        Calculation = etree.SubElement(calculationCollection, 'Calculation')
        Calculation.attrib['cagridId'] = cagridId
        Calculation.attrib['codeMeaning'] = 'Z Score'
        Calculation.attrib['codeValue'] = 'ZSCORE'
        Calculation.attrib['codingSchemeDesignator'] = 'FMRI'
        Calculation.attrib['description'] = 'Z Score'
        Calculation.attrib['uid'] = ''
        # create an AIM ImageAnnotation template - referencedCalculationCollection *UNUSED*
        referencedCalculationCollection = etree.SubElement(Calculation, 'referencedCalculationCollection')
        # create an AIM ImageAnnotation template - calculationResultCollection sub-tree w/attributes
        calculationResultCollection = etree.SubElement(Calculation, 'calculationResultCollection')
        CalculationResult = etree.SubElement(calculationResultCollection, 'CalculationResult')
        CalculationResult.attrib['cagridId'] = cagridId
        CalculationResult.attrib['numberOfDimensions'] = '1'
        CalculationResult.attrib['type'] = 'Scalar'
        CalculationResult.attrib['unitOfMeasure'] = 'Zscore'
        calculationDataCollection = etree.SubElement(CalculationResult, 'calculationDataCollection')
        CalculationData = etree.SubElement(calculationDataCollection, 'CalculationData')
        CalculationData.attrib['cagridId'] = cagridId
        CalculationData.attrib['value'] = str(zScore) # fMRI z-score goes for a single coordinate
        # create an AIM ImageAnnotation template - coordinateCollection sub-tree w/attributes
        coordinateCollection = etree.SubElement(CalculationData,'coordinateCollection')
        Coordinate = etree.SubElement(coordinateCollection, 'Coordinate')
        Coordinate.attrib['cagridId'] = cagridId
        Coordinate.attrib['dimensionIndex'] = '0'
        Coordinate.attrib['position'] = '0'
        # create an AIM ImageAnnotation template - dimensionCollection sub-tree w/attributes
        dimensionCollection = etree.SubElement(CalculationResult,'dimensionCollection')
        Dimension = etree.SubElement(dimensionCollection, 'Dimension')
        Dimension.attrib['cagridId'] = cagridId
        Dimension.attrib['index'] = '0'
        Dimension.attrib['label'] = 'Value'
        Dimension.attrib['size'] = '1'
        # create an AIM ImageAnnotation template - user sub-tree w/attributes
        user = etree.SubElement(ImageAnnotation, 'user')
        User = etree.SubElement(user, 'User')
        User.attrib['cagridId'] = cagridId
        User.attrib['loginName'] = ''
        User.attrib['name'] = ''
        User.attrib['numberWithinRoleOfClinicalTrial'] = ''
        User.attrib['roleInTrial'] = 'N/A'
        # create an AIM ImageAnnotation template - equipment sub-tree w/attributes
        equipment = etree.SubElement(ImageAnnotation, 'equipment')
        Equipment = etree.SubElement(equipment, 'Equipment')
        Equipment.attrib['cagridId'] = cagridId
        Equipment.attrib['manufacturerModelName'] = 'N/A'
        Equipment.attrib['manufacturerName'] = 'University'
        Equipment.attrib['softwareVersion'] = '3.0.0.0'
        # create an AIM ImageAnnotation template - anatomicEntityCollection sub-tree w/attributes
        anatomicEntityCollection = etree.SubElement(ImageAnnotation, 'anatomicEntityCollection')
        AnatomicEntity = etree.SubElement(anatomicEntityCollection, 'AnatomicEntity')
        AnatomicEntity.attrib['annotatorConfidence'] = '1'
        AnatomicEntity.attrib['cagridId'] = cagridId
        AnatomicEntity.attrib['codeMeaning'] = fmaLabel # this info exists outside of atlas
        AnatomicEntity.attrib['codeValue'] = fmaid # changed from RID
        AnatomicEntity.attrib['codingSchemeDesignator'] = 'FMA' 
        AnatomicEntity.attrib['isPresent'] = 'True'
        AnatomicEntity.attrib['label'] = 'Pixel in %s' % fmaLabel #should be 'Pixel in '+ %s codeMeaning
        # create an AIM ImageAnnotation template - imageReferenceCollection sub-tree w/attributes
        imageReferenceCollection = etree.SubElement(ImageAnnotation, 'imageReferenceCollection')
        ImageReference = etree.SubElement(imageReferenceCollection,'ImageReference')
        ImageReference.attrib['cagridId'] = cagridId
        ImageReference.attrib[XSI_TYPE] = 'DICOMImageReference'
        # create an AIM ImageAnnotation template - imageStudy sub-tree w/attributes
        imageStudy = etree.SubElement(ImageReference, 'imageStudy')
        ImageStudy = etree.SubElement(imageStudy, 'ImageStudy')
        ImageStudy.attrib['cagridId'] = cagridId
        ImageStudy.attrib['instanceUID'] = ''
        ImageStudy.attrib['startDate'] = ''
        ImageStudy.attrib['startTime'] = ''
        # create an AIM ImageAnnotation template - imageSeries sub-tree w/attributes
        imageSeries = etree.SubElement(ImageStudy,'imageSeries')
        ImageSeries = etree.SubElement(imageSeries, 'ImageSeries')
        ImageSeries.attrib['cagridId'] = cagridId
        ImageSeries.attrib['instanceUID'] = ''
        # create an AIM ImageAnnotation template - imageCollection sub-tree w/attributes
        imageCollection = etree.SubElement(ImageSeries, 'imageCollection')
        Image = etree.SubElement(imageCollection, 'Image')
        Image.attrib['cagridId'] = cagridId
        Image.attrib['sopClassUID'] = ''
        Image.attrib['sopInstanceUID'] = ''
        # create an AIM ImageAnnotation template - geometricShapeCollection sub-tree w/attributes
        geometricShapeCollection = etree.SubElement(ImageAnnotation, 'geometricShapeCollection')
        GeometricShape = etree.SubElement(geometricShapeCollection, 'GeometricShape')
        GeometricShape.attrib['cagridId'] = cagridId
        GeometricShape.attrib['includeFlag'] = 'true'
        GeometricShape.attrib['shapeIdentifier'] = '0'
        GeometricShape.attrib[XSI_TYPE] = 'Point'
        # create an AIM ImageAnnotation template - spatialCoordinateCollection sub-tree w/attributes
        spatialCoordinateCollection = etree.SubElement(GeometricShape, 'spatialCoordinateCollection')
        SpatialCoordinate = etree.SubElement(spatialCoordinateCollection, 'SpatialCoordinate')
        SpatialCoordinate.attrib['cagridId'] = cagridId
        SpatialCoordinate.attrib['coordinateIndex'] = '0'
        SpatialCoordinate.attrib['imageReferenceUID'] = ''
        SpatialCoordinate.attrib['referencedFrameNumber'] = '1'
        SpatialCoordinate.attrib['x'] = x.__str__()
        SpatialCoordinate.attrib['y'] = y.__str__()
        SpatialCoordinate.attrib['z'] = z.__str__()
        SpatialCoordinate.attrib[XSI_TYPE] = 'ThreeDimensionSpatialCoordinate' #updated for 3D
        # create an AIM ImageAnnotation template - person sub-tree w/attributes
        person = etree.SubElement(ImageAnnotation,'person')
        Person = etree.SubElement(person, 'Person')
        Person.attrib['cagridId'] = cagridId
        Person.attrib['id'] = 'NIF'
        Person.attrib['name'] = ''
        return ImageAnnotation



#----FMA-GRAPH-------------------------------------------------------------------------------
class fmaGraph:
    def __init__(self):
       self.url = None
       self.fmaGraph = None
       self.fmaList = []
       self.aalList = []
       self.aalDict = {}

# URL to fmaGraph object
    def fmaURL(self,url):
        print "Checking URL..."
        try:
            tempGraph = Graph()
            fmaURL = urlopen(url)
            tempGraph.parse(fmaURL)
            self.fmaGraph = tempGraph
        except:
            print "Cannot open " + str(url) + ". Exiting."
            sys.exit()

# RDF file to fmaGraph object
    def fmaRDF(self,rdf):
        print "Checking RDF file..."
        if os.path.exists(rdf):
            tempGraph = Graph()
            tempGraph.parse(rdf)
            self.fmaGraph = tempGraph
        else:
            print "Cannot find rdf file " + str(rdf) + " . Exiting!"
            sys.exit()

# Parse fmaGraph Object    
    def fmaRead(self):
        print "Reading RDF..."
        self.fma2List()
        self.list2Nested()
        self.aalMap()
         
# Parse rdf into list
    def fma2List(self):

        # Each RDF entry should be organized with something like:
        # <rdf:Description rdf:nodeID="A0">
        # <result:fmaName xml:lang="en">Right supplemental motor cortex</result:fmaName>
        # <result:aalName rdf:datatype="http://www.w3.org/2001/XMLSchema#string">Supplemental_Motor_Area_RIGHT</result:aalName>
        # <result:FMAID rdf:datatype="http://www.w3.org/2001/XMLSchema#string">276680</result:FMAID>
        # <result:AALID rdf:datatype="http://www.w3.org/2001/XMLSchema#string">20</result:AALID>
        # </rdf:Description>

        fmaTemp = []
        for s,p,o in self.fmaGraph:
            fmaTemp.append((str(s),str(p),str(o)))
            self.fmaList = sorted(fmaTemp)          

            # sorting ensures correct indexing of aalIDs, fmaids, labels, etc in aalDict.
            # fmaList will have four subsequent entries for each label, in the format:
            # fmaList[0] --> ('IsOsrOtg10', 'http://sig.uw.edu/result#AALID', '19')
            # fmaList[1] --> ('IsOsrOtg10', 'http://sig.uw.edu/result#FMAID', '276677')
            # fmaList[2] --> ('IsOsrOtg10', 'http://sig.uw.edu/result#aalName', 'Supplemental_Motor_Area_LEFT')
            # fmaList[3] --> ('IsOsrOtg10', 'http://sig.uw.edu/result#fmaName', 'Left supplemental motor cortex')

# Organize fmaList into nested sets
    def list2Nested(self):
        aalList = []
        # Since each fma label has four entires (AALID, FMAID, aalname and fmaname) we add until the length is 4
        while len(self.fmaList) >= 4:
            # Clear tempList from the last group of four
            tempList = []
            for i in range(4):
                # We are popping off four members of the fmaList and putting them into tempList
                # The first pop is the AALID, 2nd: FMAID, 3rd: aalName, 4th: fmaName
                tempList.append(self.fmaList.pop())

                # Since we "pop" from the start of the list and add to the end, this reverses our data. So tempList looks like:
                
                # tempList[0] --> ('IsOsrOtg99', 'http://sig.uw.edu/result#fmaName', 'Left thalamus') fmaName
                # tempList[1] --> ('IsOsrOtg99', 'http://sig.uw.edu/result#aalName', 'Thalamus_LEFT') aalName
                # tempList[2] --> ('IsOsrOtg99', 'http://sig.uw.edu/result#FMAID', '258716')  FMAID
                # tempList[3] --> ('IsOsrOtg99', 'http://sig.uw.edu/result#AALID', '77')  AALID

            # We are then appending ((fmaName, aalName, FMAID, AALID))  or  (('Left thalamus','Thalamus_LEFT','258716','77'))
            self.aalList.append((tempList[0][2],tempList[1][2],tempList[2][2],tempList[3][2]))

    def aalMap(self):
        aalDict = {}
        for i in range(len(self.aalList)):
            # Pop each [('Left thalamus', 'Thalamus_LEFT', '258716', '77')] entry from the list
            aalTemp = self.aalList.pop()

            # Add to the AAL dictionary, indexing by the AALID, as a string
            # aalDict[AALID] = ((fmaName, aalName, FMAID))
            aalDict[aalTemp[3]] = ((aalTemp[0],aalTemp[1],aalTemp[2]))
        self.aalDict = aalDict


#-----------------------------------------------------------------------------------
def usage():
    print __doc__

# Check Directory for Output
def checkdir(userdir):
    # Make sure we don't end in a slash
    if userdir[-1] == "/":
        userdir = userdir[:-1]

    # If the directory doesn't exist, make it.
    if not os.path.exists(userdir):
    	print "Creating output directory..."
        os.makedirs(userdir)
        return userdir
    else:
	print "Directory already exists."
        return userdir	


#-----------------------------------------------------------------------------------
# SCRIPT STARTS RUNNING HERE  
#-----------------------------------------------------------------------------------
def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hi:o:nru", ["help","input=","out=","name=","rdf","url"])

    except getopt.GetoptError:
        usage()
        sys.exit(2)
    
    # First cycle through the arguments to collect user variables
    infile = None
    outname = None
    outfol = None
    url = 'http://xiphoid.biostr.washington.edu:8080/ValueSetService/ValueSet?qid=96'
    # This rdf maps all AALIDs to all FMAIDS - see http://xiphoid.biostr.washington.edu:8080/QueryManager/QueryManager.html#qid=96 for details 
    rdf = None

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()    
        if opt in ("--input","-i"):
            infile = arg
        if opt in ("-o", "--out"):
            outfol = arg
	if opt in ("--name","-n"):
            outname = arg
        if opt in ("--rdf","r"):
            rdf = arg
        if opt in ("--url","u"):
            url = arg

      
    aimFile = checkdir(outfol)     # check output directory

    # If no outname provided, use file name
    if not outname:
        outname,ext = os.path.splitext(infile)

    # Create rdf graph object that maps AAL IDs and FMAIDs
    FMA = fmaGraph()                              
    if rdf: FMA.fmaRDF(rdf)       
    else: FMA.fmaURL(url)                             
    
    # Extract aalIDs, fmaids, labels, etc and put into FMA objects aal dictionary
    FMA.fmaRead()    

    # Create AIMTemplate Object
    print "Creating AIM Template..."
    AIM = AIMTemplate(infile)

    # Get aal dictionary from rdf object - this is the dict to look up FMAID by aalID
    aalDict = FMA.aalDict    
    
    # Generate AIM Template with all aalIDs from atlas found with activation in input image
    aimTree = AIM.aimGen(aalDict)

    # Write to file
    aimTree.write('%s/AIM-%s.xml' % (aimFile, str(outname))) # write out the aim files
    

if __name__ == "__main__":
    main(sys.argv[1:])
