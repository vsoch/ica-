#!/bin/sh

# Some sort of component filtering algorithm needs to go here.  For now, 
# we will do them all, or just one.

# This section will extract max voxel coordinates of all signicant clusters
# We will maintain the numbering schema of the components for now, but in the
# future it might make sense to link this with an ontology to provide a more
# meaningful name for the functional network.  We will then write the coordinates
# to text file in the format:
# num x y z 
#
# We want one text file per dual regression significant result.  Each
# text file will be used to extract coordinates, and query subject raw data.
#
# We can extract the timeseries for each coordinate for each subject,
# and we want one text file per max voxel, OR one text file per component,
# containing multiple max voxels

# Each group of max voxel timeseries will be used as the group of predictors
# (for one subject) and the subject's ADHD diagnosis used as the decision for
# the training set of the data

# The goal at the end of the day is to predict ADHD diagnosis based on raw
# resting data (although it is motion corrected, registered to MNI)

# Submit on command line: melodic_ts.sh /path/to/ica/experiment name.gica dual_regression_name

TSOUT="negts"
EXPERIMENT=$1	  # This is the top level, experimental directory, no slash on the end
GICA=$2           # The name of the gica run (NAME.gica)
DRRUN=$3          # The name of the dual regression run
COMPS=( {0..24} ) # List of components to extract timeseries for, number of thresh_zstat#.nii.gz under stats
                  # For a range, do ( {1..30) )   for a custom list do ( {1 3 5} ) and for just one, do ( 4 )

# Make sure that we have group ICA folder
if [ ! -e "$EXPERIMENT/gica/$GICA/groupmelodic.ica" ]; then
    echo "Cannot find groupmelodic.ica directory " $GICA " under " $EXPERIMENT"/gica. Exiting."
    exit 32 
else
    ICADIR=$EXPERIMENT/gica/$GICA
fi

# Make sure that we have dr run specified...
if [ ! -d "$EXPERIMENT/dr/$DRRUN/result" ]; then
    echo "Cannot find dual regression directory " $DRRUN " under " $EXPERIMENT"/dr Exiting."
    exit 32
else
    DRRUN=$EXPERIMENT/dr/$DRRUN/result
fi
 

# Read in filelist from group run, check that we have all single subject input
FILELIST=`cat $ICADIR/.filelist`

echo "Filelist is " $FILELIST

INPUTS=( $FILELIST )
echo "Inputs are " $INPUTS
for i in ${INPUTS[@]} ; do
    if [ ! -e "$i.nii.gz" ]; then
	    echo "Cannot find single subject input " $i " Exiting."
	    exit 32 
    fi
done

# Create output ts folder
if [ ! -e "$EXPERIMENT/ts" ]; then
    echo "Creating ts (timeseries) directory under " $EXPERIMENT "..."
    mkdir -p $EXPERIMENT/ts
fi

if [ ! -e "$EXPERIMENT/ts/$TSOUT" ]; then
    echo "Creating " $TSOUT " under ts..."
    mkdir -p $EXPERIMENT/ts/$TSOUT
else
    echo $TSOUT " already exists! Choose a different name.  Exiting."
    exit 32
fi

OUTPUT=$EXPERIMENT/ts/$TSOUT

# Create output log folder
if [ ! -e "$EXPERIMENT/ts/$TSOUT/log" ]; then
    echo "Creating ts/$TSOUT/log output folder..."
    mkdir -p $OUTPUT/log
fi

# FILTER OUT BAD ONES SOMEHOW...
# ALSO NEED TO CHECK FOR EMPTY IMAGES and REMOVE
# FOR NOW WE WILL USE COMPONENTS SPECIFIED BY USER (UNDER COMPS)

# Get list of output dr images, add to list if specified
DRRES=''
COUNT=0
for temp in `ls $DRRUN/dr_stage3_ic*_tfce_corrp_tstat1.nii.gz`; do
    for n in "${COMPS[@]}"; do
        if [ "$n" -eq "$COUNT" ] ; then
	    echo "Adding " $temp " to extract timeseries from..."
	    DRRES[$COUNT]=$temp;
        fi
    done 
    let COUNT=$COUNT+1; 
done

# GET MAX VOXEL COORDINATES FOR EACH

# Go to output directory
cd $OUTPUT
ZVALSORT=''

COUNT=0	
for IMAGE in ${DRRES[@]}; do
    echo "Getting local maximum voxel coordinates for component " ${COMPS[$COUNT]} 
    echo "Image is: " $IMAGE"..."

    # This is the text file to save the local maxima list, with a header plus (cluster_index) zval xvox yvox zvox
    OLFILE=lomax${COMPS[$COUNT]}.txt
    echo "OLFILE is " $OLFILE

    # We run cluster, with the dual regression corrp, thresholded image as input
    cluster --in=$IMAGE.nii.gz --thresh=0 --olmax=$OLFILE --no_table
    
    # We read in the output table from file
    lines=( $(cat $OLFILE) )
    let uprange=${#lines[@]}-5
    echo "Uprange is " $uprange

    # Iterate through subjects, and extract timeseries for each one, print to file
    for SUB in ${INPUTS[@]} ; do
	
    # WE SHOULD KEEP TRACK OF ZVAL - TO PRESENT LIST OF MOST SIGNIFICANT DIFFERENCES TO USER, IN ORDER?
	
	for (( cnt=6; cnt <= $uprange; cnt+=5 )); do
	    CLUSINDEX=${lines[$cnt]}
	    ZVAL=${lines[$cnt+1]}
	    XVOX=${lines[$cnt+2]}
	    YVOX=${lines[$cnt+3]}
	    ZVOX=${lines[$cnt+4]}
	    echo "Extracting for XYZ: " $XVOX " " $YVOX " " $ZVOX

            # Go through the individual subject ICA directories, extract each set for each subject, into one long row, print to file
	    # one row per subject!		
	    TS=`fslmeants -i $SUB.nii.gz -c $XVOX $YVOX $ZVOX --transpose`
            echo -ne $TS" " >>  ${COMPS[$COUNT]}_ts.txt

        done
            echo "Finished with subject " $SUB
	    echo "" >> ${COMPS[$COUNT]}_ts.txt
    done
    
    # PREPARE LIST ZVALUES - GREATEST TO LEAST
    # Here we look through the list of the top Z value for each component and print a file with these ordered from largest --> smallest

    # Cycle again through the lines and grab the ZVAL and CLUSTINDEX
    for (( cnt=6; cnt <=$uprange; cnt+=5 )); do
        CLUSINDEX=${lines[$cnt]}
        ZVAL=${lines[$cnt+1]}
        echo "Adding ZVAL " $ZVAL " to list for later sorting"
        ZVALSORT=$ZVALSORT" "$ZVAL-${COMPS[$COUNT]}-$CLUSINDEX
    done
    let COUNT=$COUNT+1
done

# Now print the string of Zvalues to file
echo "zval-component-cluster#" >> zval_sort.txt
for zv in `echo $ZVALSORT`; do
    echo $zv
done | sort -r >> zval_sort.txt

# Here is an example of the training data.  There is one row per subject with a list of the timeseries
# for each max voxel.  To start, this will be done for the component that is clearly the attention
# network.  In the future, I'd like to figure out a way to first filter out bad components, and then
# choose out the "most meaningful" based on the disorder in question (or perhaps the most deviant from
# the norm?) This method would be used to come up with an algorithm that can make an ADHD diagnosis based
# on the extracted raw data from these timepoints during a resting bold scan.

# QUESTION: As is shown below, each point in a timeseries is used as a parameter, and there are 
# (number of TRs) X (number of max voxels) predictors.  Is it possible to put an entire timeseries
# in as a predictor, and then add other variables such as the index and size of the cluster?  

# These would also need to come with labels for the coordinates, so you don't forget where the data
# has come from, and can reproduce the values for testing data.

# TRAINING DATA EXAMPLE
#	  ADHD_rX	TIMESERIES (start with for just one component - the attention one?)	
# Subject :   		[_TP, _TP, ,TP..._TP, _TP, _TP...n_TP, n_TP, n_TP]
# Subject :   0		[_TP, _TP, ,TP..._TP, _TP, _TP...n_TP, n_TP, n_TP]
# Subject :   		[_TP, _TP, ,TP..._TP, _TP, _TP...n_TP, n_TP, n_TP]
