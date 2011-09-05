#!/bin/sh

# Some sort of component filtering algorithm needs to go here.  For now, 
# we will do them all, and only use output for "good" components.

# This section will extract max voxel coordinates of all signicant clusters
# We will maintain the numbering schema of the components for now, but in the
# future it might make sense to link this with an ontology to provide a more
# meaningful name for the functional network.  We will then write the coordinates
# to text file in the format (coordinate space):
# index_num Z_val x y z 
#
# We want one text file per dual regression significant result.  Each
# text file will be used to extract timeseries and means from subject raw data.
#
# The output will be either extracted timeseries (melodic_ts) or the mean
# across time (this script).  There will be one output file per subject
# in the format #_mean.txt with a list of mean activation values for each
# voxel coordinates, and one subject per row.  These values can be
# read into matlab to calculate correlations, etc.

# Application to machine learning - each group of mean values will 
# be used as the group of predictors (for one subject) and the 
# subject's ADHD diagnosis used as the decision for the training set

# The goal at the end of the day is to predict ADHD diagnosis based on raw
# resting data (preprocessed)

# Submit on command line: melodic_ts.sh /path/to/ica/experiment name.gica dual_regression_name

TSOUT="tsneg"   # Name for output directory under $EXPERIMENT/ts
EXPERIMENT=$1	# This is the top level, experimental directory, no slash on the end
GICA=$2         # The name of the gica run (NAME.gica)
DRRUN=$3        # The name of the dual regression run
COMPS=( {0..24} ) # List of dual regression results to extract timeseries for, number of dr_stage3_ic00#_tfce_corrp*.nii.gz under dr.result
                # For a range, do ( {0..29) )   for a custom list do ( {1 3 5} ) and for just one, do ( 4 )
		# Keep in mind that dual regression results start counting at 0

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
    echo "Creating $TSOUT under ts..."
    mkdir -p $EXPERIMENT/ts/$TSOUT
else
    echo $TSOUT " already exists under ts.  Choose different name.  Exiting"
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
    echo "Getting local maximum voxel coordinates for dual regression result " ${COMPS[$COUNT]} 
    echo "Image is: " $IMAGE"..."

    # This is the text file to save the local maxima list, with a header plus (cluster_index) zval xvox yvox zvox
    OLFILE=lomaxdr${COMPS[$COUNT]}.txt
    echo "OLFILE is " $OLFILE

    # We run cluster, with the dual regression corrp, thresholded image as input
    cluster --in=$IMAGE.nii.gz --thresh=0 --olmax=$OLFILE --no_table
    
    # We read in the output table from file
    lines=( $(cat $OLFILE) )
    let uprange=${#lines[@]}-5
    echo "Uprange is " $uprange

    # Print voxel coordinates as headers to file
    for (( cnt=6; cnt <= $uprange; cnt+=5 )); do
	    CLUSINDEX=${lines[$cnt]}
	    ZVAL=${lines[$cnt+1]}
	    XVOX=${lines[$cnt+2]}
	    YVOX=${lines[$cnt+3]}
	    ZVOX=${lines[$cnt+4]}
   	    echo -ne $XVOX"_"$YVOX"_"$ZVOX" " >>  dr${COMPS[$COUNT]}_mean.txt
    done
            echo "" >> dr${COMPS[$COUNT]}_mean.txt

    # Iterate through subjects, and extract timeseries for each one, print to file
    for SUB in ${INPUTS[@]} ; do
	
    # WE KEEP TRACK OF ZVAL - TO PRESENT LIST OF MOST SIGNIFICANT DIFFERENCES TO USER, IN ORDER
	
	for (( cnt=6; cnt <= $uprange; cnt+=5 )); do
	    CLUSINDEX=${lines[$cnt]}
	    ZVAL=${lines[$cnt+1]}
	    XVOX=${lines[$cnt+2]}
	    YVOX=${lines[$cnt+3]}
	    ZVOX=${lines[$cnt+4]}
	    echo "Extracting for XYZ: " $XVOX " " $YVOX " " $ZVOX
	
	    # Extract the timeseries and save to variable TS
	    TS=`fslmeants -i $SUB.nii.gz -c $XVOX $YVOX $ZVOX --transpose`

	    # Calculate the mean
	    meanTS=`echo $TS | awk '{ tot=0; for (i=1; i<=NF; i++) tot += $i; print tot/NF; }'`

            # Print to file
	    # one row per subject!		
	    echo -ne $meanTS" " >>  dr${COMPS[$COUNT]}_mean.txt

        done
            echo "Finished with subject " $SUB
	    echo "" >> dr${COMPS[$COUNT]}_mean.txt
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
echo "zval-component-cluster#" >> zval.txt
for zv in `echo $ZVALSORT`; do
    echo $zv
done | sort -r >> zval.txt

# Here is an example of the training data.  There is one row per subject with a list of the mean
# for each max voxel.  To start, this will be done for the component that is clearly the attention
# network.  In the future, I'd like to figure out a way to first filter out bad components, and then
# choose out the "most meaningful" based on the disorder in question (or perhaps the most deviant from
# the norm?) This method would be used to come up with an algorithm that can make an ADHD diagnosis based
# on the extracted raw data from these timepoints during a resting bold scan.

# These would also need to come with labels for the coordinates, so you don't forget where the data
# has come from, and can reproduce the values for testing data.

# TRAINING DATA EXAMPLE
#	  ADHD_rX	MEAN VALUES / LOCAL MAX	
# Subject :   		[ MV1 MV2 MV3 MV4 MV5 MV6 ]
# Subject :   		[ MV1 MV2 MV3 MV4 MV5 MV6 ]
# Subject :   		[ MV1 MV2 MV3 MV4 MV5 MV6 ]
