#!/bin/sh

# This bash script does a group ICA for use with a multisession temporal concatenation.

# USER SPECIFIED VARIABLES
TR=2.5		# TR
BGTHS=10	# Brain background threshold
COMPS=30        # Number of components for GICA

# VARIABLES SET AT RUN TIME
OUTPUT=$1       
PYEXEC=$2       # Python executable to use to run filtering, will be same used to run ica+.py
FILTERSCRIPT=$2 # Path to melodic_hp.py, should be in same directory as MRTools.py
SUBJECTS=( $3 ) # List of all subject ICA directories

# make sure output directory was made by submission script.
if [ ! -d "$OUTPUT" ]; then
    mkdir -p $OUTPUT
fi

# Make sure all paths exist
for ICADIR in ${SUBJECTS[@]}; do
  if [ ! -d $ICADIR ]; then
    echo "Cannot find " $ICADIR ". Exiting"
    exit 32
  fi
done

mask_list=();
bg_list=();

# Higher-level MELODIC
for ICADIR in ${SUBJECTS[@]}; do

	# Go to the ICADIR
	cd $ICADIR

	# check n timepoints
	total_volumes=`fslnvols filtered_func_data`
	echo "Counting total volumes are " $total_volumes
	    
	# check for reg/example_func2standard.mat to make sure that registration has happened
	if [ ! -f "$ICADIR/reg/example_func2standard.mat" ]; then
	    echo "Cannot find example_func2standard.mat for " $ICADIR
	fi

	# add the mask to the list if the subject has it
	if [ -f "$ICADIR/reg_standard/mask.nii.gz" ]; then	
		masklength=`echo ${#mask_list[*]}`
		mask_list[$masklength]=$ICADIR/reg_standard/mask.nii.gz
	fi

	# add the bg_image to the list if the subject has it
	if [ -f "$ICADIR/reg_standard/bg_image.nii.gz" ]; then	
		bglength=`echo ${#bg_list[*]}`
		bg_list[$bglength]=$ICADIR/reg_standard/bg_image.nii.gz
	fi

done


# GROUP MASK PREPARATION

# GO TO THE GROUP OUTPUT DIRECTORY
cd $OUTPUT

# Merge all single subject bg_images to make a "master" image
echo "Creating group background image..."
fslmerge -t bg_image $bg_list

# Intensity normalize the bg_image so mean is 1000, and take mean across time 
fslmaths bg_image -inm 1000 -Tmean bg_image -odt float

# Merge all the subject mask images to make a master mask -t specifies to merge across time, as opposed to any dimension (x,y,z)
echo "Creating group mask..."
fslmerge -t mask $mask_list

# Make the inputreg directory for...
mkdir -p inputreg
cd inputreg

# Multiply the mask by the number of subjects, so the value at each voxel is equivalent to the number of subjects.  Then take the mean across time so we are sure to have a single image with each voxel equal to the number of subjects.  The output is masksum  
fslmaths ../mask -mul ${#mask_list[*]} -Tmean masksum -odt short

# Threshold it at the number of subjects, then add itself to itself, so you have each voxel being double the number of subjects. This seems repetitive, but there must be some rationale for it?
fslmaths masksum -thr ${#mask_list[*]} -add masksum masksum

# Calculate double the number of subjects
let dubsubs=${#mask_list[*]}*2

# Using the bg_image as a background image, overlay masksum (range between 0.9 and double the number of subjects) with a solid color type (0) and floating point output type (0), and use -a for automatic estimation of background display range.
# The output is masksum_overlay
overlay 0 0 -c ../bg_image -a masksum 0.9 $dubsubs masksum_overlay

# User slicer to print an image of the masksum_overlay, -S 2 means sample every 2 slices, and output width is 750, filename is masksum_overlay.png
slicer masksum_overlay -S 2 750 masksum_overlay.png

# multiply masksum by 0 to get maskunique - where we will add all voxels for each subject that are missing.
fslmaths masksum -mul 0 maskunique



# CREATION OF UNIQUE MASK - with voxel value identifying subject number that is missing coverage for a particular area

# Cycle again through the list of subject masks, and add each one to the mask image
for (( mstart = 0; mstart < ${#mask_list[*]}; mstart++ )); do

	# Put the mask name into a variable, and add one to the count variable
	# to get the subject number, since we start iterating with 0
	MASK=${mask_list[$mstart]}
	let subnum=$mstart+1	

	# Multiply the mask by -1 and add 1 to get it's inverse
	# Values that were 1 are now zero, and values that were zero are 1
	# Multiply this area by the subject number, so the intensity of the region
	# equals that, and we can use intensities to ID subjects.
	# Then add to maskunique, overriding the old one.
	fslmaths $MASK -mul -1 -add 1 -mul $subnum -add maskunique maskunique

done

let subs_minus=${#mask_list[*]}-1

fslmaths masksum -thr $subs_minus -uthr $subs_minus -bin -mul maskunique maskunique

# do the same process with overlay and slicer to get the png image
overlay 0 0 ../bg_image -a maskunique 0.9 $dubsubs maskunique_overlay
slicer maskunique_overlay -S 2 750 maskunique_overlay.png

# Return to output directory
cd $OUTPUT
fslmaths mask -Tmin -bin mask -odt char



# RUN GROUP MELODIC

# copy file for GUI use
cp ${FSLDIR}/etc/luts/ramp.gif .ramp.gif

# Prepare list of individual subject ica directories
inputsubs=${SUBJECTS[0]}"/reg_standard/filtered_func_data.nii.gz"

for (( sstart = 1; sstart < ${#SUBJECTS[*]}; sstart++ )); do
    inputsubs=$inputsubs","${SUBJECTS[$sstart]}"/reg_standard/filtered_func_data.nii.gz"
done

echo "inputsubs are " $inputsubs

melodic -i $inputsubs -o $OUTPUT/groupmelodic.ica -v --nobet --bgthreshold=$BGTHS --tr=$TR --report --bgimage=$OUTPUT/bg_image -d $COMPS --Ostats -a concat

# This is the standard melodic command - the above is customized
# melodic -i .filelist -o groupmelodic.ica -v --nobet --bgthreshold=10 --tr=2.5 --report --guireport=../../report.html --bgimage=bg_image -d 0 --mmthresh=0.5 --Ostats -a concat

# AUTOMATICALLY FILTER COMPONENTS
# Produce a list of IC components and (corresponding potential DR future result images) that pass a highpass filter
# This list will be used by pyMatch to match ONLY good components to a template of interest!
echo "Performing highpass filtering of output networks..."

mkdir -p $OUTPUT/filter

echo "Command is bsub -J gica_filter -o $OUTPUT/log/filter.out -e $OUTPUT/log/filter.err $PYEXEC $FILTERSCRIPT -o $OUTPUT/filter --name=gica --ts=$OUTPUT/groupmelodic.ica/report --gica=$OUTPUT/groupmelodic.ica/stats "

# Submit group filter script to produce files with lists of good ICA and (potential future) dual regression results
bsub -J "gica_filter" -o $OUTPUT/log/filter.out -e $OUTPUT/log/filter.err $PYEXEC $FILTERSCRIPT -o $OUTPUT/filter --name=gica --ts=$OUTPUT/groupmelodic.ica/report --gica=$OUTPUT/groupmelodic.ica/stats

echo "Results gica_IC-hpfilter-good.txt and gica_DR-hpfilter-good.txt will be in " $OUTPUT"/filter, for use with Match functionality of ica+ package."
