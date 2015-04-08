#!/bin/sh

# This bash script does an Individual ICA for use with a multisession temporal concatenation. It takes the data through the application of the registration, to be ready for group ICA.

# USER DEFINED VARIABLES
TR=2.5            # TR
HPFC=150          # highpass frequency filter cutoff
SKERN=6.0         # Smoothing kernel (FWSE)
BBTHRESH=10       # Brain background threshold
MATLABRUNTIME=/home/vsochat/software/matlab/Matlab64bitRuntime_v714/v714
LFILT=0.008       # Lower filter threshold
UFILT=0.1         # Upper filter threshold
                  # MATLABRUNTIME is path to Matlab Runtime

# VARIABLES DEFINED AT RUNTIME
OUTPUT=$1         # Output folder (not yet created)
FUNCDATA=$2       # Full path to raw rest data
ANATDATA=$3       # Full path to raw anatomical data
SCRIPTDIR=$PWD

# SETUP

# make sure output directory was made by submission script.
if [ ! -d "$OUTPUT" ]; then
    mkdir -p $OUTPUT
fi

# Check once more for functional data...
if [ ! -f "$FUNCDATA" ]; then
   echo "Cannot find functional data file " $FUNCDATA ". Exiting."
   exit 32
fi

# Check once more for anatomical data...
if [ ! -f "$ANATDATA" ]; then
   echo "Cannot find anatomical data file " $ANATDATA ". Exiting."
   exit 32
fi

# go to output directory
cd $OUTPUT 

# Use fslmaths to copy the raw anatomical data into the outdir
fslmaths $ANATDATA mprage



# ANATOMICAL PREPROCESSING

# Perform BET (brain extraction) on the raw anatomical data
bet mprage mprage_bet -S -f .225

# Delete the raw file
rm mprage.nii.gz



# FUNCTIONAL PREPROCESSING

# take the raw functional data and use fslmaths to convert it to float
fslmaths $FUNCDATA prefiltered_func_data -odt float

# Get the number of slices
func_slices=`fslval prefiltered_func_data dim4`
middle_slice=`echo "(($func_slices/2))" | bc`

# Create "example_func" based on the middle slice, one 3D image
fslroi prefiltered_func_data example_func $middle_slice 1

# run MCFLIRT motion correction on the prefiltered_functional_data
mcflirt -in prefiltered_func_data -out prefiltered_func_data_mcf -mats -plots -refvol $middle_slice -rmsrel -rmsabs

# make mc directory and copy new mcflirt files into it
mkdir -pv mc
mv -f prefiltered_func_data_mcf.mat prefiltered_func_data_mcf.par prefiltered_func_data_mcf_abs.rms prefiltered_func_data_mcf_abs_mean.rms prefiltered_func_data_mcf_rel.rms prefiltered_func_data_mcf_rel_mean.rms mc

# Plot the Mcflirt rotational data
fsl_tsplot -i mc/prefiltered_func_data_mcf.par -t 'MCFLIRT estimated rotations (radians)' -u 1 --start=1 --finish=3 -a x,y,z -w 640 -h 144 -o mc/rot.png

# Plot the Mcflirt translational data
fsl_tsplot -i mc/prefiltered_func_data_mcf.par -t 'MCFLIRT estimated translations (mm)' -u 1 --start=4 --finish=6 -a x,y,z -w 640 -h 144 -o mc/trans.png 

# Plot the Mcflirt mean displacement data
fsl_tsplot -i mc/prefiltered_func_data_mcf_abs.rms,mc/prefiltered_func_data_mcf_rel.rms -t 'MCFLIRT estimated mean displacement (mm)' -u 1 -w 640 -h 144 -a absolute,relative -o mc/disp.png 

# take the motion corrected functional data and calculate the mean across time - mean_func
fslmaths prefiltered_func_data_mcf -Tmean mean_func

# Perform bet2 on the mean_func data, use a threshold of .3 (.225 used for anatomical)
bet2 mean_func mask -f 0.3 -n -m

# This command moves the first input (mask_mask) to the new image (mask)
immv mask_mask mask

# Mask the motion corrected functional data with the mask to create the masked (bet) motion corrected functional data
fslmaths prefiltered_func_data_mcf -mas mask prefiltered_func_data_bet

# Calculate the difference between the 98th and 2nd percentile (the region between the tails) and use that range as a threshold minimum on the prefiltered, motion corrected, masked functional data - so we are eliminating 
# intensities outside that are below 2nd percentile, and above 98th percentile.

# Use fslstats to output the 2nd and 98th percentile
lowerp=`fslstats prefiltered_func_data_bet -p 2`
upperp=`fslstats prefiltered_func_data_bet -p 98`
# The brain/background threshold (to distinguish between brain and background is 10% - so we divide by 10)
# Anything above this value, then, is activation between the 2nd and 98th percentile that is likely to be 
# brain activation and not noise or something else!
thresholdp=`echo "scale=6; ($upperp-$lowerp)/$BBTHRESH" | bc`

# Use fslmaths to threshold the brain extracted data based on the highpass filter above
# use "mask" as a binary mask, and Tmin to specify we want the minimum across time
fslmaths prefiltered_func_data_bet -thr $thresholdp -Tmin -bin mask -odt char

# Take the motion corrected functional data, and using "mask" as a mask (the -k option)
# output the 50th percentile (the mean?)
# We will need this later to calculate the intensity scaling factor
meanintensity=`fslstats prefiltered_func_data_mcf -k mask -p 50`

# IM NOT SURE WHAT WE DO WITH THIS?
# difF is a spatial filtering option that specifies maximum filtering of all voxels
# I don't completely understand why we would filter one image with itself...
fslmaths mask -dilF mask

# We are now masking the motion corrected functional data with the mask to produce
# functional data that is motion corrected and thresholded based on the highpass filter
fslmaths prefiltered_func_data_mcf -mas mask prefiltered_func_data_thresh

# We now take this functional data that is motion corrected, high pass filtered, and
# create a "mean_func" image that is the mean across time (Tmean)
fslmaths prefiltered_func_data_thresh -Tmean mean_func

# To run susan, FSLs tool for noise reduction, we need a brightness threshold.  Here is how to calculate:
# After thresholding, the values in the image are between $upperp-lowerp and $thresholdp
# If we set the expected noise level to .66, then anything below (($upperp-$lowerp)-$thresholdp)/0.66 should be noise.
# This is saying that we want the brightness threshold to be 66% of the median value.
# Note that the FSL "standard" is 75% (.75)
# This is the value that we use for bt, the "brightness threshold" in susan
uppert=`echo "scale=6; (($upperp-$lowerp))" | bc`
difft=`echo "scale=8; (($uppert-$thresholdp))" | bc`

# We also need to calculate the spatial size based on the smoothing.
# FWHM = 2.355*spatial size. So if desired FWHM = 6mm, spatial size = 2.54...
ssize=`echo "scale=11; (($SKERN/2.355))" | bc`
# susan uses nonlinear filtering to reduce noise 
# by only averaging a voxel with local voxels which have similar intensity
susan prefiltered_func_data_thresh $difft $ssize 3 1 1 mean_func $difft prefiltered_func_data_smooth

# 3 means 3D smoothing
# 1 says to use a local median filter 
# 1 says that we determine the smoothing area from 1 secondary image, "mean_func" and then we use the same brightness threshold for the secondary image.
# prefiltered_func_data_smooth is the output image

# Now we mask the smoothed functional data with the mask image, and overwrite the smoothed image.
fslmaths prefiltered_func_data_smooth -mas mask prefiltered_func_data_smooth

# We now need to calculate the intensity scaling factor applied to the whole  4D dataset so that it's mean is 10000
inscalefactor=`echo "scale=6; ((10000/$meanintensity))" | bc`

# We now multiply the smoothed data by the intensity scaling factor to get
# the intensity normalized data
fslmaths prefiltered_func_data_smooth -mul $inscalefactor prefiltered_func_data_intnorm

#####FSL BANDPASS METHOD  NOT IN USE##################################################################
# Now bandpass temporal filter the intensity normalized data.
# We need to calculate $hp_sigma_vol before continuing:
# $HPFC is the highpass filter cutoff, and we set the second 
# argument to -1, since we don't want lowpass filtering
# hp_sigma_sec=`echo "scale=6; (($HPFC/2.0))" | bc`
# hp_sigma_vol=`echo "scale=6; (($hp_sigma_sec/$TR))" | bc`

# fslmaths prefiltered_func_data_intnorm -bptf $hp_sigma_vol -1 prefiltered_func_data_tempfilt
######################################################################################################

# Bandpass filter the data using matlab executable, bandpass

# First convert the .nii.gz to .nii
fslchfiletype NIFTI prefiltered_func_data_intnorm.nii.gz

# Launch matlab script (<scriptname> <matlabruntime> <input .nii> <name for output .nii (no extension)> <TR> <lower filter> <upper filter>
$SCRIPTDIR/run_Bandpass.sh $MATLABRUNTIME prefiltered_func_data_intnorm.nii prefiltered_func_data_tempfilt $TR $LFILT $UFILT

# The output will be .nii, which we need to change back to .nii.gz
fslchfiletype NIFTI_GZ prefiltered_func_data_tempfilt.nii

# Use fslmaths to copy the file with a new name (filtered_func_data)
fslmaths prefiltered_func_data_tempfilt filtered_func_data

# Read the header info from filtered_func_data, and write it into tmpHeader
# and change dt from whatever it is to 2.0
fslhd -x filtered_func_data | sed 's/  dt = .*/  dt = '2.0'/g' > tmpHeader

# Create the header and apply to filtered_func_data
fslcreatehd tmpHeader filtered_func_data

# delete the temporary header
rm tmpHeader

# Calculate the mean across time for the filtered_func_data, replace mean_func
fslmaths filtered_func_data -Tmean mean_func

# Remove all prefiltered_func_data files - we don't need them!
rm -rf prefiltered_func_data*

# Make a directory for registration
mkdir -p reg



# REGISTRATION TO STANDARD AND HIGHRES

# use fslmaths to copy the MNI152_T1 to the current directory
# the brain extracted anatomical is already here (mprage_bet)
fslmaths ${FSLDIR}/data/standard/MNI152_T1_2mm_brain standard

# Use flirt to register the example_func data to the highres anatomical
flirt -ref mprage_bet -in example_func -out example_func2mprage_bet -omat example_func2mprage_bet.mat -cost corratio -dof 7 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp trilinear  
 
convert_xfm -inverse -omat mprage_bet2example_func.mat example_func2mprage_bet.mat

# Use slicer to print pictures
slicer example_func2mprage_bet mprage_bet -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png 

pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png example_func2mprage_bet1.png 

slicer mprage_bet example_func2mprage_bet -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png

pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png example_func2mprage_bet2.png

pngappend example_func2mprage_bet1.png - example_func2mprage_bet2.png example_func2mprage_bet.png

rm -f sl?.png

# Now register the highres anatomical to the standard space
flirt -ref standard -in mprage_bet -out mprage_bet2standard -omat mprage_bet2standard.mat -cost corratio -dof 12 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp trilinear  

convert_xfm -inverse -omat standard2mprage_bet.mat mprage_bet2standard.mat

slicer mprage_bet2standard standard -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png

pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png mprage_bet2standard1.png 

slicer standard mprage_bet2standard -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png 

pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png mprage_bet2standard2.png

pngappend mprage_bet2standard1.png - mprage_bet2standard2.png mprage_bet2standard.png

rm -f sl?.png

# Now combine the bet2standard and example_func2mprage_bet to get example_func2standard
convert_xfm -omat example_func2standard.mat -concat mprage_bet2standard.mat example_func2mprage_bet.mat

flirt -ref standard -in example_func -out example_func2standard -applyxfm -init example_func2standard.mat -interp trilinear

convert_xfm -inverse -omat standard2example_func.mat example_func2standard.mat

slicer example_func2standard standard -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png 

pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png example_func2standard1.png

slicer standard example_func2standard -s 2 -x 0.35 sla.png -x 0.45 slb.png -x 0.55 slc.png -x 0.65 sld.png -y 0.35 sle.png -y 0.45 slf.png -y 0.55 slg.png -y 0.65 slh.png -z 0.35 sli.png -z 0.45 slj.png -z 0.55 slk.png -z 0.65 sll.png 

pngappend sla.png + slb.png + slc.png + sld.png + sle.png + slf.png + slg.png + slh.png + sli.png + slj.png + slk.png + sll.png example_func2standard2.png 

pngappend example_func2standard1.png - example_func2standard2.png example_func2standard.png

rm -f sl?.png

# Make reg_standard folder
mkdir -p reg_standard

# Copy some files into the reg folder
cp -f example_func.nii.gz reg

# Delete intermediate image files
rm -f example_func2standard1.png example_func2standard2.png example_func2mprage_bet1.png example_func2mprage_bet2.png mprage_bet2standard1.png mprage_bet2standard2.png

# Move registration files into reg folder
mv -f standard.nii.gz example_func2* mprage_bet2* standard2* reg 


# FEATREGAPPPLY

# This first section accomplishes the same functionality of running "featregapply $ICADIR" however it doesn't require the .fsf file
	

# Create a temporary file that will be used for the output two steps down!
TMPNAME=`tmpnam frgrot -n`

# Take the standard template and apply isotropic resampling resolution (4mm)
flirt -ref reg/standard -in reg/standard -out reg_standard/standard -applyisoxfm 4 

# Now register the filtered_func_data to the standard space (with the highres anatomical as an intermediate) using the example_func2standard mat
flirt -ref reg_standard/standard -in filtered_func_data -out $TMPNAME -applyxfm -init reg/example_func2standard.mat -interp trilinear -datatype float

# Copy this newly created file (that is finally the functional data registered to standard space based on anatomical as intermediate) into reg_standard folder, and get rid of intermediate file.
immv $TMPNAME reg_standard/filtered_func_data
rm $TMPNAME
	
# Take the mean of the filtered_func_data across time, and then binarize it to make the reg_standard mask (1 timepoint as opposed to many)
fslmaths reg_standard/filtered_func_data -Tstd -bin reg_standard/mask -odt char

# Copy the standard to use as the background image
imcp reg_standard/standard reg_standard/bg_image

# This would be the command for single subject ica - but we are doing gica
# fsl:exec "${FSLDIR}/bin/melodic -i filtered_func_data -o filtered_func_data.ica -# v --nobet --bgthreshold=1 --tr=$fmri(tr) -d 0 --mmthresh=\"0.5\" --report --# guireport=../../report.html "
