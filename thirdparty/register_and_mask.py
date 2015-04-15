#!/usr/bin/env python
""" 
Modified from openFMRI pipeline, Poldracklab
https://github.com/poldrack/openfmri
register image to whole-head MNI template and apply mask,
then register again to brain-only template
"""


>>> os.system("fslswapdim t1 RL PA IS t12")
0
>>> imgfile="t12.nii.gz"
>>> os.system("fslswapdim t1 LR AP SI t13")
Cannot perform requested swap (NEUROLOGICAL/RADIOLOGICAL storage inverted)
Try the following command instead:
fslswapdim t1 RL AP SI t13
256
>>> os.system("fslswapdim t1 RL AP SI t13")


import os,sys

ANTSPATH=os.environ['ANTSPATH']
FSLDIR=os.environ['FSLDIR']
test=False

def main():
    # first run N4 bias correction
    imgfile=sys.argv[1]
    subdir,infile=os.path.split(os.path.abspath(imgfile))
    if not os.path.exists(subdir):
        print '%s does not exist!'%subdir
        sys.exit(0)
      
    image_prefix = image_prefix = infile.split(".")[0]
    cmd='%s/N4BiasFieldCorrection -i %s/%s -d 3 -o %s/%s_bfc.nii.gz -c [50,0.0001]'%(ANTSPATH,subdir,infile,subdir,image_prefix)

    print cmd
    if not test:
        os.system(cmd)

 
    # then align bias-corrected whole-head image to template

    PARAMS="-r Gauss[3,0] -t SyN[0.25] -i 30x90x20 --use-Histogram-Matching --number-of-affine-iterations 10000x10000x10000x10000x10000 --rigid-affine false"

    template_wholehead='%s/data/standard/MNI152_T1_2mm.nii.gz'%FSLDIR
    template_brain='%s/data/standard/MNI152_T1_2mm_brain.nii.gz'%FSLDIR
    template_mask='%s/data/standard/MNI152_T1_2mm_brain_mask.nii.gz'%FSLDIR


    cmd='%s/ANTS 3 -m CC[%s,%s/%s_bfc.nii.gz,1,4] -o %s_ANTS.nii.gz %s'%(ANTSPATH,template_wholehead,subdir,
                                                                       image_prefix,image_prefix,PARAMS)
    print cmd
    if not test:
        os.system(cmd)
    
    # create a version of the MNI mask aligned to subject space
    cmd='WarpImageMultiTransform 3 %s %s/%s_brain_mask.nii.gz -R %s/%s.nii.gz -i %s/%s_ANTSAffine.txt %s/%s_ANTSInverseWarp.nii.gz --use-NN'%(template_mask,subdir,image_prefix,subdir,image_prefix,subdir,image_prefix,subdir,image_prefix)
    print cmd
    if not test:
        os.system(cmd)

    cmd='fslmaths %s/%s.nii.gz -mul %s/%s_brain_mask.nii.gz %s/%s_brain.nii.gz'%(subdir,image_prefix,subdir,image_prefix,subdir,image_prefix)
    print cmd
    if not test:
        os.system(cmd)

    # rerun warp from stripped highres to stripped template

    cmd='%s/ANTS 3 -m CC[%s,%s/%s_brain.nii.gz,1,4] -o %s_ANTSstd.nii.gz %s'%(ANTSPATH,template_brain,subdir,image_prefix,image_prefix,PARAMS)
    print cmd
    if not test:
        os.system(cmd)

    cmd='WarpImageMultiTransform 3 %s/%s_brain.nii.gz %s/%s_reg2std_ANTS.nii.gz -R %s %s/%s_ANTSWarp.nii.gz %s/%s_ANTSAffine.txt'%(subdir,image_prefix,subdir,image_prefix,template_brain,subdir,image_prefix,subdir,image_prefix)

    print cmd
    if not test:
        os.system(cmd)

if __name__ == '__main__':
    main()
