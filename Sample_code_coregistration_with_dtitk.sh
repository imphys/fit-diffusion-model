### Sample code to co-register the orientation prior to subject space ###
# The code below shows how to co-register the orientation atlas to subject space for a single subject. Execute all commands from the main directory.
# Joor Arkesteijn (joorarkesteijn@gmail.com), Quantitative Imaging Group, Delft University of Technology, 2018.
	
### FSL and DTI-TK need to be properly installed ###
which fsl # Quick check if FSL is installed. If not, check: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/.
which dti_rigid_reg # Quick check if DTI-TK is installed If not, check: http://dti-tk.sourceforge.net/pmwiki/pmwiki.php.

### Define frequently-used locations and files ###
subject_dir=est_data/registration
atlas=atlas_data/commonspace_singletensor_reg.nii.gz
atlas_mask=atlas_data/commonspace_mask_reg.nii.gz

### First estimate a single tensor from your (already preprocessed) data ###
mkdir -p $subject_dir
bet orig_data/data_ecc orig_data/data_ecc_brain  -f 0.35 -g 0 -n -m
dtifit -k orig_data/data_ecc -b orig_data/bval -r orig_data/bvec -m orig_data/data_ecc_brain_mask -o $subject_dir/fsl

### Use fsl_to_dtitk to transform estimated tensor into dtitk-format ###
fsl_to_dtitk $subject_dir/fsl

### Always check if the tensors are orientated correctly! ###
# TVglyphView -in $subject_dir/fsl_dtitk.nii.gz -view axial
# TVglyphView -in $subject_dir/fsl_dtitk.nii.gz -view coronal

### Register subject to common space using recommended dtitk settings and compute inverse transformation ###
dti_rigid_reg $atlas $subject_dir/fsl_dtitk.nii.gz EDS 4 4 4 0.001
dti_affine_reg $atlas $subject_dir/fsl_dtitk.nii.gz EDS 4 4 4 0.001 1
dti_diffeomorphic_reg $atlas $subject_dir/fsl_dtitk_aff.nii.gz $atlas_mask 1 6 0.002
dfRightComposeAffine -aff $subject_dir/fsl_dtitk.aff -df $subject_dir/fsl_dtitk_aff_diffeo.df.nii.gz -out $subject_dir/fsl_dtitk_combined.df.nii.gz
dfToInverse -in $subject_dir/fsl_dtitk_combined.df.nii.gz

### Transform orientation atlas into dti-tk format ###
fslsplit atlas_data/commonspace_D1_mean500_vector atlas_data/D1_split_ -t
fslmaths atlas_data/D1_split_0000 -mul atlas_data/D1_split_0000 -add 0.2 atlas_data/D1_Dxx
fslmaths atlas_data/D1_split_0001 -mul atlas_data/D1_split_0000 -mul -1 atlas_data/D1_Dyx
fslmaths atlas_data/D1_split_0001 -mul atlas_data/D1_split_0001 -add 0.2 atlas_data/D1_Dyy
fslmaths atlas_data/D1_split_0002 -mul atlas_data/D1_split_0000 -mul -1 atlas_data/D1_Dzx
fslmaths atlas_data/D1_split_0002 -mul atlas_data/D1_split_0001 atlas_data/D1_Dzy
fslmaths atlas_data/D1_split_0002 -mul atlas_data/D1_split_0002 -add 0.2 atlas_data/D1_Dzz
fslmerge -t atlas_data/commonspace_D1_mean500_dtitk atlas_data/D1_D*
rm atlas_data/D1*

fslsplit atlas_data/commonspace_D2_mean500_vector atlas_data/D2_split_ -t
fslmaths atlas_data/D2_split_0000 -mul atlas_data/D2_split_0000 -add 0.2 atlas_data/D2_Dxx
fslmaths atlas_data/D2_split_0001 -mul atlas_data/D2_split_0000 -mul -1 atlas_data/D2_Dyx
fslmaths atlas_data/D2_split_0001 -mul atlas_data/D2_split_0001 -add 0.2 atlas_data/D2_Dyy
fslmaths atlas_data/D2_split_0002 -mul atlas_data/D2_split_0000 -mul -1 atlas_data/D2_Dzx
fslmaths atlas_data/D2_split_0002 -mul atlas_data/D2_split_0001 atlas_data/D2_Dzy
fslmaths atlas_data/D2_split_0002 -mul atlas_data/D2_split_0002 -add 0.2 atlas_data/D2_Dzz
fslmerge -t atlas_data/commonspace_D2_mean500_dtitk atlas_data/D2_D*
rm atlas_data/D2*

### Use the inverse transformation field to transform atlas orientation D1 and D2 to subject space ###
deformationSymTensor3DVolume -in atlas_data/commonspace_D1_mean500_dtitk.nii.gz -trans $subject_dir/fsl_dtitk_combined.df_inv.nii.gz -target $subject_dir/fsl_dtitk.nii.gz -out $subject_dir/D1_warped.nii.gz
deformationSymTensor3DVolume -in atlas_data/commonspace_D2_mean500_dtitk.nii.gz -trans $subject_dir/fsl_dtitk_combined.df_inv.nii.gz -target $subject_dir/fsl_dtitk.nii.gz -out $subject_dir/D2_warped.nii.gz

### Convert dtitk-format to tensor format ###
TVEigenSystem -in $subject_dir/D1_warped.nii.gz -type FSL
TVEigenSystem -in $subject_dir/D2_warped.nii.gz -type FSL

### Quite likely the data the warped atlas orientations need to be flipped in x (orientation issue between FSL/DTITK) ###
fslsplit $subject_dir/D1_warped_V1 $subject_dir/D1_warped_V1_ -t
fslmaths $subject_dir/D1_warped_V1_0000 -mul -1 $subject_dir/D1_warped_V1_0000_xflipped
fslsplit $subject_dir/D2_warped_V1 $subject_dir/D2_warped_V1_ -t
fslmaths $subject_dir/D2_warped_V1_0000 -mul -1 $subject_dir/D2_warped_V1_0000_xflipped
fslmerge -t $subject_dir/D1_warped_V1_xflipped $subject_dir/D1_warped_V1_0000_xflipped $subject_dir/D1_warped_V1_0001 $subject_dir/D1_warped_V1_0002
fslmerge -t $subject_dir/D2_warped_V1_xflipped $subject_dir/D2_warped_V1_0000_xflipped $subject_dir/D2_warped_V1_0001 $subject_dir/D2_warped_V1_0002

### Always check if the warped orientations look plausible (orientations correct in cc and cst? spatially aligned?)! ###
# fslview $subject_dir/fsl_FA $subject_dir/D1_warped_V1_xflipped # In FSLview visualize orientations as lines
# fslview $subject_dir/fsl_FA $subject_dir/D2_warped_V1_xflipped # In FSLview visualize orientations as lines

### If orientations are correct, merge both orientations into 1 file ###
fslmerge -t $subject_dir/D1D2_merged $subject_dir/D1_warped_V1_0000_xflipped $subject_dir/D1_warped_V1_0001 $subject_dir/D1_warped_V1_0002 $subject_dir/D2_warped_V1_0000_xflipped $subject_dir/D2_warped_V1_0001 $subject_dir/D2_warped_V1_0002
flirt -in $subject_dir/D1D2_merged -ref orig_data/data_ecc -out orig_data/D1D2_merged -applyxfm # Use the same orientation as data_ecc