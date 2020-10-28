# Evaluate Inverse Methods

This repository contains codes to run the analyses presented in: 

*Towards optimal source reconstruction of resting MEG of the human brain: performance, precision, and parcellation*, Under review. [Preprint available on bioRxiv](https://doi.org/10.1101/2020.01.12.903302)

The main function is **compare_algorithms.m**, which runs all metrics. This function does the following: 
1. Source reconstructs the input data using the function **f_source_reconstruction.m** also available in this repository. Note that this function uses a fixed regularization parameter for all algorithms. In the later revised versions of the manuscript (under review), we altered this to calculate regularization parameter based on SNR. 
1. If any of the input measures are 'rsq' or 'rsqcv' (or input measures is left empty): Calculates sensor variance explained by the source solution, either with the full solution (if 'rsq') or by using leave-one-out cross validation ('rsqcv') using the function **f_variance_analysis.m**.
1. If any of the input measures are 'le' or 'sect' (or input measures is left empty): Calculates localization error ('le') and spatial extent of cross talk ('sect') using the function **f_resolution_analysis.m**. 

At present, all analyses are those presented in the preprint version. The version under review has some small modifications (namely 10-fold cross validation instead of leave-one-out, normalization of variance explained against variance explained in empty room data, and the inclusion of the spatial extent of point spread), and upon acceptance these updates and any further changes will be uploaded. 

The functions included in this repository are dependent on the [Fieldtrip toolbox for M/EEG Analysis](www.fieldtriptoolbox.org), which must be installed for the codes to run.

### Reduced Atlas
The sub-repository **evaluate_inverse_methods/reduced_atlas** contains files for a 230 ROI reduction of the Human Connectome Project's multimodal parcellation of the human brain based on MEG resolution properties. This 230 ROI atlas is presented in the version of the manuscript currently under review, and follows the methodology described in the preprint linked above to derive a 250 ROI atlas (available in **evaluate_inverse_methods/reduced_atlas/preprint_250**). 

The atlas (**HCP230.mat**) is uploaded as a 1mm volumetric in the same format as the Fieldtrip template atlases (i.e. the output of a Fieldtrip function ft_read_atlas), and coregistered with the [Fieldtrip template anatomy](https://www.fieldtriptoolbox.org/template/), which is the Colin27 MRI in the MNI coordinate system. Therefore the atlas should be ready to integrate with the template head models, source models, and electrode models included with Fieldtrip. Two of the Fieldtrip template sourcemodels (4mm volumetric grid and 8196 dipole canonical cortical mesh) are included already coregistered to the atlas. 

For non-template participants (e.g. when anatomical MRI is available), we have included a function **align_individual2atlas.m** which includes a range of methods to align the atlas to an individual forward model. Example scripts demonstrating the useage of this function are given in **evaluate_inverse_methods/reduced_atlas/examples**.
