%% Example scripts for aligning individual sourcemodels to the atlas
clear, clc, close all

% Make sure fieldtrip is on the path
ft_defaults ; 

% Make sure the codes and atlas are on the path
addpath ..

% Load the atlas
atlas = load('HCP230.mat') ; 

% To load the individual MRI, you require the Fieldtrip tutorial data, which
% can be downloaded from 
%    ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/tutorial/Subject01.zip
% You need to unzip this file and ensure that Subject01.mri is in this 
% directory. 
% 
% Below we slightly preprocess the MRI data, and then use it to create a
% Fieldtrip headmodel structure (vol), generated from the MRI using the
% code below (3 layer BEM model) following steps described in Fieldtrip
% tutorial (http://www.fieldtriptoolbox.org/tutorial/headmodel_eeg_bem/). 


% Load mri
mri = ft_read_mri('Subject01.mri') ; 

% Process mri
mri = ft_convert_coordsys(mri,'acpc') ; % conver to acpc coordinates
cfg = struct ; 
cfg.resolution = 1 ; 
mri = ft_volumereslice(cfg,mri) ; % reslice to 1mm
 
% Segment
cfg = struct ; 
cfg.output = {'brain','skull','scalp'} ; % for 3 layer BEM
seg = ft_volumesegment(cfg,mri) ; 
 
% Create mesh
cfg = struct ; 
cfg.tissue = {'brain','skull','scalp'} ; 
cfg.numvertices = [3000 2000 1000] ; 
bnd = ft_prepare_mesh(cfg,seg) ; 

% Make headmodel
cfg = struct ; 
cfg.method = 'bemcp' ; 
vol = ft_prepare_headmodel(cfg,seg) ;


%% Method 1: Default settings
% By default, the function does the following: 
% 1) Reslices the atlas to 4mm resolution, and creates a 3d regular grid in
%    template space based on voxels with a lael in the atlas
% 2) Nonlinear warps individual MRI onto template MRI
% 3) Applies inverse warp to grid to map template grid to individual space

sourcemodel_defaults = align_individual2atlas(mri) ; 
% This is equivalent to: 
% sourcemodel_defaults = align_individual2atlas([],mri,atlas) ; 

% plot
figure(1); plot_sourcemodel2atlas(sourcemodel_defaults,vol) ; 


%% Method 2: Specify a template grid 
% Instead of using the atlas to create a template grid, you can use one of
% Fieldtrip's templates (or make your own grid in template space) by
% specifying cfg.grid (this is the same format as cfg.grid in
% ft_prepare_sourcemodel) Fieldtrip has a range of template grids which are
% described at:
%    http://www.fieldtriptoolbox.org/template/sourcemodel/
% These include 3d volumes at 4,5,6,7.5,8,10mm resolution, and canonical
% cortical meshes with 5124,8196,20484 vertices. Here, we will show an
% example for the 5124 voxel cortical mesh.
% 
% In this case, the steps are: 
% 1) Interpolate tissue labels onto grid, since both are assumed to be in
%    the same template space.
% 2) Nonlinear warps individual MRI onto template MRI
% 3) Applies inverse warp to grid to map template grid to individual space

% get the fieldtrip template
ft_path = fileparts(which('ft_defaults')) ; % find the path to fieldtrip
cfg = struct ; 
cfg.grid = ft_read_headshape(sprintf('%s/template/sourcemodel/cortex_5124.surf.gii',ft_path)) ; % load example surface sourcemodel
cfg.method = 'template2individual' ; % we need to specify that we wish to warp the template grid to the individual

sourcemodel_template = align_individual2atlas(cfg,mri,atlas) ; 

% plot
figure(2); plot_sourcemodel2atlas(sourcemodel_template,vol) ; 


%% Method 3: Use a pre-labelled grid
% Two of the Fieldtrip template grids have been aligned to the atlas by
% default, meaning there is no need for the interpolation to the atlas in
% template space. The function align_individual2atlas will not interpolate
% if tissue labels are present, so you can alternatively inverse warp those
% pre-aligned source models. The steps are essentially identical to the
% previous section, except step 1 is not required as it has already been
% performed. 

cfg = struct ; 
cfg.grid = load('cortex_8196_HCP230.mat') ; % load example surface sourcemodel, which has tissue included
cfg.method = 'template2individual' ; % we need to specify that we wish to warp the template grid to the individual

sourcemodel_prealigned = align_individual2atlas(cfg,mri,atlas) ; 

% plot
figure(3); plot_sourcemodel2atlas(sourcemodel_prealigned,vol) ; 


%% Method 4: Define the grid in individual space
% All methods so far have defined a grid in template space, aligned this
% grid to the atlas, and then inverse warped to individual space. This
% method has the advantage that dipoles are consistent between
% participants; dipole i will have the same MNI coordinates and the same
% atlas label in all participants. However, the downside is that the
% accuracy of the individual space grid depends on the accuracy of the
% inverse warp. Furthermore, since the warp is nonlinear, the grid
% is no longer a regular 3d grid in individual space, and heads of
% different sizes will reduce in different spacing between dipoles (e.g.
% after warping, a template regular 4mm grid will not be regular and not
% 4mm spaced). Cortical surfaces may not be accurate using a template
% either, particularly with regards to orientations of dipoles. 
% 
% An alternative option is to define your grid individually for
% participants. In this case, dipoles no longer have the same MNI
% coordinates or labels, or are even necessarily defined at the same
% locations. This means groups level statistics need to be interpolated
% onto a template.  However more accurate cortical models can be defined.
% 
% To make a grid in individual space, you can use ft_prepare_sourcemodel
% to generate a 3d grid (or e.g. use Freesurfer to get a cortical mesh),
% and supply this as cfg.grid in the same manner as example 2 - the only
% difference is you must specify cfg.method = 'individual2template'.
% Alternatively, you can supply a cfg.headmodel and cfg.method =
% 'individual2template' (and optionally cfg.resolution) to generate a grid.
% 
% Steps are: 
% 1) (If no grid is supplied) Generate an individual grid from the supplied
%    individual headmodel.
% 2) Nonlinear warp individual MRI onto template MRI
% 3) Applies warp to grid to map individual grid to template space
% 4) Interpolate tissue labels onto grid, since both are assumed to be in
%    the same template space. You do not need to warp back to individual
%    space; you can just use the original (unwarped) grid with the labels
%    given to the warped grid.

cfg = struct ; 
cfg.headmodel = vol ; % supply headmodel 
cfg.method = 'individual2template' ; % required to let code know to make grid in individual space
sourcemodel_individual = align_individual2atlas(cfg,mri,atlas) ; 

% plot
figure(4); plot_sourcemodel2atlas(sourcemodel_individual,vol) ; 


%% Method 5: Using Freesurfer
% In our original paper on the atlas (doi: 10.1101/2020.01.12.903302), we
% used Freesurfer to generate surface atlases. The Freesurfer average
% HCPMMP.annot files were downloaded from
% https://figshare.com/articles/HCP-MMP1_0_projected_on_fsaverage/3498446
% and Freesurfer was used to extract cortical surfaces and labels according
% to this atlas. We used Brainstorm to import these anatomies and their
% labels, and create a 10,000 vertex cortical mesh of the mid points
% between white and pial surfaces. We then merged ROIs according to the
% description in our paper. Here, we have included the code to merge these
% ROIs and make new labels. "mesh" is the output of a Freesurfer/Brainstorm
% analysis for an example (colin27) participant, containing an example
% 10,000 vertex midpoint cortical mesh following approximately the
% methodology from our paper.
% 
% No warping is required for this method, as the output mesh is already in
% individual space and all warping to a template and labelling was
% performed by Freesurfer. Note that the vertices are in Brainstorm ctf
% coordinates not MNI coordinates, and the files need reformatting to be
% placed in Fieldtrip format. However, this script is just to demonstrate
% useage of f_reduce_atlas230. 

load('mesh.mat') % Brainstorm output mesh

HCPatlas = mid_10003V.Atlas(2).Scouts ; % We know it is index 2 because mid_10003V.Atlas(2).Name is 'template_HCPMMP1'
atlas = f_reduce_atlas230(HCPatlas) ; % downsample to 230 ROIs

% plot
figure(5)
tissue = nan(length(mid_10003V.Vertices),1) ; 
for i = 1:length(dsatlas)
    ind = dsatlas(i).Vertices ; 
    tissue(ind) = i ; 
end
patch('Vertices',mid_10003V.Vertices,'Faces',mid_10003V.Faces,'FaceColor','flat','FaceVertexCData',tissue)
rng('default') ; colormap(repmat(rand(115,3),2,1)) ;
