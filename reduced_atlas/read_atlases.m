% Read volume and surface based HCP250 atlas into Matlab R2019b
clear, clc, close all

% Read volumetric atlas
volHCP = niftiread('vol_HCP250.nii') ; % atlas
T = dlmread('vol_transform.txt') ; % transform

% Plot volumetric atlas
grey = find(volHCP>0) ; % voxels which are labelled
P = ones(4,length(grey)) ; % initialize coordinates
[P(1,:),P(2,:),P(3,:)] = ind2sub(size(volHCP),grey) ; % get coordinates 
P = (T*P)' ; % do the affine transformation
figure(1), clf 
scatter3(P(:,1),P(:,2),P(:,3),5,volHCP(grey),'filled') % plot the voxels, with color corresponding to atlas label
cmap = dlmread('ROI_colors.txt') ; % get original colours from atlas printed in Glasser et al for reproducibility
colormap(cmap)

% Read surface atlas
surf = struct ; 
surf.vertices = dlmread('surf_vertices.txt') ; 
surf.faces = dlmread('surf_faces.txt') ; 
surfHCP = dlmread('surf_HCP250.txt') ; 
figure(2), clf
patch(surf,'facecolor','flat','facevertexCdata',surfHCP,'linestyle','none')
camlight headlight; lighting gouraud; material dull ; 
colormap([0.8,0.8,0.8 ; cmap]) % include the [0.8,0.8,0.8] for the ROIs with label 0 (i.e. not included in atlas)
