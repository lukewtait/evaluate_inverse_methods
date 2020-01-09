function [le,sect] = f_resolution_analysis(LF,K,pos) ; 
% Function to calculate resolution metrics LE and SECT. 
% Inputs: 
% - LF: (Scalar) leadfield matrix
% - K: Spatial filter, i.e. inverse of LF (e.g. wLCMV, eLORETA, etc)
% - pos: Nx3 matrix with the x,y,z coordinates of each dipole, where N is
%        number of dipoles
% 
% Outputs: 
% - le: Nx1 vector of localization errors for each dipole
% - sect: Nx1 vector of spatial extent of cross talk for each dipole

%% Make resolution matrix
% Eqns 5-6 in manuscript
R = K*LF ; 

%% Localization error
% Eqn 7 in manuscript
[~,idxpeak] = max(abs(R)) ; % for each dipole, find the dipole with peak of the PSF
for i = 1:length(idxpeak) % loop over dipoles
    le(i,1) = sqrt(sum((pos(i,:)-pos(idxpeak(i),:)).^2)) ; % calculate difference in position between PSF peak and origin
end

%% Spatial extent of cross talk
% Eqn 8 in manuscript
R = R.^2 ; % formula calls only for Rij^2, so square each element of R

% get distances between each pair of dipoles
for i = 1:3 % x,y,z
    dx(:,:,i) = (pos(:,i) - pos(:,i)').^2 ; 
end
dx = sum(dx,3) ; % distance squared is used in formula, so don't worry about square root

sect = sqrt(sum(dx.*R,2)./sum(R,2)) ; % calculate SECT
