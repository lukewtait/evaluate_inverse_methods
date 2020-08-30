function sourcemodel = align_individual2atlas(cfg,mri,atlas) ; 
% Function to align an individual MRI to an atlas. 
% 
% Useage: 
% 
% sourcemodel=ALIGN_INDIVIDUAL2ATLAS(mri) will create a Fieldtrip 
% sourcemodel structure aligned to the individual MRI supplied. Input mri
% should be a Fieldtrip mri structure containing an anatomical MRI. In
% brief, the function generates a regular 3d grid with resolution of 4mm in
% template space by downsampling the atlas. The MRI is warped to the
% template MRI (fieldtrip-yyyymmdd/template/anatomy/single_subj_T1_1mm.nii)
% and then the inverse warp is applied to the 3d grid to map from template
% to individual space. The atlas HCP230.mat must be on the path.
% 
% sourcemodel=ALIGN_INDIVIDUAL2ATLAS([],mri,atlas) allows the user to
% specify the atlas, and uses the same steps as in the previous
% usage. Input atlas must be a Fieldtrip atlas structure. If
% atlas=load('HCP230.mat'), then this call is identical to the previous
% useage. 
% 
% sourcemodel=ALIGN_INDIVIDUAL2ATLAS(cfg,mri,atlas) allows the user to
% specify configurations. The input structure cfg is described below. If
% cfg.grid and cfg.grid.tissue are specified, no atlas is required. 
% 
% Options: 
% 
% cfg.method: Can be either 'template2individual' or 'individual2template'.
% If cfg.grid (see below) is not specified, 'template2individual' makes the
% code acts as above, i.e. downsample the atlas to generate a regular grid
% in template space, and inverse warp this grid to individual space.
% Conversely, if 'individual2template' is chosen and cfg.grid is not
% specified, a regular 3d grid will be generated in individual space,
% warped to template space, and the atlas will be interpolated onto the
% grid. In this case, cfg.headmodel must be specified for generating the
% individual grid. 
%    If cfg.grid is specified, this grid will be aligned to the atlas 
% instead of a regular 3d grid generated within this code. In this case,
% cfg.method specifies whether cfg.grid is in template space (cfg.method =
% 'template2individual'), in which case atlas labels will be interpolated
% onto the grid and then the grid will be inverse warped to individual
% space, or if cfg.grid is in individual space (cfg.method =
% 'individual2template'), in which case the grid will be warped to template
% space and the atlas labels interpolated. 
% 
% cfg.grid: Allows you to specify the grid on which the atlas labels should
% be interpolated. This can be in either individual space or the same space
% as the Fieldtrip template MRI, and the space should be described by
% cfg.method. It can be specified in the same manner as the input cfg.grid
% for the function ft_prepare_sourcemodel, as this function is used to
% prepare the sourcemodel from cfg.grid. In particular, you can specify
% cfg.grid.pos as the locations of the dipoles, and (if the grid is a
% cortical mesh as opposed to a 3d volume) cfg.grid.tri are the faces of
% the mesh. 
% 
% cfg.headmodel: A Fieldtrip headmodel structure, used to generate an
% individual grid if cfg.grid is not specified and cfg.method =
% 'individual2template'. 
% 
% cfg.resolution: The resolution of regular 3d grids generated if cfg.grid
% is not specified. By default, 4mm. 

%% Check inputs

% Deal with case of one input
if nargin == 1
    mri = cfg ; 
    cfg = [] ; 
end


% Check cfg is a structure
if ~isstruct(cfg)
    if isempty(cfg)
        cfg = struct ; 
    else
        error('cfg must be specified as a structure')
    end
end

% Check mri
mri = ft_datatype_volume(mri) ; 

% Check atlas
if nargin == 3
    atlas = ft_datatype_volume(atlas) ; 
    if ~isfield(atlas,'tissue')
        error('atlas must be a valid atlas structure with datatype volume and containing the field "tissue"')
    end
    tissuelabel = atlas.tissuelabel ;
    hasatlas = true ; 
elseif nargin == 1
    atlas = load('HCP230.mat') ; 
    tissuelabel = atlas.tissuelabel ; 
    hasatlas = true ; 
else
    hasatlas = false ; 
end

% Get path to fieldtrip
ft_path = fileparts(which('ft_defaults')) ; % find the path to fieldtrip
if isempty(ft_path)
    warning('Fieldtrip must be on the path')
end
ft_defaults ; % ensure fieldtrip is added

%% Set defaults

% Check method
method = ft_getopt(cfg,'method','template2individual') ; 
if ~any(strcmp(method,{'individual2template','template2individual'})) ; 
    error('cfg.method must be either individual2template or template2individual')
end

% Check for grid
grid = ft_getopt(cfg,'grid',[]) ; if isempty(grid) ; hasgrid = false ; else ; hasgrid = true ; end

% Check whether an atlas was supplied with the grid
hastissue = isfield(grid,'tissue') ; if hastissue ; tissue = grid.tissue ; end
if ~hasatlas && ~hastissue
    error('atlas must be supplied unless cfg.grid.tissue is specified')
end
if hastissue
    if isfield(grid,'tissuelabel')
        tissuelabel = grid.tissuelabel ; 
    elseif hasatlas
        warning('No tissuelabels supplied, assuming labels match supplied atlas')
    else
        warning('No tissuelabels supplied,')
        tissuelabel = 'not supplied' ; 
    end
end
if hasatlas && hastissue
    warning('Ignoring atlas, as cfg.grid.tissue is supplied')
    clear atlas; hasatlas = false ; 
end
if strcmp(method,'individual2template') && hastissue % If labels have already been supplied and the grid is already in individual space, return
    fprintf('Tissue labels already aligned to atlas in individual space, no alignment necessary')
    sourcemodel = ft_datatype_volume(grid) ;
    return
end

% Resolution of template grid
resolution = ft_getopt(cfg,'resolution',4) ; % not used if grid is specified

% Check if headmodel is specified
vol = ft_getopt(cfg,'headmodel',[]) ; if isempty(vol) ; hasvol = false ; else ; hasvol = true ; end



%% Get grid of source points

if ~hasgrid & strcmp(method,'template2individual') % make a template grid
    
    fprintf('No grid supplied, prepraring grid in template space...\n')
    
    %% downsample atlas
    tmpcfg = struct ; 
    tmpcfg.resolution = resolution ; 
    tmpcfg.method = 'nearest' ; 
    atlas_grid = ft_volumereslice(tmpcfg,atlas) ; 
    
    % get positions
    ind = find(atlas_grid.tissue>0) ; 
    [x,y,z] = ind2sub(atlas_grid.dim,ind) ; 
    pos = [x,y,z,ones(length(x),1)] ; 
    pos = (atlas_grid.transform*pos')' ; pos = pos(:,1:3) ; clear x y z
    
    % make grid
    grid = struct ; 
    grid.pos = pos ; 
    grid.unit = atlas_grid.unit ; 
    
    % ensure inside brain
    tmpcfg = struct ; 
    tmpcfg.grid = grid ;
    tmpcfg.headmodel = vol ;
    grid = ft_prepare_sourcemodel(tmpcfg) ; 
    
    % get ROI labels
    tissue = atlas_grid.tissue(ind) ; 
    hastissue = true ; 
    
elseif ~hasgrid & strcmp(method,'individual2template') % make an individual grid
    
    fprintf('No grid supplied, prepraring grid in individual space...\n')
    
    if ~hasvol % can only make individual grid if headmodel is supplied
        error('To make a grid in individual space, an indivdiaual headmoel must be supplied as cfg.headmodel')
    else
        
        vol = ft_datatype_headmodel(vol) ; % check input is headmodel
        
        % Make a template grid from the headmodel 
        tmpcfg = struct ; 
        tmpcfg.headmodel = vol ; 
        tmpcfg.grid.resolution = resolution ; 
        tmpcfg.unit = 'mm'  ; 
        grid = ft_prepare_sourcemodel(tmpcfg) ; 
        
    end
    
    
    
else % has grid
    
    fprintf('Preparing user supplied grid...\n')
    
    % Make the grid using user specified inputs
    tmpcfg = struct ; 
    tmpcfg.grid = grid ; 
    if hasvol & strcmp(method,'individual2template')
        tmpcfg.headmodel = vol ;
    end
    grid = ft_prepare_sourcemodel(tmpcfg) ; 
    
end
    

%% Warp individual MRI, make both template and individual space grids, and get tissue labels

% warp mri
cfg = struct ; 
cfg.parameter = 'all ' ; 
cfg.template = sprintf('%s/template/anatomy/single_subj_T1_1mm.nii',ft_path) ; 
cfg.nonlinear = 'yes' ; 
mri_warped = ft_volumenormalise(cfg,mri) ; 

% define grid in individual and template spaces
switch method
    case 'template2individual' 
        grid_template = grid ; % if already defined in template space
        grid.pos = ft_warp_apply(inv(mri_warped.initial), ...
            ft_warp_apply(mri_warped.params, grid_template.pos, 'sn2individual')); % inverse warp
    case 'individual2template' 
        grid_template = grid ; % initialize
        grid_template.pos = ft_warp_apply(mri_warped.params,...
            ft_warp_apply(mri_warped.initial,grid.pos), 'individual2sn'); % warp
end

% Make mask against atlas (not required for atlas or surface)
if hasatlas & (~hastissue & ~isfield(grid_template,'tri'))
    tmpcfg = struct ; 
    tmpcfg.atlas = atlas ; 
    tmpcfg.inputcoord = 'mni' ; 
    tmpcfg.roi = tissuelabel ; 
    mask = ft_volumelookup(tmpcfg,grid_template) ; 
    
    grid_template.inside = grid_template.inside & mask(:) ; 
end

% get tissue labels
if ~hastissue 
    tissue = interp_tissue(grid_template,atlas) ;
end

% make output
sourcemodel = grid ; 
sourcemodel.inside = grid_template.inside ; 
sourcemodel.tissue = tissue ; 
sourcemodel.tissuelabel = tissuelabel ; 

% Clean up output
if isfield(sourcemodel,'xgrid')
    sourcemodel = rmfield(sourcemodel,{'xgrid','ygrid','zgrid','dim'}) ; 
end
mask = sourcemodel.inside & sourcemodel.tissue > 0 ; 
sourcemodel.pos = sourcemodel.pos(mask,:) ; % remove dipoles outside of brain
sourcemodel.tissue = sourcemodel.tissue(mask) ; 
sourcemodel.inside = sourcemodel.inside(mask) ; 


end % end main function

function tissue = interp_tissue(grid_template,atlas)

        % get tissue labels
        if isfield(grid_template,'tri') % mesh
            
            % for a mesh, we assume all points are in the cortex, so make a
            % grid of scattered tissue points (i.e. cannot interpolate to
            % 0)
            ind = find(atlas.tissue>0) ; 
            [x,y,z] = ind2sub(size(atlas.tissue),ind) ; 
            pos = [x,y,z,ones(length(x),1)] ; 
            pos = (atlas.transform*pos')' ; pos = pos(:,1:3) ; clear x y z
            
            % interpolate labels
            atlgrid.pos = pos ; 
            atlgrid.tissue = atlas.tissue(ind) ; 
            tmpcfg = struct ; 
            tmpcfg.parameter = 'tissue' ;
            tmpcfg.interpmethod = 'nearest' ; 
            interp = ft_sourceinterpolate(tmpcfg,atlgrid,grid_template) ; 
            
        else % 3d grid
            
            % interpolate labels
            tmpcfg = struct ; 
            tmpcfg.parameter = 'tissue' ;
            tmpcfg.interpmethod = 'nearest' ; 
            interp = ft_sourceinterpolate(tmpcfg,atlas,grid_template) ; 
            
        end
        
        tissue = interp.tissue(:) ; 
end
    
    
