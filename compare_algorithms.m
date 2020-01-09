function measures_out = compare_algorithms(data,sourcemodel,measures)
% Main function to compute metrics described in the Tait et al (2020)
% manuscript for comparing source reconstruction algorithms. 
% 
% Inputs: 
% - data: A fieldtrip data structure. This should either be a fieldtrip raw
%         data structure (ft_datatype_raw), in which case a timelock
%         analysis will be performed to average over trials (if present)
%         and calculate covariance. Alternatively it can be a fieldtrip
%         timelock data structure (ft_datatype_timelock) including a
%         covariance output. 
% - sourcemodel: A fieldtrip leadfield structure, the output of
%         ft_prepare_leadfield. Importantly, the headmodel and grad fields
%         should not have been cleared from the input cfg structure (i.e.
%         sourcemodel.cfg.headmodel and sourcemodel.cfg.grad should not be
%         empty). 
% - measures: optional input. A cell array consisting of any combination of
%         'rsq' (r^2), 'rsqcv' (r^2_CV), 'le' (LE), or 'sect' (SECT). By
%         default all measures are calculated. 
% 
% Output: 
% - measures_out: A structure consisting of all measures that have been
%         calculated. For each measure, each column represents an
%         algorithm, and the order of algorithms is shown in
%         measures_out.methods. 
  
%% Check input measures
if nargin < 3 % default - do all
    measures = {'rsq','rsqcv','le','sect'} ; 
end
% defaults
do_rsq = false ; 
do_rsqcv = false ; 
do_le = false ; 
do_sect = false ; 
% make true if you wish to calculate each measure
if any(strcmp(measures,'rsq')) 
    do_rsq = true ; 
end
if any(strcmp(measures,'rsqcv')) 
    do_rsqcv = true ; 
end
if any(strcmp(measures,'le')) 
    do_le = true ; 
end
if any(strcmp(measures,'sect')) 
    do_sect = true ; 
end
% check at least one measure is calculated
measures_out = struct ; % initialize output
if ~(do_rsq||do_rsqcv||do_le||do_sect)
    return
end

%% Get covariance of data

if ~isfield(data,'cov')
    msg = 'Calculating data covariance...' ; fprintf(msg)
    cfg = struct ; 
    cfg.covariance = 'yes' ; 
    [~,data] = evalc('ft_timelockanalysis(cfg,data)') ; 
end

%% Choose methods 

methods = {'lcmv','wlcmv','mne','wmne','sloreta','eloreta'} ; 
measures_out.methods = methods ; 

%% Calculate r^2
if do_rsq || do_le || do_sect % we need to do this section if we are doing resolution metrics to calculate filter
    % display output
    fprintf(repmat('\b',1,length(msg))) ; 
    msg = 'Calculating variance explained...' ; fprintf(msg)

    % initialize
    rsq = zeros(size(methods)) ; 
    filt = cell(size(methods)) ; 

    % do each algorithm
    for i = 1:numel(methods)
        % update output
        msg2 = sprintf(' %s: %d of %d',methods{i},i,numel(methods)) ; fprintf(msg2) ; 

        % source reconstruct
        source = f_source_reconstruction(methods{i},data,sourcemodel) ; 

        % calculate rsq
        rsq(i) = f_variance_analysis(data.avg,source,sourcemodel) ;

        % get filter for resolution metrics
        filt{i} = cell2mat(source.avg.filter) ; 
        % account for fact that sometimes source.avg.filter is a column, sometimes a row
        if size(filt{i},1) == 1 || size(filt{i},1) == numel(filt{i})
            filt{i} = cell2mat(source.avg.filter') ; 
        end

        % update output
        fprintf(repmat('\b',1,length(msg2)))
    end
    
    % add rsq to the measures_out structure
    if do_rsq
        measures_out.rsq = rsq ; 
    end
    clear rsq
end
    
%% Calculate resolution metrics

if do_le || do_sect % we calculate both metrics either way
    % display output
    fprintf(repmat('\b',1,length(msg))) ; 
    msg = 'Calculating resolution metrics...' ; fprintf(msg)

    % initialize 
    le = zeros(size(filt{1},1),numel(methods)) ; 
    sect = zeros(size(filt{1},1),numel(methods)) ; 
    
    % get leadfield in matrix form
    lf = cell2mat(sourcemodel.leadfield) ; 

    % do each algorithm
    for i = 1:numel(methods)
        % update output
        msg2 = sprintf(' %s: %d of %d',methods{i},i,numel(methods)) ; fprintf(msg2) ; 

        % calculate resolution metrics
        [le(:,i),sect(:,i)] = f_resolution_analysis(lf,filt{i},sourcemodel.pos) ;

        % update output
        fprintf(repmat('\b',1,length(msg2)))
    end
    
    % add le to the measures_out structure
    if do_le
        measures_out.le = le ; 
    end
    clear le
    
    % add sect to the measures_out structure
    if do_sect
        measures_out.sect = sect ; 
    end
    clear sect
end


%% Loop over electrodes

if do_rsqcv
    
    % save memory - can be expensive
    clearvars -except data measures measures_out methods sourcemodel
    
    % set up parallel pool
    delete(gcp('nocreate'))
    pool = parpool(3) ; % I use 3 nodes here as I run out of memory with more 
   
    % initialize
    rsqcvi = nan(length(data.label),numel(methods)) ; 

    % Need to make sure variable names aren't ambiguous, can get a
    % transparency error otherwise
    sourcemodel = sourcemodel ; 
    headmodel = sourcemodel.cfg.headmodel ; 
    sens = sourcemodel.cfg.grad ; 
    data = data ;

    % make progress bar
    fprintf('\nCalculating cross validated rsq, this may take a while - Progress:\n');
    fprintf([repmat('.',1,length(data.label)) '\n\n']);

    parfor elc = 1:length(data.label) ; 
        % initialize cross validated variance explained for this sensor
        pp = nan(1,length(measures)) ; 

        % make training set - remove sensor elc
        sourcemodel_ds = sourcemodel ; 
        sourcemodel_ds.label(elc) = [] ;
        for i = 1:length(sourcemodel.leadfield)
            sourcemodel_ds.leadfield{i}(elc,:) = [] ;
        end
        data_ds = data ; 
        data_ds.label(elc) = [] ; 
        data_ds.avg(elc,:) = [] ; 
        data_ds.var(elc,:) = [] ; 
        data_ds.dof(elc,:) = [] ; 
        data_ds.cov(elc,:) = [] ; 
        data_ds.cov(:,elc) = [] ; 
        
        % loop over algorithms
        for i = 1:numel(methods)
            source = f_source_reconstruction(methods{i},data_ds,sourcemodel_ds) ; % do Fieldtrip source analysis
            pp(i) = f_variance_analysis(data.avg(elc,:),source,sourcemodel,elc) ; % do source data compare
        end
    
        rsqcvi(elc,:) = pp ; 
    
        fprintf('\b|\n');
    end
    
    measures_out.rsqcv = rsqcvi ; 
end


end % end function

