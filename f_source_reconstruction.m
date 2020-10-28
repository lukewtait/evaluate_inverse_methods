function source = f_source_reconstruction(method,data,sourcemodel)

% Note: Here we set regularization parameter cfg.lambda=0.05 for all algorithms. 
% This does not allow for fair comparison between metrics. We will update a later
% version upon acceptance of the manuscript (under review), which allows for 
% estimating cfg.lambda based on the predicted SNR of the data, giving consistent
% estimates between algorithms.

% General options
cfg = struct ; 
cfg.method = method ; % 'lcmv','wlcmv','mne','wmne','sloreta', or 'eloreta'
cfg.sourcemodel = sourcemodel ; % from output of ft_prepare_sourcemodel
cfg.headmodel = sourcemodel.cfg.headmodel ; % from output of ft_prepare_headmodel
cfg.grad = sourcemodel.cfg.grad ; % from output of ft_read_sens and aligned to headmodel
cfg.senstype = 'meg' ;
cfg.keepfilter = 'yes' ; % keep filters for resolution analysis

% Method specific options
switch method
    case 'lcmv'
        cfg.lambda = 0.05 ; % set regulatization parameter
        [~,source] = evalc('ft_sourceanalysis(cfg,data)') ; % do Fieldtrip source analysis
        
    case 'wlcmv'
        cfg.method = 'lcmv' ; % first construct the LCMV solution
        cfg.lambda = 0.05 ; % set regularization parameter
        [~,source] = evalc('ft_sourceanalysis(cfg,data)') ; % do Fieldtrip source analysis
        filt = cell2mat(source.avg.filter) ; % get filter
        filt = filt./vecnorm(filt,2,2) ; % normalize by vector norm of filter
        S = filt*data.avg ; % construct source data
        for i = 1:size(filt,1) 
            source.avg.mom{i} = S(i,:) ; % add source data back to the source structure
            source.avg.filter{i} = filt(i,:) ; % add filter back to the source structure
        end
        
    case {'mne','sloreta','eloreta'}
        cfg.(method).lambda = 0.05 ; % set regularization parameter
        [~,source] = evalc('ft_sourceanalysis(cfg,data)') ; % do Fieldtrip source analysis
        
    case 'wmne'
        cfg.method = 'mne' ; 
        cfg.mne.lambda = 0.05 ; % set regularization parameter
        lf = cell2mat(sourcemodel.leadfield) ; % get leadfield matrix
        cfg.mne.sourcecov = sparse(diag(1./sqrt(sum(lf.^2)))) ; % set weights to 1 over leadfield norm
        cfg.mne.scalesourcecov = true ; % scale to uniformity for comparison with MNE/eLORETA
        [~,source] = evalc('ft_sourceanalysis(cfg,data)') ; % do Fieldtrip source analysis 
        
    otherwise
        error('Input method not known')

end
