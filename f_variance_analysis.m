function rsq = f_variance_analysis(X,source,sourcemodel,elc) ; 
% Function to calculate variance explained metrics.
% 
% Usage: 
% No cross validation
%    rsq = f_variance_analysis(data.avg,source,sourcemodel),
%    where data is a fieldtrip timelock structure (ft_timelockanalysis),
%    source is a fieldtrip source structure (ft_sourceanalysis), and
%    sourcemodel is a fieldtrip leadfield structure (ft_prepareleadfield).
% 
% Cross validation
%    rsqcv = f_variance_analysis(data.avg(elc,:),source_training,sourcemodel,elc),
%    where elc is the set of test electrodes, data is a fieldtrip timelock
%    structure, source_training is a fieldtrip source structure computed
%    using the training electrodes, and sourcemodel is a fieldtrip
%    leadfield structure.

% Inputs: 
% - X: Time course of data to be predicted (i.e. test electrodes)
% - source: output of Fieldtrip's ft_sourceanalysis function
% - sourcemodel: output of Fieldtrip's ft_prepareleadfield function (note:
%   you should check in sourcemodel.cfg that fields headmodel and grad have
%   not been cleared, as these are required). 
% - elc: the set of test electrodes, in vector form corresponding to their
%   location in the leadfield matrix (sourcemodel). 
% 
% Outputs: 
% - rsq/rsqcv: variance explained

if nargin == 3 % no cross validation
    
    % unpack the leadfield into an array
    LF = cell2mat(sourcemodel.leadfield) ; 

    % unpack the source data into an array and map into sensor space - note
    % sometimes source.avg.mom is a column, sometimes it is a row. Here, we
    % try column first, and if the dimensions are incorrect we try row.
    try % column
        S = cell2mat(source.avg.mom) ; % unpack source data
        S = LF*S ; % map back into sensor space
    catch % row
        S = cell2mat(source.avg.mom') ; % unpack source data
        S = LF*S ; % map back into sensor space
    end
    
    % correlate each electrode with prediction
    r = zeros(size(S,1),1) ; 
    for i = 1:size(S,1)
        r(i) = corr(S(i,:)',X(i,:)') ; 
    end
    rsq = mean(r.^2) ; % average over all electrodes
    
    
elseif nargin == 4 % cross validation
    
    % initialize predicted data to zeros
    S = zeros(size(X)) ; 
    
    % for each dipole in the estimated source data, map to the training
    % electrodes. Since forward mapping is linear, we can simply sum the
    % contributions of each dipole. 
    for i = 1:length(sourcemodel.leadfield)
        S = S + (sourcemodel.leadfield{i}(elc,:))*source.avg.mom{i} ; 
    end
    
    % correlate each test electrode with prediction
    r = zeros(size(X,1),1) ; 
    for i = 1:size(X,1)
        r(i) = corr(S(i,:)',X(i,:)') ; 
    end
    rsq = mean(r.^2) ; % average over all electrodes

else
    error('f_variance_analysis requires at least 3 inputs')
end