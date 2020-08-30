function dsatlas = f_reduce_atlas(atlas)
% Useage: 
% dsatlas = F_REDUCE_ATLAS(atlas)
% 
% Inputs: 
% atlas: An atlas in the form of a Brainstorm output. It is a 1x360
% structure, where each entry atlas(i) corresponds to an ROI in the HCP
% atlas. Each element should contain: 
% - atlas(i).Vertices: a list of vertices in ROI i
% - atlas(i).Label: the name of ROI i (in format 'L_X_ROI L', where L (or
%     swap for R) is hemisphere and X is the ROI label, e.g. L_V1_ROI L is
%     the left V1. 
% 
% Outputs: 
% dsatlas: This atlas downsampled to 230 ROIs. 

rois = template_rois() ; % get the correct order for the ROIs
atlas = reorder_rois(atlas,rois) ; % ensure atlas is in the same order as the template
dsatlas = downsample_atlas(atlas) ; 


function rois = template_rois()

% Parcellation 
rois = {{'V1'} ; 
        {'V2','V3','V4'} ; 
        {'V6','V3A','V7','IPS1','V3B','V6A'} ; 
        {'V8','FFC','PIT','VMV1','VMV3','VMV2','VVC'} ; 
        {'MST','LO1','LO2','MT','PH','V4t','FST','V3CD','LO3'} ;
        {'4','3b','1','2','3a'} ; 
        {'5m','5mv','5L','24dd','24dv','SCEF','6ma','6mp'} ; 
        {'FEF','PEF','55b','6d','6v','6r','6a'} ; 
        {'43','OP4','OP1','OP2-3','FOP1','PFcm'} ;
        {'A1','RI','PBelt','MBelt','LBelt'} ; 
        {'TA2','STGa','A5','STSda','STSdp','STSvp','A4','STSva'} ;
        {'52','PoI2','FOP4','MI','Pir','AVI','AAIC','FOP3','FOP2','PoI1','Ig','FOP5','PI'} ; 
        {'EC','PreS','PeEc','PHA1','PHA2','PHA3'} ;
        {'TF','TGd','TE1a','TE1p','TE2a','TE2p','PHT','TGv','TE1m'} ; 
        {'PSL','STV','TPOJ1','TPOJ2','TPOJ3'} ; 
        {'7Pm','7AL','7Am','7PL','7PC','LIPv','VIP','MIP','LIPd','AIP'} ; 
        {'PFt','PGp','IP2','IP1','IP0','PFop','PF','PFm','PGi','PGs'} ; 
        {'23c','RSC','POS2','PCV','7m','POS1','23d','v23ab','d23ab','31pv','ProS','DVT','31pd','31a'} ;  
        {'p24pr','33pr','a24pr','p32pr','a24','d32','8BM','p32','10r','9m','10v','25','s32','a32pr','p24'} ; 
        {'pOFC','47m','10d','a10p','10pp','11l','13l','OFC','47s','p10p','a47r'} ; 
        {'44','45','47l','IFJa','IFJp','IFSp','IFSa','p47r'} ; 
        {'SFL','8Av','8Ad','8BL','9p','8C','p9-46v','46','a9-46v','9-46d','9a','i6-8','s6-8'}
        } ; 
    
end

function atlas = reorder_rois(atlas,rois)
    ind = nan(length(atlas),1) ; 
    n = 0 ; 
    for i = 1:length(atlas)
        lbls{i} = atlas(i).Label ; 
    end
    
    default_lbl = [1;2] ; 
    
    % Left hemisphere
    for i = 1:length(rois)
        for j = 1:length(rois{i})
            n = n+1 ; 
            lbl = {['L_' rois{i}{j} '_ROI L'],['L_' rois{i}{j} '_ROI']} ; 
            ind0 = find(strcmp(lbls,lbl{default_lbl(1)})) ; 
            if isempty(ind0)
                ind0 = find(strcmp(lbls,lbl{default_lbl(2)})) ; 
                default_lbl = flipud(default_lbl) ; 
            end
            ind(n) = ind0 ; 
            atlas(ind(n)).Cluster = i ; 
        end
    end
    % Right hemisphers
    for i = 1:length(rois)
        for j = 1:length(rois{i})
            n = n+1 ; 
            lbl = {['R_' rois{i}{j} '_ROI R'],['R_' rois{i}{j} '_ROI']} ;
            ind0 = find(strcmp(lbls,lbl{default_lbl(1)})) ; 
            if isempty(ind0)
                ind0 = find(strcmp(lbls,lbl{default_lbl(2)})) ; 
                default_lbl = flipud(default_lbl) ; 
            end
            ind(n) = ind0 ; 
            atlas(ind(n)).Cluster = i ; 
        end
    end
    ind(isnan(ind)) = [] ; 
    atlas = atlas(ind) ;     
end 

function dsatlas = downsample_atlas(atlas) ; 

    dsatlas = atlas ; 
    if strcmp(atlas(1).Label(end-1:end),' L')
        suffix = {' L',' R'} ; 
    else
        suffix = {'',''} ; 
    end
    
    % Rois which need to be merged
    merg = {{'V3B','V7'} ; 
            {'V6','V6A'} ; 
            {'V8','VVC'} ; 
            {'VMV1','VMV2','VMV3'} ; 
            {'FFC','PIT'} ; 
            {'MST','MT','V4t'} ; 
            {'LO1','LO2','LO3'} ; 
            {'24dd','24dv'} ; 
            {'5m','5mv'} ; 
            {'PEF','6v'} ; 
            {'43','FOP1'} ; 
            {'OP1','OP2-3','PFcm'} ; 
            {'A1','LBelt','MBelt'} ; 
            {'PBelt','RI'} ; 
            {'STGa','TA2'} ; 
            {'STSda','STSdp'} ; 
            {'STSva','STSvp'} ; 
            {'Ig','52','PoI1','PoI2','PI','MI'} ; 
            {'FOP2','FOP3','FOP4','FOP5'} ; 
            {'AAIC','AVI','Pir'} ; 
            {'PHA1','PHA2','PHA3','PreS'} ; 
            {'EC','PeEc'} ; 
            {'TPOJ2','TPOJ3'} ; 
            {'7PL','7Pm'} ; 
            {'LIPd','LIPv'} ; 
            {'DVT','ProS','POS1'} ; 
            {'d23ab','v23ab','23d','RSC'} ; 
            {'31a','31pd','31pv','7m'} ; 
            {'23c','PCV'} ; 
            {'a24pr','p24pr','33pr'} ; 
            {'a32pr','p32pr','p24'} ; 
            {'8BM','d32'} ; 
            {'p32','s32','a24','25'} ; 
            {'10r','10v'} ; 
            {'47m','47s'} ; 
            {'11l','13l'} ; 
            {'a10p','p10p'} ; 
            {'10d','10pp'} ; 
            {'IFJa','IFJp','IFSp'} ; 
            {'47l','p47r'}} ;

    % names of rois
    clear names
    for i = 1:length(atlas)
        names{i} = atlas(i).Label ; 
    end
    % Numerical rois
    for i = 1:length(merg)
        for k = 1:length(merg{i,1})
            % Left hemisphere
            idx = find(strcmp(sprintf('L_%s_ROI%s',merg{i}{k},suffix{1}),names)) ; 
            if isempty(idx)
                error(sprintf('Cannot find %s L',merg{i}{k}))
            end
            merg{i,2}{k} = idx ;
            
            idx = find(strcmp(sprintf('R_%s_ROI%s',merg{i}{k},suffix{2}),names)) ; 
            if isempty(idx)
                error(sprintf('Cannot find %s R',merg{i}{k}))
            end
            merg{i,3}{k} = idx ;
        end
    end

    % Next, parcellate
    keep = true(length(dsatlas),1) ; 
    for i = 1:length(merg)
        for hemi = 2:3 % left/right hemisphere
            idx = [] ; 
            lbl = [];  
            clst = [] ; 
            for k = 1:length(merg{i,hemi})
                idx = [idx;atlas(merg{i,hemi}{k}).Vertices(:)] ; 
                lbl = [lbl,atlas(merg{i,hemi}{k}).Label,' + '] ; 
                clst = [clst,atlas(merg{i,hemi}{k}).Cluster] ; 
            end
            clst = unique(clst) ; 
            lbl = lbl(1:end-3) ; 
            dsatlas(merg{i,hemi}{1}).Vertices = idx ; 
            dsatlas(merg{i,hemi}{1}).Label = lbl ;
            dsatlas(merg{i,hemi}{1}).Cluster = clst ; 
            for k = 2:length(merg{i,hemi})
                keep(merg{i,hemi}{k}) = false ; 
            end
        end
    end
    dsatlas = dsatlas(keep) ; 
end


% 
% function dsatlas = downsample_atlas(atlas) ; 
% 
%     dsatlas = atlas ; 
%     
%     % Parcellation      
%     rois = {{'V3B','V7'} ; 
%             {'V6','V6A'} ; 
%             {'V8','VVC'} ; 
%             {'VMV1','VMV2','VMV3'} ; 
%             {'FFC','PIT'} ; 
%             {'MST','MT','V4t'} ; 
%             {'LO1','LO2','LO3'} ; 
%             {'24dd','24dv'} ; 
%             {'5m','5mv'} ; 
%             {'PEF','6v'} ; 
%             {'43','FOP1'} ; 
%             {'OP1','OP2-3','PFcm'} ; 
%             {'A1','LBelt','MBelt'} ; 
%             {'PBelt','RI'} ; 
%             {'STGa','TA2'} ; 
%             {'STSda','STSdp'} ; 
%             {'STSva','STSvp'} ; 
%             {'Ig','52','PoI1','PoI2','PI','MI'} ; 
%             {'FOP2','FOP3','FOP4','FOP5'} ; 
%             {'AAIC','AVI','Pir'} ; 
%             {'PHA1','PHA2','PHA3','PreS'} ; 
%             {'EC','PeEc'} ; 
%             {'TPOJ2','TPOJ3'} ; 
%             {'7PL','7Pm'} ; 
%             {'LIPd','LIPv'} ; 
%             {'DVT','ProS','POS1'} ; 
%             {'d23ab','v23ab','23d','RSC'} ; 
%             {'31a','31pd','31pv','7m'} ; 
%             {'23c','PCV'} ; 
%             {'a24pr','p24pr','33pr'} ; 
%             {'a32pr','p32pr','p24'} ; 
%             {'8BM','d32'} ; 
%             {'p32','s32','a24','25'} ; 
%             {'10r','10v'} ; 
%             {'47m','47s'} ; 
%             {'11l','13l'} ; 
%             {'a10p','p10p'} ; 
%             {'10d','10pp'} ; 
%             {'IFJa','IFJp','IFSp'} ; 
%             {'47l','p47r'}} ; 
% 
%     % names of rois
%     clear names
%     for i = 1:length(atlas)
%         names{i} = atlas(i).Label ; 
%     end
%     % Numerical rois
%     for i = 1:length(rois)
%         for k = 1:length(rois{i,1})
%             idx = find(strcmp(sprintf('R_%s_ROI',rois{i}{k}),names)) ; 
%             if isempty(idx)
%                 error(sprintf('Cannot find %s R',rois{i}{k}))
%             end
%             rois{i,2}{k} = idx ;
%         end
%     end
% 
%     % Next, parcellate
%     keep = true(length(dsatlas),1) ; 
%     for i = 1:length(rois)
%         for hemi = 2% :3 % left/right hemisphere
%             idx = [] ; 
%             lbl = [];  
%             clst = [] ; 
%             for k = 1:length(rois{i,hemi})
%                 idx = [idx,atlas(rois{i,hemi}{k}).Vertices] ; 
%                 lbl = [lbl,atlas(rois{i,hemi}{k}).Label,' + '] ; 
%                 clst = [clst,atlas(rois{i,hemi}{k}).Cluster] ; 
%             end
%             clst = unique(clst) ; 
%             if length(clst) > 1
%                 error('what????')
%             end
%             lbl = lbl(1:end-3) ; 
%             dsatlas(rois{i,hemi}{1}).Vertices = idx ; 
%             dsatlas(rois{i,hemi}{1}).Label = lbl ;
%             dsatlas(rois{i,hemi}{1}).Cluster = clst ; 
%             for k = 2:length(rois{i,hemi})
%                 keep(rois{i,hemi}{k}) = false ; 
%             end
%         end
%     end
%     dsatlas = dsatlas(keep) ; 
% end




end