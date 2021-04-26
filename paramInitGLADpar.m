function glacfg = paramInitGLADpar(cfg,type)
% Created by Xinan Chen on 04/26/2021.
% This function returns parameters corresponding to the input 'data_name' and analysis type
%   Input:  type  := 's' if run for speed map; 'v' if run for flux vectors
%   Output: 'glad_par' := structure containing parameters for GLAD analysis

%%

glacfg.do_spMsk = 0; %1 if use alternative mask for selecting pathline start points
glacfg.spMsk_name = 'olf_nomain'; %name of mask to use for selecting start points

% find index of mask with spMsk_name
% get path of spMsk
% check if spMsk is the same as the msk, if so -> set do_spMsk = 0
if glacfg.do_spMsk && isfield(cfg,'sp_mask_opts')
    sp_msk_path_opts = {cfg.sp_mask_opts(:).name};
    glacfg.spMsk_ind = find(strcmp(glacfg.spMsk_name,sp_msk_path_opts));
    if isempty(glacfg.spMsk_ind)
        warning(fprintf('getLPPA.m: unable to find path coresponding to mask named %s\n>>> Default mask being used and ''do_spMsk'' set to 0\n',glacfg.spMsk_name));
        glacfg.do_spMsk = 0; 
        clear glacfg.spMsk_name glacfg.spMsk_ind
    end        
else
    glacfg.do_spMsk = 0; 
    clear glacfg.spMsk_name glacfg.spMsk_ind
end
if glacfg.do_spMsk
    glacfg.spMSK_str = sprintf('altSPmsk_%s_opt%d',glacfg.spMsk_name,glacfg.spMsk_ind);
else
    glacfg.spMSK_str = 'altSPmsk0';
end
    
glacfg.do_spErodeMsk = 0; %0 if don't erode SPmask to find start points, N if erode by N points
glacfg.do_spDilateMsk = 0; %0 if don't dilate SPmask for start points

if glacfg.do_spErodeMsk>0 && glacfg.do_spDilateMsk>0
    warning('getGLAcfg.m: spMask cannot be both eroded and dilated. spMask will be eroded (no dilation) by default');
    glacfg.do_spDilateMsk = 0;
end

glacfg.do_masked = 1;%1 if ensure everything is inside mask before saving to .nii

%% set parameters, start points and directory
switch cfg.dataset_name
    case 'CAA'
        switch type
            case 's'
                glacfg.do_sp = 1; %1 if use sp from max(dpsnrv)>sp_thresh, 0 o/w
                glacfg.sp_thresh = 12;
                
                glacfg.mdt = 10;
                glacfg.sl_tol = 2; %threshold for minimum Euclidean length between initial and final streamline points                
                glacfg.spType = 'ordered'; %'uniform' (default: 'ordered'); %'ordered': every fs-th possible sp is used (default)sp distribution; 'uniform': sp chosen according to uniform distribution
                switch glacfg.spType
                    case 'ordered'
                        glacfg.fs = 5;
                        glacfg.spSTR = sprintf('spDISTordered_fs%d',glacfg.fs);
                        glacfg.spSTRtitle = sprintf('spDISTordered-fs%d',glacfg.fs);
                    case 'uniform'
                        glacfg.spPerc = 40;
                        glacfg.spSTR = sprintf('spDISTunif_spPerc%d',glacfg.spPerc);
                        glacfg.spSTRtitle = sprintf('spDISTunif-spPerc%d',glacfg.spPerc);
                end
                glacfg.nEulStep = 1; %number of Eulerian steps to be taken

                % QuickBundle parameters 
                glacfg.do_centroid = 0;
                glacfg.qbthresh = 4;%5;%QuickBundles threshold distance
                glacfg.clus_tol = 12; %threshold for min # of sl's needed in a cluster
                glacfg.clus_cutoff = 600; %pick up to # of largest clusters
                glacfg.nbp = 124; %# of points that QuickBundles will resample streamline to
                glacfg.metric_str = 'AveragePointwiseEuclideanMetric';
                switch glacfg.metric_str
                    case 'AveragePointwiseEuclideanMetric'
                        glacfg.ms = 'AvgPwEuc';
                        glacfg.feature_str = sprintf('ResampleFeature(nb_points=%d)',glacfg.nbp);
                    case 'CosineMetric'
                        glacfg.ms = 'Cos';
                        glacfg.feature_str = 'VectorOfEndpointsFeature()';
                end
                
            case 'v'
                glacfg.do_sp = 0; %1 if use sp from max(dpsnrv)>sp_thresh, 0 o/w
                glacfg.sp_thresh = 12;
                
                glacfg.mdt = 10;
                glacfg.sl_tol = 2; %threshold for minimum Euclidean length between initial and final streamline points                
                glacfg.spType = 'ordered'; %'uniform' (default: 'ordered'); %'ordered': every fs-th possible sp is used (default)sp distribution; 'uniform': sp chosen according to uniform distribution
                switch glacfg.spType
                    case 'ordered'
                        glacfg.fs = 1;
                        glacfg.spSTR = sprintf('spDISTordered_fs%d',glacfg.fs);
                        glacfg.spSTRtitle = sprintf('spDISTordered-fs%d',glacfg.fs);
                    case 'uniform'
                        glacfg.spPerc = 40;
                        glacfg.spSTR = sprintf('spDISTunif_spPerc%d',glacfg.spPerc);
                        glacfg.spSTRtitle = sprintf('spDISTunif-spPerc%d',glacfg.spPerc);
                end
                glacfg.nEulStep = 1; %number of Eulerian steps to be taken
                
                % NCA parameters
                glacfg.radius = 2; % redius to dialate
                glacfg.NNnum_tol = 20; 
                glacfg.stdcut = 0.5; 
                glacfg.meancut = [0.5,1]; 
                glacfg.Npcut = 1; 
                glacfg.Avepathlcut = 8; 
                %glacfg.AveMaxpathscut = 0.02;
                glacfg.threshstr = sprintf('NNnum_tol_%d_||_stdcut_%.2f_wmeancut_[%.2f,%.2f]_Npcut_%d_Avepathlcut_%.1f',glacfg.NNnum_tol,glacfg.stdcut,glacfg.meancut(1),glacfg.meancut(2),glacfg.Npcut,glacfg.Avepathlcut);
                
                glacfg.maskdilate = 3;
                glacfg.radius2 = 3; % redius to dialate
                glacfg.NNnum_tol2 = 10; 
                glacfg.stdcut2 = 0.4; 
                glacfg.meancut2 = [0.6,1]; 
                glacfg.Npcut2 = 1; 
                glacfg.Avepathlcut2 = 10; 
                glacfg.AveMaxpathscut2 = 0.02; 
                glacfg.pathlcut2 = 10;%20;
                glacfg.threshstr2 = sprintf('NNnum_tol_%d_||_stdcut_%.2f_wmeancut_[%.2f,%.2f]_Npcut_%d_Avepathlcut_%.1f|pathlcut_%.1f',glacfg.NNnum_tol2,glacfg.stdcut2,glacfg.meancut2(1),glacfg.meancut2(2),glacfg.Npcut2,glacfg.Avepathlcut2,glacfg.pathlcut2);
                
                
                % vis parameters
                glacfg.jp = 5;
                
                
                
            otherwise
                fprintf('paramInitGLADpar: non-applicable type!');
                return
        end
    otherwise
        fprintf('paramInitGLADpar: non-applicable dataset_name!');
        return

end


glacfg.flw_type = 'vel';%'flw';
glacfg.pln = 2; %minimum number of unique points needed to be considered a pathline
glacfg.RD = 'R'; %'D' if use data; 'R' if use interpolated img
glacfg.XT = 'T'; %'X' if interp fixed spatial dist, 'T' if use time step;
glacfg.minIm0 = 0; %1 if set img = img - min(img), 0 otherwise

glacfg.cutoff_str = 'min';
switch glacfg.cutoff_str
    case 'min'
        glacfg.thresholds.conc = 0.0001;
        glacfg.thresholds.front = 0;
        glacfg.thresholds.flow = 0;
        glacfg.thresholds.speed = 0.0001;
    case 'max'
        glacfg.thresholds.conc = .0001;
        glacfg.thresholds.front = 0;
        glacfg.thresholds.flow = 0;
        glacfg.thresholds.speed = 0;
    case 'mean'
        glacfg.thresholds.conc = .00001;
        glacfg.thresholds.front = 0;
        glacfg.thresholds.flow = 0.0001;
        glacfg.thresholds.speed = 0.0001;
end



end