function glacfg = paramInitGLADpar(cfg)
% Created by Xinan Chen on 04/26/2021.
% This function returns parameters corresponding to the input 'data_name'
%%

glacfg.do_masked = 1;%1 if ensure everything is inside mask before saving to .nii
glacfg.spMSK_str = 'ROImsk';
%% set parameters, start points and directory
switch cfg.dataset_name
    case 'CAA'
        glacfg.do_sp = 1; %1 if use sp from max(data)>sp_thresh, 0 o/w
        glacfg.sp_thresh = 12;

        glacfg.mdt = 1;%3;%5;%10;
        glacfg.sl_tol = 3;%2; %threshold for minimum Euclidean length between initial and final streamline points                
        glacfg.spType = 'ordered'; %'uniform' (default: 'ordered'); %'ordered': every fs-th possible sp is used (default)sp distribution; 'uniform': sp chosen according to uniform distribution
        switch glacfg.spType
            case 'ordered'
                glacfg.spfs = 1;
                glacfg.spSTR = sprintf('spDISTordered_fs%d',glacfg.spfs);
                glacfg.spSTRtitle = sprintf('spDISTordered-fs%d',glacfg.spfs);
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
        glacfg.sfs = 5;
        glacfg.vfs = 1;

        %glacfg.do_sp = 0; %1 if use sp from max(dpsnrv)>sp_thresh, 0 o/w
        %glacfg.sp_thresh = 12;


        % NCA parameters
        glacfg.radius = 2; % redius to dialate
        glacfg.NNnum_tol = 20; 
        glacfg.stdcut = 0.5; 
        glacfg.meancut = [0.5,1]; 
        glacfg.Npcut = 1; 
        glacfg.Avepathlcut = 8; 
        glacfg.threshstr = sprintf('NNnum_tol_%d_||_stdcut_%.2f_wmeancut_[%.2f,%.2f]_Npcut_%d_Avepathlcut_%.1f',glacfg.NNnum_tol,glacfg.stdcut,glacfg.meancut(1),glacfg.meancut(2),glacfg.Npcut,glacfg.Avepathlcut);

        glacfg.maskdilate = 3;
        glacfg.radius2 = 3; % redius to dialate
        glacfg.NNnum_tol2 = 10; 
        glacfg.stdcut2 = 0.4; 
        glacfg.meancut2 = [0.6,1]; 
        glacfg.Npcut2 = 1; 
        glacfg.Avepathlcut2 = 10; 
        glacfg.AveMaxpathscut2 = 0.02; 
        glacfg.pathlcut2 = 20;%10;%20;
        glacfg.threshstr2 = sprintf('NNnum_tol_%d_||_stdcut_%.2f_wmeancut_[%.2f,%.2f]_Npcut_%d_Avepathlcut_%.1f|pathlcut_%.1f',glacfg.NNnum_tol2,glacfg.stdcut2,glacfg.meancut2(1),glacfg.meancut2(2),glacfg.Npcut2,glacfg.Avepathlcut2,glacfg.pathlcut2);


        % vis parameters
        glacfg.jp = 5;       
                

    otherwise
        %fprintf('paramInitGLADpar: non-applicable dataset_name!');
        %return
        glacfg.do_sp = 1;%0; %1 if use sp from max(data)>sp_thresh, 0 o/w
        glacfg.sp_thresh = 0.2;%1;%0;

        glacfg.mdt = 1;
        glacfg.sl_tol = 3; %threshold for minimum Euclidean length between initial and final streamline points                
        glacfg.spType = 'ordered'; %'uniform' (default: 'ordered'); %'ordered': every fs-th possible sp is used (default)sp distribution; 'uniform': sp chosen according to uniform distribution
        switch glacfg.spType
            case 'ordered'
                glacfg.spfs = 1;
                glacfg.spSTR = sprintf('spDISTordered_fs%d',glacfg.spfs);
                glacfg.spSTRtitle = sprintf('spDISTordered-fs%d',glacfg.spfs);
            case 'uniform'
                glacfg.spPerc = 40;
                glacfg.spSTR = sprintf('spDISTunif_spPerc%d',glacfg.spPerc);
                glacfg.spSTRtitle = sprintf('spDISTunif-spPerc%d',glacfg.spPerc);
        end
        glacfg.nEulStep = 1; %number of Eulerian steps to be taken
        
        glacfg.smoothv = 0; % smooth velocity field
        if glacfg.smoothv
            glacfg.Svt = 5;%10; % smoothing w.r.t time
            glacfg.Svs = 1;%5; % smoothing w.r.t space
            glacfg.smoothvSTR = sprintf('1_Tolt%d_Tols%d',glacfg.Svt,glacfg.Svs);
        else
            glacfg.smoothvSTR = '0';
        end
        
        glacfg.smoothp = 0;%1; % smooth pathline
        if glacfg.smoothp
            glacfg.smpTol = 100; % the higher, the smoother the pathlines
            glacfg.smoothpSTR = sprintf('1_Tol%d',glacfg.smpTol);
        else
            glacfg.smoothpSTR = '0';
        end
        
        glacfg.sfs = 5;
        glacfg.vfs = 1;


        % NCA parameters
        glacfg.radius = 2; % redius to dialate
        glacfg.NNnum_tol = 20; 
        glacfg.stdcut = 0.5; 
        glacfg.meancut = [0.5,1]; 
        glacfg.Npcut = 1; 
        glacfg.Avepathlcut = 8; 
        glacfg.threshstr = sprintf('NNnum_tol_%d_||_stdcut_%.2f_wmeancut_[%.2f,%.2f]_Npcut_%d_Avepathlcut_%.1f',glacfg.NNnum_tol,glacfg.stdcut,glacfg.meancut(1),glacfg.meancut(2),glacfg.Npcut,glacfg.Avepathlcut);

        glacfg.maskdilate = 3;
        glacfg.radius2 = 3; % redius to dialate
        glacfg.NNnum_tol2 = 10; 
        glacfg.stdcut2 = 0.4; 
        glacfg.meancut2 = [0.6,1]; 
        glacfg.Npcut2 = 1; 
        glacfg.Avepathlcut2 = 10; 
        glacfg.pathlcut2 = 20;%10;%20;
        glacfg.threshstr2 = sprintf('NNnum_tol_%d_||_stdcut_%.2f_wmeancut_[%.2f,%.2f]_Npcut_%d_Avepathlcut_%.1f|pathlcut_%.1f',glacfg.NNnum_tol2,glacfg.stdcut2,glacfg.meancut2(1),glacfg.meancut2(2),glacfg.Npcut2,glacfg.Avepathlcut2,glacfg.pathlcut2);


        % vis parameters
        glacfg.jp = 5;       

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
