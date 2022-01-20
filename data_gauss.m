% This functions contains the creation of a sample dataset named 'gauss2' which is from a Gaussian sphere dataset.

addpath('./Inverse','./Sensitivities','./analyzeFlows',genpath('./utilities'))

%% set directories and parameters

clear cfg

cfg.tag                     = 'gauss2'; % name of this case
cfg.dataset_name            = 'Gaussian';  % name of dataset


% optional
%cfg.anato                   = './data/12_MONTH_DATA/psnrv/WT/C294/pbase_snrv_C294_031318A_DOTA37_30ul_E53.nii'; %
%cfg.sp_mask_opts(1).name    = 'brain';
%cfg.sp_mask_opts(1).path    = '../../data/12_MONTH_DATA/12months_mask_brainCSF/C294.nii';
cfg.do_ROI_msk              = 0;

% set rOMT parameters
cfg.do_resize               = 0;
cfg.size_factor             = 1;
cfg.smooth                  = 0;
cfg.reinitR                 = 1;%1;%0; %0 if do consecutively and 1 if reinitialize rho
cfg.dilate                  = 0;

cfg.first_time              = 1;
cfg.time_jump               = 1;
cfg.last_time               = 4;

cfg.sigma                   = 0.002;%0;%0.2;%0;%2e-3;
cfg.dt                      = 0.4;
cfg.nt                      = 10;
cfg.gamma                   = 0.008;
cfg.beta                    = 0.0001;
cfg.niter_pcg               = 60;%20;
cfg.dTri                    = 1;%1 := 'closed', 3:= 'open' for boundary condition
cfg.add_source              = 0;
cfg.reInitializeU           = 1;
%% create gaussian sphere
a = 6; mv = 0.8;
Nx = 50; Ny = Nx; Nz = Nx;

cfg.x_range                 = 1:Nx;
cfg.y_range                 = 1:Ny;
cfg.z_range                 = 1:Nz;

% load ROI
if cfg.do_ROI_msk
    tmp = nii2mat(cfg.ROI_msk_path,cfg.x_range,cfg.y_range,cfg.z_range);
    cfg.msk = zeros(size(tmp));
    cfg.msk(tmp>0) = 1; 
else
    cfg.msk = ones(length(cfg.x_range),length(cfg.y_range),length(cfg.z_range));
end
if cfg.do_resize
   cfg.msk = resizeMatrix(cfg.msk,round(cfg.size_factor.*size(cfg.msk)),'linear');
   cfg.msk(cfg.msk~=1) = 0;
end

x = linspace(-a,a,Ny) - 2*mv;
y = linspace(-a,a,Nx) - 2*mv;
z = linspace(-a,a,Nz) - 2*mv;
for i = 1:(cfg.last_time-cfg.first_time)/cfg.time_jump+2
    [X,Y,Z] = meshgrid(x,y,z);
    tmp = (100/sqrt(2*pi).*exp(-(X.^2/2)-(Y.^2/2)-(Z.^2/2)));
    if cfg.do_resize
       tmp = resizeMatrix(tmp,round(cfg.size_factor.*size(tmp)),'linear');
    end
    if cfg.smooth>0
    tmp = affine_diffusion_3d(tmp,cfg.smooth,0.1,1,1);
    end
    
    % add diffusion
    if i~=1
    tmp = imgaussfilt3(tmp,sqrt(0.2)*i);
    end
    
    %figure, montageArray(tmp), axis image, caxis([0,50]), title(sprintf('rho%d',i))
    tmp(~cfg.msk) = 0;
    cfg.vol(i).data = tmp;
    x = x + mv; y = y + mv; z = z + mv;
    %fprintf('i = %d, sum = %.5f\n', i, sum(tmp(:)))
end
%%

cfg.domain_size             = size(cfg.vol(1).data);
cfg.sig_str                 = erase(num2str(cfg.sigma,'%.0e'),{'-0','+0'});
cfg.true_size               = round(cfg.size_factor*[length(cfg.x_range),length(cfg.y_range),length(cfg.z_range)]);
cfg.version                 = sprintf('diff_%s_tj_%d_dt_%2.1f_nt_%d_ti_%d_tf_%d_uini_0_beta_%5.4f_R_gamma_%4.3f_dtri%d_rsmooth%d_rreinit%d_source%d_pcg%d',...
                                cfg.sig_str,cfg.time_jump,cfg.dt,cfg.nt,cfg.first_time,cfg.last_time,cfg.beta,cfg.gamma,cfg.dTri,cfg.smooth,cfg.reinitR,cfg.add_source,cfg.niter_pcg);
cfg.out_dir                 = sprintf('./test_results/%s/%s',cfg.tag,cfg.version);

