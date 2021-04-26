% This functions contains the rOMT algorithm and post-processing of the
% results on a sample dataset named 'C294' which is from a CAA dataset containsing about 48 rat brain cases w/w.o CAA disease.

addpath('./Inverse','./Sensitivities','./analyzeFlows',genpath('./utilities'))

%% set directories and parameters

clear cfg

cfg.tag                     = 'C294'; % name of this case
cfg.dataset_name            = 'CAA';  % name of dataset
% set directories
cfg.data_dir                = './data/12_MONTH_DATA/psnrv/WT/C294/C294_031318A/psnrv_C294_031318A_DOTA37_30ul_E';
cfg.extension               = '.nii'; % format of data
cfg.max_dpsnrv              = './data/12_MONTH_DATA/MAXpsnrv/C294_031318A_psnrv_max.nii';
cfg.anato                   = './data/12_MONTH_DATA/psnrv/WT/C294/pbase_snrv_C294_031318A_DOTA37_30ul_E53.nii';
cfg.sp_mask_opts(1).name    = 'brain';
cfg.sp_mask_opts(1).path    = '../../data/12_MONTH_DATA/12months_mask_brainCSF/C294.nii';
cfg.data_mask_path          = '../../data/12_MONTH_DATA/masksforOMT/WT/C294_MASK.nii';
cfg.domain_size             = [100,106,100];
cfg.x_range                 = 20:80;
cfg.y_range                 = 1:106;
cfg.z_range                 = 39:85;

% set rOMT parameters
cfg.do_resize               = 1;
cfg.size_factor             = 0.5;
cfg.data_index_E            = 19:53;
cfg.smooth                  = 1;
cfg.reinitR                 = 0; %0 if do consecutively and 1 if reinitialize rho
cfg.dilate                  = 3;

cfg.first_time              = cfg.data_index_E(13);
cfg.time_jump               = 3;
cfg.last_time               = cfg.data_index_E(31);

cfg.sigma                   = 2e-3;
cfg.sig_str                 = '2e3';
cfg.dt                      = 0.4;
cfg.nt                      = 10;
cfg.gamma                   = 0.008;
cfg.beta                    = 0.0001;
cfg.niter_pcg               = 20;

cfg.dTri                    = 1;%1 := 'closed', 3:= 'open' for boundary condition
cfg.add_source              = 0;


%% Run rOMT
runROMT(cfg);

%% Run post-processing
runGLAD(cfg);










