% This functions contains the rOMT algorithm and post-processing of the
% results on a sample dataset named 'C294' which is from a CAA dataset containsing about 48 rat brain cases w/w.o CAA disease.

addpath('./Inverse','./Sensitivities','./analyzeFlows',genpath('./utilities'))

%% set directories and parameters

clear cfg

cfg.tag                     = 'C294'; % name of this case
cfg.dataset_name            = 'CAA';  % name of dataset
% set directories
cfg.data_dir                = './data/12_MONTH_DATA/psnrv/WT/C294/C294_031318A/psnrv_C294_031318A_DOTA37_30ul_E'; % data dir
cfg.extension               = '.nii'; % format of data
cfg.do_ROI_msk              = 1;
cfg.ROI_msk_path            = '../../data/12_MONTH_DATA/masksforOMT/WT/C294_MASK.nii';
cfg.x_range                 = 20:80;
cfg.y_range                 = 1:106;
cfg.z_range                 = 39:85;

% optional
%cfg.max_dpsnrv              = './data/12_MONTH_DATA/MAXpsnrv/C294_031318A_psnrv_max.nii'; %
cfg.anato                   = './data/12_MONTH_DATA/psnrv/WT/C294/pbase_snrv_C294_031318A_DOTA37_30ul_E53.nii'; %
cfg.sp_mask_opts(1).name    = 'brain';
cfg.sp_mask_opts(1).path    = '../../data/12_MONTH_DATA/12months_mask_brainCSF/C294.nii';

% set rOMT parameters
cfg.do_resize               = 0;%1;
cfg.size_factor             = 1;%0.5;
%cfg.data_index_E            = 19:53;
cfg.smooth                  = 1;
cfg.reinitR                 = 1;%1;%0; %0 if do consecutively and 1 if reinitialize rho
cfg.dilate                  = 3;

cfg.first_time              = 31;%cfg.data_index_E(13);
cfg.time_jump               = 2;%3;
cfg.last_time               = 51;%cfg.data_index_E(33);%;cfg.data_index_E(31);

cfg.sigma                   = 2e-3;
cfg.dt                      = 0.4;
cfg.nt                      = 10;
cfg.gamma                   = 0.008;
cfg.beta                    = 0.0001;
cfg.niter_pcg               = 60;%20;
cfg.dTri                    = 1;%1 := 'closed', 3:= 'open' for boundary condition
cfg.add_source              = 0;
%
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
if cfg.dilate>0
    [xr,yr,zr] = meshgrid(-cfg.dilate:cfg.dilate,-cfg.dilate:cfg.dilate,-cfg.dilate:cfg.dilate);
    strel = (xr/cfg.dilate).^2 + (yr/cfg.dilate).^2 + (zr/cfg.dilate).^2 <= 1;
    cfg.msk = imdilate(cfg.msk,strel);
end

% load vol
for i = 1:(cfg.last_time-cfg.first_time)/cfg.time_jump+2
    tmp = nii2mat(sprintf('%s%02d%s',cfg.data_dir,cfg.first_time+(i-1)*cfg.time_jump,cfg.extension),cfg.x_range,cfg.y_range,cfg.z_range);
    if cfg.do_resize
       tmp = resizeMatrix(tmp,round(cfg.size_factor.*size(tmp)),'linear');
    end
    if cfg.smooth>0
    tmp = affine_diffusion_3d(tmp,cfg.smooth,0.1,1,1);
    end
    tmp(~cfg.msk) = 0;
    cfg.vol(i).data = tmp;
end

cfg.domain_size             = size(cfg.vol(1).data);
cfg.sig_str                 = erase(num2str(cfg.sigma,'%.0e'),'-0');
cfg.true_size               = round(cfg.size_factor*[length(cfg.x_range),length(cfg.y_range),length(cfg.z_range)]);
cfg.version                 = sprintf('diff_%s_tj_%d_dt_%2.1f_nt_%d_ti_%d_tf_%d_uini_0_beta_%5.4f_R_gamma_%4.3f_dtri%d_rsmooth%d_rreinit%d_source%d_dilate%d_pcg%d',...
                                cfg.sig_str,cfg.time_jump,cfg.dt,cfg.nt,cfg.first_time,cfg.last_time,cfg.beta,cfg.gamma,cfg.dTri,cfg.smooth,cfg.reinitR,cfg.add_source,cfg.dilate,cfg.niter_pcg);
cfg.out_dir                 = sprintf('./test_results/%s/%s',cfg.tag,cfg.version);

%% Run rOMT
[cfg, flag] = runROMT(cfg);
%runROMT_par(cfg);

%% Run post-processing

%[cfg, map, SL, stream, PATH] = runGLAD(cfg);
[cfg, s, SL, PATH] = runGLAD2(cfg);


%% Visualization of speed map
x = 1:cfg.true_size(1);
y = 1:cfg.true_size(2);
z = 1:cfg.true_size(3);
figure,
hs=slice(y,x,z,s,x,y,z); 
set(hs,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.04);
alpha('color'),alphamap(linspace(0,1,100))
title(sprintf('Test: tag = %s, Speed Map',cfg.tag),'Fontsize',20)
grid off, box off, axis image
xlabel('x-axis','FontSize',20),ylabel('y-axis','FontSize',20),zlabel('z-axis','FontSize',20)
colormap(jet)
caxis([0,0.5])
view([-188.3500   13.7463])
set(gca,'Color',[0.85,0.85,0.93]), set(gcf,'unit','normalized','position',[0.1,1,0.4,0.5],'Color',[0.85,0.85,0.93]), set(gcf, 'InvertHardcopy', 'off')
colorbar, grid on,

saveas(gcf, sprintf('%s/%s/%s_LagSpeed_E%02d_%02d.png',cfg.out_dir,cfg.outdir_s,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 
%% Visualization of flux vectors
figure,
strid = 10;
magnify = .5;%1;
[x, y, z] = meshgrid(1:cfg.true_size(2), 1:cfg.true_size(1), 1:cfg.true_size(3));
mskfv = isosurface(x,y,z,cfg.msk,0.5);
mskp = patch(mskfv);
mskp.FaceColor = [.17,.17,.17];
mskp.FaceAlpha= 0.031;
mskp.EdgeColor = [.17,.17,.17];
mskp.EdgeAlpha= 0;
mskp.DisplayName = 'mask';

grid on, axis image
hold on,
q = quiver3(PATH.startp(PATH.ind_brain(1:strid:end),2),PATH.startp(PATH.ind_brain(1:strid:end),1),PATH.startp(PATH.ind_brain(1:strid:end),3),PATH.disp(PATH.ind_brain(1:strid:end),2)*magnify,PATH.disp(PATH.ind_brain(1:strid:end),1)*magnify,PATH.disp(PATH.ind_brain(1:strid:end),3)*magnify,'color','r','LineWidth',2,'MaxHeadSize',0.3,'AutoScale','off','DisplayName','flux vectors');
title(sprintf('%s: velocity flux vectors',cfg.tag),'FontSize',20, 'Interpreter', 'none'), set(gcf, 'Position', [376 49 1256 719])
xlabel('x-axis','FontSize',20),ylabel('y-axis','FontSize',20),zlabel('z-axis','FontSize',20)
%// Compute the magnitude of the vectors
mags = sqrt(sum(cat(2, q.UData(:), q.VData(:), ...
            reshape(q.WData, numel(q.UData), [])).^2, 2));

%// Get the current colormap
currentColormap = colormap(jet);
[~, ~, ind] = histcounts(mags, size(currentColormap, 1));
cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
cmap(:,:,4) = 255;
cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);
set(q.Head, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:3,:,:), [], 4).');   
set(q.Tail, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:2,:,:), [], 4).');
legend('Location','best','Fontsize',12); colorbar;
view([-188.3500   13.7463]);
set(gca,'Color',[0.85,0.85,0.93]), set(gcf,'unit','normalized','position',[0.1,1,0.4,0.5],'Color',[0.85,0.85,0.93]), set(gcf, 'InvertHardcopy', 'off')

saveas(gcf, sprintf('%s/%s/%s_LagFluxVector_E%02d_%02d.png',cfg.out_dir,cfg.outdir_v,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 

%{
%% Visualization of advective & diffusive vectors
[InDadv,~,~] = intersect(PATH.ADVind,PATH.ind_brain);
[InDdiff,~,~] = intersect(PATH.DIFFind,PATH.ind_brain);

figure,
strid = 10;
magnify = .5;%1;
[x, y, z] = meshgrid(1:cfg.true_size(2), 1:cfg.true_size(1), 1:cfg.true_size(3));
mskfv = isosurface(x,y,z,cfg.msk,0.5);
mskp = patch(mskfv);
mskp.FaceColor = [.17,.17,.17];
mskp.FaceAlpha= 0.031;
mskp.EdgeColor = [.17,.17,.17];
mskp.EdgeAlpha= 0;
mskp.DisplayName = 'mask';

grid on, axis image
hold on,
q = quiver3(PATH.startp(InDadv(1:strid:end),2),PATH.startp(InDadv(1:strid:end),1),PATH.startp(InDadv(1:strid:end),3),PATH.disp(InDadv(1:strid:end),2)*magnify,PATH.disp(InDadv(1:strid:end),1)*magnify,PATH.disp(InDadv(1:strid:end),3)*magnify,'color','m','LineWidth',2,'MaxHeadSize',0.3,'AutoScale','off','DisplayName','Advecitve Vectors');
title(sprintf('%s: Advective and Diffusive flux vectors',cfg.tag),'FontSize',20, 'Interpreter', 'none'), set(gcf, 'Position', [376 49 1256 719])
xlabel('x-axis','FontSize',20),ylabel('y-axis','FontSize',20),zlabel('z-axis','FontSize',20)

hold on,
l = quiver3(PATH.startp(InDdiff(1:strid:end),2),PATH.startp(InDdiff(1:strid:end),1),PATH.startp(InDdiff(1:strid:end),3),PATH.disp(InDdiff(1:strid:end),2)*magnify,PATH.disp(InDdiff(1:strid:end),1)*magnify,PATH.disp(InDdiff(1:strid:end),3)*magnify,'color','g','LineWidth',2,'MaxHeadSize',0.3,'AutoScale','off','DisplayName','Diffusive Vectors');

view([-188.3500   13.7463]); legend('Location','best','Fontsize',12);
set(gca,'Color',[0.85,0.85,0.93]), set(gcf,'unit','normalized','position',[0.1,1,0.4,0.5],'Color',[0.85,0.85,0.93]), set(gcf, 'InvertHardcopy', 'off')

saveas(gcf, sprintf('%s/%s/%s_LagAdvDiffVector_E%02d_%02d.png',cfg.out_dir,cfg.outdir_v,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 
%}
%% Visualization of pathlines
figure,
SL2 = SL(PATH.ind_brain);
nSL = length(SL2);
%colors = jet(round(max(PATH.displen)));

for ind = 1:10:nSL
    SL_tmp = SL2{ind};
    colors = jet(size(SL_tmp,1));
    hlines = patch([SL_tmp(:,2);NaN],[SL_tmp(:,1);NaN],[SL_tmp(:,3);NaN],[1,1,1]);
    %set(hlines,'EdgeColor',colors(round(PATH.displen(ind)),:));
    set(hlines,'FaceVertexCData',[colors;colors(end,:)],'EdgeColor','flat','FaceColor','none');
end
%
hold on;
[x, y, z] = meshgrid(1:cfg.true_size(2), 1:cfg.true_size(1), 1:cfg.true_size(3));
mskfv = isosurface(x,y,z,cfg.msk,0.5);
mskp = patch(mskfv);
mskp.FaceColor = [.37,.37,.37];
mskp.FaceAlpha= 0.1;
mskp.EdgeColor = [.37,.37,.37];
mskp.EdgeAlpha= 0;
view([-188.3500   13.7463])
axis image
set(gca,'Color',[0.85,0.85,0.93]), set(gcf,'unit','normalized','position',[0.1,1,0.4,0.5],'Color',[0.85,0.85,0.93]), set(gcf, 'InvertHardcopy', 'off')
colormap('jet'); colorbar, grid on,
title(sprintf('%s: Pathlines Starting with brain',cfg.tag), 'Interpreter', 'none', 'FontSize',18);
saveas(gcf, sprintf('%s/%s/%s_LagPathlines_E%02d_%02d.png',cfg.out_dir,cfg.outdir_v,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 









