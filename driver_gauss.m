% This functions contains the rOMT algorithm and post-processing of the
% results on a sample dataset named 'gauss2' which is from a Gaussian sphere dataset.

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

%% Run rOMT
if cfg.reinitR
    [cfg, flag] = runROMT_par(cfg);
else
    [cfg, flag] = runROMT(cfg);
end

%% Run post-processing

[cfg, map, SL, stream, PATH] = runGLAD(cfg);



%% Visualization of speed map (full)
x = 1:cfg.true_size(1);
y = 1:cfg.true_size(2);
z = 1:cfg.true_size(3);
figure,
hs=slice(y,x,z,map.s_full,y,x,z); 
set(hs,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.04);
alpha('color'),alphamap(linspace(0,1,100))
%title(sprintf('Speed Map'),'Fontsize',20)
grid off, box off, axis image
xticks(0:10:cfg.true_size(1)); yticks(0:5:cfg.true_size(2)); zticks(0:5:cfg.true_size(3));
ax = gca; ax.FontSize = 16; 
xlabel('x-axis','FontSize',25),ylabel('y-axis','FontSize',25),zlabel('z-axis','FontSize',25)
colormap(jet)
caxis([0,1])%caxis([0,1])
view([242.1011   14.4475])
set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.5195 0.2292 0.2508 0.3667],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
%colorbar('FontSize',18), 
grid on,
xlim([0 cfg.true_size(1)]); ylim([0 cfg.true_size(2)]); zlim([0 cfg.true_size(3)])

cb = colorbar;
cb.Ticks = linspace(0, 1, 6);
cb.TickLabels = num2cell(0:0.2:1);
cb.FontSize = 20;
text(1,-3,20,'speed (a.u.)','Rotation',90,'FontSize',20);

saveas(gcf, sprintf('%s/%s/%s_LagSpeed_E%02d_%02d.png',cfg.out_dir,cfg.outdir_s,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 
%
%% Visualization of Pe map (full)
x = 1:cfg.true_size(1);
y = 1:cfg.true_size(2);
z = 1:cfg.true_size(3);
figure,
hs=slice(y,x,z,map.Pe_full,y,x,z); 
set(hs,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.04);
alpha('color'),alphamap(linspace(0,1,100))
%title(sprintf('Speed Map'),'Fontsize',20)
grid off, box off, axis image
xticks(0:10:cfg.true_size(1)); yticks(0:5:cfg.true_size(2)); zticks(0:5:cfg.true_size(3));
ax = gca; ax.FontSize = 16; 
xlabel('x-axis','FontSize',25),ylabel('y-axis','FontSize',25),zlabel('z-axis','FontSize',25)
colormap(jet)
caxis([0,800])%caxis([0,20])%caxis([0,10])%caxis([0,800])
view([242.1011   14.4475])
set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.5195 0.2292 0.2508 0.3667],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
%colorbar('FontSize',18), 
grid on,
xlim([0 cfg.true_size(1)]); ylim([0 cfg.true_size(2)]); zlim([0 cfg.true_size(3)])

cb = colorbar;
cb.Ticks = linspace(0, 800, 6);
cb.TickLabels = num2cell(0:160:800);
cb.FontSize = 20;
text(1,-3,30,'\it{Pe}','Rotation',90,'FontSize',20);

saveas(gcf, sprintf('%s/%s/%s_LagPe_E%02d_%02d.png',cfg.out_dir,cfg.outdir_s,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 
%
%% Visualization of flux vectors
figure,
strid = 40;%40;
magnify = 1;%1;%.5;%1;
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
q = quiver3(PATH.startp(PATH.ind_brain(1:strid:end),2),PATH.startp(PATH.ind_brain(1:strid:end),1),PATH.startp(PATH.ind_brain(1:strid:end),3),PATH.disp(PATH.ind_brain(1:strid:end),2)*magnify,PATH.disp(PATH.ind_brain(1:strid:end),1)*magnify,PATH.disp(PATH.ind_brain(1:strid:end),3)*magnify,...
    'color','r','LineWidth',1.5,'MaxHeadSize',0.5,'AutoScale','off','DisplayName','flux vectors');
%title(sprintf('Velocity Flux Vectors'),'FontSize',20, 'Interpreter', 'none'), 
set(gcf, 'Position', [376 49 1256 719])
xticks(0:10:cfg.true_size(1)); yticks(0:5:cfg.true_size(2)); zticks(0:5:cfg.true_size(3));
ax = gca; ax.FontSize = 16; 
xlabel('x-axis','FontSize',25),ylabel('y-axis','FontSize',25),zlabel('z-axis','FontSize',25)
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

view([242.1011   14.4475])
set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.5195 0.2292 0.2508 0.3667],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
%colorbar('FontSize',18), 
grid on,
xlim([0 cfg.true_size(1)]); ylim([0 cfg.true_size(2)]); zlim([0 cfg.true_size(3)])

cb = colorbar;
cb.Ticks = linspace(0, 800, 6);
cb.TickLabels = num2cell(0:160:800);
cb.FontSize = 20;
text(1,-3,10,'distance of movement (a.u.)','Rotation',90,'FontSize',20);

saveas(gcf, sprintf('%s/%s/%s_LagFluxVector_E%02d_%02d.png',cfg.out_dir,cfg.outdir_v,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 

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
axis image
colormap('jet');
view([242.1011   14.4475])
xticks(0:10:cfg.true_size(1)); yticks(0:5:cfg.true_size(2)); zticks(0:5:cfg.true_size(3));
ax = gca; ax.FontSize = 16; 
xlabel('x-axis','FontSize',25),ylabel('y-axis','FontSize',25),zlabel('z-axis','FontSize',25)
set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.5195 0.2292 0.2508 0.3667],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
%colorbar('FontSize',18), 
grid on,
xlim([0 cfg.true_size(1)]); ylim([0 cfg.true_size(2)]); zlim([0 cfg.true_size(3)])

%title(sprintf('Pathlines'), 'Interpreter', 'none', 'FontSize',18);

cb = colorbar;
set(cb,'YTick',[])
text(1,-3,53,'end','Rotation',90,'FontSize',20);
text(1,-3,-1,'start','Rotation',90,'FontSize',20);

saveas(gcf, sprintf('%s/%s/%s_LagPathlines_E%02d_%02d.png',cfg.out_dir,cfg.outdir_v,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 

%% Visualization of speedlines
figure,
spdmax = 1;%3;%1;%0.3;
Num = 100; % the larger, the more intervals in colors

SL2 = SL(PATH.ind_brain);
SL_spd = stream.sstream(PATH.ind_brain);
nSL = length(SL2);
colors = jet(Num);

for ind = 1:10:nSL
    SL_tmp = SL2{ind};
    SL_spd_tmp = SL_spd{ind};
    SL_spd_tmp(SL_spd_tmp>spdmax) = spdmax;
    SL_spd_rk = round(SL_spd_tmp/spdmax*Num);
    SL_spd_rk(SL_spd_rk<1) = 1;
    hlines = patch([SL_tmp(:,2);NaN],[SL_tmp(:,1);NaN],[SL_tmp(:,3);NaN],[1,1,1]);
    set(hlines,'FaceVertexCData',[colors(SL_spd_rk,:);colors(SL_spd_rk(end),:)],'EdgeColor','flat','FaceColor','none');
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
view([242.1011   14.4475])
xticks(0:10:cfg.true_size(1)); yticks(0:5:cfg.true_size(2)); zticks(0:5:cfg.true_size(3));
ax = gca; ax.FontSize = 16; 
xlabel('x-axis','FontSize',25),ylabel('y-axis','FontSize',25),zlabel('z-axis','FontSize',25)
set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.5195 0.2292 0.2508 0.3667],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
axis image
colormap('jet'); grid on;
%cbh = colorbar ; 
%cbh.Ticks = linspace(0, 1, 11) ; 
%cbh.TickLabels = num2cell(linspace(0, spdmax, 11));
%colorbar('FontSize',18)
%title(sprintf('Speedlines'), 'Interpreter', 'none', 'FontSize',18);

cb = colorbar;
cb.Ticks = linspace(0, 1, 6);
cb.TickLabels = num2cell(0:0.2:1);
cb.FontSize = 20;
text(1,-3,20,'speed (a.u.)','Rotation',90,'FontSize',20);

xlim([0 cfg.true_size(1)]); ylim([0 cfg.true_size(2)]); zlim([0 cfg.true_size(3)])
saveas(gcf, sprintf('%s/%s/%s_LagSpdlines_E%02d_%02d.png',cfg.out_dir,cfg.outdir_v,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 


%% Visualization of masslines
%{
figure,
massmax = 30;
Num = 100; % the larger, the more intervals in colors

SL2 = SL(PATH.ind_brain);
SL_mass = stream.rstream(PATH.ind_brain);
nSL = length(SL2);
colors = jet(Num);

for ind = 1:10:nSL
    SL_tmp = SL2{ind};
    SL_mass_tmp = SL_mass{ind};
    SL_mass_tmp(SL_mass_tmp>massmax) = massmax;
    SL_mass_rk = round(SL_mass_tmp/massmax*Num);
    SL_mass_rk(SL_mass_rk<1) = 1;
    hlines = patch([SL_tmp(:,2);NaN],[SL_tmp(:,1);NaN],[SL_tmp(:,3);NaN],[1,1,1]);
    set(hlines,'FaceVertexCData',[colors(SL_mass_rk,:);colors(SL_mass_rk(end),:)],'EdgeColor','flat','FaceColor','none');
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
view([242.1011   14.4475])
xticks(0:10:cfg.true_size(1)); yticks(0:5:cfg.true_size(2)); zticks(0:5:cfg.true_size(3));
ax = gca; ax.FontSize = 16; 
xlabel('x-axis','FontSize',25),ylabel('y-axis','FontSize',25),zlabel('z-axis','FontSize',25)
set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.5195 0.2292 0.2508 0.3667],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
axis image
colormap('jet'); grid on;
%cbh = colorbar ; 
%cbh.Ticks = linspace(0, 1, 11) ; 
%cbh.TickLabels = num2cell(linspace(0, massmax, 11));
%colorbar('FontSize',18)
%title(sprintf('Speedlines'), 'Interpreter', 'none', 'FontSize',18);

xlim([0 cfg.true_size(1)]); ylim([0 cfg.true_size(2)]); zlim([0 cfg.true_size(3)])
saveas(gcf, sprintf('%s/%s/%s_LagMasslines_E%02d_%02d.png',cfg.out_dir,cfg.outdir_v,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 
%}
%% Visualization of Pe-lines
figure,
Pemax = 800;%20;%10;%800;
Num = 100; % the larger, the more intervals in colors

SL2 = SL(PATH.ind_brain);
SL_Pe = stream.pestream(PATH.ind_brain);
nSL = length(SL2);
colors = jet(Num);

for ind = 1:10:nSL
    SL_tmp = SL2{ind};
    SL_Pe_tmp = SL_Pe{ind};
    SL_Pe_tmp(SL_Pe_tmp>Pemax) = Pemax;
    SL_Pe_rk = round(SL_Pe_tmp/Pemax*Num);
    SL_Pe_rk(SL_Pe_rk<1) = 1;
    hlines = patch([SL_tmp(:,2);NaN],[SL_tmp(:,1);NaN],[SL_tmp(:,3);NaN],[1,1,1]);
    set(hlines,'FaceVertexCData',[colors(SL_Pe_rk,:);colors(SL_Pe_rk(end),:)],'EdgeColor','flat','FaceColor','none');
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
view([242.1011   14.4475])
xticks(0:10:cfg.true_size(1)); yticks(0:5:cfg.true_size(2)); zticks(0:5:cfg.true_size(3));
ax = gca; ax.FontSize = 16; 
xlabel('x-axis','FontSize',25),ylabel('y-axis','FontSize',25),zlabel('z-axis','FontSize',25)
set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.5195 0.2292 0.2508 0.3667],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
axis image
colormap('jet'); grid on;
%cbh = colorbar ; 
%cbh.Ticks = linspace(0, 1, 11) ; 
%cbh.TickLabels = num2cell(linspace(0, massmax, 11));
%colorbar('FontSize',18)
%title(sprintf('Speedlines'), 'Interpreter', 'none', 'FontSize',18);

cb = colorbar;
cb.Ticks = linspace(0, 1, 6);
cb.TickLabels = num2cell(0:160:800);
cb.FontSize = 20;
text(1,-3,30,'\it{Pe}','Rotation',90,'FontSize',20);

xlim([0 cfg.true_size(1)]); ylim([0 cfg.true_size(2)]); zlim([0 cfg.true_size(3)])
saveas(gcf, sprintf('%s/%s/%s_LagPelines_E%02d_%02d.png',cfg.out_dir,cfg.outdir_v,cfg.tag,cfg.first_time,cfg.last_time+cfg.time_jump)); 

%% raw data
x = 1:Nx;
y = 1:Ny;
z = 1:Nz;

for i = 1:(cfg.last_time-cfg.first_time)/cfg.time_jump+2

    figure,
    hs=slice(y,x,z,cfg.vol(i).data,y,x,z); 
    set(hs,'EdgeColor','none','FaceColor','interp','FaceAlpha',0.04);
    alpha('color'),alphamap(linspace(0,1,100))
    %title(sprintf('Gaussian Data t = %d/%d',i,(cfg.last_time-cfg.first_time)/cfg.time_jump+2),'Fontsize',20)
    grid off, box off, axis image
    xticks(0:10:Nx); yticks(0:5:Ny); zticks(0:5:Nz);
    ax = gca; ax.FontSize = 16; 
    xlabel('x-axis','FontSize',25),ylabel('y-axis','FontSize',25),zlabel('z-axis','FontSize',25)
    xlim([0 Nx]); ylim([0 Ny]); zlim([0 Nz])
    colormap(jet)
    caxis([0,30])
    view([242.1011   14.4475])
    set(gca,'Color',[1,1,1]), set(gcf,'unit','normalized','position',[0.5195 0.2292 0.2508 0.3667],'Color',[1,1,1]), set(gcf, 'InvertHardcopy', 'off')
    %colorbar('FontSize',18), 
    grid on,
    saveas(gcf, sprintf('test_results/%s/%s_data_E%02d.png',cfg.tag,cfg.tag,i)); 

end
