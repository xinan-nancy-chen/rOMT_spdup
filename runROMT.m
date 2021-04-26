function [] = runROMT(cfg)


cfg.true_size = round(cfg.size_factor*[length(cfg.x_range),length(cfg.y_range),length(cfg.z_range)]);
cfg.version = sprintf('diff_%s_tj_%d_dt_%2.1f_nt_%d_ti_%d_tf_%d_uini_0_beta_%5.4f_R_gamma_%4.3f_dtri%d_rsmooth%d_rreinit%d_source%d_dilate%d_pcg%d',...
            cfg.sig_str,cfg.time_jump,cfg.dt,cfg.nt,cfg.first_time,cfg.last_time,cfg.beta,cfg.gamma,cfg.dTri,cfg.smooth,cfg.reinitR,cfg.add_source,cfg.dilate,cfg.niter_pcg);
cfg.out_dir = sprintf('./test_results/%s/%s',cfg.tag,cfg.version);

reInitializeU = 1; %1 if reinitialize u to 0 before each time step; 0 if not, unless first time step

if ~exist(cfg.out_dir,'dir')
    mkdir(cfg.out_dir)
end

fname = sprintf('%s/record.txt',cfg.out_dir);
if ~exist(sprintf('%s/record.txt',cfg.out_dir),'file')
    csvwrite_with_headers(fname,[0 0 0 0 0 0 0 0 0],{'time-ind','ti','tf','phi','mk','Ru','phiN','max(u)','toc'});
end

mask = nii2mat(cfg.data_mask_path,cfg.x_range,cfg.y_range,cfg.z_range);
msk = zeros(size(mask));
msk(mask>0) = 1;

if cfg.do_resize
   msk = resizeMatrix(msk,round(cfg.size_factor.*size(msk)),'linear');
   msk(msk~=1) = 0;
end

if cfg.dilate>0
    [xr,yr,zr] = meshgrid(-cfg.dilate:cfg.dilate,-cfg.dilate:cfg.dilate,-cfg.dilate:cfg.dilate);
    strel = (xr/cfg.dilate).^2 + (yr/cfg.dilate).^2 + (zr/cfg.dilate).^2 <= 1;
    msk = imdilate(msk,strel);
end

rho_n = nii2mat(sprintf('%s%02d%s',cfg.data_dir,cfg.first_time,cfg.extension),cfg.x_range,cfg.y_range,cfg.z_range);
if cfg.do_resize
   rho_n = resizeMatrix(rho_n,round(cfg.size_factor.*size(rho_n)),'linear');
end
if cfg.smooth>0
    rho_n = affine_diffusion_3d(rho_n,cfg.smooth,0.1,1,1);
end
m=min(rho_n(:));

rho_n(~msk) = m;
rho_n = rho_n(:);

if ~exist(sprintf('%s/rho_%s_%d_t_0.mat',cfg.out_dir,cfg.tag,cfg.first_time),'file')
    save(sprintf('%s/rho_%s_%d_t_0.mat',cfg.out_dir,cfg.tag,cfg.first_time),'rho_n');
end

fprintf('\n =============== rOMT Starts ===============\n')
fprintf('______________________________________________\n\n')
fprintf(' tag:\t\t%s\n dataset:\t%s\n sigma:\t\t%.4f\n gamma:\t\t%.4f\n beta:\t\t%.4f\n nt:\t\t%d\n dt:\t\t%.2f\n pcg:\t\t%d\n',cfg.tag,cfg.dataset_name,cfg.sigma,cfg.gamma,cfg.beta,cfg.nt,cfg.dt,cfg.niter_pcg)
fprintf(' size:\t\t[%d,%d]\n do_resize:\t%d\n resize_factor:\t%.2f\n start frame:\t%d\n end frame:\t%d\n frame jump:\t%d\n\n\n',cfg.true_size(1),cfg.true_size(2),cfg.do_resize,cfg.size_factor,cfg.first_time,cfg.last_time+cfg.time_jump,cfg.time_jump)
%%
profile on
clear T; T = 0;
for tind = 1:length(cfg.first_time:cfg.time_jump:cfg.last_time)
    fprintf('tind = %d\n',tind)
    tic
    %{
    if exist(sprintf('%s/rhoNe_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.tag,cfg.first_time+(tind-1)*cfg.time_jump,cfg.first_time+(tind-1)*cfg.time_jump+cfg.time_jump,tind),'file')==2 && exist(sprintf('%s/u0_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.tag,cfg.first_time+(tind-1)*cfg.time_jump,cfg.first_time+(tind-1)*cfg.time_jump+cfg.time_jump,tind),'file') == 2        
        rho_n = load(sprintf('%s/rhoNe_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.tag,cfg.first_time+(tind-1)*cfg.time_jump,cfg.first_time+(tind-1)*cfg.time_jump+cfg.time_jump,tind));
        rho_n = rho_n.rho_n;
        u = load(sprintf('%s/u0_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.tag,cfg.first_time+(tind-1)*cfg.time_jump,cfg.first_time+(tind-1)*cfg.time_jump+cfg.time_jump,tind));
        u = reshape(u.u,[],cfg.nt);
        continue
    end
    %}
    if cfg.reinitR
        rho_0 = nii2mat(sprintf('%s%02d%s',cfg.data_dir,cfg.first_time+(tind-1)*cfg.time_jump,cfg.extension),cfg.x_range,cfg.y_range,cfg.z_range);
        if cfg.do_resize
           rho_0 = resizeMatrix(rho_0,round(cfg.size_factor.*size(rho_0)),'linear');
        end
        if cfg.smooth>0
	    rho_0 = affine_diffusion_3d(rho_0,cfg.smooth,0.1,1,1);
        end
        rho_0(~msk) = m;
        rho_0 = rho_0(:);
    else
        rho_0 = rho_n(:);
    end

    %true final density
    rho_N = nii2mat(sprintf('%s%02d%s',cfg.data_dir,cfg.first_time+tind*cfg.time_jump,cfg.extension),cfg.x_range,cfg.y_range,cfg.z_range);
    if cfg.do_resize
       rho_N = resizeMatrix(rho_N,round(cfg.size_factor.*size(rho_N)),'linear');
    end
    if cfg.smooth>0
        rho_N = affine_diffusion_3d(rho_N,cfg.smooth,0.1,1,1);
    end
    rho_N(~msk) = m;
    
    par = paramInitFunc(cfg.true_size',cfg.nt,cfg.dt,cfg.sigma,cfg.add_source,cfg.gamma,cfg.beta,cfg.niter_pcg,cfg.dTri);
    par.drhoN     = rho_N(:);
    
    % initial guess for u:
    if tind == 1 || reInitializeU
    u = zeros(par.dim*prod(par.n),par.nt);%zeros(2*prod(par.n),par.nt);
    end
    
    %% Descent for u
    fprintf('\n =============== Descent on u ===============\n')
    fprintf('______________________________________________\n\n')
    fprintf('i.lsiter\tphi    \t      descent output\n')
    fprintf('________    ___________     __________________\n')
    [u,phi,dphi] = GNblock_u(rho_0,u,par.nt,par.dt,par);
    
    
    [phi,mk,phiN,Rho,Ru]  = get_phi(rho_0,u,par.nt,par.dt,par);
    rho_n = Rho(:,end);
    btoc = toc;
    T = T + btoc;
    
    dlmwrite(fname,[tind,cfg.first_time+(tind-1)*cfg.time_jump,cfg.first_time+(tind-1)*cfg.time_jump+cfg.time_jump,phi,mk,Ru,phiN,max(u(:)),btoc],'-append');
    
    save(sprintf('%s/u0_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.tag,cfg.first_time+(tind-1)*cfg.time_jump,cfg.first_time+(tind-1)*cfg.time_jump+cfg.time_jump,tind),'u');
    save(sprintf('%s/rhoNe_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.tag,cfg.first_time+(tind-1)*cfg.time_jump,cfg.first_time+(tind-1)*cfg.time_jump+cfg.time_jump,tind),'rho_n');
    
    fprintf('tind = %d, max(u) = %5.4f\n',tind,max(u));
end
fprintf('\n =============== rOMT Ends ===============\n')
fprintf('\n Elapsed Time: %s\n',datestr(seconds(T),'HH:MM:SS'))

profile viewer
profile off

end