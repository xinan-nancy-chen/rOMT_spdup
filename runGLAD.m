function [] = runGLAD(cfg)


fprintf('\n =============== Post-processing Starts ===============\n')
fprintf('_________________________________________________________\n\n')

types = {'s','v'};

tag = cfg.tag;
cfg.true_size = round(cfg.size_factor*[length(cfg.x_range),length(cfg.y_range),length(cfg.z_range)]);
cfg.version = sprintf('diff_%s_tj_%d_dt_%2.1f_nt_%d_ti_%d_tf_%d_uini_0_beta_%5.4f_R_gamma_%4.3f_dtri%d_rsmooth%d_rreinit%d_source%d_dilate%d_pcg%d',...
            cfg.sig_str,cfg.time_jump,cfg.dt,cfg.nt,cfg.first_time,cfg.last_time,cfg.beta,cfg.gamma,cfg.dTri,cfg.smooth,cfg.reinitR,cfg.add_source,cfg.dilate,cfg.niter_pcg);
cfg.out_dir = sprintf('./test_results/%s/%s',cfg.tag,cfg.version);

for l = 1:length(types)
type = types{l};
switch type 
    case 's'
        paper_fig_str = sprintf('set0%02d',1);
    case 'v'
        paper_fig_str = sprintf('set0%02d',2);
    otherwise
        fprintf('getGLAD: non-applicable type!');
        return
end
formatOut = 'mmddyy';
date_str = datestr(now,formatOut);

% get GLAD parameters
glacfg = paramInitGLADpar(cfg,type);

% set usefule variables:
nt = cfg.nt;
n = cfg.true_size;
ti = cfg.first_time;
tf = cfg.last_time;
tj = cfg.time_jump;

%% load mask:
model_mask = load_untouch_nii(cfg.data_mask_path);
model_mask.hdr.dime.datatype = 16;
model_mask.hdr.dime.bitpix = 32;
if cfg.do_resize
    S = round(cfg.domain_size.*cfg.size_factor);
    model_mask.img = zeros(S);
    model_mask.hdr.dime.dim(2:4) = S;
else
    S = cfg.domain_size;
    model_mask.img = zeros(S);
end

mask2 = nii2mat(cfg.data_mask_path,cfg.x_range,cfg.y_range,cfg.z_range);
msk = zeros(size(mask2));
msk(mask2>0) = 1;
clear mask2

if cfg.do_resize
   msk = resizeMatrix(msk,round(cfg.size_factor.*size(msk)),'linear');
   msk(msk~=1) = 0;
end

if glacfg.do_spMsk
    mskSP = nii2mat(cfg.sp_mask_opts(glacfg.spMsk_ind).path,cfg.x_range,cfg.y_range,cfg.z_range);
    if cfg.do_resize
       mskSP = resizeMatrix(mskSP,round(cfg.size_factor.*size(mskSP)),'linear');
       mskSP(mskSP~=1) = 0;
    end
    mskSP(~msk) = 0; %only include part of spMSK that is inside the whole msk
else
    mskSP = msk;
end


if glacfg.do_spErodeMsk > 0
    [x_er,y_er,z_er] = meshgrid(-glacfg.do_spErodeMsk:glacfg.do_spErodeMsk,-glacfg.do_spErodeMsk:glacfg.do_spErodeMsk,-1:1);
    strel_e = (x_er/glacfg.do_spErodeMsk).^2 + (y_er/glacfg.do_spErodeMsk).^2 + (z_er).^2 <= 1;    
    mskSP = imerode(mskSP, strel_e);
elseif glacfg.do_spDilateMsk > 0    
    %first find if original mask was dilated for rOMT run
    mskDIL = msk;
    if cfg.dilate>0
        [xr,yr,zr] = meshgrid(-cfg.dilate:cfg.dilate,-cfg.dilate:cfg.dilate,-cfg.dilate:cfg.dilate);
        strel = (xr/cfg.dilate).^2 + (yr/cfg.dilate).^2 + (zr/cfg.dilate).^2 <= 1;
        mskDIL = imdilate(msk,strel);
    end
    %now dilate spMsk and remove any points outside of mask used for rOMT
    [x_dil,y_dil,z_dil] = meshgrid(-glacfg.do_spDilateMsk:glacfg.do_spDilateMsk,-glacfg.do_spDilateMsk:glacfg.do_spDilateMsk,-glacfg.do_spDilateMsk:glacfg.do_spDilateMsk);
    strel_d = (x_dil/glacfg.do_spDilateMsk).^2 + (y_dil/glacfg.do_spDilateMsk).^2 + (z_dil/glacfg.do_spDilateMsk).^2 <= 1;
    mskSP = imdilate(mskSP, strel_d);
    mskSP(~mskDIL) = 0;    
end
%mskSP(~msk) = 0;  

if glacfg.do_sp
    dpsnrv_max = nii2mat(cfg.max_dpsnrv,cfg.x_range,cfg.y_range,cfg.z_range);
    if cfg.do_resize
       dpsnrv_max = resizeMatrix(dpsnrv_max,round(cfg.size_factor.*size(dpsnrv_max)),'linear');
    end
    mind = find((mskSP>0) & (dpsnrv_max>glacfg.sp_thresh));
    glacfg.do_sp_str = sprintf('_dpsnrv_min_%d',glacfg.sp_thresh);
else
    mind = find(mskSP>0);
    glacfg.do_sp_str = '';
end

sig_str = cfg.sig_str;
outdir = sprintf('LPPA_%s_%s',paper_fig_str,date_str);
switch type
    case 's'
        outdir_long = sprintf('LPPA_%s_type_%s_%s_img%s_mdt%d%s_%s_erode%d_dilate%d_nEul%d_cutoffStr_%s_concThresh%5.4f_spdThresh%5.4f_minIm0_%d_%s_slTol%d_clusters_centroid%d_nbp%d_metric%s_qbthresh%d_clusTol%d_clusCutoff%d_diff%s_tj%d_nt%d%s_%s_%s',...
            cfg.tag,'speedmap',glacfg.flw_type,glacfg.RD,glacfg.mdt,glacfg.XT,glacfg.spMSK_str,glacfg.do_spErodeMsk,glacfg.do_spDilateMsk,glacfg.nEulStep,glacfg.cutoff_str,glacfg.thresholds.conc,glacfg.thresholds.speed,glacfg.minIm0,glacfg.spSTR,glacfg.sl_tol,glacfg.do_centroid,glacfg.nbp,glacfg.ms,glacfg.qbthresh,glacfg.clus_tol,glacfg.clus_cutoff,sig_str,tj,nt,glacfg.do_sp_str,paper_fig_str,date_str);
    case 'v'
        outdir_long = sprintf('LPPA_%s_type_%s_%s_img%s_mdt%d%s_%s_erode%d_dilate%d_nEul%d_cutoffStr_%s_concThresh%5.4f_spdThresh%5.4f_minIm0_%d_%s_slTol%d_%s_%s_diff%s_tj%d_nt%d%s_%s_%s',...
            cfg.tag,'vectors',glacfg.flw_type,glacfg.RD,glacfg.mdt,glacfg.XT,glacfg.spMSK_str,glacfg.do_spErodeMsk,glacfg.do_spDilateMsk,glacfg.nEulStep,glacfg.cutoff_str,glacfg.thresholds.conc,glacfg.thresholds.speed,glacfg.minIm0,glacfg.spSTR,glacfg.sl_tol,glacfg.threshstr,glacfg.threshstr2,sig_str,tj,nt,glacfg.do_sp_str,paper_fig_str,date_str);
    otherwise
        fprintf('getGLAD: non-applicable dataset_name!');
        return
end
if ~exist(sprintf('%s/%s',cfg.out_dir,outdir),'dir')
    mkdir(sprintf('%s/%s',cfg.out_dir,outdir))
end

%%
%
fid = fopen(sprintf('%s/%s/%s_record_%s_%s.txt',cfg.out_dir,outdir,cfg.tag,paper_fig_str,date_str),'a+');
fprintf(fid,'%s/%s/%s_record_%s_%s\noutdir-long: %s',cfg.out_dir,outdir,cfg.tag,paper_fig_str,date_str,outdir_long);
switch type
    case 's'
        title_str = sprintf('============= initiating...\n\nLagrangian-Pathline (%s data, affSmooth = %d, dilate = %d), \nanalysis type = %s\n \nflw_type = %s, img = %s, mdt = %d(%s), %s, spErode = %d, spDilate = %d, nEulStep= %d, \ncutoffStr = %s, concThresh = %5.4f, spdThresh = %5.4f, minIm0 = %d, %s, slTol = %d, \nclusters, centroid = %d, nbp = %d, metric = %s, qbthresh = %d, clusTol = %d, clusCutoff = %d, \ndiff = %s, tj = %d, nt= %d%s_%s_%s\n\n',...
            cfg.tag,cfg.smooth,cfg.dilate,'speedmap',glacfg.flw_type,glacfg.RD,glacfg.mdt,glacfg.XT,glacfg.spMSK_str,glacfg.do_spErodeMsk,glacfg.do_spDilateMsk,glacfg.nEulStep,glacfg.cutoff_str,glacfg.thresholds.conc,glacfg.thresholds.speed,glacfg.minIm0,glacfg.spSTRtitle,glacfg.sl_tol,glacfg.do_centroid,glacfg.nbp,glacfg.ms,glacfg.qbthresh,glacfg.clus_tol,glacfg.clus_cutoff,sig_str,tj,nt,glacfg.do_sp_str,paper_fig_str,date_str);
    case 'v'
        title_str = sprintf('============= initiating...\n\nLagrangian-Pathline (%s data, affSmooth = %d, dilate = %d), \nanalysis type = %s\n \nflw_type = %s, img = %s, mdt = %d(%s), %s, spErode = %d, spDilate = %d, nEulStep= %d, \ncutoffStr = %s, concThresh = %5.4f, spdThresh = %5.4f, minIm0 = %d, %s, slTol = %d, \n%s\n%s\ndiff = %s, tj = %d, nt = %d%s_%s_%s\n\n',...
            cfg.tag,cfg.smooth,cfg.dilate,'vectors',glacfg.flw_type,glacfg.RD,glacfg.mdt,glacfg.XT,glacfg.spMSK_str,glacfg.do_spErodeMsk,glacfg.do_spDilateMsk,glacfg.nEulStep,glacfg.cutoff_str,glacfg.thresholds.conc,glacfg.thresholds.speed,glacfg.minIm0,glacfg.spSTRtitle,glacfg.sl_tol,glacfg.threshstr,glacfg.threshstr2,sig_str,tj,nt,glacfg.do_sp_str,paper_fig_str,date_str);
    otherwise
        fprintf('getGLAD: non-applicable dataset_name!');
        return
end
fprintf(title_str)
%%

[x, y, z] = meshgrid(1:n(2), 1:n(1), 1:n(3));
[syind,sxind,szind] = ind2sub(n,mind); %find indices of all voxels inside the ROI-SP
switch glacfg.spType
    case 'ordered'
        sy = syind(1:glacfg.fs:end);
        sx = sxind(1:glacfg.fs:end);
        sz = szind(1:glacfg.fs:end);
    case 'uniform'
        mskSPvol = sum(mskSP(:)); %volume of mask used to select start points
        NSP = round(glacfg.spPerc*mskSPvol/100);
        if ~exist(sprintf('%s_%s%d_spPerc%d_nsp%d.mat',tag,glacfg.spMsk_name,glacfg.spMsk_ind,glacfg.spPerc,NSP),'file')
            [spIND,spINDid] = datasample(mind,NSP,'Replace',false);
            %spIND = mind(spINDid);
            save(sprintf('%s_%s%d_spPerc%d_nsp%d.mat',tag,glacfg.spMsk_name,glacfg.spMsk_ind,glacfg.spPerc,NSP),'spIND');
        else
            load(sprintf('%s_%s%d_spPerc%d_nsp%d.mat',tag,glacfg.spMsk_name,glacfg.spMsk_ind,glacfg.spPerc,NSP));
        end
        [sy,sx,sz] = ind2sub(n,sort(spIND,'ascend'));
end
 
%% Lagrangian pathlines
%variables if nt > 1
cfg.n = n';
h1 = 1; h2 = 1; h3 = 1;
cfg.h1 = h1.*ones(n(1),1);
cfg.h2 = h2.*ones(n(2),1);
cfg.h3 = h3.*ones(n(3),1);
cfg.dim = length(cfg.n);
switch cfg.dTri
    case 1
        cfg.bc = 'closed';
    case 3
        cfg.bc = 'open';
end
[cfg.Xc,cfg.Yc,cfg.Zc] = getCellCenteredGrid(cfg.h1,cfg.h2,cfg.h3);
cfg.Grad = getCellCenteredGradMatrix({'ccn' 'ccn' 'ccn'},cfg.h1,cfg.h2,cfg.h3);
Mdis = -cfg.sigma.*cfg.Grad'*cfg.Grad;
%initialize streamlines:
nsp = length(sx);
%convert from matlab grid to cell-centered grid:
s1 = (sy-0.5).*h1; %i/y-axis
s2 = (sx-0.5).*h2; %j/x-axis
s3 = (sz-0.5).*h3; %k/z-axis
sp_123 = [s1,s2,s3];


pcur = sp_123; %current point i.e. list of current location in each streamline that hasn't been terminated
npoints = length(pcur); %keep track of the # of streamlines that have not yet been terminated 

SL = cell(1,nsp);
RHO_SL = cell(1,nsp);
AUGSPD_SL = cell(1,nsp);
SPD_SL = cell(1,nsp);
T_SL = cell(1,nsp); %store timestep
FLX_SL = cell(1,nsp); %store flux
dR_SL = cell(1,nsp); %store drho/dt

% Initialize masks
maskCluster = model_mask;
maskRho = model_mask;
maskSpeed = model_mask;
maskAugSpeed = model_mask;
maskTime = model_mask;
maskFlux = model_mask;
maskdRho = model_mask;

% make sure don't use sl-info from previous timestep
File1 = fullfile(cd, 'pl_cur.mat');
File2 = fullfile(cd, 'pli_array.mat');
File3 = fullfile(cd, 'pl_centroid_array.mat');
if exist(File1,'file')
    delete(File1);
end
if exist(File2,'file')
    delete(File2);
end
if exist(File3,'file')
    delete(File3);
end

%initiate pathlines
pl = NaN(npoints,3,length(ti:tj:tf)*nt*glacfg.nEulStep+1);
plrho = NaN(npoints,length(ti:tj:tf)*nt*glacfg.nEulStep+1);
plspd = NaN(npoints,length(ti:tj:tf)*nt*glacfg.nEulStep+1);
plaugspd = NaN(npoints,length(ti:tj:tf)*nt*glacfg.nEulStep+1);
plt = NaN(npoints,length(ti:tj:tf)*nt*glacfg.nEulStep+1);
plflx = NaN(npoints,length(ti:tj:tf)*nt*glacfg.nEulStep+1);
pldr =  NaN(npoints,length(ti:tj:tf)*nt*glacfg.nEulStep+1);

step = 1;
xt = pcur; %current position in rcl orientation
xt = max([h1*0.5,h2*0.5,h3*0.5],min(xt,[h1*(n(1)-.5001),h2*(n(2)-.5001),h3*(n(3)-.5001)]));

%% running pathlines

for t1 = ti:tj:tf
    U = load(sprintf('%s/u0_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.tag,t1,t1+tj,(t1-ti)/tj + 1), 'u');
    U = reshape(U.u,[],nt);
    
    if strcmp(glacfg.RD,'R')
        if t1 == ti
            RHO = load(sprintf('%s/rho_%s_%d_t_0.mat',cfg.out_dir,cfg.tag,ti));
        else
            RHO = load(sprintf('%s/rhoNe_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.tag,t1-tj,t1,(t1-ti)./tj));
        end
        RHO = RHO.rho_n;
        if nt > 1
            RHO_t = [RHO, advecDiff(RHO,U(:),nt,cfg.dt,cfg)];
        else
            RHO_t = RHO;
        end
    end
    
    for t2 = 1:nt
        TIND = ((t1 - ti)/tj)*nt + t2;
        T = t1+(t2-1)*(tj/nt);
        fprintf('t = %d (t1 = %d, t2 = %d -> T = %.3f)\n',TIND,t1,t2,T);
        
        switch glacfg.RD
            case 'D'
                d = getData(tag,round(T),'none');
                if cfg.smooth>0
                    d = affine_diffusion_3d(d,cfg.smooth,0.1,1,1);
                end
            case 'R'
                d = reshape(RHO_t(:,t2),n);
        end
        
        if glacfg.minIm0
            d = d-min(d(:));
        end
        
        %compute streamlines
        u = reshape(U(:,t2),[],3);
        v1 = reshape(u(:,1),n);
        v2 = reshape(u(:,2),n);
        v3 = reshape(u(:,3),n);
        
        %add eps to d so can take log(d) and not log(0)
        [w2,w1,w3] = gradient(log(d+2*eps));
        u1 = v1 - cfg.sigma.*w1;
        u2 = v2 - cfg.sigma.*w2;
        u3 = v3 - cfg.sigma.*w3;
        
        switch glacfg.flw_type
            case 'vel'
                a1=u1;
                a2=u2;
                a3=u3;
            case 'flw'
                a1=u1.*d;
                a2=u2.*d;
                a3=u3.*d;
        end
           
        [Gx, Gy, Gz] = gradient(d);
        speed = reshape(sqrt(sum(u.^2,2)),n);
        speedAug = sqrt(u1.^2 + u2.^2 + u3.^2);
        img_flow = speed.*d;
        drdt_dif = Mdis*d(:);
        drdt_ad1 = Gx.*reshape(u(:,2),n) + Gy.*reshape(u(:,1),n) + Gz.*reshape(u(:,3),n);
        DIVu = divergence(x,y,z,reshape(u(:,2),n),reshape(u(:,1),n),reshape(u(:,3),n));
        drdt_ad2 = d(:).*DIVu(:);
        drdt = reshape(drdt_dif - (drdt_ad1(:)+drdt_ad2),n);
        
        %update first step of density and speed
        if step == 1

            pl(:,:,1) = xt;
            plrho(:,1) = d(sub2ind(n,sy,sx,sz));
            plspd(:,1) = speed(sub2ind(n,sy,sx,sz));
            plaugspd(:,1) = speedAug(sub2ind(n,sy,sx,sz));
            plt(:,1) = step;%tstep;
            plflx(:,1) = img_flow(sub2ind(n,sy,sx,sz));
            pldr(:,1) = drdt(sub2ind(n,sy,sx,sz));
        end
        
        switch glacfg.cutoff_str
            case 'min'
                conf.conc = 1;
                conf.speed = 1;
            case 'max'
                conf.conc = mean(d(d>0)) + std(d(d>0));
                conf.speed = mean(speed(speed>0)) + std(speed(speed>0));              
            case 'mean'
                conf.conc = mean(d(d>0));
                conf.speed = mean(speed(speed>0));
        end
        
        %vector field to be integrated in order to compute streamlines
        V = [a1(:),a2(:),a3(:)];
        V(msk==0,:) = 0; %don't want to move outside of the masked ROI
        u(msk==0,:) = 0;
        
        for Estep = 1:glacfg.nEulStep
            step = step + 1;
            V_interp = interp_vel(V,n,xt(:,1),xt(:,2),xt(:,3),[h1,h2,h3]);
            u_interp = interp_vel(u,n,xt(:,1),xt(:,2),xt(:,3),[h1,h2,h3]);
            conc_interp = interpF(d,n,xt(:,1),xt(:,2),xt(:,3),[h1,h2,h3]);
            spdAug_interp = sqrt(sum(V_interp.^2,2));
            spd_interp = sqrt(sum(u_interp.^2,2));
            flx_interp = interpF(img_flow,n,xt(:,1),xt(:,2),xt(:,3),[h1,h2,h3]);
            dr_interp = interpF(drdt,n,xt(:,1),xt(:,2),xt(:,3),[h1,h2,h3]);
            
            conc_conf_interp = conc_interp./conf.conc;
            spd_conf_interp = spd_interp./conf.speed;
            
            thrsh_ind = find(conc_conf_interp > glacfg.thresholds.conc & spd_conf_interp > glacfg.thresholds.speed);

            if isempty(thrsh_ind)
                break
            else
                switch glacfg.XT
                    case'T'
                        xt(thrsh_ind,:) = xt(thrsh_ind,:) + glacfg.mdt.*V_interp(thrsh_ind,:);%so keep all current locations in case can take a step later when conc gets there
                    case 'X'
                        xt(thrsh_ind,:) = xt(thrsh_ind,:) + glacfg.mdt.*(V_interp(thrsh_ind,:)./spd_interp);%so keep all current locations in case can take a step later when conc gets there
                end
                %make sure it stays in bounds:
                xt = max([h1*0.5,h2*0.5,h3*0.5],min(xt,[h1*(n(1)-.5001),h2*(n(2)-.5001),h3*(n(3)-.5001)]));
                pl(thrsh_ind,:,step) = xt(thrsh_ind,:);
                plrho(thrsh_ind,step)  = conc_interp(thrsh_ind);
                plspd(thrsh_ind,step) = spd_interp(thrsh_ind);
                plaugspd(thrsh_ind,step) = spdAug_interp(thrsh_ind);
                plt(thrsh_ind,step) =  step;
                plflx(thrsh_ind,step) = flx_interp(thrsh_ind);
                pldr(thrsh_ind,step) = dr_interp(thrsh_ind);
            end
        end
    end
end


%% select qualified pathlines
pli = 0; %counter for added pathlines
for sli = 1:npoints
    pl_cur = squeeze(pl(sli,:,:))';
    aind = any(~isnan(pl_cur),2); %1 if row is not NaN, 0 if row is NaN
    pl_cur = pl_cur(aind,:);
    %check that unique sline has more than 1 point (remove not-move coordinates)
    [pl_cur,ia,ic] = unique(pl_cur,'rows','stable');
    if size(pl_cur,1)>glacfg.pln
        pli = pli + 1;
        SL{pli} = pl_cur;
        
        plr_cur = plrho(sli,aind)';
        pls_cur = plspd(sli,aind)';
        plas_cur = plaugspd(sli,aind)';
        plt_cur = plt(sli,aind)';
        plflx_cur = plflx(sli,aind)';
        pldr_cur = pldr(sli,aind)';
        
        RHO_SL{pli} = plr_cur(ia);
        SPD_SL{pli} = pls_cur(ia);
        AUGSPD_SL{pli} = plas_cur(ia);
        T_SL{pli} = plt_cur(ia);
        FLX_SL{pli} = plflx_cur(ia);
        dR_SL{pli} = pldr_cur(ia);
    end
end

%remove empty cell spaces from pathlines that were removed:
SL = SL(1:pli);
%RHO_SL = RHO_SL(1:pli);
SPD_SL = SPD_SL(1:pli);
%AUGSPD_SL = AUGSPD_SL(1:pli);
%T_SL = T_SL(1:pli);
%FLX_SL = FLX_SL(1:pli);
%dR_SL = dR_SL(1:pli);

sl_euc = cellfun(@(x) sqrt(sum((x(end,:)-x(1,:)).^2)),SL);% getcell array with euclidean length of sl between first and last point
%figure, histogram(sl_euc),title('sl-euc'),axis tight;
SL2 = SL(sl_euc>glacfg.sl_tol);%only keep streamlines whose Euclidean length between first and last points is larger than the threshold
%rstream2 = RHO_SL(sl_euc>glacfg.sl_tol);
sstream2 = SPD_SL(sl_euc>glacfg.sl_tol);
%asstream2 = AUGSPD_SL(sl_euc>glacfg.sl_tol);
%tstream2 = T_SL(sl_euc>glacfg.sl_tol);
%fstream2 = FLX_SL(sl_euc>glacfg.sl_tol);
%dstream2 = dR_SL(sl_euc>glacfg.sl_tol);

fprintf(' # of start points = %d\n # of effective pathlines after pathline-number (pln) threshold = %d \n # of effective pathlines after Euclidean dist (sl_tol) threshold = %d\n',npoints,pli,length(SL2))
pl_cur = cellfun(@(x) x(:,[2,1,3]),SL2,'UniformOutput',false);


switch type
    case 's'
        fprintf('getting speed map WITHOUT running QuickBundle.py...\n')
        % initialize temporary masks:
        s = zeros(n);%speed

        %getting clustered pathline start points
        spTMP = zeros(n);
        for k = 1:length(pl_cur)
            % streamline cluster mask:
            slines_tmp = pl_cur(k);
            spdlines_tmp = sstream2(k);
            %getting start points in world coords
            spx = cellfun(@(x) x(1,1),pl_cur(k));
            spy = cellfun(@(x) x(1,2),pl_cur(k));
            spz = cellfun(@(x) x(1,3),pl_cur(k));
            %getting start points in matlab coords
            spm1 = round(spy./h1 + 0.5);
            spm2 = round(spx./h2 + 0.5);
            spm3 = round(spz./h3 + 0.5);
            %save to spTMP
            indtmp = sub2ind(n,spm1,spm2,spm3);
            spTMP(indtmp) = 1;

            n_slines = size(slines_tmp,2);
            for ind_line = 1:n_slines
                %convert back to MATLAB grid
                sline = round((slines_tmp{ind_line}./[h1,h2,h3])+0.5);
                [sline,ia,~] = unique(sline, 'rows', 'stable');
                ssl = spdlines_tmp{ind_line}(ia);


                subs_1 = sline(:,2);
                subs_2 = sline(:,1);
                subs_3 = sline(:,3);
                subs_1(subs_1 < 1) = 1;
                subs_2(subs_2 < 1) = 1;
                subs_3(subs_3 < 1) = 1;
                subs_1(subs_1 > n(1)) = n(1);
                subs_2(subs_2 > n(2)) = n(2);
                subs_3(subs_3 > n(3)) = n(3);
                inds = sub2ind(n, subs_1, subs_2, subs_3);

                s(inds) = ssl;
            end
        end
        % QuickBundle.py version
        %{
        save('pl_cur.mat','pl_cur');
        %% === run QB.py Here ===
        
        fprintf('\n======= waiting to run run_dipyQB_pl.py =======\n')
        fprintf('\nInstructions:\n\npl_cur.mat has been saved at the current directory.\nDirectly run run_dipyQB_pl.py with Python also in this directory.\nResults will be saved automatically.\nThen come back to Matlab and press any key to continue.\n\n')
        pause
        fprintf('...Matlab code sucessfully continues...\n')

        %% add start points of pathlines considered for clustering analysis
        load(fullfile(cd,'pli_array.mat'));
        load(fullfile(cd,'pl_centroid_array.mat'));
        % only want clusters with more than tol # of sls:
        clus_length = cellfun('size',pli_array,2)'; 
        fprintf(' # of original clusters = %d\n',length(clus_length))
        
        sli = pli_array(clus_length>glacfg.clus_tol);
        sl_centroid = pl_centroid_array(clus_length>glacfg.clus_tol,:,:);
        nclus = size(sli,2);%# of clusters that made the cutoff
        clus_length = cellfun('size',sli,2)'; 
        fprintf(' # of clusters after cluster-length (clus_tol) threshold = %d\n',length(clus_length))
        [~,clus_size_rank] = sort(clus_length);

        %only keep largest N clusters where N = clus_cutoff:
        if glacfg.clus_cutoff>0
            if nclus > glacfg.clus_cutoff
                ind_tmp = zeros(nclus,1);
                ind_tmp(clus_size_rank(end-glacfg.clus_cutoff+1:end)) = 1;
                %% to keep clusters in same order that they were returned from QB alg:
                sli = sli(ind_tmp==1);
                sl_centroid = sl_centroid(ind_tmp==1,:,:);
                nclus = size(sli,2);%# of clusters that made the cutoff
                fprintf('nclus=%d\n',nclus)
                clus_length = cellfun('size',sli,2)';
                [~,clus_size_rank] = sort(clus_length);
            end
        end
        fprintf(' # of clusters after max-cluster-number (clus_cutoff) threshold = %d\n',nclus)

        sl_centroid = double(sl_centroid);
        slc = cell(1,nclus);
        for ind = 1:nclus
            slc{1,ind}=squeeze(sl_centroid(ind,:,:));
        end
        [B,I] = sort(squeeze(sl_centroid(:,1,1)));%sl_centroid is nclus x npoints x 3

        % initialize temporary masks:
        c = zeros(n);%clusters
        r = zeros(n);%density
        s = zeros(n);%speed
        as = zeros(n);%augspeed
        tt = zeros(n);%time step
        ff = zeros(n);%flux
        drt = zeros(n);%dr/dt

        % make cluster masks
        %getting clustered pathline start points
        spTMP = zeros(n);
        for ind_clus = 1:nclus
            slines_tmp = pl_cur([sli{I(ind_clus)}]+1);
            rlines_tmp = rstream2([sli{I(ind_clus)}]+1);
            spdlines_tmp = sstream2([sli{I(ind_clus)}]+1);
            aspdlines_tmp = asstream2([sli{I(ind_clus)}]+1);
            tlines_tmp = tstream2([sli{I(ind_clus)}]+1);
            flines_tmp = fstream2([sli{I(ind_clus)}]+1);
            dlines_tmp = dstream2([sli{I(ind_clus)}]+1);

            %getting start points in world coords
            spx = cellfun(@(x) x(1,1),pl_cur([sli{I(ind_clus)}]+1));
            spy = cellfun(@(x) x(1,2),pl_cur([sli{I(ind_clus)}]+1));
            spz = cellfun(@(x) x(1,3),pl_cur([sli{I(ind_clus)}]+1));
            %getting start points in matlab coords
            spm1 = round(spy./h1 + 0.5);
            spm2 = round(spx./h2 + 0.5);
            spm3 = round(spz./h3 + 0.5);
            %save to spTMP
            indtmp = sub2ind(n,spm1,spm2,spm3);
            spTMP(indtmp) = 1;

            n_slines = size(slines_tmp,2);
            for ind_line = 1:n_slines
                %convert back to MATLAB grid
                sline = round((slines_tmp{ind_line}./[h1,h2,h3])+0.5);
                [sline,ia,~] = unique(sline, 'rows', 'stable');
                rsl = rlines_tmp{ind_line}(ia);
                ssl = spdlines_tmp{ind_line}(ia);
                assl = aspdlines_tmp{ind_line}(ia);
                tsl = tlines_tmp{ind_line}(ia);
                fsl = flines_tmp{ind_line}(ia);
                dsl = dlines_tmp{ind_line}(ia);

                subs_1 = sline(:,2);
                subs_2 = sline(:,1);
                subs_3 = sline(:,3);
                subs_1(subs_1 < 1) = 1;
                subs_2(subs_2 < 1) = 1;
                subs_3(subs_3 < 1) = 1;
                subs_1(subs_1 > n(1)) = n(1);
                subs_2(subs_2 > n(2)) = n(2);
                subs_3(subs_3 > n(3)) = n(3);
                inds = sub2ind(n, subs_1, subs_2, subs_3);

                c(inds) = ind_clus*1;
                r(inds) = rsl;
                s(inds) = ssl;
                as(inds) = assl;
                tt(inds) = tsl;
                ff(inds) = fsl;
                drt(inds) = dsl;
            end
            %centroid mask:
            if glacfg.do_centroid
                %centroid_tmp = round(slc{clus_size_rank(ind_clus)});
                centroid_tmp = round(slc{I(ind_clus)});
                centroid_tmp = unique(centroid_tmp,'rows', 'stable');
                subs_1 = centroid_tmp(:,2);
                subs_2 = centroid_tmp(:,1);
                subs_3 = centroid_tmp(:,3);
                subs_1(subs_1 < 1) = 1;
                subs_2(subs_2 < 1) = 1;
                subs_3(subs_3 < 1) = 1;
                subs_1(subs_1 > n(1)) = n(1);
                subs_2(subs_2 > n(2)) = n(2);
                subs_3(subs_3 > n(3)) = n(3);
                inds = sub2ind(n, subs_1, subs_2, subs_3);
                cent_tmp(inds) = ind_clus*1;
            end
        end
        
        if glacfg.do_masked
            %Added to ensure nothing outside of masked region
            c(~msk) = 0;
            r(~msk) = 0;
            s(~msk) = 0;
            as(~msk) = 0;
            tt(~msk) = 0;
            ff(~msk) = 0;
            drt(~msk) = 0;
        end
        % only save within brain mask
        sp_ind = find(strcmp('brain',{cfg.sp_mask_opts(:).name}));
        mskROI = nii2mat(cfg.sp_mask_opts(sp_ind).path,cfg.x_range,cfg.y_range,cfg.z_range);
        msk_brain = mskROI>1;
        if cfg.do_resize
            msk_brain = resizeMatrix(double(msk_brain),round(cfg.size_factor.*size(msk_brain)),'linear');
            msk_brain(msk_brain~=1) = 0;
            c(msk_brain==0) = 0;
            r(msk_brain==0) = 0;
            s(msk_brain==0) = 0;
            as(msk_brain==0) = 0;
            tt(msk_brain==0) = 0;
            ff(msk_brain==0) = 0;
            drt(msk_brain==0) = 0;
            
            maskCluster.img(round(cfg.x_range(1)*cfg.size_factor):round(cfg.x_range(end)*cfg.size_factor),round(cfg.y_range(1)*cfg.size_factor):round(cfg.y_range(end)*cfg.size_factor),round(cfg.z_range(1)*cfg.size_factor):round(cfg.z_range(end)*cfg.size_factor)) = c;
            maskRho.img(round(cfg.x_range(1)*cfg.size_factor):round(cfg.x_range(end)*cfg.size_factor),round(cfg.y_range(1)*cfg.size_factor):round(cfg.y_range(end)*cfg.size_factor),round(cfg.z_range(1)*cfg.size_factor):round(cfg.z_range(end)*cfg.size_factor)) = r;
            maskSpeed.img(round(cfg.x_range(1)*cfg.size_factor):round(cfg.x_range(end)*cfg.size_factor),round(cfg.y_range(1)*cfg.size_factor):round(cfg.y_range(end)*cfg.size_factor),round(cfg.z_range(1)*cfg.size_factor):round(cfg.z_range(end)*cfg.size_factor)) = s;
            maskAugSpeed.img(round(cfg.x_range(1)*cfg.size_factor):round(cfg.x_range(end)*cfg.size_factor),round(cfg.y_range(1)*cfg.size_factor):round(cfg.y_range(end)*cfg.size_factor),round(cfg.z_range(1)*cfg.size_factor):round(cfg.z_range(end)*cfg.size_factor)) = as;
            maskTime.img(round(cfg.x_range(1)*cfg.size_factor):round(cfg.x_range(end)*cfg.size_factor),round(cfg.y_range(1)*cfg.size_factor):round(cfg.y_range(end)*cfg.size_factor),round(cfg.z_range(1)*cfg.size_factor):round(cfg.z_range(end)*cfg.size_factor)) = tt;
            maskFlux.img(round(cfg.x_range(1)*cfg.size_factor):round(cfg.x_range(end)*cfg.size_factor),round(cfg.y_range(1)*cfg.size_factor):round(cfg.y_range(end)*cfg.size_factor),round(cfg.z_range(1)*cfg.size_factor):round(cfg.z_range(end)*cfg.size_factor)) = ff;
            maskdRho.img(round(cfg.x_range(1)*cfg.size_factor):round(cfg.x_range(end)*cfg.size_factor),round(cfg.y_range(1)*cfg.size_factor):round(cfg.y_range(end)*cfg.size_factor),round(cfg.z_range(1)*cfg.size_factor):round(cfg.z_range(end)*cfg.size_factor)) = drt;
        else
            c(msk_brain==0) = 0;
            r(msk_brain==0) = 0;
            s(msk_brain==0) = 0;
            as(msk_brain==0) = 0;
            tt(msk_brain==0) = 0;
            ff(msk_brain==0) = 0;
            drt(msk_brain==0) = 0;
            
            maskCluster.img(cfg.x_range,cfg.y_range,cfg.z_range) = c;
            maskRho.img(cfg.x_range,cfg.y_range,cfg.z_range) = r;
            maskSpeed.img(cfg.x_range,cfg.y_range,cfg.z_range) = s;
            maskAugSpeed.img(cfg.x_range,cfg.y_range,cfg.z_range) = as;
            maskTime.img(cfg.x_range,cfg.y_range,cfg.z_range) = tt;
            maskFlux.img(cfg.x_range,cfg.y_range,cfg.z_range) = ff;
            maskdRho.img(cfg.x_range,cfg.y_range,cfg.z_range) = drt;
        end
        %}
        
        if glacfg.do_masked
            %Added to ensure nothing outside of masked region
            s(~msk) = 0;
        end
        % only save within brain mask
        sp_ind = find(strcmp('brain',{cfg.sp_mask_opts(:).name}));
        mskROI = nii2mat(cfg.sp_mask_opts(sp_ind).path,cfg.x_range,cfg.y_range,cfg.z_range);
        msk_brain = mskROI>1;
        if cfg.do_resize
            msk_brain = resizeMatrix(double(msk_brain),round(cfg.size_factor.*size(msk_brain)),'linear');
            msk_brain(msk_brain~=1) = 0;
            s(msk_brain==0) = 0;
            maskSpeed.img(round(cfg.x_range(1)*cfg.size_factor):round(cfg.x_range(end)*cfg.size_factor),round(cfg.y_range(1)*cfg.size_factor):round(cfg.y_range(end)*cfg.size_factor),round(cfg.z_range(1)*cfg.size_factor):round(cfg.z_range(end)*cfg.size_factor)) = s;
       
        else

            s(msk_brain==0) = 0;

            maskSpeed.img(cfg.x_range,cfg.y_range,cfg.z_range) = s;
        end
        
        % save to nifty (view in Amira later)
        
        save_untouch_nii(maskSpeed,sprintf('%s/%s/%s_LagSpeed_E%02d_%02d_%s_%s.nii',cfg.out_dir,outdir,cfg.tag,ti,tf+tj,paper_fig_str,date_str));
        
        %fprintf('For Speed Map, the next step is to visualize the speed nifty file in Amira.\nHowever, we can still plot in Matlab\n')
        fprintf('Speed Map in nifty format saved in %s/%s\n\n',cfg.out_dir,outdir)
        
        %% Visualization
        x = 1:n(1);
        y = 1:n(2);
        z = 1:n(3);
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

        saveas(gcf, sprintf('%s/%s/%s_LagSpeed_E%02d_%02d.png',cfg.out_dir,outdir,cfg.tag,ti,tf+tj)); 
    case 'v'

        outversion = sprintf('%s_%s',paper_fig_str,date_str);
        outdir = sprintf('LPPA_%s',outversion);

        if ~exist(sprintf('%s/%s',cfg.out_dir,outdir),'dir')
            fprintf('%s/%s\n Directory does not exist :(\n',cfg.out_dir,outdir);
            return
        else
            fprintf('%s: %s Directory exists :)\n',tag,outdir);
        end

        %sstream = sstream2;
        fprintf('Total original %d pathlines\n',length(SL2)); 

        %convert from cell-centered grid to matlab grid:
        SL = cellfun(@(x) [x(:,1)./h1+0.5,x(:,2)./h2+0.5,x(:,3)./h3+0.5],SL2,'UniformOutput',false); 
        clear SL2;
        
        dispnor = cellfun(@(x) (x(end,:)-x(1,:))/norm(x(end,:)-x(1,:)),SL,'UniformOutput',false); 
        disp = cellfun(@(x) x(end,:)-x(1,:),SL,'UniformOutput',false); 
        startp = cellfun(@(x) x(1,:),SL,'UniformOutput',false); 
        endp = cellfun(@(x) x(end,:),SL,'UniformOutput',false); 
        
        PATH.NPtsInPath = cellfun(@(x) size(x,1),SL); % number of points in each pathline
        PATH.LengthInPath = cellfun(@(x) sum(sqrt(sum(diff(x).^2,2))),SL); % total length of path in each pathline
        PATH.disp = reshape([disp{:}]',3,[])'; % displacement field
        PATH.dispnor = reshape([dispnor{:}]',3,[])'; % normalized displacement field
        PATH.startp = reshape([startp{:}]',3,[])'; % startp points
        PATH.endp = reshape([endp{:}]',3,[])'; % endp points
        PATH.displen = sqrt(PATH.disp(:,1).^2+PATH.disp(:,2).^2+PATH.disp(:,3).^2); % displacement length in each pathline

        anato = load_untouch_nii(cfg.anato); 
        anato = anato.img(cfg.x_range,cfg.y_range,cfg.z_range);
        sp_ind = find(strcmp('brain',{cfg.sp_mask_opts(:).name}));
        mskROI = nii2mat(cfg.sp_mask_opts(sp_ind).path,cfg.x_range,cfg.y_range,cfg.z_range);
        msk_brain = mskROI>1;
        
        if cfg.do_resize
           anato = resizeMatrix(anato,round(cfg.size_factor.*size(anato)),'linear');
           msk_brain = resizeMatrix(double(msk_brain),round(cfg.size_factor.*size(msk_brain)),'linear');
           msk_brain(msk_brain~=1) = 0;
        end
        anato(~msk) = 0;
        anato(msk_brain==0) = 0;
        indb = find(msk_brain(sub2ind(n,PATH.startp(:,1),round(PATH.startp(:,2)),PATH.startp(:,3)))==1); % index of those starting within brain
        
        % save to vtk
        SL_tmp = SL(indb);
        vtkwrite_pathlines(sprintf('%s/%s/%s_pathlines_lentol_%.1f_jp_%d_%s.vtk',cfg.out_dir,outdir,tag,glacfg.sl_tol,glacfg.jp,outversion),'polydata','lines',SL_tmp(1:glacfg.jp:end));
        vtkwrite(sprintf('%s/%s/%s_disp_lentol_%.2f_%s.vtk',cfg.out_dir,outdir,tag,glacfg.sl_tol,outversion), 'structured_grid', PATH.startp(indb,1), PATH.startp(indb,2), PATH.startp(indb,3), ... 
            'vectors', 'vector_field', PATH.disp(indb,1), PATH.disp(indb,2), PATH.disp(indb,3));
        vtkwrite(sprintf('%s/%s/%s_anato_%s.vtk',cfg.out_dir,outdir,tag,outversion), 'structured_points', 'mask', anato);
        
        %% NCA
        [xx, yy, zz] = meshgrid(1:n(2),1:n(1),1:n(3));

        startpind = sub2ind(n,round(PATH.startp(:,1)),round(PATH.startp(:,2)),round(PATH.startp(:,3)));

        Mean = zeros(length(PATH.displen),1); % mean in the neighbor
        WMean = Mean; % weighted mean in the neighbor
        STD = Mean; % L_2^2 distance to 1
        Np = Mean; % number of all points in neighbored paths
        Avepathlen = Mean; % mean path length in neighborhood
        %AveMaxpathspd = Mean; % average max speed  
        NNnum = zeros(length(PATH.displen),1); % number of neighbors

        for sp = 1:length(PATH.displen)
            ctr = round(PATH.startp(sp,:));
            binaryMap = (yy-ctr(1)).^2+(xx-ctr(2)).^2+(zz-ctr(3)).^2 <=glacfg.radius^2;
            indd = find(binaryMap); % startpoint ind in n1*n2*n3 coord.
            % pick neighboring points
            [~,Locb] = ismember(indd,startpind);
            Lia = Locb(Locb~=0);

            NNnum(sp) = length(Lia);
            costheta = dot(repmat(PATH.disp(sp,:),length(Lia),1),PATH.disp(Lia,:),2)/norm(PATH.disp(sp,:))./vecnorm(PATH.disp(Lia,:),2,2);
            WMean(sp) = sum(costheta.*PATH.displen(Lia))/sum(PATH.displen(Lia));
            Mean(sp) = mean(costheta);
            STD(sp) = sum((costheta-1).^2)/length(costheta);%std(costheta);
            Np(sp) = sum(PATH.NPtsInPath(Lia));
            Avepathlen(sp) = mean(PATH.LengthInPath(Lia));
            %AveMaxpathspd(sp) = mean(cellfun(@(x) max(x),sstream(Lia)));  

        end

        %%
        INDD = find(NNnum > glacfg.NNnum_tol & STD<glacfg.stdcut & Np>glacfg.Npcut & WMean>=glacfg.meancut(1) & Avepathlen>glacfg.Avepathlcut);
        INDD2 = find(NNnum > glacfg.NNnum_tol);

        %tmpind = INDD;
        %tmpind2 = setdiff(INDD2,INDD);
        %PATH.ADVind = tmpind;
        %PATH.DIFFind = tmpind2;
        %PATH.threshstr = glacfg.threshstr;

        %% further choose potential adv vectors
        maskadv = zeros(n);
        maskadv(sub2ind(n,round(PATH.startp(INDD,1)),round(PATH.startp(INDD,2)),round(PATH.startp(INDD,3))))=1;

        [xr,yr,zr] = meshgrid(-glacfg.maskdilate:glacfg.maskdilate,-glacfg.maskdilate:glacfg.maskdilate,-glacfg.maskdilate:glacfg.maskdilate);
        strel = (xr/glacfg.maskdilate).^2 + (yr/glacfg.maskdilate).^2 + (zr/glacfg.maskdilate).^2 <= 1;
        maskadvdia = imdilate(maskadv,strel);

        % find potential start points
        indp = find(maskadvdia(sub2ind(n,PATH.startp(:,1),round(PATH.startp(:,2)),round(PATH.startp(:,3))))==1 ...
            & maskadv(sub2ind(n,PATH.startp(:,1),round(PATH.startp(:,2)),round(PATH.startp(:,3))))==0);
        indpADV = find(maskadv(sub2ind(n,PATH.startp(:,1),round(PATH.startp(:,2)),round(PATH.startp(:,3))))==1);

        % set parameters
        %%
        Mean2 = -1*ones(length(PATH.displen),1); % mean in the neighbor
        WMean2 = Mean2; % weighted mean in the neighbor
        STD2 = Mean2; % L_2^2 distance to 1
        Np2 = Mean2; % number of all points in neighbored paths
        Avepathlen2 = Mean2; % mean path length in neighborhood
        %AveMaxpathspd2 = Mean2; % average max speed 
        pathl2 = Mean2; % path length in sp
        NNnum2 = -1*ones(length(PATH.displen),1); % number of neighbors

        for sp = 1:length(PATH.displen)
            if ~ismember(sp,indp)%sp is not in indp
                continue
            end
            ctr = round(PATH.startp(sp,:));
            binaryMap = (yy-ctr(1)).^2+(xx-ctr(2)).^2+(zz-ctr(3)).^2 <=glacfg.radius2^2;
            indd = find(binaryMap==1); % startpoint ind in n1*n2*n3 coord. 
            [~,Locb] = ismember(indd,indpADV); % only find neighbor in maskadv
            Lia = Locb(Locb~=0);

            NNnum2(sp) = length(Lia);
            costheta = dot(repmat(PATH.disp(sp,:),length(Lia),1),PATH.disp(Lia,:),2)/norm(PATH.disp(sp,:))./vecnorm(PATH.disp(Lia,:),2,2);
            %costheta(costheta<=-0.75) = 1; % remove those are too opposite trend
            WMean2(sp) = sum(costheta.*PATH.displen(Lia))/sum(PATH.displen(Lia));
            Mean2(sp) = mean(costheta);
            STD2(sp) = sum((costheta-1).^2)/length(costheta);%std(costheta);
            Np2(sp) = sum(PATH.NPtsInPath(Lia));
            Avepathlen2(sp) = mean(PATH.LengthInPath(Lia));
            pathl2(sp) = PATH.LengthInPath(sp);
            %AveMaxpathspd2(sp) = mean(cellfun(@(x) max(x),sstream(Lia)));
        end

        %%
        INDD_dia = find(NNnum2>glacfg.NNnum_tol2 & STD2<glacfg.stdcut2 & Np2>glacfg.Npcut2 & WMean2>=glacfg.meancut2(1) & Avepathlen2>glacfg.Avepathlcut2 | pathl2>glacfg.pathlcut2);
        fprintf('After further dilate to add more ADV, %d vectors are added among %d candidates to %d already ADV vectors\n',length(INDD_dia),length(indp),length(INDD))

        tmpind = [INDD;INDD_dia];
        tmpind2 = setdiff(INDD2,tmpind);
        PATH.ADVind = tmpind;
        PATH.DIFFind = tmpind2;
        PATH.threshstr = sprintf('Overall separation: \n%s, \nandfurther dilate to select more ADV: \n%s',glacfg.threshstr,glacfg.threshstr2);
        
        [InDadv,~,~] = intersect(PATH.ADVind,indb);
        [InDdiff,~,~] = intersect(PATH.DIFFind,indb);
        % save
        vtkwrite(sprintf('%s/%s/%s_ADVdisp_%s.vtk',cfg.out_dir,outdir,tag,outversion), 'structured_grid', PATH.startp(InDadv,1), PATH.startp(InDadv,2), PATH.startp(InDadv,3), ... 
            'vectors', 'vector_field', PATH.disp(InDadv,1), PATH.disp(InDadv,2), PATH.disp(InDadv,3));

        vtkwrite(sprintf('%s/%s/%s_DIFFdisp_%s.vtk',cfg.out_dir,outdir,tag,outversion), 'structured_grid', PATH.startp(InDdiff,1), PATH.startp(InDdiff,2), PATH.startp(InDdiff,3), ... 
            'vectors', 'vector_field', PATH.disp(InDdiff,1), PATH.disp(InDdiff,2), PATH.disp(InDdiff,3));
        
        fprintf('Flux vectors in vtk format saved in %s/%s\n\n',cfg.out_dir,outdir)
        
        %% Visualization
        figure,
        strid = 5;
        magnify = .5;%1;
        [x, y, z] = meshgrid(1:n(2), 1:n(1), 1:n(3));
        mskfv = isosurface(x,y,z,msk,0.5);
        mskp = patch(mskfv);
        mskp.FaceColor = [.17,.17,.17];
        mskp.FaceAlpha= 0.031;
        mskp.EdgeColor = [.17,.17,.17];
        mskp.EdgeAlpha= 0;
        mskp.DisplayName = 'mask';
        view([-178.3000   -2.0000]);

        grid on, axis image
        hold on,
        q = quiver3(PATH.startp(indb(1:strid:end),2),PATH.startp(indb(1:strid:end),1),PATH.startp(indb(1:strid:end),3),PATH.disp(indb(1:strid:end),2)*magnify,PATH.disp(indb(1:strid:end),1)*magnify,PATH.disp(indb(1:strid:end),3)*magnify,'color','r','LineWidth',2,'MaxHeadSize',0.3,'AutoScale','off','DisplayName','flux vectors');
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
            'ColorData', reshape(cmap(1:3,:,:), [], 4).');   %'
        set(q.Tail, ...
            'ColorBinding', 'interpolated', ...
            'ColorData', reshape(cmap(1:2,:,:), [], 4).');
        %hold on,
        %scatter3(PATH.startp(indb(1:strid:end),2),PATH.startp(indb(1:strid:end),1),PATH.startp(indb(1:strid:end),3),4,'MarkerFaceColor',[0 .75 .75],'DisplayName','start point')
        legend; colorbar;
        view([-188.3500   13.7463])
        set(gca,'Color',[0.85,0.85,0.93]), set(gcf,'unit','normalized','position',[0.1,1,0.4,0.5],'Color',[0.85,0.85,0.93]), set(gcf, 'InvertHardcopy', 'off')
        
        saveas(gcf, sprintf('%s/%s/%s_LagFluxVector_E%02d_%02d.png',cfg.out_dir,outdir,cfg.tag,ti,tf+tj)); 
end

end

end