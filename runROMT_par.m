function [cfg,flag] = runROMT_par(cfg)

reInitializeU = 1; %1 if reinitialize u to 0 before each time step; 0 if not, unless first time step

if ~exist(cfg.out_dir,'dir')
    mkdir(cfg.out_dir)
end

fname = sprintf('%s/record.txt',cfg.out_dir);
if ~exist(sprintf('%s/record.txt',cfg.out_dir),'file')
    csvwrite_with_headers(fname,[0 0 0 0 0 0 0 0 0],{'time-ind','ti','tf','phi','mk','Ru','phiN','max(u)','toc'});
end

rho_n = cfg.vol(1).data(:);

if ~exist(sprintf('%s/rho_%s_%d_t_0.mat',cfg.out_dir,cfg.tag,cfg.first_time),'file')
    save(sprintf('%s/rho_%s_%d_t_0.mat',cfg.out_dir,cfg.tag,cfg.first_time),'rho_n');
end

fprintf('\n =============== rOMT Starts ===============\n')
fprintf('______________________________________________\n\n')
fprintf(' tag:\t\t%s\n dataset:\t%s\n sigma:\t\t%.4f\n gamma:\t\t%.4f\n beta:\t\t%.4f\n nt:\t\t%d\n dt:\t\t%.2f\n pcg:\t\t%d\n',cfg.tag,cfg.dataset_name,cfg.sigma,cfg.gamma,cfg.beta,cfg.nt,cfg.dt,cfg.niter_pcg)
fprintf(' size:\t\t%s\n do_resize:\t%d\n resize_factor:\t%.2f\n start frame:\t%d\n end frame:\t%d\n frame jump:\t%d\n\n\n',sprintf('%d ', cfg.true_size),cfg.do_resize,cfg.size_factor,cfg.first_time,cfg.last_time+cfg.time_jump,cfg.time_jump)
%%
%{
tag = cfg.tag;
data_dir = cfg.data_dir;
first_time = cfg.first_time;
time_jump = cfg.time_jump;
extension = cfg.extension;
x_range = cfg.x_range;
y_range = cfg.y_range;
z_range = cfg.z_range;
do_resize = cfg.do_resize;
size_factor = cfg.size_factor;
smooth = cfg.smooth;
true_size = cfg.true_size;
nt = cfg.nt;
dt = cfg.dt;
sigma = cfg.sigma;
add_source = cfg.add_source;
gamma = cfg.gamma;
beta = cfg.beta;
niter_pcg = cfg.niter_pcg;
dTri = cfg.dTri;
out_dir = cfg.out_dir;
%}
%%
profile on
clear T; T = 0;
parfor (tind = 1:length(cfg.first_time:cfg.time_jump:cfg.last_time),2)
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
    rho_0 = cfg.vol(tind).data(:);

    %true final density
    par = paramInitFunc(cfg.true_size',cfg.nt,cfg.dt,cfg.sigma,cfg.add_source,cfg.gamma,cfg.beta,cfg.niter_pcg,cfg.dTri);
    par.drhoN     = cfg.vol(tind+1).data(:);
    
    %par = paramInitFunc(true_size',nt,dt,sigma,add_source,gamma,beta,niter_pcg,dTri);
    %par.drhoN     = rho_N(:);
    
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
    
    save_un(sprintf('%s/u0_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.tag,cfg.first_time+(tind-1)*cfg.time_jump,cfg.first_time+(tind-1)*cfg.time_jump+cfg.time_jump,tind),u);
    save_rhon(sprintf('%s/rhoNe_%s_%d_%d_t_%d.mat',cfg.out_dir,cfg.tag,cfg.first_time+(tind-1)*cfg.time_jump,cfg.first_time+(tind-1)*cfg.time_jump+cfg.time_jump,tind),rho_n);
    
    fprintf('tind = %d, max(u) = %5.4f\n',tind,max(u));
end
fprintf('\n =============== rOMT Ends ===============\n')
fprintf('\n Elapsed Time: %s\n',datestr(seconds(T),'HH:MM:SS'))

profile viewer
profile off
flag = 1;
end