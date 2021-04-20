function [u,phi,g] = GNblock_u(Rho_i,u,nt,dt,par,tag_str)
%%
if par.dim==2
%% 2D verison
%% Note: Changes made by Xinan Chen in Apr 21: (modified 3D version in this project too)
% (1) in GNblock_u.m: move "U = reshape(u,2*prod(par.n),[]);" and "dmk2 = zeros(2*prod(par.n),nt);" into the loop
% (2) in get_dRudu.m: remove negative sign computing G
% (3) in get_drNduT.m: remove unnecessary variable sensi to avoid computing pcg twice
% (4) add new function get_drNduTdrNdu.m to combine
% get_drNduT(Rho_i,u,nt,dt,par,get_drNdu(Rho_i,u,nt,dt,par,x)) togeter in
% computing H
% (5) add new version in dTrilinear.m with equivalent calculation in math,
% even though no obvious running improvement has shown
% (6) add Rho as input in get_drNduTdrNdu.m and get_drNduT2.m, to avoid repeated computing of
% Rho = advecDiff.m
% (7) add new function get_drNduTdrNdu2.m and get_drNduT3.m to store S,Tx,Ty,Tz
% (8) combine 2d and 3d together. Modified functions include: GNblock_u.m,
% paramInitFunc.m, advecDiff.m, dTrilinear2d.m (new),  dTrilinear3d.m
% (new), get_drNduT3.m, get_drNduTdrNdu2.m
% (9) add new input dTri in paramInitFunc.m to pass type of boundary condition to
% dTrilinears2d.m and dTrilinears3d.m
% (10) in GNblock_u.m, precompute part of g and H
% (11) in paramInitFunc.m, precompute par.B = I-dt*Mdis and save to par,
% resulting in modifications in get_drNduT3.m, get_drNduTdrNdu2.m and
% advecDiff.m

if nargin < 6
    tag_str = '';
end

phi          = get_phi(Rho_i,u,nt,dt,par);
A            = kron(ones(1,2),speye(prod(par.n)));
Abig         = kron(speye(nt),A);
flag         = 0;
%U            = reshape(u,2*prod(par.n),[]);
%dmk2         = zeros(2*prod(par.n),nt);

%% Loop
for i = 1:par.maxUiter
    Rho      = advecDiff(Rho_i,u,nt,dt,par);
    U            = reshape(u,2*prod(par.n),[]);
    dmk2         = zeros(2*prod(par.n),nt);
    
    RHO0 = [Rho_i,Rho(:,1:end-1)];
    for k = 1:nt
        U1 = reshape(U(1:prod(par.n),k),par.n');
        U2 = reshape(U(prod(par.n)+1:end,k),par.n');

        [M.S{k},M.Tx{k},M.Ty{k}]  = dTrilinears2d(RHO0(:,k),par.Xc + dt*U1, par.Yc + dt*U2,...
                         par.h1(1),par.h2(1),par.bc);
    end
            
    for      j = 1:nt
             dmk2(:,1:j) = dmk2(:,1:j) + reshape(par.hd*dt*get_drNduT3(M,j,dt,par,A*(U(:,j).*U(:,j))),2*prod(par.n),j);
    end
    
    g        = (par.beta*2*par.hd*dt*Rho(:)'*Abig*sdiag(u(:)))' + par.beta*dmk2(:) + ...
                get_drNduT3(M,nt,dt,par,Rho(:,end) - par.drhoN) + par.gamma*dt*par.hd*get_dRudu(u,nt,par)';

    fprintf('%3d.%d\t      %3.2e \t     ||g|| = %3.2e\n',i,0,phi,norm(g));

    H13       = par.beta*2*dt*par.hd*sdiag(Rho(:)'*Abig) + par.gamma*dt*par.hd.*kron(speye(nt*2),(par.Grad)'*par.Grad);
    H        = @(x) H13*x + get_drNduTdrNdu2(M,nt,dt,par,x);    
                       
    %H        = @(x) par.beta*2*dt*par.hd*(Rho(:)'*Abig*sdiag(x))' + ...
    %           get_drNduTdrNdu2(M,nt,dt,par,x)+...%get_drNduT(Rho_i,u,nt,dt,par,get_drNdu(Rho_i,u,nt,dt,par,x))+...
    %           par.gamma*dt*par.hd.*kron(speye(nt*2),par.Grad'*par.Grad)*x; 
 
    
    [s,pcgflag,relres,iter]    = pcg(H,-g,0.01,par.niter_pcg);
    
    if pcgflag ~= 0
      warning('MATLAB:pcgExitFlag','Warning: GNblock_u.m >>> iter %d, while finding s, pcg exit flag = %d \nrelres = %3.2e, iter = %d, %s',i,pcgflag,relres,iter,tag_str)
    end
    
    muls     = 0.7; 
    lsiter = 1;
    while 1
        ut   = u(:) + muls*s;
        
        phit = get_phi(Rho_i,ut,nt,dt,par);
        
        fprintf('%3d.%d\t      %3.2e \t     phit  = %3.2e        %s\n',i,lsiter,phi,phit,tag_str);
        
        % test for line search termination
        if phit < phi + 1e-8*s'*g%1e-8*s'*g
            break;                      %breaks while loop entirely (and goes to next statement after end of while loop)
        end
        muls = muls/2; lsiter = lsiter+1;
        
        % fail if lsiter is too large
        if lsiter > 4
            fprintf('LSB\n');
            ut = u;     
            flag = 1;    
            return;                     % returns and exits function
        end
    end
    
    if flag
        return; 
    end                % returns and exits function
    u   = ut;
    phi = phit;
    
end

elseif par.dim==3
%% 3D verison
%%
if nargin < 6
    tag_str = '';
end

phi          = get_phi(Rho_i,u,nt,dt,par);
A            = kron(ones(1,3),speye(prod(par.n)));
Abig         = kron(speye(nt),A);
flag         = 0;
%U            = reshape(u,3*prod(par.n),[]);
%dmk2         = zeros(3*prod(par.n),nt);

%% Loop
for i = 1:par.maxUiter
    Rho      = advecDiff(Rho_i,u,nt,dt,par);
    U            = reshape(u,3*prod(par.n),[]);
    dmk2         = zeros(3*prod(par.n),nt);
    
    RHO0 = [Rho_i,Rho(:,1:end-1)];
    for k = 1:nt
        U1 = reshape(U(1:prod(par.n),k),par.n');
        U2 = reshape(U(prod(par.n)+1:2*prod(par.n),k),par.n');
        U3 = reshape(U(2*prod(par.n)+1:end,k),par.n');

        [M.S{k},M.Tx{k},M.Ty{k},M.Tz{k}]  = dTrilinears3d(RHO0(:,k),par.Xc + dt*U1, par.Yc + dt*U2, par.Zc + dt*U3,...
                         par.h1(1),par.h2(1),par.h3(1),par.bc);
    end
    
    for      j = 1:nt
             dmk2(:,1:j) = dmk2(:,1:j) + reshape(par.hd*dt*get_drNduT3(M,j,dt,par,A*(U(:,j).*U(:,j))),3*prod(par.n),j);
    end
    
    g        = (par.beta*2*par.hd*dt*Rho(:)'*Abig*sdiag(u(:)))' + par.beta*dmk2(:) + ...
                get_drNduT3(M,nt,dt,par,Rho(:,end) - par.drhoN) + par.gamma*dt*par.hd*get_dRudu(u,nt,par)';
    
    fprintf('%3d.%d\t      %3.2e \t     ||g|| = %3.2e       %s\n',i,0,phi,norm(g),tag_str);
    
    H13       = par.beta*2*dt*par.hd*sdiag(Rho(:)'*Abig) + par.gamma*dt*par.hd.*kron(speye(nt*3),(par.Grad)'*par.Grad);
    H        = @(x) H13*x + get_drNduTdrNdu2(M,nt,dt,par,x);  
    
    %H        = @(x) par.beta*2*dt*par.hd*(Rho(:)'*Abig*sdiag(x))' + ...
    %           get_drNduTdrNdu2(M,nt,dt,par,x)+...%get_drNduT(Rho_i,u,nt,dt,par,get_drNdu(Rho_i,u,nt,dt,par,x))+...
    %           par.gamma*dt*par.hd.*kron(speye(nt*3),par.Grad'*par.Grad)*x;      
 
    
    [s,pcgflag,relres,iter]    = pcg(H,-g,0.01,par.niter_pcg);
    
    if pcgflag ~= 0
      warning('MATLAB:pcgExitFlag','Warning: GNblock_u.m >>> iter %d, while finding s, pcg exit flag = %d \nrelres = %3.2e, iter = %d, %s',i,pcgflag,relres,iter,tag_str)
    end
    
    muls     = 0.7; lsiter = 1;
    while 1
        ut   = u(:) + muls*s;
        
        phit = get_phi(Rho_i,ut,nt,dt,par);
        
        fprintf('%3d.%d\t      %3.2e \t     phit  = %3.2e        %s\n',i,lsiter,phi,phit,tag_str);
        
        % test for line search termination
        if phit < phi + 1e-8*s'*g
            break;                      %breaks while loop entirely (and goes to next statement after end of while loop)
        end
        muls = muls/2; lsiter = lsiter+1;
        
        % fail if lsiter is too large
        if lsiter > 4
            fprintf('LSB\n');
            ut = u;     
            flag = 1;    
            return;                     % returns and exits function
        end
    end
    
    if flag
        return; 
    end                % returns and exits function
    u   = ut;
    phi = phit;
    
end
    
else
    warning('In GNblock_u.m: dimension of data should be either 2 or 3')
end
end


