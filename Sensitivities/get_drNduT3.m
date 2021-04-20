function [drNduT] = get_drNduT3(M,nt,dt,par,y)
%% Sensitivity of rho(:,end) w.r.t 'u' transpose times a vector x
%  y... vector length(prod(n),1)

% rho(:,0)  -----------> rho(:,1)  --------------> rho(:,2)
%             u(:,1)                   u(:,2)
%               S1                      S2
if par.dim==2
%% 2D verison
n                  = par.n;
%u                  = reshape(u,2*prod(n),nt);
%rho                = advecDiff(rho_0,u,nt,dt, par);
%Mdis               = - par.sigma*par.Grad'*par.Grad;
%I                  = speye(prod(n));
%B                  = (I - dt*Mdis);

sensTx             = zeros(2*prod(n),nt);
sens               = y;
%Rho                = zeros(prod(n),nt+1);
%Rho(:,1)           = rho_0; Rho(:,2:end) = rho;

for i = nt:-1:1
    %U1 = reshape(u(1:prod(n),i),n');
    %U2 = reshape(u(prod(n)+1:end,i),n');
    
    %[S,Tx,Ty]  = dTrilinears(Rho(:,i),par.Xc + dt*U1, par.Yc + dt*U2,...
    %                 par.h1(1),par.h2(1));
    
    % Sensitivity:          
    [sensI,pcgflag1]   = pcg(par.B',sens);
    if pcgflag1 ~= 0
        warning('MATLAB:pcgExitFlag','Warning: get_drNduT3.m >>> while finding drho_%dT/du, pcg exit flag1 = %d',i,pcgflag1)
    end
    sensTx(:,i)              = dt*[M.Tx{i},M.Ty{i}]'*sensI;
    
    if  i>1 
        sens = M.S{i}'*sensI;
    end
end

drNduT = sensTx(:);

elseif par.dim==3
%% 3D verison
%%
n                  = par.n;
%u                  = reshape(u,3*prod(n),nt);
%rho                = advecDiff(rho_0,u,nt,dt, par);
%Mdis               = - par.sigma*par.Grad'*par.Grad;
%I                  = speye(prod(n));
%B                  = (I - dt*Mdis);

sensTx             = zeros(3*prod(n),nt);
%sens               = repmat(y,1,nt);
sens               = y;
%Rho                = zeros(prod(n),nt+1);
%Rho(:,1)           = rho_0; Rho(:,2:end) = rho;

for i = nt:-1:1
    %U1 = reshape(u(1:prod(n),i),n');
    %U2 = reshape(u(prod(n)+1:2*prod(n),i),n');
    %U3 = reshape(u(2*prod(n)+1:end,i),n');
    
    %[S,Tx,Ty,Tz]  = dTrilinears(Rho(:,i),par.Xc + dt*U1, par.Yc + dt*U2, par.Zc + dt*U3,...
    %                 par.h1(1),par.h2(1),par.h3(1));
    
    % Sensitivity:          
    [sensI,pcgflag1]   = pcg(par.B',sens);
    if pcgflag1 ~= 0
        warning('MATLAB:pcgExitFlag','Warning: get_drNduT3.m >>> while finding drho_%dT/du, pcg exit flag1 = %d',i,pcgflag1)
    end
    sensTx(:,i)              = dt*[M.Tx{i},M.Ty{i},M.Tz{i}]'*sensI;
    
    if  i>1 
        sens = M.S{i}'*sensI;%S'*sensi;
    end
end

drNduT = sensTx(:);
    
    
else
    warning('In get_drNduT3.m: dimension of data should be either 2 or 3')
end
end

