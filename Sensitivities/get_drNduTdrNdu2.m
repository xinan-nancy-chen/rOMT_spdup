function [drNduTdrNdu] = get_drNduTdrNdu2(M,nt,dt,par,x)
%% Created by Xinan Chen on 04/19/21
%% This function combines get_drNduT.m and get_drNdu.m to compute the second term in Hessian
%% Sensitivity of rho(:,end) transpose w.r.t 'u' * rho(:,end) w.r.t 'u' * full vector

% rho(:,0)  -----------> rho(:,1)  --------------> rho(:,2)
%             u(:,1)                   u(:,2)

if par.dim==2
%% 2D verison

n                  = par.n;
%u                  = reshape(u,2*prod(n),nt);
%rho                = zeros(prod(n),nt+1);
%rho(:,1)           = rho_0;
%rho(:,2:end)       = Rho;%advecDiff(rho_0,u,nt,dt,par);
%Mdis               = -par.sigma*par.Grad'*par.Grad;
%I                  = speye(prod(n));
%B                  = (I - dt*Mdis);

X                  = reshape(x,prod(n)*2,nt);
sensx              = zeros(prod(n),nt);

% first part:  rho(:,end) w.r.t 'u' * full vector

for i = 1:nt
    %U1 = reshape(u(1:prod(n),i),n');
    %U2 = reshape(u(prod(n)+1:end,i),n');
    
    %[S,Tx,Ty]  = dTrilinears(rho(:,i),par.Xc + dt*U1, par.Yc + dt*U2,...
    %    par.h1(1),par.h2(1));
    
    
%     rho(:,i+1)  = S*rho(:,i);
%     [rho(:,i+1),pcgflag1]  = pcg((I - dt*Mdis),rho(:,i+1));
%     if pcgflag1 ~= 0
%         warning('MATLAB:pcgExitFlag','Warning: get_drNdu.m >>> while finding drho(:,%d)/du, pcg exit pcgflag1 = %d',i,pcgflag1)
%     end
    
    % Sensitivity:
    if  i>1
        for j = 1:i-1
            [sensx(:,j),pcgflag2] = pcg(par.B,(M.S{i}*sensx(:,j)));
            if pcgflag2 ~= 0
                warning('MATLAB:pcgExitFlag','Warning: get_drNduTdrNdu.m >>> while finding drho(:,n)/du_%d, pcg exit pcgflag2 = %d',j,pcgflag2)
            end
        end
    end
    [sensx(:,i),pcgflag3] = pcg(par.B,(dt*[M.Tx{i},M.Ty{i}]*X(:,i)));
    if pcgflag3 ~= 0
        warning('MATLAB:pcgExitFlag','Warning: get_drNduTdrNdu.m >>> while finding drho(:,%d)/du, pcg exit pcgflag3 = %d',i,pcgflag3)
    end
end

sens = sum(sensx,2);


% second part: rho(:,end) transpose w.r.t 'u' * full vector

sensTx             = zeros(2*prod(n),nt);


for i = nt:-1:1
    %U1 = reshape(u(1:prod(n),i),n');
    %U2 = reshape(u(prod(n)+1:end,i),n');
    
    %[S,Tx,Ty]  = dTrilinears(rho(:,i),par.Xc + dt*U1, par.Yc + dt*U2,...
    %    par.h1(1),par.h2(1));
    
    % Sensitivity:          
    [sensI,pcgflag1]   = pcg(par.B',sens);
    if pcgflag1 ~= 0
        warning('MATLAB:pcgExitFlag','Warning: get_drNduTdrNdu.m >>> while finding drho_%dT/du, pcg exit flag1 = %d',i,pcgflag1)
    end
    sensTx(:,i)              = dt*[M.Tx{i},M.Ty{i}]'*sensI;
    
    if  i>1 
    sens = M.S{i}'*sensI;
    end
end

drNduTdrNdu = sensTx(:);

elseif par.dim==3
%% 3D verison
%%
n                  = par.n;
%u                  = reshape(u,3*prod(n),nt);
%rho                = zeros(prod(n),nt+1);
%rho(:,1)           = rho_0;
%rho(:,2:end)       = Rho;%advecDiff(rho_0,u,nt,dt,par);
%Mdis               = -par.sigma*par.Grad'*par.Grad;
%I                  = speye(prod(n));
%B                  = (I - dt*Mdis);

X                  = reshape(x,prod(n)*3,nt);
sensx              = zeros(prod(n),nt);

% first part:  rho(:,end) w.r.t 'u' * full vector

for i = 1:nt
    %U1 = reshape(u(1:prod(n),i),n');
    %U2 = reshape(u(prod(n)+1:2*prod(n),i),n');
    %U3 = reshape(u(2*prod(n)+1:end,i),n');
    
    %[S,Tx,Ty,Tz]  = dTrilinears(rho(:,i),par.Xc + dt*U1, par.Yc + dt*U2, par.Zc + dt*U3,...
    %    par.h1(1),par.h2(1),par.h3(1));
    
    
%     rho(:,i+1)  = S*rho(:,i);
%     [rho(:,i+1),pcgflag1]  = pcg((I - dt*Mdis),rho(:,i+1));
%     if pcgflag1 ~= 0
%         warning('MATLAB:pcgExitFlag','Warning: get_drNdu.m >>> while finding drho(:,%d)/du, pcg exit pcgflag1 = %d',i,pcgflag1)
%     end
    
    % Sensitivity:
    if  i>1
        for j = 1:i-1
            [sensx(:,j),pcgflag2] = pcg(par.B,(M.S{i}*sensx(:,j)));
            if pcgflag2 ~= 0
                warning('MATLAB:pcgExitFlag','Warning: get_drNdu.m >>> while finding drho(:,n)/du_%d, pcg exit pcgflag2 = %d',j,pcgflag2)
            end
        end
    end
    [sensx(:,i),pcgflag3] = pcg(par.B,(dt*[M.Tx{i},M.Ty{i},M.Tz{i}]*X(:,i)));
    if pcgflag3 ~= 0
        warning('MATLAB:pcgExitFlag','Warning: get_drNdu.m >>> while finding drho(:,%d)/du, pcg exit pcgflag3 = %d',i,pcgflag3)
    end
end

sens = sum(sensx,2);


% second part: rho(:,end) transpose w.r.t 'u' * full vector

sensTx             = zeros(3*prod(n),nt);


for i = nt:-1:1
    %U1 = reshape(u(1:prod(n),i),n');
    %U2 = reshape(u(prod(n)+1:2*prod(n),i),n');
    %U3 = reshape(u(2*prod(n)+1:end,i),n');
    
    %[S,Tx,Ty,Tz]  = dTrilinears(rho(:,i),par.Xc + dt*U1, par.Yc + dt*U2, par.Zc + dt*U3,...
    %                 par.h1(1),par.h2(1),par.h3(1));
    
    % Sensitivity:          
    [sensI,pcgflag1]   = pcg(par.B',sens);
    if pcgflag1 ~= 0
        warning('MATLAB:pcgExitFlag','Warning: get_drNduT.m >>> while finding drho_%dT/du, pcg exit flag1 = %d',i,pcgflag1)
    end
    sensTx(:,i)              = dt*[M.Tx{i},M.Ty{i},M.Tz{i}]'*sensI;
    
    if  i>1 
    sens = M.S{i}'*sensI;
    end
end

drNduTdrNdu = sensTx(:);

    
else
    warning('In get_drNduTdrNdu2.m: dimension of data should be either 2 or 3')
end
end

