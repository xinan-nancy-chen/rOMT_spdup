function [drNduT] = get_drNduT4(M,nt,dt,par,y)
%% Sensitivity of rho(:,end) w.r.t 'u' transpose times a vector x
%  y... vector length(prod(n),1)

% rho(:,0)  -----------> rho(:,1)  --------------> rho(:,2)
%             u(:,1)                   u(:,2)
%               S1                      S2

%% 
n                  = par.n;
sensTx             = zeros(par.dim*prod(n),nt);
sens               = y;

for i = nt:-1:1
    % Sensitivity:          
    [sensI,pcgflag1]   = pcg(par.B',sens);
    if pcgflag1 ~= 0
        warning('MATLAB:pcgExitFlag','Warning: get_drNduT4.m >>> while finding drho_%dT/du, pcg exit flag1 = %d',i,pcgflag1)
    end
    if par.dim==2
        y2 = [M.Tx{i},M.Ty{i}]'*(dt*sensI);
    elseif par.dim==3
        y2 = [M.Tx{i},M.Ty{i},M.Tz{i}]'*(dt*sensI);
    end
    sensTx(:,i)              = y2;
    
    if  i>1 
        sens = M.S{i}'*sensI;
    end
end

drNduT = sensTx(:);
end

