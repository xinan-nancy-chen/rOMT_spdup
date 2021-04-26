%function[phi,mk,phi0,phiN] =  get_phi(rho0,u,nt,dt,par)
function[phi,mk,phiN,rho,Ru] =  get_phi(rho0,u,nt,dt,par)

%% 
rho   = advecDiff(rho0,u,nt,dt,par);

mk    = MKdist(u,rho,nt,dt, par); %=hd*dt*rho'*||v||^2

phiN  = 0.5*norm(rho(:,end) - par.drhoN)^2;

%% smoothing deformation field
%
Ru = 0;

if par.gamma~=0
    uvec=vec2mat(u(:),par.dim*nt);
    %GTG=par.Grad'*par.Grad;
    for ii = 1:par.dim*nt
        %Ru=Ru+0.5*par.hd*dt*(uvec(:,ii)'*GTG*uvec(:,ii));
        Ru=Ru+0.5*par.hd*dt*(norm(par.Grad*uvec(:,ii))^2); 
    end
    %}
end

phi   = par.beta*mk + phiN + par.gamma.*Ru;
end


