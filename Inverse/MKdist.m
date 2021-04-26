function[mk] = MKdist(u,rho,nt,dt, par)
%[I] = MKdist(u,rho,n,nt,h,dt)
%

%% 
% build averaging matrix
n   = par.n;

A   = kron(ones(1,par.dim),speye(prod(n)));

U   = reshape(u,par.dim*prod(n),nt);
R   = reshape(rho,prod(n),nt);

mk = 0;


for i=1:nt
    mk = mk + par.hd*dt*R(:,i)'*A*(U(:,i).*U(:,i));
end
end
