function [Jmat_num]=numJacobi(mesh,mode,r)
%
nn=mesh.nn;
pp=2*nn;
% 
uvet = mode.u_num;
rvet = r;
%
delta=1e-4;
uvet_inc = repmat(uvet,1,length(uvet));
Jmat_num=sparse(length(uvet),length(uvet));

for ind = 1:length(uvet)
    uvet_inc(ind,ind)=uvet_inc(ind,ind)*(1+delta);
end

% parfor indJ = 1:length(uvet)-2
for indJ = 1:length(uvet)-2

    uvet = uvet_inc(:,indJ);
    [rvet_inc]=evalMatrices(mesh,mode,uvet);
    Jmat_num(:,indJ)=(rvet_inc-rvet)./(delta*mode.u_num(indJ));
end

% Boundary conditions: Poisson equation
Jmat_num(1,:) = 0;
Jmat_num(nn,:) = 0;
Jmat_num(1,1) = 1; % line contact on the left
Jmat_num(nn,nn) = 1; % line contact on the left
%
% Electron equation boundary condFitions: zero charge at contacts
Jmat_num(nn+1,:) = 0;
Jmat_num(nn+nn,:) = 0;
Jmat_num(nn+1,nn+1) = 1; % line contact on the left!
Jmat_num(nn+nn,nn+nn) = 1; % ground contact on the right!
%
% Hole equation boundary conditions: zero charge at contacts
Jmat_num(pp+1,:) = 0;
Jmat_num(pp+nn,:) = 0;
Jmat_num(pp+1,pp+1) = 1; % line contact on the left!
Jmat_num(pp+nn,pp+nn) = 1; % ground contact on the right!
%
% Additional BCs coming from voltage line and impedences
Jmat_num(nn,nn+pp+1) = -1;
Jmat_num(nn+pp+2,nn+pp+2) = -1;
Jmat_num(nn+pp+1,nn+pp+1) = 1;
%

end