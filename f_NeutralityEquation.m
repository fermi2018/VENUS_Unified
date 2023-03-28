function [rho,drho] = f_NeutralityEquation(phi,mode,mesh,indContact)
%
qel=1.6021766208e-019; % Elementary charge (C)
kB=1.3806488e-23; % Boltzmann constant (J/K)
Vt = mesh.T(indContact).*kB./qel;
qelNorm=qel*mode.CarrierNorm;
%
phi_r=mesh.phi_r(indContact);
Eg=mesh.Eg(indContact);
Nc=mesh.Nc(indContact);
Nv=mesh.Nv(indContact);
dop_a=mesh.dop_a(indContact);
dop_d=mesh.dop_d(indContact);
%
Ec = - phi + Eg/2 - Vt/2*log(Nv/Nc) + phi_r;
Ev = Ec - Eg;
%
if((isfield(mode,'stats'))&&(strcmp(mode.stats,'Fermi')))
    elec = Nc.*ferdr(-Ec./Vt,1/2); % 1/cm^3
    delec = Nc./Vt.*ferdr(-Ec./Vt,-1/2); % 1/cm^3
    elecB = Nc.*exp(-Ec./Vt); % Boltzmann
    gamma_n = elec./elecB;
    dgamma_n = (1-elec./(delec.*Vt))./elecB; % here, delec is w.r.t. phi
    %
    hole = Nv.*ferdr(Ev./Vt,1/2); % 1/cm^3;
    dhole = - Nv./Vt.*ferdr(Ev./Vt,-1/2); % 1/cm^3;
    holeB = Nv.*exp(Ev./Vt); % Boltzmann
    gamma_p = hole./holeB;
    dgamma_p = (1-hole./(-dhole.*Vt))./holeB; % here, dhole is w.r.t. phi
else % Boltzmann statistics
    elec = Nc.*exp(-Ec./Vt); % 1/cm^3
    delec = elec./Vt;
    gamma_n=1*ones(size(elec));
    dgamma_n=0*ones(size(elec));
    %
    hole = Nv.*exp(Ev./Vt); % 1/cm^3;
    dhole = -hole./Vt;
    gamma_p=ones(size(hole));
    dgamma_p=0*ones(size(hole));
end
%
if((isfield(mode,'ionization'))&&(strcmp(mode.ionization,'Incomplete')))
    DeltaEa=mesh.DeltaEa(indContact);
    DeltaEd=mesh.DeltaEd(indContact);
    gD = 2;
    n1 = Nc.*exp(-DeltaEd/Vt);
    gA = 4;
    p1 = Nv.*exp(-DeltaEa/Vt);
    dop_dp = dop_d./(1+gD*elec./(gamma_n.*n1));
    dop_am = dop_a./(1+gA*hole./(gamma_p.*p1));
    %
    ddop_dp = - dop_d./((1 + gD.*elec./(n1.*gamma_n)).^2) ...
        .*gD.*((n1.*gamma_n - elec.*n1.*dgamma_n)./((n1.*gamma_n).^2));
    ddop_am = - dop_a./((1 + gA.*hole./(p1.*gamma_p)).^2) ...
        .*gA.*((p1.*gamma_p - hole.*p1.*dgamma_p)./((p1.*gamma_p).^2));
else % Full ionization
    dop_dp=dop_d;
    dop_am=dop_a;
    ddop_dp=0;
    ddop_am=0;
end

rho = qelNorm.*(hole - elec + dop_dp - dop_am);
drho = qelNorm.*(dhole - delec + ddop_dp.*delec - ddop_am.*dhole);

return