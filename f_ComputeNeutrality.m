
function [mode] = f_ComputeNeutrality(geom,mesh,mode)

Vt=mode.Vt;

% Update neutrality conditions for Fermi-Dirac and/or Incomplete
% ionization
for ic=1:geom.nc % cycle on number of contacts
    iter=0;
    deltaphi=inf; rho=inf;
    ii=(mesh.contact==ic); jj=find(ii&mesh.iq);
    phi=mode.phi_neutr(jj);
    [phiU,ia,ib]=unique(phi);
    while(abs(deltaphi)>mode.tolconv_neutr)
        [rho,drho]=f_NeutralityEquation(phiU,mode,mesh,jj(ia));
        deltaphi=-rho./drho;
        phiU=phiU+deltaphi;
        iter=iter+1;
    end
    phi=phiU(ib);
    mode.phi_neutr(jj)=phi;
    mode.phi(jj)=phi;
    ecb = - mode.phi_neutr(jj) + mesh.Eg(jj)/2 - ...
        Vt(jj)./2.*log(mesh.Nv(jj)/mesh.Nc(jj)) + mesh.phi_r(jj);
    evb = ecb - mesh.Eg(jj);
    mode.ecb=ecb; mode.evb=evb;
    if((isfield(mode,'stats'))&&(strcmp(mode.stats,'Fermi')))
        mode.elec(jj) = mesh.Nc(jj).*ferdr((+ 0 - ecb)./Vt(jj),1/2); % 1/cm^3;
        mode.hole(jj) = mesh.Nv(jj).*ferdr((- 0 + evb)./Vt(jj),1/2); % 1/cm^3;
    else
        mode.elec(jj) = mesh.Nc(jj).*exp((+ 0 - ecb)./Vt(jj)); % 1/cm^3;
        mode.hole(jj) = mesh.Nv(jj).*exp((- 0 + evb)./Vt(jj)); % 1/cm^3;
    end
end