phi=uvet(1:nn);
iq = mesh.iq;

    mode.ecb(iq) = - phi(iq) + mesh.Eg(iq)/2 - mode.Vt(iq)./2.*log(mesh.Nv(iq)./mesh.Nc(iq)) + mesh.phi_r(iq);
    mode.evb(iq) = mode.ecb(iq) - mesh.Eg(iq);

elec=abs(uvet(nn+(1:nn)))*mode.CarrierNorm;
if((isfield(mode,'stats'))&&(strcmp(mode.stats,'Fermi')))
    tmp = elec(iq)./mesh.Nc(iq)*mode.CarrierNorm;
    nferinv = invferdr(tmp,mode.tolinvferdr);
    EFn(iq) = mode.ecb(iq) + Vt(iq).*nferinv; 
end
%
hole=abs(uvet(pp+(1:nn)))*mode.CarrierNorm;
if((isfield(mode,'stats'))&&(strcmp(mode.stats,'Fermi')))
    tmp = hole(iq)./mesh.Nv(iq)*mode.CarrierNorm;
    nferinv = invferdr(tmp,mode.tolinvferdr);
    EFp(iq) = mode.evb(iq) - Vt(iq).*nferinv; 
end
%

% plotvenus