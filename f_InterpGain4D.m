function [g,dgE,dgH,rsp,drspE,drspH,Dn]=f_InterpGain4D(vportE,vportH,lambda_i,T)
% function [g,dg,rsp,drsp,Rsp,dRsp]=f_InterpGain(vportE,vportH,lambda_i,T)
%
% function [g,dg,rsp,drsp]=f_InterpGain(vportE,vportH,lambda_i,T)
%
% This function interpolates the look-up tables from Debernardi (also
% published in Michalzik, book chapter 3) and provides the stimulated
% emission gain g, the spontaneous emission rsp, and the corresponding
% derivatives to be inserted in the generalized Newton method jacobian
% matrix.
%
%--------------------------------------------------------------------------
% Inputs
% vportE,vportH: carrier concentration N=(n+p)/2 (1/cm^3)
% lambda_i: optical mode wavelength (nm)
% T: temperature (K)
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Outputs
% g: stimulated emission gain (1/s)
% dg: derivative of g with respect to n and/or p (mind the "/2"!)
% rsp: spontaneous emission gain (1/s)
% drsp: derivative of rsp with respect to n and/or p (mind the "/2"!)
% Rsp: spontaneous emission recombination rate (1/(s*cm^3))
% dRsp: derivative of Rsp with respect to n and/or p (mind the "/2"!)
%--------------------------------------------------------------------------
%
% Alberto Tibaldi, v1.0 30/11/2016
%                  v1.1 02/12/2016

% %-- Rsp from book Michalzik, (3.29)

global Gtot Rtotm DeltaN
global LAV PO_E PO_H TV
global dGdH dRdH dRspdH 
global dGdE dRdE dRspdE

mip=(min(squeeze(PO_E(1,:,1,1))));
fi=find(vportE<=mip);
vportE(fi)=mip*1.001;
fi=find(vportH<=mip);
vportH(fi)=mip*1.001;

map=(max(squeeze(PO_E(1,:,1,1))));
fi=find(vportE>=map);
if length(fi>0) 
 '!!!! Warning: elettroni fuori range '
end
vportE(fi)=map*.95;
map=(max(squeeze(PO_H(:,1,1,1))));
fi=find(vportH>=map);
if length(fi>0) 
 '!!!! Warning: lacune fuori range '
end
vportH(fi)=map*.95;

map=(max(squeeze(TV(1,1,1,:))));
fi=find(T>=map);
if length(fi>0) 
 '!!!! Warning: Temperatura fuori range '
end
T(fi)=map*.95;


%'fine caricamento: ora calcolo', keyboard
%g=interpn(PO_H,PO_E,LAV,TV,real(Gtot),vportH,vportE,lambda_i,T,'linear');
tipo='spline';
%tipo='linear';
g=interpn(PO_H,PO_E,LAV,TV,Gtot,vportH,vportE,lambda_i,T,tipo);
rsp=interpn(PO_H,PO_E,LAV,TV,Rtotm,vportH,vportE,lambda_i,T,tipo);
Dn=interpn(PO_H,PO_E,LAV,TV,DeltaN,vportH,vportE,lambda_i,T,tipo);
%figure, plot(g),
%'dopo', keyboard


tipo='linear';
tipon='spline';
LamGrid=squeeze(LAV(1,1,:,1));
if max(LamGrid)<max(lambda_i) | min(LamGrid)>max(lambda_i)
    'LUT fuori Range lambda'
end    
dgE=interpn(PO_H,PO_E,LAV,TV,dGdE,vportH,vportE,lambda_i,T,tipo);
fin=find(isnan(dgE)==1);
if length(fin)>0
 dgE(fin)=interpn(PO_H,PO_E,LAV,TV,dGdE,vportH(fin),vportE(fin),lambda_i(fin),T(fin),tipon);
% 'dgw',  keyboard
end

dgH=interpn(PO_H,PO_E,LAV,TV,dGdH,vportH,vportE,lambda_i,T,tipo);
fin=find(isnan(dgH)==1);
if length(fin)>0
 dgH(fin)=interpn(PO_H,PO_E,LAV,TV,dGdH,vportH(fin),vportE(fin),lambda_i(fin),T(fin),tipon);
% 'dgw',  keyboard
end


drspE=interpn(PO_H,PO_E,LAV,TV,dRdE,vportH,vportE,lambda_i,T,tipo);
fin=find(isnan(drspE)==1);
if length(fin)>0
 drspE(fin)=interpn(PO_H,PO_E,LAV,TV,dRdE,vportH(fin),vportE(fin),lambda_i(fin),T(fin),tipon);
end
drspH=interpn(PO_H,PO_E,LAV,TV,dRdH,vportH,vportE,lambda_i,T,tipo);
fin=find(isnan(drspH)==1);
if length(fin)>0
 drspH(fin)=interpn(PO_H,PO_E,LAV,TV,dRdH,vportH(fin),vportE(fin),lambda_i(fin),T(fin),tipon);
end
%'dopo Inter', keyboard