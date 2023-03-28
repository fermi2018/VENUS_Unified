function [Rsp_0,dRsp_E,dRsp_H]=f_InterpRspEH(vportE,vportH,T)
% function [g,dg,rsp,drsp,Rsp,dRsp]=f_InterpGain(vport,lambda_i,T)
%
% function [g,dg,rsp,drsp]=f_InterpGain(vport,lambda_i,T)
%
% This function interpolates the look-up tables from Debernardi (also
% published in Michalzik, book chapter 3) and provides the stimulated
% emission gain g, the spontaneous emission rsp, and the corresponding
% derivatives to be inserted in the generalized Newton method jacobian
% matrix.
%
%--------------------------------------------------------------------------
% Inputs
% vport: carrier concentration N=(n+p)/2 (1/cm^3)
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
%load(fileName)
global Rsp TV2 PO_E2 PO_H2
global dRspdH 
global dRspdE 

%global dGdH dEsdH dRspdH 
%global dGdE dEsdE dRspdE dDndE 
% load Gres_jo2_AT
% load Gres_jos
% load Gres_gus

% Legend of saved variables
% Gtot: stimulated emission gain (1/ns)
% Rtotm: spontaneous emission gain (1/ns)
% Rsp: spontaneous emission integrated with respect to lambda (??????)
% lav: wavelength grid (um)
% port: carrier concentrations grid ((1/cm^3)/1e18)
% Tv: temperature grid (K)

%-- Normalization factor introduced by Debernardi

%
%--------------------------------------------------------------------------
% Spontaneous emission recombination rate
%--------------------------------------------------------------------------
%ConvRsp=10^18*10^9; % Conversion factor to 1/(s*cm^3)

%'in Rsp', keyboard

mip=(min(squeeze(PO_H2(:,1,1))));
fi=find(vportE<=mip);
vportE(fi)=mip*1.001;
fi=find(vportH<=mip);
vportH(fi)=mip*1.001;


map=(max(squeeze(PO_E2(1,:,1))));
fi=find(vportE>=map);
if length(fi>0) 
 '!!!! Warning: elettroni fuori range '
end
vportE(fi)=map*.95;
map=(max(squeeze(PO_H2(:,1,1))));
fi=find(vportH>=map);
if length(fi>0) 
 '!!!! Warning: lacune fuori range '
end
vportH(fi)=map*.95;

map=(max(squeeze(TV2(1,1,:))));
fi=find(T>=map);
if length(fi>0) 
 '!!!! Warning: Temperatura fuori range '
end
T(fi)=map*.95;

%' Rsp', keyboard

Rsp_0=interpn(PO_H2,PO_E2,TV2,Rsp,vportH,vportE,T,'spline');
%figure, plot(g),

tipo='spline';
tipo='linear';
dRsp_E=interpn(PO_H2,PO_E2,TV2,dRspdE,vportH,vportE,T,tipo);
dRsp_H=interpn(PO_H2,PO_E2,TV2,dRspdH,vportH,vportE,T,tipo);

%'rsp', keyboard