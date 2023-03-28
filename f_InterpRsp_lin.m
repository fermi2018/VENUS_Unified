function [Rsp_0,dRsp_E,dRsp_H]=f_InterpRsp_lin(vportE,vportH,indQW)


global N2lin P2lin indGan
global Rsp0 dRspE0 dRspH0 

%
%--------------------------------------------------------------------------
% Spontaneous emission recombination rate
%--------------------------------------------------------------------------
%ConvRsp=10^18*10^9; % Conversion factor to 1/(s*cm^3)



%'in Rsp', keyboard

ind=indGan{indQW};
dE=vportE-N2lin(1,ind);
dH=vportH-P2lin(1,ind);

dRsp_E=dRspE0(1,ind);
dRsp_H=dRspH0(1,ind);
Rsp_0=Rsp0(1,ind)+dRsp_E.*dE+dRsp_H.*dH;


