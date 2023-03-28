function Glut4Dinc(mode)

global Gtot Rtotm Rsp DeltaN
global LAV PO_E PO_H TV TV2 PO_E2 PO_H2
global dGdH dRdH dRspdH 
global dGdE dRdE dRspdE
% global TFasano EFasano xmolFasano alphaFasano
% 
% load(mode.Fasano);
load(mode.GLUT);

if isfield(mode,'CarrierNorm2D')==0
    mode.CarrierNorm2D=1;
end
% subscript "2": for spontaneous emission
por_E2=porE_Rsp/mode.CarrierNorm2D;
por_H2=porH_Rsp/mode.CarrierNorm2D;

% no subscript: for gain
por_E=por_E/mode.CarrierNorm2D;
por_H=por_H/mode.CarrierNorm2D;

% Gain
Gtot=G;
dGdH=dGdH*mode.CarrierNorm2D;
dGdE=dGdE*mode.CarrierNorm2D;

% Total spontaneous emission
Rsp=Ric/mode.CarrierNorm2D;
dRspdH=dRicdH;
dRspdE=dRicdE;

% Spontaneous emission in the lasing mode
Rtotm=Es;
dRdE=dEdE*mode.CarrierNorm2D;
dRdH=dEdH*mode.CarrierNorm2D;


lavnm=lav*1000;
[PO_H,PO_E,LAV,TV]=ndgrid(por_H,por_E,lavnm,Tv);
[PO_H2,PO_E2,TV2]=ndgrid(por_H2,por_E2,Tv);