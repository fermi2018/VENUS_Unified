%************************************************************************76
function [Dn,Dp,dDndphi1,dDndphi2,dDndphi3,dDpdphi1,dDpdphi2,dDpdphi3] = ...
   mobility(mode,mobn0,mobp0,Se,r1x,r2x,r3x,r1y,r2y,r3y)
%
%========================================================================76
switch(mode.mobility)
case{'none'}, % low-field diffusivity
nt=length(Se);
dDndphi1=zeros(1,nt); dDndphi2=zeros(1,nt); dDndphi3=zeros(1,nt);
dDpdphi1=zeros(1,nt); dDpdphi2=zeros(1,nt); dDpdphi3=zeros(1,nt);
Dn=mobn0.*mode.Vt_tr; Dp=mobp0.*mode.Vt_tr;
%========================================================================76
case{'gaas_x','gaas_y'}
%
if(strcmp(mode.mobility,'gaas_x')),  
E=mode.efieldx; % driving force: electric field along x
% electric field derivatives
dEdphi1=-1./(2*Se).*r1y.*sign(E)/mode.scl;
dEdphi2=-1./(2*Se).*r2y.*sign(E)/mode.scl;
dEdphi3=-1./(2*Se).*r3y.*sign(E)/mode.scl;
elseif(strcmp(mode.mobility,'gaas_y')), 
E=mode.efieldy; % driving force: electric field along y
% electric field derivatives
dEdphi1=-1./(2*Se).*r1x.*sign(E)/mode.scl;
dEdphi2=-1./(2*Se).*r2x.*sign(E)/mode.scl;
dEdphi3=-1./(2*Se).*r3x.*sign(E)/mode.scl; end 
%
E=abs(E); % we have already included the sign in dEdphi
%vsatn=9.7565e6; % GaAs elec saturation velocity, cm/s, Dessis Table 11.77
vsatn=1.0099e7; % cm/s, Perez
vsatp=8.37e6; % GaAs hole saturation velocity, cm/s, Dessis Table 11.77
%E0=4.6461e3; % V/cm
E0=4.3653e3; % V/cm, Perez
mobn = (mobn0 + vsatn.*(E.^3)/(E0.^4))./(1 + (E/E0).^4);
mobp = (mobp0 + vsatp.*(E.^3)/(E0.^4))./(1 + (E/E0).^4);
dmobndE = -E.^2.*(-3.*vsatn.*E0.^4 + vsatn.*E.^4 + 4.*E.*mobn0.*E0.^4)./ ...
          (E0.^4 + E.^4).^2;
dmobpdE = -E.^2.*(-3.*vsatp.*E0.^4 + vsatp.*E.^4 + 4.*E.*mobp0.*E0.^4)./ ...
          (E0.^4 + E.^4).^2;
dDndphi1=mode.Vt_tr.*dmobndE.*dEdphi1; 
dDndphi2=mode.Vt_tr.*dmobndE.*dEdphi2; 
dDndphi3=mode.Vt_tr.*dmobndE.*dEdphi3;
dDpdphi1=mode.Vt_tr.*dmobpdE.*dEdphi1; 
dDpdphi2=mode.Vt_tr.*dmobpdE.*dEdphi2; 
dDpdphi3=mode.Vt_tr.*dmobpdE.*dEdphi3;
Dn=mode.Vt_tr.*mobn; Dp=mode.Vt_tr.*mobp;
%========================================================================76
case{'caughey_thomas_x','caughey_thomas_y'}
%
if(strcmp(mode.mobility,'caughey_thomas_x')),  
E=mode.efieldx; % driving force: electric field along x
% electric field derivatives
dEdphi1=-1./(2*Se).*r1y.*sign(E)/mode.scl;
dEdphi2=-1./(2*Se).*r2y.*sign(E)/mode.scl;
dEdphi3=-1./(2*Se).*r3y.*sign(E)/mode.scl;
elseif(strcmp(mode.mobility,'caughey_thomas_y')), 
E=mode.efieldy; % driving force: electric field along y
% electric field derivatives
dEdphi1=-1./(2*Se).*r1x.*sign(E)/mode.scl;
dEdphi2=-1./(2*Se).*r2x.*sign(E)/mode.scl;
dEdphi3=-1./(2*Se).*r3x.*sign(E)/mode.scl; end 
%
E=abs(E); % we have already included the sign in dEdphi
betan=2; betap=1;
%vsatn=1.07e7; % Si saturation velocity, cm/s, Dessis Table 11.76
vsatn=1.4e7; % GaAs saturation velocity, cm/s, Dessis Table 11.76
vsatp=8.37e6; % Si saturation velocity, cm/s, Dessis Table 11.76
mobn=mobn0./((1+(mobn0.*E/vsatn).^betan).^(1/betan));
mobp=mobp0./((1+(mobp0.*E/vsatp).^betap).^(1/betap));
dmobndE=-mobn0.*(1+(mobn0.*E./vsatn).^betan).^((-1-betan)./betan).* ...
   (mobn0./vsatn).^betan.*E.^(betan-1);
dmobpdE=-mobp0.*(1+(mobp0.*E./vsatp).^betap).^((-1-betap)./betap).* ...
   (mobp0./vsatp).^betap.*E.^(betap-1);
dDndphi1=mode.Vt_tr.*dmobndE.*dEdphi1; 
dDndphi2=mode.Vt_tr.*dmobndE.*dEdphi2; 
dDndphi3=mode.Vt_tr.*dmobndE.*dEdphi3;
dDpdphi1=mode.Vt_tr.*dmobpdE.*dEdphi1; 
dDpdphi2=mode.Vt_tr.*dmobpdE.*dEdphi2; 
dDpdphi3=mode.Vt_tr.*dmobpdE.*dEdphi3;
Dn=mode.Vt_tr.*mobn; Dp=mode.Vt_tr.*mobp;
%========================================================================76
otherwise, error('which mobility model!?'), end
%************************************************************************76
% 
% E=sym('E','real'); 
% vsatn=sym('vsatn','real'); E0=sym('E0','real'); mobn0=sym('mobn0','real');
% mobn = (mobn0 + (vsatn./E)*(E/E0)^4)/(1 + (E/E0)^4);
% dmobndE=(simplify(diff(mobn,'E')));
% vectorize(dmobndE)


