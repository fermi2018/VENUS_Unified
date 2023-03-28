clc
clear
% close all
colordef white
dbstop if error

addpath('OtticoBar','-end')
addpath('OtticoBar\new17Optica','-end')

load VELM
ifp_Opt=-10;

[St_wave,camzT,cam2_0,gain_0,lamod_0,modc,xro,fian,fPdif_0,fPES_0,fPGA,...
tyPmod_0,omP0_0,eta_eff_0,Tper,Pa,Ppol,Plot,Eo_x,Eo_y]=...
 caloptdd(MODC0,ro_campo,ro_in,fil_str,iLP,ifp_Opt,calop);
