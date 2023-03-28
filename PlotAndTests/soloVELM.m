clear all
close all
 dbstop if error
colordef black
dbstop if error



addpathVENUS
rmpath OtticoBar/new22Optica
addpath OtticoBar/new23OpticaGR
% rmpath OtticoBar/new23OpticaGR

% load saVELM
% load dati/VELM/VELM_MarkusN_TJ_oxAbove_LgDBR_fixed
load('VELM_GRmetal.mat')

fil_str=[nomeSR,mode.strName,'.str'];

VelmOptions.ivett=0;
%VelmOptions.ipolar=1;      %1 TE
%VelmOptions.ipolar=-1;      %-1 TM, 1 TE
VelmOptions.ipolar=1;      %1 TM, -1 TE; 2 entrambe
%VelmOptions.ipolar=2;      %1 TM, -1 TE; 2 entrambe
%VelmOptions.iBWnew=1

mode.verbVELM=1;
% VelmOptions.Dlam=[-30 30 151 0 .4];
% VelmOptions.Dlam=[-10 0 51 0 .4];
% VelmOptions.Dlam=[-0 1 3 0 .4];
% VelmOptions.NP_k=100;

% VelmOptions.krel_max=.3;

% % VelmOptions.dissfun='diss_SimaZoom';
if VelmOptions.ivett==0
 VelmOptions.iany=0;
end 
%fil_str='C:\Users\Utente\Desktop\VENUSJulian\dati\Stephan_GR0.str'
%fil_str='C:\Users\Utente\Desktop\VENUSJulianOLD\dati\Julian_907047_850nmPD.str'
%ParVet=ParVet(1:3);
%VelmOptions.Pf.iBEL=102;
VelmOptions.Pf.ipolar=VelmOptions.ipolar;
% ParVet=[3.5 9.3/2 22/2 3];
VelmOptions.Pf.nmasce=-3;
% VelmOptions.num_azim=2;

% mode.DT0=200;

pausak
[velm] = f_CallVELM(mesh,mode,mode1,ParVet,VelmOptions,fil_str);



