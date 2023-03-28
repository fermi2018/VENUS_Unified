close all
clear
colordef black
% dbstop if error

addpathVENUS

% rmpath('OtticoBar\new22Optica')

% load('out\LW_MarkusN_FINALE_LgDBR_Thesis.mat')
% load('out\LW_MarkusN_TJ_oxAbove_LgDBR_Thesis.mat')
% load('out\LW_MarkusN_TJ_oxBelow_LgDBR_Thesis.mat')
% load('out\LW_MarkusN_BTJetch_DD_ROSSO_LgDBR_singleMode_Isize=2_nmodes=2.mat')

% load LW_Stephan_GRok_.mat
% load LW_MarkusN_FINALE_LgDBR_saveVELM_newVELM.mat
% iold=input('  old VELM? [ENTER]: Yes; [Any key]: No  ');
% if isempty(iold)
%     rmpath('OtticoBar\new23OpticaGR')
%     load LW_MarkusN_FINALE_LgDBR_oldVELM.mat
% %     load LW_MarkusN_FINALEdisFitto_oldVELM.mat
% else
%     rmpath('OtticoBar\new22Optica')
%     load LW_MarkusN_FINALE_LgDBR_newVELM.mat
% %     load LW_MarkusN_FINALEdisFitto_newVELM.mat
% end
load LW_MarkusN_FINALEdisFitto_.mat

ParVet=MODEplot{1}.ParVet;
VelmOptions=MODEplot{1}.VelmOptions;

structureName=[nomeSR,mode.strName];
fis= strfind(structureName,'\');
strName=structureName(fis(end)+1:end);
DirName=structureName(1:fis(end));

fil_str=[structureName,'.str'];

verVE=1;    % verbVELM
ian=1;      % guiding
MulT=1;     % Temperature pre-factor

nqw=3;      % Number of QW

% mode.DT0=200;
% NP_k=20;

% load velmSa

if nqw==1
    fat_gain=2.6;
else
    fat_gain=1;
end

s_LoadConstants
vph=Clight./mode.nindexQW;


vind=[];    % DD indexes of VELM called (indVoltage)
veind=[];   % Progressive vector of VELM calling (indVELM)
for k=1:length(VELMInput)-1
    vv=VELMInput(k).indVoltage;
    vev=VELMInput(k).indVELM;
    if length(vv)==1
        vind=[vind vv];
        veind=[veind vev];
        fPold(k,:)=VELMInfo(k).fPdif;
        Laold(k,:)=VELMInfo(k).vlambda;
        Ga(k,:)=VELMInfo(k).Lm/vph;
        Tvet(k)=VELMInfo(k).DeltaTmax;
        E2v(:,:,k)=VELMInfo(k).E2';
    end
end
modeold=mode;

%'corrente', keyboard
Pcor=input(' Corrente = ')

Cor=1000*modeold.ii_dd(vind);
% Cor=1000*modeold.ii_dd;         % Computed current from DD  
[~,fim]=min(abs(Pcor-Cor));     % find closest current index to Cor (DD)
pu=fim;                         % for plot purposes

% ves= find(vind==fim);           % index of VELM calling (1st, 2nd,...)      
ves=fim;           % index of VELM calling (1st, 2nd,...)      
vindm=vind(1:end-1);

figure, hold on
plot(1000*mode.ii_dd,modeold.Gmod,'g','LineWidth',2)
plot(1000*mode.ii_dd,modeold.Lmod,'r--','LineWidth',2)
if(mode.nmodes>1)
    plot(1000*modeold.ii_dd(vind),modeold.Gmod(:,vind)','go','LineWidth',2)
    plot(1000*modeold.ii_dd(vind),modeold.Lmod(:,vind)','ro','LineWidth',2)
else
    plot(1000*modeold.ii_dd(vindm),modeold.Gmod(vindm),'go','LineWidth',2)
    plot(1000*modeold.ii_dd(vind),modeold.Lmod(vind),'ro','LineWidth',2)
end
xlabel('Current, mA')
ylabel('Modal gain vs losses, cm^{-1}')
ylim([0 100])
pausak

I=1000*modeold.ii_dd(vind);     % Current when VELM is called
I(ves)
I0=Pcor;                        % target current

figure, 
plot(I,Tvet(veind),'LineWidth',2)   % VELM current vs VELM Temperature
hold on, plot(I(ves),Tvet(ves),'mo','LineWidth',2),
xlabel('Corrente')
ylabel(' Tmax ')
pausak

%figure, plot(mesh.xQW,E2v,'LineWidth',2),
figure, hold on
plot(1e4*mesh.xgrid,squeeze(E2v(:,1,:)),'LineWidth',1.5),
plot(1e4*mesh.xgrid,squeeze(E2v(:,1,ves)),'r','LineWidth',2),
if size(E2v,2)>1
    chold, plot(1e4*mesh.xgrid,squeeze(E2v(:,2,:)),'--','LineWidth',1.5),
    plot(1e4*mesh.xgrid,squeeze(E2v(:,2,ves)),'r--','LineWidth',2),
end
xlabel(' \rho, \mum ')
ylabel(' Field intensity ')
xlim([0 10])
pausak


figure, hold on
plot(I,squeeze(E2v(10,:,veind)),'LineWidth',2),
if size(E2v,2)>1
    hold on, plot(I(ves),squeeze(E2v(10,:,ves)),'mo','LineWidth',2),
end
xlabel('Corrente, mA')
ylabel(['Fields intensity at x(10)=',num2str(mesh.xgrid(10)*1e4),' \mum'])
pausak

figure, plot(I,fPold(veind,:),'LineWidth',2),
hold on, plot(I(ves),fPold(ves,:),'mo','LineWidth',2),
xlabel('Corrente')
ylabel(' fPdif ')
pausak

figure, plot(I,Laold(veind,:),'LineWidth',2),
hold on, plot(I(ves),Laold(ves,:),'mo','LineWidth',2),
xlabel('Corrente')
ylabel(' Lambda ')
pausak

%return

%'corrente', keyboard
%Pcor=input(' Corrente = ')

Cor=1000*modeold.ii_dd(vind);
[~,fim]=min(abs(Pcor-Cor));
veind=fim;

ico=0;
for indVELM=veind
    ico=ico+1;
    
    mode.verbVELM=verVE;
%     VelmOptions.ianti_gui=ian;  % 0 per LP
    DL=VelmOptions.Dlam;  % 0 per LP
    %DL(5)=.8;
    %VelmOptions.Dlam=DL;  % 0 per LP
%     VelmOptions.gain_gui=ian;  % 0 per LP
%     VelmOptions.NP_k=NP_k;  % 0 per LP
    VelmOptions.itutmir=0; % 1: caso termico completo; 0: caso ridotto (+ veloce)
%     VelmOptions.Pf.nmasce=-3;
    
    mode.matgain=VELMInput(indVELM).matgain;
    mode.DeltaN=VELMInput(indVELM).DeltaN;
    mode.Deltan=VELMInput(indVELM).Deltan;
    mode.efield_y=VELMInput(indVELM).efield_y;
    mode.efield_x=VELMInput(indVELM).efield_x;
    mode.vlambda=VELMInput(indVELM).vlambda;
    mode.alpha=VELMInput(indVELM).alpha;
    mode.Lm=VELMInput(indVELM).Lm;
    mode.NQW=VELMInput(indVELM).NQW;
    mode.Gamma_z=VELMInput(indVELM).Gamma_z;
    mode.nindexQW=VELMInput(indVELM).nindexQW;
    mode.fPES=VELMInput(indVELM).fPES;
    mode.fPdif=VELMInput(indVELM).fPdif;
    mode.E2=VELMInput(indVELM).E2;
    mode.Gmod=VELMInput(indVELM).Gmod;
    if indVELM>1
        mode1.TmVelm=VELMInput(indVELM).TmVelm;
        mode1.LamVelm=VELMInput(indVELM).LamVelm;
    else
        mode1.a=0;
    end
    mesh.DeltaTvelm=MulT*VELMInput(indVELM).DeltaTvelm;
    mesh.ygrid=VELMInput(indVELM).ygrid;
    mesh.xgrid=VELMInput(indVELM).xgrid;
    mesh.nnx=VELMInput(indVELM).nnx;
    mesh.nny=VELMInput(indVELM).nny;
    
    'ver prima di call', keyboard
    
    [velm] = f_CallVELM(mesh,mode,mode1,ParVet,VelmOptions,fil_str);
    indv=vind(ico);
    
    mode.Lmod(:,indv)=velm.Lm/vph;
    mode.lambda(:,indv)=velm.vlambda;
    fPvet(:,indv)=velm.fPdif;
    E2v(:,:,indv)=velm.E2'; %dà errore qui (multimodo)
    Nv(:,indv)=mode.matgain;
    Nrv(:,indv)=mode.DeltaN;
end

% colordef white
vindm=vind(1:end-1);
figure, hold on
plot(mode.ii_dd*1e3,modeold.Gmod,'k','LineWidth',2)
plot(mode.ii_dd*1e3,modeold.Lmod,'r--','LineWidth',2)
if(mode.nmodes>1)
    plot(modeold.ii_dd(vind)*1e3,modeold.Gmod(:,vind)','ko','LineWidth',2)
    plot(modeold.ii_dd(vind)*1e3,modeold.Lmod(:,vind)','ro','LineWidth',2)
    %       plot(mode.ii_dd(vind),mode.Gmod(:,vind)','k+','LineWidth',2)
    %     plot(mode.ii_dd(vind),mode.Lmod(:,vind)','r+','LineWidth',2)
else
    plot(modeold.ii_dd(vindm)*1e3,modeold.Gmod(vindm),'ko','LineWidth',2)
    plot(modeold.ii_dd(vind)*1e3,modeold.Lmod(vind),'ro','LineWidth',2)
    plot(mode.ii_dd(vindm)*1e3,mode.Gmod(vindm)','k+','LineWidth',2)
    plot(mode.ii_dd(vind)*1e3,mode.Lmod(vind)','r+','LineWidth',2)
end
xlabel('Current, mA')
ylabel('Modal gain vs losses, cm^{-1}')
% 
indQW=2;
inQW = mesh.inMQW{indQW};
% xQW = mesh.node(1,inQW);
xQW=mesh.xgrid(1:mesh.nnxQW{1})*1e4;
% figure, plot(xQW,squeeze(E2v(1:mesh.nnxQW{1},:,pu)),'LineWidth',2),
figure, plot(xQW,velm.E2(:,1:mesh.nnxQW{1}),'LineWidth',2),
xlabel(' \rho, \mum ')
ylabel('Optical Field intensity, a.u.')
pausak
%
% figure, plot(mode.ii_dd(vindm),fPvet(vindm),'LineWidth',2),
% xlabel('Voltage, V')
% ylabel(' fPdif ')

