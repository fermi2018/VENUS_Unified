clear
clear global
close all
colordef black
dbstop if error

addpath('Termico')
addpath('Ottico')
addpath('Dati')
addpath('out')


inewpat=1;
if inewpat==1
 rmpath('generageom')
 rmpath('Ottico')
 rmpath('Ottico\new17Optica')
 addpath('generageomBar')
 addpath('OtticoBar') 
 addpath('OtticoBar\new17Optica') 
else 
 rmpath('generageomBar')
 rmpath('OtticoBar') 
 rmpath('OtticoBar\new17Optica') 
 addpath('generageom')
 addpath('Ottico')
 addpath('Ottico\new17Optica')
end 


verVE=1;
ian=1;           %antiguiding
MulT=1         % fattore moltiplicativo Temperatura

nqw=3



NP_k=20;

%iExte=input(' Este ? [0/1] ');

%   load velmS1
  load velmSa

 if nqw==1
  fat_gain=2.6;
 else
  fat_gain=1;
end

mode.ABS_Texp=0;

%'ver', keyboard

s_LoadConstants
vph=Clight./mode.nindexQW;


vind=[];
veind=[];
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

%Cor=1000*modeold.ii_dd(vind);
Cor=1000*modeold.ii_dd;
[du,fim]=min(abs(Pcor-Cor));

ves= find(vind==fim);

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
%figure, plot(
I=1000*modeold.ii_dd(vind);
I0=Pcor;

figure, plot(I,Tvet(veind),'LineWidth',2),
hold on, plot(I0,Tvet(ves),'mo','LineWidth',2),
        xlabel('Corrente')
        ylabel(' Tmax ')
        pausak

%figure, plot(mesh.xQW,E2v,'LineWidth',2),
figure, plot(1e4*mesh.xgrid,squeeze(E2v(:,1,:)),'LineWidth',1.5),
if size(E2v,2)>1
hold on, plot(1e4*mesh.xgrid,squeeze(E2v(:,2,:)),'--','LineWidth',1.5),
end
%        xlabel('Voltage, V')
        ylabel(' Field intensity ')
pausak

figure, plot(1000*modeold.ii_dd(veind),squeeze(E2v(10,:,veind)),'LineWidth',2),
if size(E2v,2)>1
hold on, plot(I0,squeeze(E2v(10,:,ves)),'mo','LineWidth',2),
end
        xlabel('Corrente')
        ylabel(' Field intensity at xqw(10) ')
pausak

figure, plot(I,fPold(veind,:),'LineWidth',2),
hold on, plot(I0,fPold(ves,:),'mo','LineWidth',2),
        xlabel('Corrente')
        ylabel(' fPdif ')
        pausak
        
figure, plot(I,Laold(veind,:),'LineWidth',2),
hold on, plot(I0,Laold(ves,:),'mo','LineWidth',2),
        xlabel('Corrente')
        ylabel(' Lambda ')        
pausak

%return

%'corrente', keyboard
%Pcor=input(' Corrente = ')

Cor=1000*modeold.ii_dd(vind);
[du,fim]=min(abs(Pcor-Cor));
veind=fim

ico=0;
for indVELM=veind
ico=ico+1;

 mode.verbVELM=verVE;
 %VelmOptions.ianti_gui=ian;  % 0 per LP
 %DL=VelmOptions.Dlam;  % 0 per LP
 %DL(5)=-1;
 %VelmOptions.Dlam=DL;  % 0 per LP
 %VelmOptions.gain_gui=ian;  % 0 per LP 
 %VelmOptions.ianti_gui=ian;  % 0 per LP 
 VelmOptions.NP_k=NP_k;  % 0 per LP 
 VelmOptions.itutmir=0; % 1: caso termico completo; 0: caso ridotto (+ veloce) 

            mode.matgain=VELMInput(indVELM).matgain;
            mode.DeltaN=VELMInput(indVELM).DeltaN;
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
            mode1.TmVelm=VELMInput(indVELM).TmVelm;
            mode1.LamVelm=VELMInput(indVELM).LamVelm;
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
         fPvet(indv)=velm.fPdif;
         E2v(:,indv)=velm.E2';
         Nv(:,indv)=mode.matgain;
         Nrv(:,indv)=mode.DeltaN;
end

vindm=vind(1:end-1);
figure, hold on
        plot(mode.vv0_dd,modeold.Gmod,'k','LineWidth',2)
        plot(mode.vv0_dd,modeold.Lmod,'r--','LineWidth',2)
        if(mode.nmodes>1)
            plot(modeold.vv0_dd(vind),modeold.Gmod(:,vind)','ko','LineWidth',2)
            plot(modeold.vv0_dd(vind),modeold.Lmod(:,vind)','ro','LineWidth',2)
            plot(mode.vv0_dd(vind),mode.Gmod(:,vind)','k+','LineWidth',2)
            plot(mode.vv0_dd(vind),mode.Lmod(:,vind)','r+','LineWidth',2)            
        else
            plot(modeold.vv0_dd(vindm),modeold.Gmod(vindm),'ko','LineWidth',2)
            plot(modeold.vv0_dd(vind),modeold.Lmod(vind),'ro','LineWidth',2)
            plot(mode.vv0_dd(vindm),mode.Gmod(vindm)','k+','LineWidth',2)
            plot(mode.vv0_dd(vind),mode.Lmod(vind)','r+','LineWidth',2)              
        end
        xlabel('Voltage, V')
        ylabel('Modal gain vs losses, cm^{-1}')
        xlim([1.9 2.8])
pausak

figure, plot(mesh.xQW,E2v(:,vind),'LineWidth',2),
        xlabel('Voltage, V')
        ylabel(' Field intensity ')
pausak

figure, plot(mode.vv0_dd(vindm),fPvet(vindm),'LineWidth',2), 
        xlabel('Voltage, V')
        ylabel(' fPdif ')

