
% clear all
close all
%clc
warning off
colordef black

addpath('generageom')
addpath('rumenta')
addpath('rumenta\new17Optica')
addpath('termico')


verVE=1;
ian=1;

nqw=3



%  if nqw==3
% %  load VELMinput_3qw
%   load VELMinput
%  else
%   load VELMinput_1qw
%  end
 
 if nqw==1
  fat_gain=2.6;
 else
  fat_gain=1;
end




s_LoadConstants
vph=Clight./mode.nindexQW;


vind=[];
veind=[];
kk=0;
for k=2:length(VELMInput)-1
 kk=kk+1;
 vv=VELMInput(k).indVoltage;
 vev=VELMInput(k).indVELM;
 if length(vv)==1
  vind=[vind vv];
  veind=[veind vev];
  fPold(:,kk)=VELMInfo(k).fPdif;
  Laold(:,kk)=VELMInfo(k).vlambda;
  Ga(:,kk)=VELMInfo(k).Lm/vph;
  Tvet(kk)=VELMInfo(k).DeltaTmax;
  E2v(:,:,kk)=VELMInfo(k).E2';
  MatG(:,kk)=VELMInput(k).matgain;
  MatD(:,kk)=VELMInput(k).DeltaN;
 end 
end
modeold=mode;
x=VELMInput(2).xgrid*1e4;

Cor=1000*modeold.ii_dd(vind);
%'corrente', keyboard
Pcor=input(' Corrente = ')
if length(Pcor)==0
 Pcor=Cor(end);
end

[du,fim]=min(abs(Pcor-Cor));
ves=fim


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
%        xlim([1.9 2.8])
pausak
%figure, plot(
I=1000*modeold.ii_dd(vind);
I0=Pcor;
pu=1:length(I);
figure, plot(I,Tvet(pu),'LineWidth',2),
hold on, plot(I0,Tvet(ves),'mo','LineWidth',2),
        xlabel('Corrente')
        ylabel(' Tmax ')
        pausak

figure, 
nmod=size(E2v,2);
X=mesh.xQW*1e4;
for nm=1:nmod
plot(X,squeeze(E2v(:,nm,:)),'LineWidth',2),
hold on, plot(X,squeeze(E2v(:,nm,ves)),'mo','LineWidth',2),
%        xlabel('Voltage, V')
        ylabel(' Field intensity ')
end        
pausak
puvi=1:length(vind);
figure, plot(1000*modeold.ii_dd(vind),squeeze(E2v(1,:,puvi)),'LineWidth',2),
hold on, plot(I0,E2v(1,:,ves),'mo','LineWidth',2),
        xlabel('Corrente')
        ylabel(' Field intensity at center ')
pausak
puviL=1:length(vind);
figure, plot(1000*modeold.ii_dd(vind(1:end)),fPold(:,puviL),'LineWidth',2),
hold on, plot(I0,fPold(:,ves),'mo','LineWidth',2),
        xlabel('Corrente')
        ylabel(' fPdif ')
        pausak
        
figure, plot(1000*modeold.ii_dd(vind(1:end)),Laold(:,puviL),'LineWidth',2),
hold on, plot(I0,Laold(:,ves),'mo','LineWidth',2),
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
 VelmOptions.ianti_gui=ian;  % 0 per LP
 VelmOptions.gain_gui=ian;  % 0 per LP 
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
            mesh.DeltaT=VELMInput(indVELM).DeltaT;
            mesh.ygrid=VELMInput(indVELM).ygrid;
            mesh.xgrid=VELMInput(indVELM).xgrid;
            mesh.nnx=VELMInput(indVELM).nnx;
            mesh.nny=VELMInput(indVELM).nny;



[velm] = f_CallVELM(mesh,mode,mode1,ParVet,VelmOptions,fil_str);
         indv=vind(ico);
      
         mode.Lmod(:,indv)=velm.Lm/vph;
         mode.lambda(:,indv)=velm.vlambda;    
         fPvet(:,indv)=velm.fPdif';
         E2v(:,:,indv)=velm.E2';
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

