%clear all
%clear global
close all
colordef white
dbstop if error

addpath('Termico')
addpath('OtticoBar')
addpath('Dati')
%addpath('out')
icarica=input(' Carica dati  [enter for NO]');
%keyboard
if length(icarica)==0
 icarica=0;
end 
wqw=8e-7;

iold=2;  %0 _AT,  2 _PD

 %NOMELUT='LUT4D_Feb_nMark_40';
%  NOMELUT='LUT4D_May_nMark_40';
 NOMELUT='LUT4D_Jun_Markus_nMark_40';
 
 
 mode.Fasano=['Fasano'];
 
mode.GLUT=[NOMELUT,'_Der.mat'];
mode.iderGLUT=1;   % = 0 fa vecchio fitN2


mode.iderGLUTnew=iold;   % = 0 fa vecchio fit
geom.GLUTm=[NOMELUT,'_more.mat'];

%carica LUT e la mette nei common
if icarica==1
 Glut4Dinc(mode)
end

global LAV

% Scelgo i portatori che voglio far vedere in logspace
N2lin=logspace(-2,1,31)*1e12;
P2lin=N2lin;
N2pl=N2lin/wqw;
mode.Deltalam=0;
%
% Facciamo 3 tagli in temperatura, con lambda=860 nm;
% T2=linspace(300,350,length(N2lin));
T2=[300 450];
L2=850;

% [Rsp0,dRspE0,dRspH0]=f_InterpRsp4D(np_E,np_H,T2);
ggT=zeros(length(T2),length(N2lin));
rspT=zeros(length(T2),length(N2lin));
for indT=1:length(T2)

    [G0,dgE0,dgH0,rsp0,drspE0,drspH0,DeltaN0] = f_InterpGain4D(N2lin,P2lin,mode.Deltalam+L2,T2(indT));
    for ind_diag=1:length(N2lin)
        ggT(indT,ind_diag)=G0(ind_diag,ind_diag,1,1);
        rspT(indT,ind_diag)=rsp0(ind_diag,ind_diag,1,1);
    end
    Tlegend{indT}=['T=',num2str(T2(indT)),' K'];
    
end

% figure,plot(N2lin,G0)
figure(1)
set(gcf,'Position',[101 550 560 420])
hold on
grid on
plot(N2pl,rspT,'LineWidth',2)
xlabel('n2D, p2D, cm^{-2}')
xlabel('n2D, p2D, cm^{-3}')
ylabel('r_{sp}, s^{-1}')
legend(Tlegend,'location','best')
title('r_{sp} vs ehDensity for different temperatures')

%
% Facciamo 3 tagli in lambda, con T=300K;
% T2=linspace(300,350,length(N2lin));
T2=[300];
L2=[840 850 860];

% [Rsp0,dRspE0,dRspH0]=f_InterpRsp4D(np_E,np_H,T2);
ggL=zeros(length(L2),length(N2lin));
rspL=zeros(length(L2),length(N2lin));
for indlam=1:length(L2)

    [G0,dgE0,dgH0,rsp0,drspE0,drspH0,DeltaN0] = f_InterpGain4D(N2lin,P2lin,mode.Deltalam+L2(indlam),T2);
    for ind_diag=1:length(N2lin)
        ggL(indlam,ind_diag)=G0(ind_diag,ind_diag,1,1);
        rspL(indlam,ind_diag)=rsp0(ind_diag,ind_diag,1,1);
    end
    Llegend{indlam}=['\lambda=',num2str(L2(indlam)),' nm'];

end

% figure,plot(N2lin,G0)
figure(2)
set(gcf,'Position',[695 550 560 420])
hold on
grid on
plot(N2pl,rspL,'LineWidth',2)
xlabel('n2D, p2D, cm^{-2}')
xlabel('n2D, p2D, cm^{-3}')
ylabel('r_{sp}, s^{-1}')
legend(Llegend,'location','best')
title('r_{sp} vs ehDensity for different wavelengths')


N2lin=[3 6]*1e12;
P2lin=N2lin;
N2pl=N2lin/wqw;
L2temp=squeeze(LAV(1,1,:,1));
L2=linspace(830,L2temp(end),101);
T2=300;
T2l=T2;
for indeDensity=1:length(N2lin)
    [G0,dgE0,dgH0,rsp0,drspE0,drspH0,DeltaN0] = f_InterpGain4D(N2lin(indeDensity),P2lin(indeDensity),mode.Deltalam+L2,T2);
    ggeDens(indeDensity,:)=squeeze(G0);
    rspeDens(indeDensity,:)=squeeze(rsp0);
    eDLegend{indeDensity}=['n,p = ',num2str(N2lin(indeDensity)/1e12),'\times10^{12} cm^{-2}'];
end


% figure,plot(N2lin,G0)
figure(3)
set(gcf,'Position',[104 20 560 420])
hold on
grid on
plot(L2,rspeDens,'LineWidth',2)
xlabel('wavelength, nm')
ylabel('r_{sp}, s^{-1}')
legend(eDLegend,'location','best')
title(['r_{sp} vs wavelength for different edensities, T = ',num2str(T2),' K'])






N2lin=[6]*1e12;
P2lin=N2lin;

L2temp=squeeze(LAV(1,1,:,1));
L2=linspace(830,L2temp(end),101);
T2=[300 450];
%T2=[300 350 400 450 500];
for indT=1:length(T2)
    [G0,dgE0,dgH0,rsp0,drspE0,drspH0,DeltaN0] = f_InterpGain4D(N2lin,P2lin,mode.Deltalam+L2,T2(indT));
    ggT2(indT,:)=squeeze(G0);
    rspT2(indT,:)=squeeze(rsp0);
    T2legend{indT}=['T = ',num2str(T2(indT)),' K'];
end


% figure,plot(N2lin,G0)
figure(4)
set(gcf,'Position',[694 20 560 420])
hold on
grid on
plot(L2,rspT2,'LineWidth',2)
xlabel('wavelength, nm')
ylabel('r_{sp}, s^{-1}')
legend(T2legend,'location','best')
title(['r_{sp} for different temperatures; N=',num2str(N2lin/1e12),'e12/cm^2'])

figure(30)
set(gcf,'Position',[ 244          92        1155         870])
subplot(221)
hold on
grid on
plot(L2,rspeDens*1e-9,'LineWidth',2)
xlabel('wavelength, nm')
ylabel('r_{sp}, ns^{-1}')
legend(eDLegend,'location','best')
title([' T = ',num2str(T2l),' K'])
subplot(223)

hold on
grid on
plot(L2,rspT2*1e-9,'LineWidth',2)
xlabel('wavelength, nm')
ylabel('r_{sp},n s^{-1}')
legend(T2legend,'location','best')
title([' N=',num2str(N2lin/1e12),'e12/cm^2'])


'salvare: procedendo chiude tutto '
pausak
close(1:4)

N2lin=logspace(-2,1,31)*1e12;
P2lin=N2lin;
mode.Deltalam=0;
T2=[300 450];
L2=850;

% [Rsp0,dRspE0,dRspH0]=f_InterpRsp4D(np_E,np_H,T2);
ggT=zeros(length(T2),length(N2lin));
for indT=1:length(T2)

    [G0,dgE0,dgH0,rsp0,drspE0,drspH0,DeltaN0] = f_InterpGain4D(N2lin,P2lin,mode.Deltalam+L2,T2(indT));
    for ind_diag=1:length(N2lin)
        ggT(indT,ind_diag)=G0(ind_diag,ind_diag,1,1);
    end
    Tlegend{indT}=['T=',num2str(T2(indT)),' K'];
    
end

N2pl=N2lin/wqw;

figure(1)
set(gcf,'Position',[101 550 560 420])
hold on
grid on
plot(N2pl,ggT,'LineWidth',2)
xlabel('n2D, p2D, cm^{-2}')
xlabel('n2D, p2D, cm^{-3}')
ylabel('gain, cm^{-1}')
legend(Tlegend,'location','best')
pausak
title('Gain vs ehDensity for different temperatures')

%
% Facciamo 3 tagli in lambda, con T=300K;
% T2=linspace(300,350,length(N2lin));
T2=[300];
L2=[840 850 860];

% [Rsp0,dRspE0,dRspH0]=f_InterpRsp4D(np_E,np_H,T2);
ggL=zeros(length(L2),length(N2lin));
for indlam=1:length(L2)

    [G0,dgE0,dgH0,rsp0,drspE0,drspH0,DeltaN0] = f_InterpGain4D(N2lin,P2lin,mode.Deltalam+L2(indlam),T2);
    for ind_diag=1:length(N2lin)
        ggL(indlam,ind_diag)=G0(ind_diag,ind_diag,1,1);
    end
    Llegend{indlam}=['\lambda=',num2str(L2(indlam)),' nm'];

end

% figure,plot(N2lin,G0)
figure(2)
set(gcf,'Position',[695 550 560 420])
hold on
grid on
plot(N2pl,ggL,'LineWidth',2)
xlabel('n2D, p2D, cm^{-2}')
xlabel('n2D, p2D, cm^{-3}')
ylabel('gain, cm^{-1}')
legend(Llegend,'location','best')
title('Gain vs ehDensity for different wavelengths')


N2lin=[3 6]*1e12;
P2lin=N2lin;
L2temp=squeeze(LAV(1,1,:,1));
L2=linspace(830,L2temp(end),101);
T2=300;
T2l=T2;
for indeDensity=1:length(N2lin)
    [G0,dgE0,dgH0,rsp0,drspE0,drspH0,DeltaN0] = f_InterpGain4D(N2lin(indeDensity),P2lin(indeDensity),mode.Deltalam+L2,T2);
    ggeDens(indeDensity,:)=squeeze(G0);
    eDLegend{indeDensity}=['n,p = ',num2str(N2lin(indeDensity)/1e12),'\times10^{12} cm^{-2}'];
end


% figure,plot(N2lin,G0)
figure(3)
set(gcf,'Position',[104 20 560 420])
hold on
grid on
plot(L2,ggeDens,'LineWidth',2)
xlabel('wavelength, nm')
ylabel('gain, cm^{-1}')
legend(eDLegend,'location','best')
title(['Gain vs wavelength for different edensities, T = ',num2str(T2),' K'])






N2lin=[6]*1e12;
P2lin=N2lin;
L2temp=squeeze(LAV(1,1,:,1));
L2=linspace(830,L2temp(end),101);
T2=[300 450];
for indT=1:length(T2)
    [G0,dgE0,dgH0,rsp0,drspE0,drspH0,DeltaN0] = f_InterpGain4D(N2lin,P2lin,mode.Deltalam+L2,T2(indT));
    ggT2(indT,:)=squeeze(G0);
    T2legend{indT}=['T = ',num2str(T2(indT)),' K'];
end


% figure,plot(N2lin,G0)
figure(4)
set(gcf,'Position',[694 20 560 420])
hold on
grid on
plot(L2,ggT2,'LineWidth',2)
xlabel('wavelength, nm')
ylabel('gain, cm^{-1}')
legend(T2legend,'location','best')
title('Gain vs ehDensity for different temperatures')

'salvare: procedendo chiude tutto '
pausak
close(1:4)

figure(30)
%set(gcf,'Position',[   818    92   574   870])
subplot(222)
hold on
grid on
plot(L2,ggeDens,'LineWidth',2)
xlabel('wavelength, nm')
ylabel('gain, cm^{-1}')
legend(eDLegend,'location','best')
title([' T = ',num2str(T2l),' K'])
subplot(224)

hold on
grid on
plot(L2,ggT2,'LineWidth',2)
xlabel('wavelength, nm')
ylabel('gain, cm^{-1}')
legend(T2legend,'location','best')
title([' N=',num2str(N2lin/1e12),'e12/cm^2'])
ylim([-1000 5000])