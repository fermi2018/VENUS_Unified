clear
clear global
close all
colordef white
dbstop if error

addpath('termico')
addpath('rumenta')
addpath('generageom')
addpath('rumenta\new17Optica')
  TITE{1}='EXPer';
  TITE{2}='None';
  TITE{3}='fPES';
  TITE{4}='Gam_z';
  TITE{5}='fPdif';
  TITE{6}='E2';
  TITE{7}='lambda';
  TITE{8}='TeVelm';
  TITE{9}='anti gui';
  TITE{10}='Diff';
  TITE{11}='Strutt';
  TITE{12}='LUT';

iLoad=input(' Vuoi ultima simulazione? [0, no, Enter, si] ')
%keyboard
if length(iLoad)==0
 load ultimonome
% save nomeLoopFinito nomeSave

 eval(['load  ',nomeSave,'_pmat'])
else

 nomeSave='FE26sta';
 nomeSave='FE27rel';
 nomeSave='FE27sta';
 nomeSave='FE28sta';
 %IPvet=2:5
 %nomeSave='FE26rel';
 nomeSave='Effetti_totF';
 nomeSave='Spegni1';   % originale
 nomeSave='SpegniO';
 nomeSave='SpegniNotteTnul';
 nomeSave='Ven9d';
 nomeSave='Lu';
 nomeSave='De1';
 nomeSave='De6Lo';
 %nomeSave='Effetti_totM';
 %nomeSave='Effetti_0';
 %nomeSave='Effetti_fPdif';
 eval(['load  ',nomeSave,'_pmat'])
 IPvet=0
 IPvet=[27]
% keyboard
end 

nomeSave
pausak
nome=nomeSave;

if irel==1
    load('Meas_ReliefVCSEL.mat')
else
    load('Meas_StandardVCSEL.mat')
end

CUR=Cur;
Rmeas=diff(Vmeas)./diff(Imeas/1000);
Rmeas=[Rmeas(1:end-2)];
Imeas_res=Imeas(1:end-3);

%IPvet=9

for IPAR=IPvet

clear MODEplot

 if IPAR==4
  Fpar=1e-18;
 else
  Fpar=1;
 end

% eval(['load All',num2str(IPAR),'.mat MODE MESH VInf VInp'])
 eval(['exi=exist( ''',nome,num2str(IPAR),'.mat',''');'])
 if exi==0
  return
 end 
 eval(['load ',nome,num2str(IPAR)])
 %keyboard
colo='rgbcmrb';
%col{1}='k'; col{2}='r'; col{3}='g'; col{4}='b'; col{5}='c'; col{6}='m'; col{7}='r.'; col{8}='b.';
col{1}='r'; col{2}='g'; col{3}='b'; col{4}='c'; col{5}='m'; col{6}='r.'; col{7}='b.';
if IPAR>0
 ppar=PMAT{IPAR}
else
 ppar=[];
end
%keyboard

for kpar=1:length(MODEplot)
 kpar
 mode=MODEplot{kpar};
 modePlot=MODEplot{kpar};

 
 if IPAR>0
 
 if kpar<length(ppar)
  tita=[tit{IPAR},'  Valore attuale= ',num2str(Fpar*ppar(kpar))]; 
 else
  tita=tit{IPAR};  
 end
% [fPES Gamma_z fPdif  E2 lambda TempVELM anti_gui]

  if IPAR==30
   titeff=['EFFETTI: ',TITE{PMAT{IPAR}(kpar)+2}]
  else
   titeff=[];
  end
 else
  tita='Set 0';
  titeff='Set 0';
 end
 %mesh=MESH{kpar};
 I=mode.ii_dd*1000;
 %'pa', keyboard
 V=mode.vv_dd;
 P=sum(mode.Pst_dd,1)+mode.Psp_dd;


Rec=mode.IntRec;
Ccap=-sum(mode.IntCcapN,2);
Rnor=diag(1./Ccap)*Rec;
        
         
                               Tstart=250;    
	                       taubottom=1;
	                       DeT=100;                      
	                       tautop=1/50;
	                       tau_lea=1e-7;

         n2Di=.1;         p2Di=.1;
         x=mode.x*1e-4;
         qel=1.6e-19;
         NumQW=3;
for kv=2:length(modePlot.ii_dd)         
 n2D=modePlot.Elqw(kv,:);
 p2D=modePlot.Hoqw(kv,:);
 TQW=modePlot.Tqw(kv,:);         
 np2D=n2D.*p2D-n2Di.*p2Di;
 tauleak=tau_lea*f_RaisedCosine(TQW,Tstart,DeT,taubottom,tautop);              
 denLeakage2D = tauleak.*(n2D+n2Di)+tauleak.*(p2D+p2Di);
 Rlea = NumQW*qel*np2D./denLeakage2D;
 Nlea(kv)=2*pi*trapz(x,x.*Rlea)/Ccap(kv);
% 'ver;,', keyboard
end 
% 'ver;,', keyboard
fi=find(Nlea>1);
Nlea(fi)=1;
Cd=mode.ii_dd*1000;
figure, semilogy(Cd,Nlea), pausak


figure
         plot(mode.ii_dd*1000,abs(Rnor),'LineWidth',1.5), hold on
         plot(mode.ii_dd*1000,Nlea,'ko','LineWidth',2)
         legend('IntSRH','IntRad','IntAug','IntStim','IntLea','location','best')
         title('Recombinations/Ccap')
         xlabel('Corrente (mA)')
%         ylim([1e-3 1])
         pausak

figure
         semilogy(mode.ii_dd*1000,abs(Rnor),'LineWidth',1.5), hold on
         semilogy(mode.ii_dd*1000,Nlea,'ko','LineWidth',2)
         legend('IntSRH','IntRad','IntAug','IntStim','IntLea','location','best')
         title('Recombinations/Ccap')
         xlabel('Corrente (mA)')
         ylim([1e-3 1])
         pausak
         
F=2./(1+exp(Cd/(9-0.7)));         
figure
plot(Imeas,Lmeas,'k.','linewidth',2)
hold on
plot(Cd,1000*Rec(:,4),'linewidth',2)
grid
%plot(Cd,1000*(Rec(:,4)-Rec(:,1)+Nlea.'.*Ccap),'o','linewidth',2)
         
         return
         
figure
         semilogy(mode.ii_dd*1000,Rnor,'LineWidth',2)
         legend('IntSRH','IntRad','IntAug','IntStim','IntLea','location','best')
         title('Recombinations/Ccap')
         xlabel('Corrente (mA)')
         ylim([1e-3 1])
         pausak 

figure
         plot(mode.ii_dd*1000,Ccap*1000,'LineWidth',2)
         xlabel('Corrente (mA)')
         title('Ccap')
         pausak
         
         return

iDDPlot=input(' Vuoi DD_PLOT? [1, s1, Enter, NO] ')
if length(iDDPlot)==1
 DD_PLOT_pris
 pausak
 close(1229:1231) 
 close(1234:1238) 
end
end
 pausak
close all
end