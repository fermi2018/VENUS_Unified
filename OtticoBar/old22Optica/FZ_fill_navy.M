close all
clear all
load plotzro
%' in Fz_tot.m ', keyboard
%clear
%close all
%load toto1
%load toto2
%load GRd45
ionlyout=0;
ischema=0;     % disegna profile struttura
ipau=1;
faX=2;  % fattore allargamento asse x
%faX=1  % fattore allargamento asse x
fadK=100;


Zma=300; 
fadK=3.0;

Zma=5; 
fadK=5;

Zsotto=1;

irot=0;
%Zma=5; 
dzc=.1;

 
irighe=0;
nox=1.6;
imo=1;
ixv=1;

imo=input(' modo ? ')
if length(imo)==0
 imo=1;
end

for ix=ixv


%ix=1
iy=1
Pp=Plot;
Ppo=Ppol;
FF=Pp.FF;
uL=Ppo.uL;
faGam=uL(2)/uL(1);
vg=3e10/rr*faGam;
vg=faGam;
%' FF ', keyboard
ipoi=1;   % Pointing (1 reale, 2 imag) o campo (0)
ieh=1;    % campo E (1) o H (0)
iex=1;  % componente x (1) o y (2)
%iex=0;  % somma moduli quadri
if ipoi>0
 iex=3;
end 
%iex=1
isalf=0;    % 1 salta primi plots
ifita=1;
npadd=5;
npadd=20;
%load pro
iztm=1;

Nz=FF.Nz;
segem=FF.segem;
iredmat=FF.iredmat;
fmlstot=FF.fmlstot;
iauto=FF.iauto;
lambda=FF.lambda;
Litot=FF.Litot;
aitot=FF.aitot;
nitot=FF.nitot;
aiat=max(FF.aiat);
KKv=FF.KKv;
KKt=FF.KKt;

keyboard
ze=FF.ze;
rr=FF.rr;
besm=FF.besm;
iFFsect=FF.iFFsect;
xvero=FF.xvero;
xvero=xvero*faX;
xvero=linspace(0,max(xvero),500);
XP=xvero'*cos(fian);
YP=xvero'*sin(fian);
ama=max(xvero)/2;
assi=ama;
%YP=FF.YP*faX;
%XP=FF.XP*faX;

Pusas=FF.Pusas;
nv=FF.nv;
Gas=FF.Gas;
Gad=FF.Gad;
Znor=FF.Znor;
Pus=   FF.Pus;
iLP=   FF.iLP;
Mvefm0=FF.Mvefm0;
if iLP==0
 Mvefp0=FF.Mvefp0;
 Mvegm0=FF.Mvegm0;
 Mvegp0=FF.Mvegp0;
 Mvez0 =FF.Mvez0;
 Mvhz0 =FF.Mvhz0;
 besp=FF.besp;
 besz=FF.besz;
end
%' qui', keyboard
pimud=FF.mbv;
pimu=pimud(2:end-1);
numo=length(pimu);
nudiv=numo-numodiacc;
nudiv=1:2
dP=find(diff(KKt(Pusas))<0 | diff(Pusas')>1)+1;
dPai=[1; dP];
dPaf=[dP-1; length(Pusas)];
Psep=[dPai dPaf];
dK=max(diff(KKv))/fadK;
pes=FF.pes;
%pes=1;
Acoz=FF.Acoz;
zi=FF.zi;
fiu=find(diff([-10 zi])~=0);

zi=zi(fiu);
if iredmat==1
 Pusff=Pusas;
else
 Pusff=Pus;
end

%keyboard
vg=3e10/rr*uL(2)/uL(1);

 puuu1=1:length(Gas);
 puuu2=length(Gas)+1:2*length(Gas);

Efd=Pp.Ef{1}.^2;
lte=ceil(size(Efd,2)/2);
Efp=Efd(:,[1 lte])/Efd(1);
te=Pp.X{1}(:,[1 lte]);
[du,imi]=min(abs(Efp-exp(-2)));
figure,  figure, plot(te,Efp,'m','linewidth',2)
title([' FF width = ', num2str(2*te(imi(1),1),3)])
xlim([-20 20])
ylim([0 1.15])
xlabel('angle')
ylabel('Far-Field profile')
grid



%	Efd=Ef.^2;
%	lte=ceil(size(Efd,2)/2);
%	Ef=Efd(:,[1 lte])/max(max(Efd));
%	te=X(:,[1 lte]);
%	figure,  figure, plot(te,Ef,'r','linewidth',2)
%	xlim([-20 20])
%	xlabel('angle')
%	ylabel('Far-Field profile')

 offset=0;
 zcut=-([-(floor(zi(end)/dzc)*dzc):dzc:Zma]+offset);
 zcut=-([-(floor(zi(end)/dzc)*dzc):dzc:Zma]+offset);
 zcut=-[0 1 2 5];
 zcut=[0 5 10 15];
 %zcut=-[0 1 2 3]*10;
% zcut=(1-[.1 .5 .7 .8 .9 .95])*11.3;


%  pvet=pola;
%  PUA=pua;
iFFsect=1;


sM=size(Mvefm0);
sm2=sM(2);
if0=1;
if90=(sm2-1)/4+1;
if180=(sm2-1)/2+1;
if irot==1
if0=if90;
if180=if180+if0;
end

 dz=dzc;
 lambdam=lambda*(1+ze);
 if iFFsect==1
  rff=nv(1,1);
  rffb=(nv(end,1));
 else
  rffb=nv(1,1);
  rff=nv(end,1);
 end
 
 
if size(Acoz,2)<10
 ionlyout=1;
end 

 if iFFsect==1
  segdir=-1;
  sectu=[0:dz:Zma];

  sectb=[];
  if ionlyout==0
   sectu(2)=1e-3;
   sectb=[0:dz:Zsotto];
   sectb=sectb(2:end);
   sectu=-fliplr(sectu(2:end));
   ztot=[sectu zi sectb+zi(end)];
   Nztot=[ones(size(sectu))*rff Nz ones(size(sectb))*rffb];
  else
   sectu=-fliplr(sectu(1:end));
   ztot=sectu;
   Nztot=ones(size(sectu))*rff;
  end
 else
  segdir=1;
  sectb=[0:dz:Zsotto];
  sectb=-fliplr(sectb(2:end));
  sectu=[0:dz:Zma];
  sectu=sectu(2:end);
  ztot=[sectb zi zi(end)+sectu];
  Nztot=[ones(size(sectb))*rff Nz ones(size(sectu))*rffb];
 end
 fatEzv=(rr./Nz).^2;

 fatEz=fatEzv(1);




% fatcut=ones(size(KKt));
% fatcut=fatcut(Pusff);
fatcut=1;

 rrfa=1;
 lim=1-(rrfa*KKv).^2;
 SQ=sqrt(lim);

 if iLP==0
  Znorn=[1./SQ; SQ];
  Ynorn=1./Znorn;
  Ynorn=Ynorn(Pusff);

  fkttn=KKt./[SQ; SQ];
  fkttn=fkttn(Pusff);
 else
%  fkttn=KKt;
  fkttn=KKt(Pusff);
  Ynorn=1;
 end
%'ferma cont', keyboard

 Are0=diag(pes)*(Acoz(puuu1,:)+Acoz(puuu2,:));
 Arh0=segdir*diag(pes)*(Acoz(puuu1,:)-Acoz(puuu2,:));

 Are=Are0;
 Arez=(-2*j)*(-segem)*fatEz*diag(fkttn)*Arh0;

% Arh=diag(Ynorn)*Arh0*diag(Nz);
% Arhz=(-2*j)*diag(Ynorn.*fkttn)*Are0*diag(Nz);

 Arh=diag(Ynorn)*Arh0;
 Arhz=(-2*j)*diag(Ynorn.*fkttn)*Are0;

 isezE=1;
 Aue1=fatcut.*Are(:,isezE);
 Auez1=fatcut.*Arez(:,isezE);
 Auh1=fatcut.*Arh(:,isezE);
 Auhz1=fatcut.*Arhz(:,isezE);


if ionlyout==0
 Aue2=fatcut.*Are;
 Auez2=fatcut.*Arez;
 Auh2=fatcut.*Arh;
 Auhz2=fatcut.*Arhz;

 Aue3=fatcut.*Are(:,end);
 Auez3=fatcut.*Arez(:,end);
 Auh3=fatcut.*Arh(:,end);
 Auhz3=fatcut.*Arhz(:,end);
end
% keyboard

% dK=max(diff(KKv));
 KKtdu=KKt(Pusas);



  Mfm=[];
  Mgm=[];
  Mfp=[];
  Mgp=[];
  Mez =[];
  Mhz =[];

   Aue21=[];
   Aue22=[];
   Aue23=[];
   Auez21=[];
   Auez22=[];
   Auez23=[];
   Auh21=[];
   Auh22=[];
   Auh23=[];
   Auhz21=[];
   Auhz22=[];
   Auhz23=[];

%if ionlyout==0
   squv=[];
   sqbv=[];
   fkv=[];
   Zdv=[];
   Kfit=[];
   besp=[];
   besm=[];
   besz=[];
   Mvefm=[];
   Mvefp=[];
   Mvegm=[];
   Mvegp=[];
   Mvez =[];
   Mvhz =[];

   kref=2*pi/lambdam*rr;
%    shi=length(Aue1)/2;
    shiv=[0 length(Aue1)/2];
   
  for shi=shiv

   kin=0;
   kint=0;
  
   while kin<numo
    [shi, kin]
    ialt=0;
    kin=kin+1;
    kint=kint+1;
    if iLP==0
     ibm=pimu(kin)-1;
     ibp=pimu(kin)+1;
     ibz=pimu(kin);
    else
     ibm=pimu(kin);
    end
    puk=(Psep(kint,:));
    Ke=KKtdu(puk);
    if Ke(1)==min(KKv)
     Ke(1)=0;
    end
    Npn=ceil(diff(Ke)/dK)+1;
    Npo=diff(puk)+1;
    if Npo<Npn
     Kad=linspace(Ke(1),Ke(2),Npn)';
     Np=Npn;
    else
     Kad=linspace(Ke(1),Ke(2),Npo)';
     Np=Npo;
    end
    if fadK==1
     Kad=KKv;
    end
    Kfit=[Kfit; Kad];
    
    lim=1-(rr*Kad/rff).^2;
    SQ=sqrt(lim);
    SQV=[SQ];
    squv=[squv; SQV];

    lim=1-(rr*Kad/rffb).^2;
    SQ=sqrt(lim);
    SQV=[SQ];
    sqbv=[sqbv; SQV];
   %'prima ', keyboard
    
%    Zd=[1./SQ; SQ];
%    fkttn=[Kad; Kad]./[SQ; SQ];    
%    fkv=[fkv; fkttn];
%    Zdv=[Zdv; Zd];

     sub_fiMod
   %'dopo ', keyboard

    if length(find(nudiv-kin)==0)==1 & numo>1 & ikr==2
     ialt=0;
     kint=kint+1;
     puk=(Psep(kint,:));
     Ke=KKtdu(puk);
     Np=ceil(diff(Ke)/dK)+1;
     Kad=linspace(Ke(1),Ke(2),Np)';
     Kfit=[Kfit; Kad; Kad];
     
         lim=1-(rr*Kad/rff).^2;
         SQ=sqrt(lim);
         SQV=[SQ];
         squv=[squv; SQV];
     
         lim=1-(rr*Kad/rffb).^2;
         SQ=sqrt(lim);
         SQV=[SQ];
         sqbv=[sqbv; SQV];
     
     sub_fiMod
%   'dopo ', keyboard

    end
    
   end  %fine kin
   
 end %shi  TE TM
%end   

% figure, subplot(211), plot(Kfit,'g.'),
% title(' new K ')
% subplot(212), plot(KKt(Pusas(puP)),'r.'),
% title(' old K ')
% pausak
KKt=Kfit;
' fine subfi', keyboard

if ionlyout==0
% sotto
 beta0=2*pi/lambdam*rffb;
 sqr=sqrt(1-(KKt*rr/rffb).^2);
 rrfa=1;
 lim=1-(rrfa*Kfit).^2;
 SQ=sqrt(lim);
 if iLP==0
  Znorn=[1./SQ; SQ];
  Ynorn=1./Znorn;
  %Ynorn=Ynorn(Pusff);
  fkttn=[Kfit; Kfit]./[SQ; SQ];
%  fkttn=fkttn(Pusff);
 else
%  fkttn=KKt(Pusff);
  fkttn=Kfit;
  Ynorn=1;
 end

 
  clear Ae Aez  Ah Ahz
  sqr=sqbv;
  
  for ks=1:length(sectb)
   Lu=sectb(ks);
   fae=-j*beta0*Lu*sqr;
   fifa=find(real(fae)>0);
   fae(fifa)=-fae(fifa);
   Axp=exp(fae);
 
   Adu=(Aue23.*Axp);
   Ae(:,ks)=Adu;
   Adu=(Auh23.*Axp);
   Ah(:,ks)=Adu;
   Adu=(Auez23.*Axp);
   Aez(:,ks)=Adu;
   Adu=(Auhz23.*Axp);
   Ahz(:,ks)=Adu;
  end
end

  if length(sectb)>0
  Aeb=Ae;
  Ahb=Ah;
  Aezb=Aez;
  Ahzb=Ahz;
  else
  Aeb=[];
  Ahb=[];
  Aezb=[];
  Ahzb=[]; 
  end

%sopra
 beta0=2*pi/lambdam*rff;
% sqr=sqrt(1-(KKt*rr/rff).^2);
 sqr=squv;
 clear Ae Aez  Ah Ahz
 for ks=1:length(sectu)
  Lu=sectu(ks);
  fae=j*beta0*Lu*sqr;
  fifa=find(real(fae)>0);
  fae(fifa)=-fae(fifa);
  Axp=exp(fae);
  Adu=(Aue21.*Axp);
  Adue=(Aue21.*Axp);
  Ae(:,ks)=Adu;
  Adu=(Auh21.*Axp);
  Aduh=(Auh21.*Axp);
  Ah(:,ks)=Adu;
  Adu=(Auez21.*Axp);
  Aez(:,ks)=Adu;
  Adu=(Auhz21.*Axp);
  Ahz(:,ks)=Adu;
  
 end
'dopo z', keyboard

  if iLP==0
%   besp=[besp; besp];
%   besz=[besz; besz];
%   besm=[besm; besm];

   Azue=[Ae Aue22 Aeb ];
   Azuh=[Ah Auh22 Ahb];
   Azuhz=[Ahz Auhz22 Ahzb];
   Azuez=[Aez Auez22 Aezb];
  end

%'cont Axu', keyboard

  if iLP==1

   if exist('if90')
    bes_an=diag(Mvefm0(:,if90))*besm;
   else
    bes_an=besm;
   end
   Esr=(bes_an'*Azu);
   Es=abs(Esr/max(max(Esr))).^2;
   tit=' LP ';

  else  %vettoriale

%   Mfm=[ Mvefm;  segem*Mvefm];
%   Mgm=[ Mvegm;  segem*Mvegm];
%   Mfp=[ Mvefp; -segem*Mvefp];
%   Mgp=[-Mvegp;  segem*Mvegp];
%   Mez =[ 0*Mvez;         Mvez];
%   Mhz =[   Mvhz;       0*Mvhz];

    bes_amfs=diag(Mfm(:,if180))*fliplr(besm);
    bes_amgs=diag(Mgm(:,if180))*fliplr(besm);
    bes_apfs=diag(Mfp(:,if180))*fliplr(besp);
    bes_apgs=diag(Mgp(:,if180))*fliplr(besp);
    bes_aezs=diag(Mez(:,if180))*fliplr(besz);
    bes_ahzs=diag(Mhz(:,if180))*fliplr(besz);

    bes_amfd=diag(Mfm(:,if0))*besm;
    bes_amgd=diag(Mgm(:,if0))*besm;
    bes_apfd=diag(Mfp(:,if0))*besp;
    bes_apgd=diag(Mgp(:,if0))*besp;
    bes_aezd=diag(Mez(:,if0))*besz;
    bes_ahzd=diag(Mhz(:,if0))*besz;

    bes_amf= [  bes_amfs   bes_amfd ];
    bes_amg= [ bes_amgs  bes_amgd ];
    bes_apf= [   bes_apfs   bes_apfd ];
    bes_apg= [  bes_apgs  bes_apgd ];
    bes_aez= [  bes_aezs   bes_aezd ];
    bes_ahz= [   bes_ahzs  bes_ahzd ]; 

   Ex=Azue.'*(bes_apf+bes_amf);
   Ey=Azue.'*(bes_apg+bes_amg);
   Hy=Azuh.'*(bes_apf+bes_amf);
   Hx=-Azuh.'*(bes_apg+bes_amg);


   if iztm==1
    Ez=Azuez.'*bes_aez;
    Hz=Azuhz.'*bes_ahz;
   else
    Ez=0;
    Hz=0;
   end

' dopo campo ', keyboard

   Pvectz=(Ex.*conj(Hy)-Ey.*conj(Hx));
   Pvectx=(Ey.*conj(Hz)-Ez.*conj(Hy));
   Pvecty=-((Ex.*conj(Hz)-Ez.*conj(Hx)));
   if iex==0
    Pointingr=sqrt(real(Pvectx).^2+real(Pvecty).^2+real(Pvectz).^2);
    Pointingi=sqrt(imag(Pvectx).^2+imag(Pvecty).^2+imag(Pvectz).^2);
   elseif iex==1
    Pointingr=real(Pvectx);
    Pointingi=imag(Pvectx);
   elseif iex==2
    Pointingr=real(Pvecty);
    Pointingi=imag(Pvecty);
   elseif iex==3
    Pointingr=real(Pvectz);
    Pointingi=imag(Pvectz); 
   end 
   
%   Es=Pointing';
%   Es=Es'/max(max(Es));
   if ipoi>=1
    if ipoi==2
     Edu=Pointingi;
      tit= 'Stored energy ';
    else
     Edu=Pointingr;
      tit= 'Power flux ';
    end
%    Edu=Pvectz;
%      tit= 'Pointing ';
   else
    if iex==0
      Edu=sqrt(abs(Ex.^2)+abs(Ey).^2);
   elseif iex==1
     if ieh==1
      Edu=Ex;
      tit= '|Ex|^2 ';
     else
      Edu=Hx;
      tit= '|Hx|^2';
     end
    elseif iex==2
     if ieh==1
      Edu=Ey;
      tit= '|Ey|^2 ';
     else
      Edu=Hy;
      tit= '|Hy|^2 ';
     end
    elseif iex==3
     if ieh==1
      Edu=Ez;
      tit= '|Ez |^2';
     else
      Edu=Hz;
      tit= '|Hz |^2';
     end
    end

    Edu=abs(Edu/max(max(Edu))).^2;
   end
  
   Es=Edu.';
   Es=Es/max(max(Es));
  end
fi0=find(abs(Es)==0);
Es(fi0)=1e-10;

 fiz=1:length(ztot);
 [du,pdu1]=min(abs(xvero-aiat/2));
 [du,pdu2]=min(abs(xvero-aiat));

 [du,fiz0]=min(abs(ztot));
% figure
 if exist('soloplot')
  if soloplot==1
   return
  end
 end

if isalf==0
 isez=[1 pdu1 pdu2]+length(xvero);
% isez=1;
if ipau==1
 figure, semilogy(ztot(fiz),Es(isez,fiz),'.-'),
 pausak
 figure, plot(ztot(fiz),Es(isez,fiz),'.-'),
 xlabel(' zeta toward device (micron)')
 ylabel(' Es ')
 title([tit,' in dis-fz: 3 radial coord: 0, Ox/2, Ox '])
 pausak
end
 Nor=max(abs(Es(1,:)));

isa=1;
if isa==0
 Ex=Ex.';
 Ey=Ey.';
 Ez=Ez.';
 Hx=Hx.';
 Hy=Hy.';
 Hz=Hz.';
 fiz=[187 188];
 fiz=[267-40:268+40];
 figure, plot(ztot(fiz),real(Ex(isez,fiz)),ztot(fiz),imag(Ex(isez,fiz))),
 title(' Real, imag Ex')
 pausak
 figure, plot(ztot(fiz),real(Hy(isez,fiz)),ztot(fiz),imag(Hy(isez,fiz))),
 title(' Real, imag Hy')
 pausak
 figure, plot(ztot(fiz),real(Ex(isez,fiz)),ztot(fiz),real(Hy(isez,fiz))),
 title(' Real ')
 pausak
 figure, plot(ztot(fiz),imag(Ex(isez,fiz)),ztot(fiz),imag(Hy(isez,fiz))),
 title(' Imag ')
 pausak
 figure, plot(ztot(fiz),real(Ex(isez,fiz).*conj(Hy(isez,fiz)))),
% figure, plot(ztot(fiz),imag(Ex(isez,fiz).*conj(Hy(isez,fiz)))),
 title(' Ex*Hy')
 pausak
end

xvero_tot=[-fliplr(xvero) xvero];
mac=max(xvero);
% map(Es(:,fiz)*diag(1./max(Es(:,fiz))),ztot(fiz),xvero)

mapxzHCG_Mayer

caxis([-4 0])
pausak
caxis([-2 0])
pausak
%mapxzoo


if iex==1
Edu=Ex;
tit='Ex^2';
else
Edu=Ey;
tit='Ey^2';
end

end
% fiz0=find(ztot==0);
% figure
% plot(xvero,Es(:,fiz0)*diag(1./max(Es(:,fiz0))))
% pausak
%' dis_fzt', keyboard

 fatEzv=(rr./Nztot).^2;

 h=figure;
set(h,'pos',[ 140    400   300*length(zcut)   250])
% assi=linspace(5,27,9)

 for k=1:length(zcut)
%  subplot(1,length(zcut),length(zcut)+1-k)
  subplot(1,length(zcut),k)
  [du,ip]=min(abs(zcut(k)-ztot));


   Aue=Azue(:,ip);

  if iLP==1

   Mfme=diag(Aue)*Mvefm0;
   Esr=besm'*Mfme;
%   Edu=abs(Esr/max(max(Esr))).^2;
   Edu=abs(Esr).^2/Nor;

  else  %vettoriale
  
   Aue=Azue(:,ip);
   Auez=Azuez(:,ip);
   Auh=Azuh(:,ip);
   Auhz=Azuhz(:,ip);

   Mfme=diag(Aue)*Mfm;
   Mgme=diag(Aue)* Mgm;  
   Mfpe=diag(Aue)*Mfp;     
   Mgpe=diag(Aue)*Mgp;    

   Mfmh=diag(Auh)* Mfm;   
   Mgmh=diag(Auh)*Mgm;   
   Mfph=diag(Auh)*Mfp;     
   Mgph=diag(Auh)*Mgp;    

   Mezi =diag(Auez)*Mez;        
   Mhzi =diag(Auhz)*Mhz;        

%   Mvegp=diag(Az)*Mgp;
%   Ex=besp'*Mvefp+besm'*Mvefm;
   Ex=besp'*Mfpe+besm'*Mfme;
   Ey=besp'*Mgpe+besm'*Mgme;

   Hy=besp'*Mfph+besm'*Mfmh;
   Hx=-(besp'*Mgph+besm'*Mgmh);


   if iztm==1
    Ez=besz'*Mezi;
    Hz=besz'*Mhzi;
   else
    Ez=0;
    Hz=0;
   end
   Pvectz=(Ex.*conj(Hy)-Ey.*conj(Hx));
   Pvectx=(Ey.*conj(Hz)-Ez.*conj(Hy));
   Pvecty=-((Ex.*conj(Hz)-Ez.*conj(Hx)));
   Pointingr=sqrt(real(Pvectx).^2+real(Pvecty).^2+real(Pvectz).^2);
   Pointingi=sqrt(imag(Pvectx).^2+imag(Pvecty).^2+imag(Pvectz).^2);

   if ipoi>=1
    if iex==0
     Pointingr=sqrt(real(Pvectx).^2+real(Pvecty).^2+real(Pvectz).^2);
     Pointingi=sqrt(imag(Pvectx).^2+imag(Pvecty).^2+imag(Pvectz).^2);
    elseif iex==1
     Pointingr=real(Pvectx);
     Pointingi=imag(Pvectx);
    elseif iex==2
     Pointingr=real(Pvecty);
     Pointingi=imag(Pvecty);
    elseif iex==3
     Pointingr=real(Pvectz);
     Pointingi=imag(Pvectz); 
    end 
    if ipoi==2
     Edu=Pointingi;
      tit= 'Stored energy ';
    else
     Edu=Pointingr;
      tit= 'Dissipated/Radiated energy ';
    end    
    
   else
    if iex==1
     if ieh==1
      Edu=Ex;
     else
      Edu=Hx;
     end
    elseif iex==2
     if ieh==1
      Edu=Ey;
     else
      Edu=Hy;
     end
    elseif iex==3
     if ieh==1
      Edu=Ez;
     else
      Edu=Hz;
     end
    end

   end
%   'edu', keyboard
    if k==1
    Nor=max(max((Edu)));
    end
%    Edu=(Edu/max(max(Edu)));
Edus=Edu;
    Edu=(Edu/Nor);

  end

  Esez{k}=Edu;
%  map(log10(abs(Edu)),YP,XP,h,-1)
  map(Edu,XP,YP,h)
%  axis([-1 1 -1 1]*assi(k))
  title(' real part' )

  axis equal
  axis square
    axis([-1 1 -1 1]*assi)
%    caxis([-1 0])
% if ipoi==0 & length(zcut)<4
%  pausak
%  map(imag(Edu),YP,XP)
%  title(' imag part' )
%  axis([-1 1 -1 1]*assi)
%  axis equal
%  pausak
% end
  %if k==3
  title([' Section = ',num2str(-zcut(k)),' \mum'])
  %else
  %title([num2str(-zcut(k))])
  %end
 end

if ipau==0
 close all
end
end


dx=[0 diff(xvero_tot)];
xdx=dx.*xvero_tot;
dz=[0 diff(ztot)];
if iex==3
figure, 
subplot(211)
plot(ztot,xdx*Es)
subplot(212)
semilogy(ztot,abs(xdx*Es))
else
figure, 
subplot(211)
plot(xvero_tot,Es*dz')
subplot(212)
semilogy(xvero_tot,abs(Es*dz'))
end