clear all
close all

Temp0=300; % temperatura bagno termico (per correzione formula)


 load datitermici

 Tprec=T3D;
 Ti=(1+Tprec/Temp0).^(-1.25);
 %Ti=1;

%keyboard

 %-- inizio Tibaldi
 ro=mesh.xgrid;
 z=mesh.ygrid;
 
 ro=ro*1e4;
 zf=z*1e4;
 
 Joule=reshape(mode.Joule,mesh.nny,mesh.nnx);
 ind=find(zf>9);
 Joule(ind,:)=0; 
 Joule=Joule.';

 Rec_nr=reshape(mode.Rec_nr,mesh.nny,mesh.nnx);
 Rec_nr=Rec_nr.'; 
 
 %-- fine Tibaldi


 T_DD=max(zf);

ifig=0;     % =1 plot risultati intermedi


Qtot0=5e-5;  % caso forzato vecchio: valore massimo di Qtot

% spessori z

Tbuf=350; % spessore buffer (um)
%Tbuf=5;

fsplit=.8;	% punto di divisione nel buffer per ridurre funzioni di base

TBuf_DD=1;
Tdbr_inf=4.7; % spessore specchio dbr inferiore (um)
Tcav=.6; % spessore cavita'
Tdbr_sup=2.721;  % spessore specchio dbr superiore
Tmetallo=.2; % spessore contatto metallico
Tmetallo=0; % spessore contatto metallico
Tcentrale=T_DD-Tdbr_inf-Tmetallo-Tdbr_sup-TBuf_DD;

Tcav=Tcentrale;
%'spessori', keyboard


% dimensioni radiali

ipass=0;    % se 1 metto passivazione sui bordi del mesa in problema termico (altrimenti aria)
% variabili di ingresso: geometria radiale; mesa, passivazione (facoltativo), e ro_max
ro_pass=22; % questi ro sono tutti radiali (um)
ro_met=3;
ro_mesa=20;
ro_max=120;
Rox=2;

% parametri sorgente di calore qtot
Ri0=10; % semiasse radiale dell'ellisse
El0=0.5; % distorsione dell'ellisse
icir=1 % se 1 e` circolare
NPfi=20; %numero punti azimutali (caso non circolare); si usano 4*NPfi+1

%icir=0; 


iads=1; % aggiunge strato relativo a isplit in z
%imedia=input(' Media cond. termica? [0/1] ' )
imedia=0;  % fa media conducibilita termica
itras=0;  %media trasversale
iaria=0;   % aggiunge strato d'aria sopra la struttura
th_aria=10;  % usato solo se iaria =1



iKcost=0; % se 1 considera costante a tratti la conducibilita`, se no mette l'autoconsistenza con T
iploCond=1
if ifig==0
 iploCond=0;
end

tic



nMax=2; % armoniche azimutali

nu_num_max=nMax;  % caso circolare
if icir==1
 nu_num_max=1;  %1 caso circolare
end



ispliro=1;  % divide intervallo in rho per migliorare set funzioni


roh=[ro_met ro_mesa ro_pass ro_max];   % raggi interfacce
if ipass==0
 roh=[ro_met ro_mesa ro_max];   % raggi interfacce
end


if ispliro==1
 IPspli=2;   % intervallo in rho da dividere
% roh=sort([roh roh(IPspli)*[.8 2]]); 
 roh=sort([roh roh(IPspli)*[ 2]]); 
end

Tvet=[Tdbr_inf Tcav Tdbr_sup];

imet=1;  % includo o meno contatto nel problema termico
if Tmetallo==0
 imet=0;
end




if imet==1
Th=[Tbuf*fsplit Tbuf*(1-fsplit) Tvet Tmetallo];  %spessori in z
else
Th=[Tbuf*fsplit Tbuf*(1-fsplit) Tvet];  %spessori in z
end
 iPa=3;

Lzd=sum(Th(1:end-1));

zza=[0 cumsum(Th)];
xxa=[0 roh];

z_act=sum(Th(1:iPa))+Th(iPa+1)/2;



Nr=50;
xx1=linspace(0,Rox*3/2,Nr);
xx2=linspace(Rox*3/2,roh(1),Nr);
xx3=linspace(roh(1),roh(2),Nr);
xx4=linspace(roh(2),roh(end),Nr);
xx=[xx1 xx2(2:end-1) xx3(2:end-1) xx4];

Nr=50;
xx1=linspace(0,Rox*3/2,Nr);
xx2=linspace(Rox*3/2,roh(1),Nr);
xx3=linspace(roh(1),roh(2)*1.2,Nr);
xx4=linspace(1.2*roh(2),roh(end),Nr);

xx=[xx1 xx2(2:end-1) xx3(2:end-1) xx4];



%L=60;  %numero funzioni Bessel espansione Temperature
%J=30;  % numero funzioni trigonometriche temperatura


nomet='provaAT';



load ParT  % da settare e sistemare; per ora vecchia struttura

fian=linspace(0,2*pi,4*(NPfi)+1);

ang=fian(1:NPfi+1); % [1 * 26] ,il 26^ è pi/2


% e' settato in calopt !!!!!!
%if (i2Dyn==1 & versione_top ~=0) % ------ > 2D misto



if nu_num_max==1
 ang=0;
 NPfi=0;     %circolare
end


% e' settato in calopt !!!!!!
%if (i2Dyn==1 & versione_top ~=0) % ------ > 2D misto

nu_num=[1:nu_num_max];
nu_vett=[0 nu_num(1:end-1).*2];



if length(ang)>1
 difan=diff(ang);
else
 difan=pi/2;
end 


 zzd=zza(2:end-1);
 
 if iaria==1
  zza=[zza zza(end)+th_aria];
  zzd=zza(2:end-1);
  'aria',
%  keyboard
 end
 

 if isfield(Bcond,'Condpas')==1
  Cond_pas=Bcond.Condpas;
  pasmes=Bcond.pasmes;
%  ipass=1;
 else
  Cond_pas=Cond_air;
end
 

  
 



   CondZc=Bcond.CondZc; 
   Cond_air=Bcond.Cond_air;
   CondZb=Bcond.CondZb;
   CondZmet=Bcond.CondZmet;
   CondZm=Bcond.CondZm;
   CondZmet=Bcond.CondZmet;
   CondTm=Bcond.CondTm;


 
   condz1=flipud([Cond_air;    CondZm; CondZc; CondZm;  CondZb]);
   condz2=flipud([CondZmet;  CondZm; CondZc; CondZm;  CondZb]);
   condz3=flipud([Cond_air; Cond_air; CondZc; CondZm; CondZb]);
  
   condt1=flipud([Cond_air; CondTm; CondZc; CondTm;  CondZb]);
   condt2=flipud([CondZmet; CondTm; CondZc; CondTm; CondZb]);
   condt3=flipud([Cond_air; Cond_air; CondZc;  CondTm; CondZb]);
   

if ispliro==0
  if ipass==1
    condz23=flipud([Cond_pas; Cond_pas; CondZc;  CondZm; CondZb]);
    condz=[condz1 condz2 condz23 condz3];
    condt23=flipud([Cond_pas; Cond_pas;  CondZc;  CondTm; CondZb]);
    condt=[condt1 condt2 condt23 condt3];
   else
    condz=[condz1 condz2 condz3];
    condt=[condt1 condt2 condt3];
  end
else
    condz13=flipud([CondZmet; CondZm; CondZc;  CondZm; CondZb]);
    condz23=flipud([Cond_pas; Cond_pas; CondZc;  CondZm; CondZb]);
    condz33=flipud([Cond_air; Cond_air; CondZc;  CondZm; CondZb]);

    condt13=flipud([CondZmet; CondTm; CondZc;  CondTm; CondZb]);
    condt23=flipud([Cond_pas; Cond_pas; CondZc;  CondTm; CondZb]);
    condt33=flipud([Cond_air; Cond_air; CondZc;  CondTm; CondZb]);   

  if ipass==1
    condzT=[ condz23 condz33];
    condtT=[ condt23 condt33];    
    condzT=[ condz23 condz33];
    condtT=[ condt23 condt33];        
   else
    condtT=[condt13 condt33];    
    condzT=[condz13 condz33];
    condtT=[condt33];    
    condzT=[condz33];    
  end
    condz=[condz1 condz2 condzT condz3];
    condt=[condt1 condt2 condtT condt3];  

end

  if imet==0
   condz=condz(1:end-1,:);
   condt=condt(1:end-1,:);
  end
  
  if iaria==1
   condt=[condt; ones(size(condt(end,:)))*Cond_air ];
   condz=[condz; ones(size(condt(end,:)))*Cond_air ];
  %'aria cond', keyboard
   %condz=
  end
  
  if iads==1
   condt=[ones(size(condt(end,:)))*CondZb; condt ];
   condz=[ones(size(condt(end,:)))*CondZb; condz ];
  end 

 if itras==1
  condt=repmat(condt(:,1),size(condz(1,:)));
  condz=repmat(condz(:,1),size(condz(1,:)));
%  condt=condt(:,1);
%  condz=condz(:,1);
 end
%  'cond', keyboard
 
 if imedia==1
%     condzm=mean(mean(condz)); 
%     condtm=mean(mean(condt)); 
%     condz(:,:)=condzm;
%     condt(:,:)=condtm;
 if iads==1
  condzm=mean(condz(2:end,1))
  condtm=mean(condt(2:end,1))
 else
  condzm=mean(condz(:,1))
  condtm=mean(condt(:,1)) 
 end

  condz(:,:)=condzm;
  condt(:,:)=condtm;
%  condz=repmat(condz(:,1),size(condz(1,:)));
%  condt=repmat(condt(:,1),size(condz(1,:)));

 end



  Lv=diff(zza);
  Nv=3*ones(1,(length(zza)-1));

  Nv(2:end-1)=8;   
  %Nv(end)=12;   
  Nv=5*ones(1,(length(zza)-1));
  %Nv(1)=3;
  
   zz=[];
   dz=[];
   for ki=1:length(zza)-1
    xini=zza(ki);
    L=Lv(ki);
    Nquad=10*Nv(ki)+1; %-- conservativo    
    [nodes,weights]=quadad('legen',1,Nquad);
    xi=L/2*nodes+L/2+xini;
    wx=L/2*weights;
%    'ver',keyboard
    fiZ{ki}=length(zz)+[1:length(xi)];
    zz=[zz xi];
    dz=[dz wx];

   end
%  'cont fine', keyboard
  
  
 
 [f,fp,fs,Nfun,MI,MIs,MIp]=SEM_Nver(Lv,Nv,zz);
 
 %[f_int,fp_int]=SEM_Nver(Lv,Nv,zzd);


%' proma di iLO', keyboard

% S=sin(kdisL'*zz');
% Su=sin(kdisL'*(zeta'+thick_bu));
S=f;
s2=f;
s2s=fp;
L=Nfun;
I=L;

if iKcost==1
 for kj=1:L
  for ki=1:L
   SIm{ki,kj}=squeeze(MI(kj,ki,:))';
   SIm2{ki,kj}=squeeze(MIp(kj,ki,:))';
  end
 end 
end
 
 if ispliro==0
  Nvr=[10 10 5];
  Nvr=[4 6 16];
  of=[0 0 0];
 else
  Nvr=[4 6 6 16 16];
  
  of=[0 0 0 0 0]; 
  nv0=10;
  Nvr=nv0*[1 1 1 1];
  Nvr(end)=2;
%  Nvr=[4 4  8 2];
  of=0*[2 1  0  0];   
 end
 
 % forse conviene mettere in testa il numero di funzioni che vuoi usare su ciascun patch
 
 
 
 
 
 Rv=diff(xxa);

tic
% ioldB=input(' vecchio metodo = ');

 
%  Nvr=[30];
%  Rv=xxa(end);
% [fb,xb,wb,NfunB,MIB,MIpB]=SEM_BasM(Rv,Nvr);
 
 
 %[f_int,fp_int]=SEM_Nver(Lv,Nv,zzd);


%' proma di iLO', keyboard

% S=sin(kdisL'*zz');
% Su=sin(kdisL'*(zeta'+thick_bu));




icnu=0;
Jsu=0;
 for nua=nu_vett
  icnu=icnu+1;
  [fb,fbp,xb,wb,NfunB,MIB,MIpB,MI1B,fiX]=SEM_Bessel(Rv,Nvr,of,nua);
  fM{icnu}=fb;
  fTrr{icnu}=fbp;
  fTrf{icnu}=fb*diag(1./xb)*nua;
  J=NfunB;  
  Js(icnu)=J;
  Jsu=Jsu+J;
  if iKcost==1
   for jp=1:J
    for jj=1:J
    JIm{icnu,jp,jj}=squeeze(MIB(jp,jj,:));
    JImp{icnu,jp,jj}=squeeze(MIpB(jp,jj,:));
    JIm1{icnu,jp,jj}=squeeze(MI1B(jp,jj,:));
    end
   end  
  end
 end
 
 R=[];
% angi=ang+difan(1)/2;
 angi=ang;
 Coa=[];
  for knu=1:length(nu_vett)
    Bd0=fM{knu}';
    Bd=[];
    for kfi=1:NPfi+1
     Bd=[Bd; Bd0*cos(nu_vett(knu)*angi(kfi))];
    end
    R=[R Bd];
  end 
  Ru=R;
  ro_out=xb;
  xx=ro_out;
  
  xdx=wb*2*pi;
  xdx0=wb;
  if length(nu_vett)>1
    xdx=[wb/2 repmat(wb,1,NPfi-1) wb/2]*difan(1)*4;
  end  
  
% if iKcost==0 
  Fcondt=zeros(length(zz),length(xx));
  Fcondz=Fcondt;
  for kx=1:length(fiX)
   px=fiX{kx};
   for kz=1:length(fiZ)
    pz=fiZ{kz};
    Fcondt(pz,px)=condt(kz,kx);
    Fcondz(pz,px)=condz(kz,kx);
   end
  end 
  
  Fcondt0=Fcondt;
  Fcondz0=Fcondz;
  
  Fcondt=Fcondt0.*Ti';
  Fcondz=Fcondz0.*Ti';
  
  [xxc,zzc]=meshgrid(xx,zz);  
  
  if iploCond==1
  figure, surf(xxc,zzc,log10(Fcondt)), view(2), shading flat, 
  title(' Transverse Thermal conducibility')
  pausak
  
  figure, surf(xxc,zzc,log10(Fcondz)), view(2), shading flat, 
  title(' Longitudinal Thermal conducibility')
  pausak
  end
%   'cond', keyboard  

  icnu=0;
  for nua=nu_vett
   icnu=icnu+1;
   fTr=fTrr{icnu};
   fTf=fTrf{icnu};
   fT=fM{icnu};
   Mulr0=fTr*diag(xdx0);
   Mulf0=fTf*diag(xdx0);
   Mulx0=fT*diag(xdx0);
   Mulz0=f*diag(dz);
   Mulpz0=fp*diag(dz);
%   Mulx=repmat(Mulx0,L,1);
%   Mulz=repmat(Mulz0,Js(icnu),1);
   pJL=L*Js(icnu);
   ic=0;
   clear MLi
   for kx=1:Js(icnu)
    for kz=1:L
     ic=ic+1;    
     MLtr=Mulz0*diag(f(kz,:))*Fcondt*diag(fTr(kx,:))*Mulr0';
     MLtf=Mulz0*diag(f(kz,:))*Fcondt*diag(fTf(kx,:))*Mulf0';
     MLz=Mulpz0*diag(fp(kz,:))*Fcondz*diag(fT(kx,:))*Mulx0';
     MLi(ic,:)=reshape(MLtr+MLtf+MLz,pJL,1);     
    end
   end 
     MLm{icnu}=MLi;  
     
  %'quib', keyboard  
  end  %nua
    
 %end %iKcost
 

 
%'fine iKost', keyboard

 
 J=Jsu;
 pJL=J*L; 
 
  icnu=0;
 for nua=nu_vett
   icnu=icnu+1;
    icont=0;

    Ji=Js(icnu);

   if iKcost==1
    clear ML
    for jp=1:Ji
%     kra=kdisJ(jp);
  
     for ip=1:L
  
      icont=icont+1;
       icont_int=0;
       clear Mdu
  
       for jj=1:Ji
  

         JI=JIm{icnu,jp,jj};
         JIp=JImp{icnu,jp,jj};
         if nua~=0
          JI1=JIm1{icnu,jp,jj}*nua^2;
         else
          JI1=0;
         end
%         krb=kdisJ(jj);         
  
        for ii=1:L

         icont_int=icont_int+1;
	          SI=SIm{ip,ii};
	          SI2=SIm2{ip,ii};  
 	          Dr=SI*condt*(JIp+JI1);
	          Dz=SI2*condz*JI;
%         Dz=SI2*condz*JI;
  
         Mdu(1,icont_int)=Dr+Dz;
        end  %jj
       end   %ii
      ML(icont,:)=Mdu;
%      'qui ML', keyboard
     end  %jp
    end   %ip
   else %iKcost 
    ML=MLm{icnu};
   end

    iML=inv(ML);
%   'qui inv', keyboard    
    if nua==0
     iMLt{icnu}=iML/(2*pi);
    else
     iMLt{icnu}=iML/pi;
    end
%      'fine calcolo 1 azimut', keyboard
   end %azimut
 % save tutto1

% ' exist iMLt ?', keyboard
 
 if exist('iMLt')==0  
   iMLt=iMLt_sa;
 end


 fr0=exp(-(ro_out/10).^2);
 fi=find(ro_out>20);
 fr0(fi)=0;
 
 nVa=2
 if icir==1
  El=0;
 else
  El=El0;
 end
 fr=[];
 rt=[];
 

 for kfi=1:NPfi+1
  ri0=Ri0*(1-El*cos(nVa*fian(kfi)));
  fr0=exp(-(ro_out/ri0).^2);
  fi=find(ro_out>20);
  fr0(fi)=0;
  fr=[fr fr0];
  rt=[rt ro_out];
 end
 
 
 fz=exp(-abs(zz-z_act)/2);
 fi=find(zz>Lzd | zz<Tbuf);
 fz(fi)=0;
 
 [xfm,yfm]=meshgrid(ro_out,zz);
 
 % sorgente (magari mettere all'inizio)
 %qtot=1e-5*fr'*fz;
 qtotv=Qtot0*fr'*fz;
 

 
 T_VCSEL=sum(Tvet);
 z_dd=zf(end);
 z_ddT=zz(end)-z_dd;
 Tr_buf=zf(end)-T_VCSEL;
 
 fi_dd=find(zz>=z_ddT);
 zTdd=zz(fi_dd)-z_ddT;
 fir=find(xx<=ro(end));
 xTdd=xx(fir);
 z=z_dd-fliplr(zf);
 qq=fliplr(Joule+Rec_nr);
% DT=flipud(DeltaT);
 fiN=find(isnan(qq)==1);
 qq(fiN)=0;
 %keyboard
 %fiP=find(z<3);
 %qq(:,fiP)=0;

 %roM=repmat(ro',1,length(z));
 %zM=repmat(z,length(ro),1);
 
 [roM,zM]=meshgrid(ro,z);
 qtotV=interp2(roM,zM,qq.',xTdd,zTdd').';

 
 qtot=zeros(length(xx),length(zz));
 qtot(fir,fi_dd)=qtotV;

 


 ' qui controllo termico', keyboard

% [zM,roM]=meshgrid(z,ro);
% qtot=interp2(roM,zM,qq,xTdd,zTdd');
 % cambio grigliato
 
% 'qtot', keyboard
 
 if itras==10
 qm=mean(qtot);
 qtotv=qtot;
 qtot=repmat(qm,length(ro_out),1);
 
% 'qtot', keyboard
 end
 
  Rp=R';
  for kd=1:J
   Rdx(kd,:)=Rp(kd,:).*xdx;
  end
  
  %dz=diff(ztot);
  %dz=[dz dz(end)];
  %dz([1 end])=dz([1 end])/2;  
  
  Pqs=qtot*(diag(dz)*S.');
  P=Rdx*Pqs;
  
  
   cijL=[];
   for inu=1:length(nu_vett)
    Ji=Js(inu);
    Jsh=sum(Js(1:inu-1));
    QL=reshape(P((1:Ji)+Jsh,:)',Ji*L,1);
    cijL=[cijL; iMLt{inu}*QL];
   end
   CoL=reshape(cijL,L,J)';

   T3D1=Ru*CoL*S;
   
   if icir==0
    T3D=reshape(T3D1,length(xx),NPfi+1,length(zz));
    T3Dt=squeeze(mean(T3D,2));
   else
    T3D=T3D1;
    T3Dt=T3D1;
   end 

Ti=(1+T3Dt/Temp0).^(-1.25);
Fcondt0=Fcondt;
Fcondz0=Fcondz;

T3D0=T3Dt;
Tprec0=T3Dt;
Tprec=T3Dt;
eRmin=1e-5;

col='mcrgw';
[du,fiat]=min(abs(zz-z_act));
h=figure;
set(h,'pos',[624         474        1098         504])
subplot(131)
plot(xx,Tprec(:,fiat)), 
ylabel('Temp rise (K)')
xlabel('ro (um)')
a=axis;
a(2)=2*ro_mesa;
a(4)=a(4)*1.2;
axis(a)
hold on
subplot(132)
plot(zz,Tprec(1,:)),
a=axis;
xlabel('z (um)')

a(1)=Tbuf*.95;
a(2)=zz(end);
a(4)=a(4)*1.2;
axis(a) 
subplot(133)
plot(zz,Tprec(1,:)),
xlabel('z (um)')
pausak

if iKcost==0
hold on
'qui prima di loop', keyboard
 
for kit=1:10
 Fcondt=Fcondt0.*Ti';
 Fcondz=Fcondz0.*Ti';
 
 subTEMP
   if icir==0
    T3D=reshape(T3D1,length(xx),NPfi+1,length(zz));
    T3Dt=squeeze(mean(T3D,2));
    Tat=squeeze(T3D(:,:,fiat));
   else
    T3Dt=T3D1;
    Tat=T3Dt(:,fiat);
   end 
   
   
   Ti=(1+T3Dt/Temp0).^(-1.25);
   su=sum(sum(Tprec));
   sui=sum(sum(T3Dt));
   er0=1-su/sui
   er(kit)=er0;
   'iter =', kit
   'errore relativo Temp ', 
   Tprec=T3Dt;
   figure(h);
   subplot(121)
   plot(xx,Tprec(:,fiat),col(kit)), 
   a=axis;
   a(2)=2*ro_mesa;
   axis(a)
   subplot(122)
   plot(zz,Tprec(1,:),col(kit)), 
   a=axis;
   a(1)=Tbuf*.95;
   a(2)=zz(end);
   axis(a)
pausak
   
   
   %keyboard
   if er0<=eRmin
    break
   end

end
else
   if icir==0
    T3D=reshape(T3D1,length(xx),NPfi+1,length(zz));
    T3Dt=squeeze(mean(T3D,2));
    Tat=squeeze(T3D(:,:,fiat));
   else
    T3Dt=T3D1;
    Tat=T3Dt(:,fiat);
   end 
   
end  %iKcost

if ifig==1

[du,fiat]=min(abs(zz-z_act));
figure, 
subplot(211)
plot(xx,Tprec0(:,fiat),xx,Tprec(:,fiat)), 
a=axis;
a(2)=2*ro_mesa;
axis(a)
subplot(212)
plot(zz,Tprec0(1,:),zz,Tprec(1,:)),
a=axis;
a(1)=Tbuf*.95;
a(2)=zz(end);
axis(a)
pausak

toc   
   figure
%      surf(xfm',yfm',squeeze(T3D(:,1,:))), shading interp, view(2), colorbar,
      surf(xfm',yfm',T3Dt), shading interp, view(2), colorbar,
%      title([' Cor = ',num2str(Cor)])
      pausak




figure, plot(ro_out,Tat,'.-')
title(' Temp zona attiva ')
pausak

end

Tbuf_aggiunto=max(max(zzd))-max(z);


 [xxd,zzd]=meshgrid(xx,zz);
 
 ro(1)=min(min(xxd));
 
 Tdd=interp2(xxd,zzd,T3D',ro,Tbuf_aggiunto+z').';

