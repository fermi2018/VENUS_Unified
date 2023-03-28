%' prima di inv_mix ' ,keyboard
%if ifr==1
%' prima di inv_mix ' ,keyboard
%end
%if exist('gpla')
%if ~exist('gpla')
% gpla=gth1d*NQW_ef;
%end
 %if gpla<2000
% gplam=gpla;
 %else
 %gplam=0;
 %end
% 'qui g0', keyboard

finatt=find(real(niat)>0);
if length(finatt)>1
 ipilat=1;
 apilat=aiat;
 if ifp==-10
 ' da sistemare il case attivo con altro .... per ora ae posto infinito '
 end
%  'Moi'
%keyboard 
end


%if ifr==1
%if ifr<1000
ra=niat(1);
global del_n_ag ianti_gui del_n_agInf
 if ianti_gui==1
%  ra_centrale=niat(1)+del_n_ag*NQW*fatqw;
  ra=niat(1)+del_n_agInf*NQW*fatqw;  
 end 
% DGp=-g0*1e9;  %in ampiezza
 pqw=-imag(niat(1))*k0*1e4; 
 
 DGp=g0;  %in ampiezza

ra_centrale=ra;
%' contre g0', keyboard


if iplan~=1
%
%' prima di inv_mix ' ,keyboard

  shaat=shav.a(1); 
  if pol==1
   MK=MKospd;
    if iztm==1
     MKz=MKoszpd;
    end
  else
   MK=MKosmd;
   if iztm==1   
    MKz=MKoszmd;  
   end	
  end  
  fiat=find(aytot(:,shaat)==aiat(1));
  KA=squeeze(MK{shaat}(:,:,fiat));
 
  Kosi=KA;
  if iztm==1
   KAz=squeeze(MKz{shaat}(:,:,fiat)); 
   Koszi=KAz;
   KOt=Kosi*diag(pes);
   KOz=(rr/ra_centrale)^4*Koszi*diag(pes);


   Ksu=KOt+KOz;
   Kdi=KOt-KOz;
   P=[Ksu Kdi; -Kdi -Ksu];
  else
   KOt=Kosi*diag(pes);
   P=[KOt KOt; -KOt -KOt];
  end
  Moi0=P(Pust,Pust);
  %'Moi', keyboard
 

 if ipilat==0 

  dea=(ra^2-rr^2)/rr^2;
  deaz=1-(rr/ra)^2; 
  drc=ra_centrale-ra;
  dea_c=(2*ra+drc)*drc/rr^2;
  deaz_c=(rr/ra_centrale)^4*dea_c;    
  ck=-j*kcav*d*1e6;
  pI=ck*(dea*Idelta+deltany*Ideltaa/4);  % fattore 4 per la normalizzazione
  pIn=ck*dea_c* KA*diag(pes);


    if ~exist('KA_ag')
     KA_ag=KA;
     KAz_ag=KAz;
    end
    iplat=-10;
     if ianti_gui==1
      iplat=0;
     end 

%      ' del_eps ', keyboard
%      keyboard
      if iLP==0
       pInz=ck*deaz_c*KAz*diag(pes);
      else
       pInz=0;
      end
      pI=pI+pIn;
      pIz=ck*deaz*Ideltaz+pInz;


  Belo=-j*be*kcav*d*1e6*fatd;


  Ksu=pI+pIz*flp;
  Kdi=pI-pIz*flp;
  Pk=[Ksu Kdi; -Kdi -Ksu];
%  Pe=P+diag(([Belo; -Belo]));
%  iama=[];
%  Mon=expm_mio(Pe,iplat,iama,Pust);
  P=Pk+diag(([Belo; -Belo]));
  %'ferma Pk', keyboard
%   [du,D]=eig(P);
% Ood = du * diag(exp(diag(D))) / du;
%  [edu,eiga]=eig(P);
%  evet=diag(eiga);
%  ei=imag(evet);
%  fi=find(ei<0);
%  eim=ei(fi);
%  figure, plot(ei(fi))
%  fi=find(ei>0);
%  eip=ei(fi);
%  figure, plot(ei(fi))

  iama=[];
  Mon=expm_mio(P,iplat,iama,Pust);
%  'Mon prima', keyboard  
%  Moi0=Moi;
  Moi=Mon*Moi0;       % esce dal calcolo esatto della linearizzazione
%  Moi=Moi0*Mon;       % esce dal calcolo esatto della linearizzazione
%  'Moi, Mon', keyboard  
 else

%'questo ramo in inv_mix per ipilat (Mesa)', keyboard
  fiat0=find(aitot~=0);
  aitd0=aitot(fiat0);
  fiat=find(aitd0==apilat(1));
  istatt=0;
  ipri=fiat(1);
  istr=1;
  Li=d*1e6;
  ni=niat([1 3:end]);
%  ni(2)=ni(2)+j*imag(ni(1));
  ni(1)=ra;
  ai=aiat([2:end]);
%  ai=ai*0;
%  ni(2:end)=0;

  bi=aral.y.a';
  pai=aral.p.a';
  shai=shav.a';
  tyari=aral.ar.a';
%  'nattivo', keyboard
  bi=bi(([1 3:end]));
  pai=pai(([1 3:end]));
  shai=shai(([1 3:end]));


 if iany==0
  Kan=0;
  ani=zeros(size(Li));
 elseif iany==1
  Kan=Kosan;
  ani=aniat*fattany;
 elseif iany==2
  Kan=Kosan*diag(pes);
  ani=aniat*fattany;
 end


 if ianys~=0
  if ianys-iany==0
   ani=ani+aniats;
   ani1=0;
   Kan1=0;
  else
   if ianys==1
    Kan1=Kosan1;
    ani1=aniats;
   elseif ianys==2
    Kan1=Kosan1*diag(pes);
    ani1=aniats;
   end
  end
 else
  Kan1=0;
  ani1=0;
 end

  istatt


%  'passo emme', keyboard
  
  

%      iplat=0;
%      global del_n_ag
%      del_eps_ag=(2*ra*del_n_ag+del_n_ag.^2)/rr^2;
   DelN=del_n_ag*NQW_eq;
   del_eps_ag=real(2*DelN*ra+DelN^2)/rr^2;
   del_eps_agZ=real((2*DelN*ra+DelN)/ra^2);
      
      
  if exist('KA_ag')
   K_antiguiding=  del_eps_ag*KA_ag*diag(pes);
   K_antiguidingZ= del_eps_agZ*KAz_ag*diag(pes);
  end


  ilaymem=0;
  layskip=0;
  ilayfast=0; 

% tolgo tutte anisotropie da zona attiva 
  ani=0;
  ani1=0;
  anir=0;
  dovesono=3;  %zona attiva
  
 % 'passo emme', keyboard  
   
   eval(emme)
   
%   'fine emme', keyboard
   K_antiguiding=0;
   K_antiguidingZ=0;   

   istatt=0;


   Mon=Oo;
   Moi=Mon*Moi0;

%disp(' dopo attivo '), pausak

%  eval(emme)
   if ifp==-11
    disp(' verifica att. '), keyboard
   end 
 end

else  %iplan =1

  if idyn==1
    dea=(real(ra)^2-rr^2)/rr^2+2j*imag(ra)/rr;
    ck=-j*kcav*d*1e6;
    pI=ck*(dea*Idelta+deltany*Ideltaa/4);  % fattore 4 per la normalizzazione
    global del_n_ag ianti_gui

    if ~exist('KA_ag')
     KA_ag=KA;
     KAz_ag=KAz;
    end
    iplat=-10;
     if ianti_gui==1
%      iplat=0;
%      global del_n_ag
%      del_eps_ag=(2*ra*del_n_ag+del_n_ag.^2)/rr^2;
      del_eps_ag=2*del_n_ag/rr*NQW_eq;
      pIn=ck*del_eps_ag*KA_ag*diag(pes);
%      ' del_eps ', keyboard
      if iLP==0
       pInz=ck*del_eps_ag*KAz_ag*diag(pes);
      else
       pInz=0;
      end
      pI=pI+pIn;
     else
      pInz=0;
     end

    pIz=ck*dea*Ideltaz+pInz;
    Belo=-j*be*kcav*d*1e6*fatd;


    Ksu=pI+pIz*flp;
    Kdi=pI-pIz*flp;
    P=[Ksu Kdi; -Kdi -Ksu];
    P=P+diag(([Belo; -Belo]));
    iama=[];
    Mon=expm_mio(P,iplat,iama,Pust);

  else %idyn

   DGp=-g0*1e9;  %in ampiezza
   DGp=-g0;  %in ampiezza
   
   ra0=(ra);
   ck=-j*kcav*d*1e6;
   Belo=-j*be*kcav*d*1e6;   
   ded=(ra0^2-rr^2)/rr^2;
   dedz=1-(ra/rr)^2;
   pI=ck*ded*Idelta;
   pIz=ck*dedz*Ideltaz;
   Ksu=pI+pIz*flp;
   Kdi=pI-pIz*flp;
   Psal=[Ksu Kdi; -Kdi -Ksu];

   P=diag([Belo; -Belo])+Psal;
   
   Mon=expm(P(Pust,Pust));
' qui mon', keyboard
  end %idyn
' qui mon', keyboard
 if iztm==1
  KOt=Idelta;
  KOz=(rr/ra).^4*Ideltaz;
  Ksu=KOt+KOz;
  Kdi=KOt-KOz;
  P=[Ksu Kdi; -Kdi -Ksu];
 else
  KOt=Idelta;
  P=[KOt KOt; -KOt -KOt];
 end

 Moi=Mon*P(Pust,Pust);
 %'Altro Moi', keyboard
end


if ifiez==1
%' prima di invpro ' ,keyboard
 Mois=Moi;
 Mons=Mon;
 return
end
if ifp==-10
 % disp(' prima di invpro '), keyboard
end
  global ired_ret
%  'cont Pured', keyboard
 if length(ired_ret)==1
  if ired_ret==1
  
   Tts=Tt;
   Tbs=Tb;
   Ga1s=Ga1;
   Ga2s=Ga2;
   Gass=Gas;
   Gads=Gad;
   Mons=Mon;
   Mois=Moi;
   Moi=Moi(Pured,Pured);
   Mon=Mon(Pured,Pured);
   Tt=Tt(Pured,Pured);
   Tb=Tb(Pured,Pured);
   if pol==pvet(1)
    Ga1=Ga1(Pured1,Pured1);
    Ga2=Ga2(Pured1,Pured1);
    Gas=Gas(Pured1);
    Gad=Gad(Pured1);
    pu1=1:length(Gas);
    pu2=length(Gas)+1:2*length(Gas);
   end 
  end
 end
 
% 'Mn Optica', keyboard
 
% load sa1
% Mon=MonB;
% Moi=MoiB;
 
  Mn=Tb*Mon*Tt;
  Mi=Tb*Moi*Tt;
  
   %'dopo Mn Optica', keyboard
  
%Tn=Td*Tan*Ts;
%Ti=Td*Tai*Ts;  


%	Lt=length(Tb)/2;
%	pu1=1:Lt;
%	pu2=pu1+Lt;
%	T11=Tb(pu1,pu1);
%	T12=Tb(pu1,pu2);
%	T21=Tb(pu2,pu1);
%	T22=Tb(pu2,pu2);
%	Gi=inv(T22-Ga2*T12)*(Ga2*T11-T21);

%' verifica Tt Tb', keyboard

%  verifica conversione interfacce
%   I=eye(size(Ga1));
%   V1=[I Ga1; Ga1 I]; 
%   V2=[I -Ga2; -Ga2 I]; 
%   MnR=V2*Mn*V1;

%  'Mi, Mn', keyboard
 if isfield(Ps,'istop_fie')==1
  if Ps.istop_fie==1
%   ' Mn Mi ', keyboard
  end
 end 
  if ifr>=2
%   ' Mn Mi ', keyboard
  end

if ifp==-11
  disp(' prima di invpro '), keyboard
end
%  clear P Pe



  invpro

if ifp>=-1,  pausak, end
if ifr==1
% ' dopo di inv_mix ' ,keyboard
end