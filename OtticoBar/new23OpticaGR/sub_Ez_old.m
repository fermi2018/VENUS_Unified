global Ps

ifp=ifpsal;
%'entro', keyboard
%iraffi=1;
%  pvet=pola;

  Lizi=Ps.Lizi;    % step interno
%  Liz0=Ps.Liz0;    % step interno
  PUA=pua;

 if isav_Az>0
  if iauto(1,1)==3
%  if iauto(1,1)==0
   fiau=find(iauto(:,1)==1);
   iauto(fiau,1)=0;
   iauto(1,1)=1;
   reassign
   Litn=Lr.t;
   nitn=nr.t.';
   aitn=ar.t.';

   fmlsp=abs(frp.t);
   fmlst=abs(frp.t);
   fi1=find(fmlsp(:,1)==0);
   fmlst(fi1,1)=1;
   fi2=find(fmlsp(:,2)==1);
   fmlst(fi1,2)=0;

   bi=aral.y.t';
   pai=aral.p.t';
   shai=shav.t';
   istfie=istfi.t;
   tyari=aral.ar.t';
   fmlsi=fmlst;

  end
iFie=1;  % non calcolo ORTA
ifphold=ifp;
%ifp=-4;
if ifp==-10
 hpro=figure
end 
  ifiez=0;
ifasasa=ifast;
 if iraffi==0 
  freq=ze;
  [du,ifr]=min(abs(fre_camp-ze));
%  ifr=ifr+1;
%  ifr=1;
%'  ifast ', keyboard

  ifast=0;
  kcav=kcav0*(1+freq);
  isaMo=1;
  pvet00=pola;
% ' adesso e corretto ', keyboard
%ifp=-111
      Pf.isemp_plan=0;
      isemp_plan=0;
      ifr=1;
      clear iloo
      ifast=0;
%      'prima acc0', keyboard
%if Ps.izoom==0
   tic
   acc0
   toc
%end
%      'dopo acc0', keyboard
  
%    ' salva pol=-1 ', keyboard
     clear Tstof Tstortp Tstortm Tstorbm Tstorbp
     lesa=size(Tstort,3);
     Tstortm{1}=0;
     Tstortp{1}=0;
     if lesa>0
      ipsat=1:lesa;
      for ip=ipsat
        if pol==-1
          Tstortm{ip}=Tstort(:,:,ip);
        else
          Tstortp{ip}=Tstort(:,:,ip);
        end
      end
     end      
     lesa=size(Tstorb,3);
     Tstorbm{1}=0;
     Tstorbp{1}=0;
     if lesa>0
      ipsat=1:lesa;
      for ip=ipsat
        if pol==-1
          Tstorbm{ip}=Tstorb(:,:,ip);
        else
          Tstorbp{ip}=Tstorb(:,:,ip);
        end
      end
     end        
      %'dopo Tsto', keyboard      
%fie_new
end

ifast=1;
%ifr=2;
  pua=PUA;
  [du,ilam]=min(abs(imag(Lamv(Nd1vu))));
  iN=Nd1vu(ilam);
  Lami=Lamv(iN);
  Ar=V(:,iN);
  Ar_sav=Ar;
%  if exist('Trcrit')
%    Ar=Trcrit*V(:,iN);
%  end
%  ' verifica Ar !!!!!!!! ',  keyboard
   Mto=Mns+Lami*Mis;
   Mn11=Mto(1:l1,1:l1);
   Mn12=Mto(1:l1,l2:l3);
   Mn21=Mto(l2:l3,1:l1);
   Mn22=Mto(l2:l3,l2:l3);
   dRi=Ga2*Mn12-Mn22;
   dPi=Mn21-Ga2*Mn11;
   Geq=inv(dRi)*dPi;

   A_in=Ga1*Ar;
   A_ref=Geq*A_in;
   if ifp==-10
%    figure(hpro), hold on, plot(abs(diag(Geq)),'r.'),
%    title(' Coeff. di riflessione equivalente: rosso con parte attiva ')
%    pausak
    figure, plot(abs(Ar)), hold on, plot(abs(A_ref),'r.'),
    title(' verifica autoconsistenza: giallo A_back; punti: A_back dopo un giro ')
    if ifp==-10
     keyboard
    end 
   end 
 end  %isav_Az

%if length(Tdis)>1 | isav_Az==2
if  isav_Az==2
   [du,ilam]=min(abs(imag(Lamv(Nd1vu))));
   iN=Nd1vu(ilam);
   Lami=Lamv(iN);
   ggg1=real(Lami)*c/(rr.*d*2);
   Ar=V(:,iN);
   A0=[Gas.*Ar; Ar];
   Aqw=Tt*A0;
   %Aqw1=(Mons+Mois*Lami)*Aqw;

   Aqw1=(Mons+Mois*Lami)*Aqw;
   Aend=Tbfield*Aqw1;
   ' qui va bene ! ', keyboard

   Acoz(:,1)=[Gas.*Ar; Ar];
   Acoz(:,2)=Aqw;
   Acoz(:,3)=Aend;
   Nz(1)=nv(1,1);
   Nz(2)=rqw;
   Nz(3)=nv(end,1);
   FF.Acoz=Acoz;
   duf=fmlstot(:,2);
   fidu=find(duf==0);
   duf(fidu)=1;
   spesdu=Litot.*duf;
   puqw=1:find(iauto(:,1)==2)-1;
   puend=1:find(iauto(:,1)==3)-2;
   lqw=sum(spesdu(puqw));
   lend=sum(spesdu(puend));

   FF.zi=[0 lqw lend];
%   ' sono qui Acoz '
%   keyboard
else



  clear dmem Acoz
  ifiez=1;
 if isav_Az>0
  if ifp~=-4
   ' dopo acc0 ',  pausak
  end

  ' chiamo re_new per il calcolo di E(z) ',   %keyboard
 if ifp==-10
  keyboard
 end
  ifr=2;
  re_new_critField

 else
  ifiez=1;
  freq=ze;
  [du,ifr]=min(abs(fre_camp-ze));
  ifr=ifr+1;
  kcav=kcav0*(1+freq);
%  isaMo=1;
  pvet00=pola;
  Lami=gg0/(c/(rr.*d*2));
   tic
   acc0
   toc
   Acozp=Acoz;
   Acoz(:,1)=Acoz_sav(:,1);
   Acoz(:,end)=Acoz_sav(:,end);
%   ' rimetto a posto'
%   keyboard
   FF.Acoz=Acoz;

 end
 if ifp~=-4
  'dopo campi z', keyboard
  keyboard
 end

 ifiez=0;
end
  ifast=ifasasa;

ifp=ifphold;