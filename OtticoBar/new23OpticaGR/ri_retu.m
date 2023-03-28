global ispeed  k_ret nucomp ikr

mbdu=(nubesi:pasnu:nubesu)+.1;
if iLP==0 
 mbdu=[mbdu -mbdu];
end
if ipolar==0
 mbdu=[mbdu mbdu];
end
%' ri_retu ', keyboard
if length(nucomp)==0
 nucomp=1;
end
if iLP==0
 if mm==0 & nucomp==1
  nucomp=2;
 end
%else
% if mm==1 & nucomp==1
%  nucomp=2;
% end
end
%nucomp=0
%keyboard
if length(fir)>0
 nmv=[1:nucomp];
 if iLP==0
  if igr_app~=2 
   ianmat=1;
  end
 else
  ianmat=0;
 end
else
 nmv=[1:numodiacc+1];
end



if ifp~=-4
'% riduco la base'
end
%keyboard

 if length(irid_bas)==0
  irid_bas=1;
 end
% if ispeed==0 | ipolar==0
%  irid_bas=0;
% end

if irid_bas==1

 lk=length(KK);
 Dm=ones(lk*(numodiacc+1),1);
 pu=1:lk;
 pre=.3;
 Pus=[];
 Moe=[];
 ldap=[0];

 Dmv=Dm;
 npam=1;
% if iLP==0
%  Dm=(Dm(1:end/2)+Dm([1:end/2]+end/2))/2;
% end
%global tolk kretvero farm


iNM=0;
if exist('tolk')
 tolku=tolk(1);
 if tolku>0
  iNM=1;
 end
 if length(tolk)>1
  fashif=tolk(2);
 else
  fashif=0;
 end
end

%'iNM', keyboard
if iNM==0
 for nm=1:npam:numodi
%  if nm==1
  if length(find(nmv==nm))==1
   k_red=0;
  else
   k_red=k_ret;
  end
  pui=pu+(nm-1)*lk;
  Du=abs(Dm(pui));
  Du=Du/max(Du);
  der=diff([0; Du]);
  if length(k_red)==0 | k_red==0
   fiu=find(Du-pre>0 & der>0);
   pua=fiu(1):lk
%   keyboard
  else
%   fiu=find(KK>k_red & Du-pre>0);
   fiu=find(KK>k_red & Du-pre>0);
   if length(fiu)>5
    pua=fiu;
   else
    pua=lk-4:lk;
   end
  end

  if npam==1
   Pus=[Pus pui(pua)];
   ldap=[ldap ldap(end)+length(pua)];
   Moe=[Moe ones(1,length(pua))*mbdu(nm)];   
%   'Moe', keyboard
  else
   Pus=[Pus pui(pua) pui(pua)+lk];
   ldad=ldap(end)+length(pua);
   ldap=[ldap ldad ldad+length(pua)];
  end
 end
else %iNM


%krv=[.2 .4]
%tolku=.05;
if ikr~=0
  krv=kretvero*[1:farm]+fashif;
  kin=krv(1)-tolku;
  kfi=krv(1)+tolku;
%  'estremi k', keyboard
  puk_ret=find(KK>kin & KK<=kfi)';
else
  krv=0;
end
%  ' puk_ret', keyboard
for nm=1:npam:numodi
 for ik=1:length(krv)
  kin=krv(ik)-tolku;
  kfi=krv(ik)+tolku;
%  if nm==1
  if length(find(nmv==nm))==1
   k_red=0;
  else
   k_red=k_ret;
  end
  pui=pu+(nm-1)*lk;
  Du=abs(Dm(pui));
  Du=Du/max(Du);
  der=diff([0; Du]);
%  fius=find(KK>kin & Du-pre>0 & KK<=kfi);
  if length(k_red)==0 | k_red==0
   fiu=find(Du-pre>0 & der>0);
   pua=fiu(1):lk;
  else
%   fiu=find(KK>k_red & Du-pre>0);
   fiu=find(KK>kin & Du-pre>0 & KK<=kfi);
%   keyboard
   if length(fiu)>5
    pua=fiu;
   else
    pua=lk-4:lk;
   end
  end
% if nm~=1 | ik~=1
 if length(find(nmv==nm))==0 | ik~=1
  if npam==1
    Pus=[Pus pui(pua)];
    ldap=[ldap ldap(end)+length(pua)];
  else
   Pus=[Pus pui(pua) pui(pua)+lk];
   ldad=ldap(end)+length(pua);
   ldap=[ldap ldad ldad+length(pua)];
  end
 end

 end  %krv

end
%' enntro proar', keyboard
%' enntro proar', keyboard
%' enntro proar', keyboard
%' enntro proar', keyboard
proar

end %iNM

%' qui ri_retu 1', keyboard
if ifp~=-4
%' qui ri_retu 1', keyboard
pri=figure; plot(KKv(Pus),'.')
if exist('Pusf')
hold on, plot(KKv(Pusf),'ro')
end
%keyboard
end

 lknumodi=lk*numodi;
%' lknumodi=lk*numodi; ' , keyboard
 lknumodi1=ldap(end);

 if iLP==0
  Pus0=Pus;
  Pusas0=Pus;
  Pus=[Pus Pus+lknumodi];
  if exist('Pusf')
   Pusf=[Pusf Pusf+lknumodi];
  end
  ldapu=[ldap ldap(2:end)+lknumodi1];
  ldap=[ldap ldap+lknumodi1];
 else
  Pus0=Pus;
  ldapu=ldap;
 end
 if iredmat==1
  Pusas=Pus;
  if exist('Pusf')
   Pusasf=Pusf;
  else
   Pusasf=Pusas;
  end
  Pus=1:length(Pus);
  Pusa=Pus;
  Pusasav=Pus;
 else
  Pusa=1:si2(1);
 end

else

 lKA=si2(1);
 Pus=1:lKA;
 numoa=0:numodi-1;
% ld=[numoa numodi+numoa]*nk1max;
% ldap=[ld ld+lKA];
 ld=[numoa numodi+numoa]*nk1max;
 ldap=[ld ld(end)+nk1max];

end

lPus=length(Pus);
if iredmat==0
 si21=si2(1);
else
 si21=lPus;
 si2=[lPus lPus];
end

%Pusa=1:si21;
Pust=[Pus Pus+si21];

% set global IdeOo pMu0u pMu0u1 lKA nk1max nures sim0

lKAn=length(Pus);
 global ive
 ive=0;
 if iredmat==1
  ive=0;
 end
 if iany==0
  ive=1;
 end
 ive=0;
if iredmat==0
 if iLP==0
  lKA=nk1max*numodi*2;
  nures=numodi*4;
  sim0=nk1max*2;
 else
  lKA=nk1max*numodi;
  nures=numodi*2;
  sim0=nk1max*2;
 end
  numoa=0:numodi-1;
  ld=[numoa numodi+numoa]*nk1max;
  ldu=[ld ld+lKA];
  pMu0d=[];
  for kdu=1:nures
   pMu0d=[pMu0d ldu+(kdu-1)*2*lKA*nk1max];
  end
 if ive==1
  pMu0u=pMu0d;
 else
  for ki=1:nk1max
   pMu0u{ki}=pMu0d+ki+(ki-1)*2*lKA;
  end
 end
  pMu0u1=[0 2*lKA^2; lKA lKA+2*lKA^2];

  ld1=1:nk1max;
  if iLP==0
   ldud0=[ld1];
   ldud=[ldud0 ldud0+lKA];
  else
   ldud0=[ld1 ld1+lKA];
   ldud=ldud0;
  end
  pMc=[];
  for kdu=ldud
   pMc=[pMc ldud+(kdu-1)*lKA*2];
  end
  ldus=ldud;

  ldud=ld1;
  pMei=[];
  for kdu=ldud
   pMei=[pMei ldud+(kdu-1)*lKA];
  end

else  %riduco matrici

 if iLP==0
  lKA=lKAn;
  nures=numodi*4;
  sim0=nk1max*2;
 else
  lKA=lKAn;
  nures=numodi*2;
  sim0=nk1max*2;
 end
  nk1maxv=diff(ldapu);
%  if iLP==0
   nk1maxv=[nk1maxv diff(ldapu)];
%  end
  numoa=0:numodi-1;
  ld=ldapu;
  ldu=[ld ld(2:end-1)+lKA];

%  pMu0u=[];
%  for kdu=1:nures
%   pMu0u=[pMu0u ldu+(kdu-1)*2*lKA*nk1maxv(kdu)];
%  end
  if iLP==0
   nktot=4*nk1max*numodi;
  else
   nktot=2*nk1max*numodi;
  end
  VPU=1:nktot^2;
  MPU=reshape(VPU,nktot,nktot);
%  Pusta=[Pusa Pusa+Pusa(end)];
  Pusta=[Pusas Pusas+nktot/2];
  MPUR=MPU(Pusta,Pusta);

  numoa=0:numodi-1;
  ld=[numoa numodi+numoa]*nk1max;
  lKA1=nktot/2;
  if iLP==0
   ldu=[ld ld+lKA1];
  else
   ldu=ld;
  end
    pMu0d=[];
    for kdu=1:nures
     pMu0d=[pMu0d ldu+(kdu-1)*2*lKA1*nk1max];
    end

  if ive==1
   pMu0u=pMu0d;
  else

ivemepu=0;
    for ki=1:nk1max
     pMu0du=pMu0d+ki+(ki-1)*2*lKA1;
%     tic
     MPUdu=MPU;
     MPUdu(pMu0du)=-MPU(pMu0du);
     MPURdu=MPUdu(Pusta,Pusta);
     fib=find(MPURdu<0);
%     toc
     if ivemepu==1
     tic
     Vdu=MPU(pMu0du);
     lV=length(Vdu);
     ic=0;
     clear lVd
     for kj=1:lV
      fi=find(Vdu(kj)==MPUR);
      if length(fi)==1
       ic=ic+1;
       lVd(ic)=fi;
      end
     end
     toc
     length(find(lVd-fib'~=0))
     [ki, length(lVd)], pausak
     else
      lVd=fib;
     end
     pMu0u{ki}=lVd;
 % ' ki = ', ki
 %  keyboard
    end
   end

% ' qui a ', keyboard
%' pMu01 ', keyboard


  pMu0u1=[0 2*lKA^2; lKA lKA+2*lKA^2];

%  for
    imo=1;
    ld1p=ldapu(imo+[0 1])+[1 0];
    ld1=ld1p(1):ld1p(2);
    if iLP==0
     ldud0=[ld1];
     ldud=[ldud0 ldud0+lKA];
    else
     ldud0=[ld1 ld1+lKA];
     ldud=ldud0;
    end

    pMc=[];
    for kdu=ldud
     pMc=[pMc ldud+(kdu-1)*lKA*2];
    end

    ldus=ldud;
    ldud=ld1;
    pMei=[];
    for kdu=ldud
     pMei=[pMei ldud+(kdu-1)*lKA];
    end

end

 IdeOo=diag(ones(2*lKA,1));
 IdeOon=diag(ones(2*lKAn,1));

pure0=1:npuntik;
lad=ldapu(end)/2;
Pured=[];
for kpun=1:4
Pured=[Pured pure0+lad*(kpun-1)];
 if kpun==2
  Pured1=Pured;
 end
end



%'fine riretu', keyboard
if ifp~=-4
 disp(' RI_retu'),
 [nures/2*length(KK) length(Pus)]
 if ifp~=-4, pausak, end
%end
%keyboard
%disp(' IdeOo')
%' riretu', keyboard
% if ifp~=-4
 figure(pri), hold on, plot(KKv(Pus0),'ro')
 pausak
%' qui ri_retu 2 VoRtex', keyboard
 end

% global IdeOo IdeOon pMu0u pMu0u1 lKA nk1max nures pMc sim0 ldap

