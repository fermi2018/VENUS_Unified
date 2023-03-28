% function [Gacrit,Trcrit]=gam_critU(Tstor,Ga1,Mcri,ficri);
function [Gacritu,Trcritu,Trc,Grc,puesp]=gam_critScatt(Ga1,Mcri,ficrin,fmlsi,dovesono,ipol);

ficri=abs(ficrin);
iPv=sign(ficrin);

lc=length(ficri);

FM=fmlsi(ficri,2)';
firep=find(FM>1);

dfirep=diff(firep);
fdif=find(dfirep>1)+1;

if length(fdif)==0
 fdif=length(firep)+1;
end

% ' sono dentro  gam_crit scatt', keyboard
 
 puesp=ficri';
if length(firep)>0
 firepid=1:fdif(1)-1;
 firepi=firep(firepid);
 fire=puesp(firepi); 
 re=FM(fire);
 puesp=[puesp(1:firepi(1)-1) repmat(fire,1,re(1)) puesp(firepi(end)+1:end)];
 FM=[FM(1:firepi(1)-1) repmat(ones(size(fire)),1,re(1)) FM(firepi(end)+1:end)]; 
 firep=find(FM>1);
% 'qui', keyboard
 while length(firep)>0
  dfirep=diff(firep);
  fdif=find(dfirep>1)+1;
  if length(fdif)>0
   firepid=1:fdif(1)-1;
  else
    firepid=1:length(firep);
  end
  firepi=firep(firepid);
  re=FM(firepi);
  fire=puesp(firepi);
  puesp=[puesp(1:firepi(1)-1) repmat(fire,1,re(1)) puesp(firepi(end)+1:end)];
  FM=[FM(1:firepi(1)-1) repmat(ones(size(fire)),1,re(1)) FM(firepi(end)+1:end)];
   firep=find(FM>1);
 end 
end
 
% ' sono dentro  gam_crit scatt', keyboard 
 Ga1_sa=Ga1;
% TscaU
 
 lcri=length(puesp);
 s=size(Ga1);
 Mn=1;
 Tr1=1;
%'prima di crit', keyboard

kcrit=puesp;
if dovesono==2
 pacri=-1;
 kcrit=puesp(end:-1:1);
end
%' crit', keyboard
% for k=1:lcri
ick=0;
 for k=kcrit
  ick=ick+1;
%  iP=iPv(ick);
  iP=iPv(k);
  kt=k;
 
  if imag(iP)==0
   Mn=Mcri{kt,ipol};
    Ga1_s=Ga1;
   [Ga1,Tr1d]=GaScattHCG(Mn,Ga1,kt,iP);
   Tr1=Tr1*Tr1d;
  end
   Trd(:,:,ick)=Tr1d;
   Grc(:,:,ick)=Ga1;
 end
 Tpre=1;

% for k=lcri:-1:1
%  Tri=Trd(:,:,k);
%  Trc(:,:,lcri-k+1)=Tri;
% end 
Trc=Trd;
Trcritu=Tr1; 
Gacritu=Ga1; 
% 'fine gacritu', keyboard 
return 

% 'fine for', keyboard
  s=size(Mn);
 if s(1)>1
  s=size(Mn);
  l1=s(1)/2;
  l2=s(1)/2+1;
  l3=s(1);
  M11=Mn(1:l1,1:l1);
  M12=Mn(1:l1,l2:l3);
  M21=Mn(l2:l3,1:l1);
  M22=Mn(l2:l3,l2:l3);
%iT=inv(M22);
%s11=-iT*M21;
%s12=iT;
%s21=M11-M12*iT*M21;
%s22=M12*iT;

   Gadu=Ga1;
   
   A=M21*Gadu+M22;
   B=M11*Gadu+M12;
   inva=inv(A);
   Gacritu=B*inva;
%   Trcritu=inv(M22);
   Trcritu=Tr1*inva;
else
 Gacritu=Ga1;
 Trcritu=Tr1;
end

%map(log10(Gacritu))
%' contro ', keyboard   
   
   return
   
   Id=eye(length(Ga1));

%   Gacritu=s22+s21*inv(Id-Gadu*s11)*Gadu*s12;
   
   
   
   Gacrit=s22u+s21u * inv(Id-Gacritu*s11u) *Gacritu*s12u;

   iM=inv(M11);   
   Tcritu=M22-M21*iM*M12;   
   Gcritu=M21*iM;   
   Trcrit=Tcritu * inv(Id-s22u*Gcritu) *s21u;
   
   
'gacrit', 
keyboard