function [gamod,la,zet,nz,Ez,Hz,uL]=CMM_1D(lam0,L_i,n_i,iat,icav,fiQW,ifp,Dlam,rr,kt,ibast,par_grat,i_campi,a_i);

ifpxont=-4;

if length(Dlam)>1
 dlam=Dlam(1);
 Ndlam=Dlam(2);
else
 dlam=Dlam;
 Ndlam=60;
end

i_campi=1;
if ~exist('rr')
 rr=3;
end 
if ~exist('kt')
 kt=0;
end 
if ~exist('ibast')
 ibast=[];
 par_grat=0;
else
 if ibast>0
  if ~exist('par_grat')
   'missing grating parameters ', keyboard
  end
 end
end 

PO=[634    51   498   711];

lav=lam0+linspace(-1,1,Ndlam)*dlam*lam0;
%lav=lam0+linspace(0,2,20)*dlam*lam0;
%lav=lam0-.01;

%'in CMM', keyboard
for k=1:length(lav)
lai=lav(k);
%(L_in,n_i,iat,icav,la0,rr,kt,ibast,par_grat,g0);
%if k==1
[ei,g0]=Flam_gamLu(L_i,n_i,iat,icav,lai,rr,kt,ibast,par_grat);
%'cont 1', keyboard
%else
%[ei,g0]=Flam_gamLu(L_i,n_i,iat,icav,lai,rr,kt,ibast,par_grat,real(ei));
%'passo'
%end
gov(k)=g0;
eiv(k)=ei;
end
alf0=imag(eiv);
gai0=real(eiv);

der=diff(alf0);
der=[der der(end)];

ne=length(find(der<=0));
po=length(find(der>0));
if po>=1 | ne>=1
 if po>ne
  puf=(find(der>0)); 
 else
  puf=(find(der<=0)); 
 end
 isw=1
else
 isw=0
end
PU=[];
  adu=alf0;
  se=adu(1:end-1).*adu(2:end);
  puf0=(find(der>0)); 
  se(puf0)=median(se);
  fib=find(se<0);
  if ifpxont==-10
   sed=se/median(se);
   figure, plot(sed),
   hold on
   plot(fib,sed(fib),'ro')
   axis([1 length(se) -.1 .1])
   ' cont 1D', keyboard 
  end

   
  ipu=0;
  if  length(fib)>0
  
   for kk=1:length(fib)
    pud=fib(kk)+[-2:2];
    fi=find(pud>0&pud<=length(alf0));
    pud=pud(fi);
    ipu=ipu+1;
    PU{ipu}=pud; 
   end
  else 
   PU{1}=1:length(lav);
  end 

if length(PU)==0
  NOInt1D
end
%    ' PU', keyboard 

for kint=1:length(PU)
 pui=PU{kint};
 [du,imid] =min(abs(alf0(pui)));
 imi=pui(imid);
  if ifpxont==-10
  figure,  
  subplot(211)
  plot(lav(pui),alf0(pui),lav(imi),0,'wo'), grid
  subplot(212)
  plot(lav(pui),gai0(pui),lav(imi),gai0(imi),'wo'), grid, pausak
  end 
GSTIM=gai0(imi);
 dlau=diff(lav(1:2));
 lav1=lav(imi)+linspace(-1,1,10)*dlau;
 clear gov eiv
 for k=1:length(lav1)
  lai=lav1(k);
  [ei,g0]=Flam_gamLu(L_i,n_i,iat,icav,lai,rr,kt,ibast,par_grat,GSTIM);
  gov(k)=g0;
  eiv(k)=ei;
 end
 


%'qui CMMkk', keyboard

g00=g0;
alf=imag(eiv);
gv=real(eiv);
if length(lav1>7)
 [du,im]=min(abs(alf));
 pud=[-2:2]+im;
 puf=find(pud>0 & pud<=length(lav1));
 pu=pud(puf);
 lavf=lav1(pu);
 gvf=gv(pu);
 alvf=alf(pu);
 govf=gov(pu); 
else
 lavf=lav1;
 gvf=gv;
 alvf=alf;
 govf=gov;
end



coa=polyfit(lavf,alvf,1);
LAM=roots(coa);

if ifpxont==-10
figure, plot(lavf,alvf,LAM,0,'wo'), pausak
end
cog=polyfit(lavf,gvf,2);
GTH=polyval(cog,LAM);
G0=GTH;
cog0=polyfit(lavf,govf,2);
Gga=polyval(cog0,LAM);

if ifpxont==-10
h=figure, 
set(h,'pos',PO)
subplot(211)
plot(lav1,eiv,LAM,GTH,'wo'), title('real part'), grid
subplot(212)
plot(lav1,imag(eiv),LAM,0,'wo'), title('imag part'), grid
pausak
end
 GV(kint)=G0;
 LV(kint)=LAM;

end   %fine kint


% inserisco ulteriore zoom
[GSTIM,ig]=min(GV);
LAV=LV(ig);

for iter=1:2
 dlau=diff(lav1(1:2)/4);
 lav1=LAM+linspace(-1,1,5)*dlau;
 clear gov eiv
 for k=1:length(lav1)
  lai=lav1(k);
  [ei,g0]=Flam_gamLu(L_i,n_i,iat,icav,lai,rr,kt,ibast,par_grat,GSTIM);
  gov(k)=g0;
  eiv(k)=ei;
 end
 
alvf=imag(eiv);
gvf=real(eiv);
lavf=lav1;
coa=polyfit(lavf,alvf,1);
LAM=roots(coa);

if ifpxont==-10
figure, plot(lavf,alvf,'.-',LAM,0,'wo'), pausak
end
cog=polyfit(lavf,gvf,1);
GTH=polyval(cog,LAM);
G0=GTH;

if ifpxont==-10
h=figure, 
set(h,'pos',PO)
subplot(211)
plot(lav1,eiv,'.-',LAM,GTH,'wo'), title('real part'), grid
subplot(212)
plot(lav1,imag(eiv),'.-',LAM,0,'wo'), title('imag part'), grid
pausak
end
 GSTIM=G0;
 LAV=LAM;
 Gmon(iter)=G0;
 Lmon(iter)=LAM;
end %iter

% fine ulteriore zoom


%' FINE ITER', keyboard


ga=GSTIM;
la=LAV;


if i_campi==1

if ifp==-10
h=figure, 
set(h,'pos',PO)
subplot(211)
plot(lav1,eiv,'.-',LAM,GTH,'wo'), title('real part'), grid
subplot(212)
plot(lav1,imag(eiv),'.-',LAM,0,'wo'), title('imag part'), grid
' campi,' , pausak
end

%' campi,' , keyboard
Nx=50;
[fiez,nz,zet,Gaqw,NQW_ef,az,uL]=Flam_field(L_i,n_i,iat,fiQW,la,rr,kt,ibast,par_grat,ga,Nx);

%[fiez,nz,zet,Gaqw,NQW_ef,az]=Flam_field(L_i,n_i1,iat,fiQW,la,rr,kt,ibast,par_grat,ga,Nx,a_i);


Ez=sum(fiez);
Ez=Ez/max(Ez);
Hz=diff(fiez);
Hz=Hz/max(Hz);
lab='';
if iscell(par_grat)==1
if par_grat.itetm==1
 lab='TE';
else
 lab='TM';
end
end
if ifp==-10
figure, plot(zet,abs(Ez).^2*3, zet, nz,'r'), 
 title(['CMM: ',lab,'  lambda_{res} = ',num2str(la),' Gth pa= ',num2str(ga/NQW_ef)]), 
'fine campi', pausak
PF=log10(abs(Ez).^2)+3;
figure, plot(zet,PF, zet, nz,'r'), 
 title(['CMM: ',lab,'  lambda_{res} = ',num2str(la),' Gth pa= ',num2str(ga/NQW_ef)]), 
 a=axis;
 a(3)=PF(end)-.1;
 axis(a)
'fine campi', pausak
end

else
 zet=0;
 nz=0;
 Ez=0;
 Hz=0;
end
gamod=ga/NQW_ef;
%'fine CMM', keyboard