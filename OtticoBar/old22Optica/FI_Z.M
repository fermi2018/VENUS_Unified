if ~exist('irepl_ind')
 irepl_ind=0;
end
%save fi
%save fi0
%' sopo dave', keyboard
%' sopo dave', keyboard
for kstra=1:length(nv0)
% idu=find(xm(kstra,:)~=0);
 idu=1;
 [du,ip]=max(xm(kstra,idu));
 n1(kstra,1)=nv0(kstra,idu(ip));
 c1(kstra,1)=xm(kstra,idu(ip));
end
if exist('shavet0')
fi_re=find(shavet0(:,1)==6);
else
fi_re=find(shavet(:,1)==6);
end
if length(fi_re)>0
 if exist('nret')
  n1(fi_re)=nret;
 else 
  n1_re=nv0(fi_re,:);
  fival=find(n1_re~=0);
  n1(fi_re)=sqrt(mean(n1_re(fival).^2));
 end
end
fi_ce=find(shavet==4);
race=imag(radii.c(fi_ce));
raca=radii.a(fi_ce);
p=(race+raca).*(race-raca);
fiout=find(p>0);
if length(fiout)>0
 fi_out=fi_ce(fiout);
 n1(fi_out)=nv0(fi_out,2);
end
%' cont cce nv', keyboard
if irepl_ind==1
%' irepl=1', keyboard
 fi_re=find(shavet==4);
 if length(fi_re)>0
  nrepla=nv0(fi_re(1)-1,1);
  n1(fi_re)=nrepla;
 end
end 

xd=[zeros(size(dv)) dv];
fl=abs(fst(:,2));
flu=find(iauto(:,1)==3);
fim=find(fl>1);
%fim=[fim; flu];
fw=length(fim);

fime=find(real(n1)<0);
n1(fime)=-n1(fime);

if length(fim)>0
miup=fst(fim(1),2);
midw=fst(fim(end),2);
else
miup=0;
midw=0;
end
poup=miup-4;
podw=4;

dto=[];
nto=[];
ato=[];
cto=[];
puo=[];

piv=1;
ico=1;
leqw=0;
fsto=[];
while fw>0 
%while fw>0 | leqw==0
 pd=fim(ico);
 nl=fst(pd,1);
 nlay=fix(abs(fst(pd,2)));
% [pd piv nl nlay], pausak
 if piv-pd~=0
  pup0=piv:pd-1;
  ddu=dv(pup0);
  imz=find(ddu>0);
  pup0=pup0(imz);
  if length(find(pup0==pua))==1
   fi0=find(pup0==pua);
   leqw=length(fsto)+fi0;
  end
  ay=anyf(pup0);
  nto=[nto; n1(pup0)];
  cto=[cto; c1(pup0)];
  fsto=[fsto; fst(pup0,2)];
  dto=[dto; dv(pup0)];
  ato=[ato; dv(pup0).*ay];
  puo=[puo; pup0'.*ay];
 end
%'ICI 0 bir', keyboard
 pupe=pd:pd+nl-1;
 ay=anyf(pupe);
 for inp=1:nlay
  nto=[nto; n1(pupe)];
  cto=[cto; c1(pupe)];
  fsto=[fsto; fst(pupe,2)];
  dto=[dto; dv(pupe)];
  if nlay==miup | nlay==midw
   if (inp==poup & nlay==miup) | (inp==podw & nlay==midw)
    ato=[ato; dv(pupe).*ay];
    puo=[puo; pupe'.*ay];
   else
    ato=[ato; -dv(pupe).*ay];
    puo=[puo; pupe'*0];
   end
  else
   if inp==fix(nlay/2)
    ato=[ato; dv(pupe).*ay];
    puo=[puo; pupe'.*ay];
   else
    ato=[ato; -dv(pupe).*ay];
    puo=[puo; pupe'*0];
   end
  end
 end
 piv=pupe(length(pupe))+1;
 fw=fw-nl;
 ico=ico+nl;
end

%'ICI 1 bir', keyboard

if piv<length(dv)+1
 pup0=piv:length(dv);
 nto=[nto; n1(pup0)];
 cto=[cto; c1(pup0)];
 fsto=[fsto; fst(pup0,2)];
 dto=[dto; dv(pup0)];
 ay=anyf(pup0);
 ato=[ato; dv(pup0).*ay];
 puo=[puo; pup0'.*ay];
end

puac=puo;
puae=find(ato>0);
%puae=find(ato>=0);
%' puae ', keyboard


nd=[nto nto];
ccd=[cto cto];
%pudto=
xd=[zeros(size(dto)) dto];
xt=cumsum(reshape(xd',2*length(dto),1));
nt=real(reshape(nd',2*length(dto),1));
ct=(reshape(ccd',2*length(dto),1));
L_i=dto;
n_i=nd(:,1);
%'fi_z', keyboard
lmaxim=1e6;
%lmaxim=1e3;
fitro=find(dto>lmaxim);
if length(fitro)>0
 dto(fitro)=lmaxim;
end
%'qui dto', keyboard
global Lmax
if length(Lmax)==0
 Lmax=15000;
end
fidt=find(dto>Lmax);
dto(fidt)=Lmax;
dJ=cumsum([0; dto]);
fJ=([0; fsto]);

nto(leqw)=real(nto(leqw));
nJ=([nto(1); nto]);
xJ=([cto(1); cto]);
lJ=length(dJ);
hn=1;
hz=[dJ(1):hn:dJ(lJ)];
uFunc=hz*0;
uF0=hz*0;
perm=[];
xfit=[];
for il=1:lJ-1
 fiz=find(hz>=dJ(il) & hz<dJ(il+1));
 perm=[perm nJ(il+1)*ones(size(fiz))];
 xfit=[xfit xJ(il+1)*ones(size(fiz))];
 if fJ(il+1)==-1
  uFunc([fiz(1)-1 fiz])=1;
 end
 if fJ(il+1)==-1 & il==leqw
  uF0([fiz(1)-1 fiz])=1;
 end
end
xfit(length(hz))=xfit(length(hz)-1);
perm(length(hz))=perm(length(hz)-1);
uFunc(length(hz))=0;
uF0(length(hz))=0;
relPerm=conj(perm.^2);
%relPerm=real(relPerm);

if ifp<=-10
rp=real(relPerm);
rpc=rp;
figure, plot(hz/1000,rp)
hzc=hz/1000;
save sa hzc  rpc
%load sa,  hold on, plot(hzc,rpc,'r')
%keyboard
%pausak
end
%'qui prima di mio dopo', keyboard
%'qui prima di mio dopo', keyboard
if ~exist('lambda0')
 lambda0=lambda/1e6
end 
%lambda0=866*1e-9
%'qui prima di mio dopo', keyboard
% dorigin=0;
% x=(hz'-dorigin)*1e-3;
% rperm=real(sqrt(relPerm));
% load pri
% rpermp=real(sqrt(relPermp));
%  figure, plot(rperm,'r'), hold on, plot(rpermp,'g')
% figure, plot(
% lambda0p=lambda0; uFuncp=uFunc; uF0p=uF0; relPermp=relPerm;
% save pri lambda0p uFuncp uF0p relPermp
if nmir.a==1
 uF0=uFunc;
end
fiat=find(iauto(:,1)==2);
datt=dv(fiat)*1e-9;
 [Ksi,lambdas,Fi,uLong,uLong0,Fa]=eiglmiu(lambda0,uFunc,uF0,relPerm,hn,datt);
 uL(1)=uLong0;
 uL(2)=uLong;
 fatqw=uLong/(uLong0*nmir.a);
 lambda=lambdas*1e6;
 gpla=2e4*pi*rr/lambda*imag(Ksi)/uLong0/fatqw; 
 gpla=2e4*pi*rr/lambda*imag(Ksi)/uLong; 
 manv=max(uF0.*Fi.');
 confzv=uL;
 dorigin=0;
 x=(hz'-dorigin)*1e-3;
 rperm=real(sqrt(relPerm));
 faca=3;
 if ifp<=-10
  figure, plot(x,rperm,'r',x,faca*Fi/manv,'w')
%    title([' lambda_{res} = ',num2str(lambda),'  Ksi = ',num2str(Ksi)]), pausak
    title([' lambda_{res} = ',num2str(lambda),'  Gth = ',num2str(gpla)]), pausak
  pausak
 man=max(uF0.*Fi.')/3.5; 
 figure, plot(hz,real(perm),'w',hz,faca*Fi/man,'r',hz,-imag(perm)*6000,'c'),
 axe=axis;
 axe(4)=5;
 axis(axe)
  title([' lambda_{res} = ',num2str(lambda)]), pausak  
 end
