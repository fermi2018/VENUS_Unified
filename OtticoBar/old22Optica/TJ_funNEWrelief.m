function [dv_out,nv_out,av_out,iauto_out, shavet_out, fmlstot_out, anyf_out, xm_out, radii_out, ifield_out]=...
TJ_fun(ifp,PaTJ,L_i1,n_in1,aitot_in1,iauto1, shavet_in1, fmlstot1, anyf1, xm_in1, radii, ifield1)

itolgo=0;
if itolgo==1
[fr,fc]=find(-imag(n_in1)>.5);
%fr=[fr; fr(end)+1];
%fr=[fr(end)+1];
fnomet=[1:fr(1)-1 fr(end)+1:length(n_in1)];
fnometS=[1:fr(1)-2 fr(end):length(aitot_in1)];
anyf=anyf1(fnomet);
xm_in=xm_in1(fnomet);
ifield=ifield1(fnomet);
L_i=L_i1(fnomet);
n_in=n_in1(fnomet,:);
aitot_in=aitot_in1(fnometS,:);
iauto=iauto1(fnomet,:);
fmlstot=fmlstot1(fnomet,:);
shavet_in=shavet_in1(fnomet,:);
else
anyf=anyf1;
xm_in=xm_in1;
ifield=ifield1;
L_i=L_i1;
n_in=n_in1;
aitot_in=aitot_in1;
iauto=iauto1;
fmlstot=fmlstot1;
shavet_in=shavet_in1;
fr=[];

end


%'fermo', keyboard
n_i=n_in(:,1);

aitot=aitot_in;
xm=-10*ones(size(n_i));
%shavet=0*ones(size(n_i));
shavet=shavet_in(:,1);

puTJ=PaTJ.puTJ-length(fr);
ifi=0;


if ifi==1

figure, plot(cumsum(L_i),real(n_i))

hold on

D=sum(L_i(1:puTJ(1)-1));
plot(cumsum(L_i(puTJ))+D,real(n_i(puTJ)),'r','linewidth',2)
pausak
end


Detch=sum(L_i(puTJ));

n_lateral=PaTJ.nex;
fino1=find(n_i==1);
fino1=1;
pui=[fino1(end)+1:puTJ(1)-1];



Le=[L_i(fino1); sum(L_i(puTJ)); L_i(pui)];

nex=[n_i(fino1,:); n_lateral*ones(size(n_i(1,:))); n_i(pui,:)];
Ltop=flipud(Le);
Les=cumsum(Ltop);
ntopE=flipud(nex);

LeInter=[[0;Les(1:end-1)+1e-2] Les];

neInter=[ntopE ntopE];
Lint_ext=reshape(LeInter',2*length(Les),1);
nint_ext=reshape(neInter.',2*length(Les),1);


puiTJ=[1:puTJ(end)];
LtopTJ=flipud(L_i(puiTJ));
ntopTJ=flipud(n_i(puiTJ));
ragTJ=flipud(aitot(puiTJ));
Les=cumsum(LtopTJ);
ntopI=ntopTJ;

fmlst= fmlstot;

fst1=flipud(fmlst(puiTJ,1));
fst2=flipud(fmlst(puiTJ,2));

LeInter=[[0;Les(1:end-1)+1e-3] Les];
LeInterF=[[0;Les(1:end-1)] Les-1e-3];
%LeInter=[[0;Les(1:end-1)+1e-5] Les];
%LeInter=[[0;Les(1:end-1)] Les];
neInter=[ntopI ntopI];



f1Inter=[fst1 fst1];
f2Inter=[fst2 fst2];
Lint_int=reshape(LeInter',2*length(Les),1);
Lint_intF=reshape(LeInterF',2*length(Les),1);
nint_int=reshape(neInter.',2*length(Les),1);

f1_int=reshape(f1Inter.',2*length(Les),1);
f2_int=reshape(f2Inter.',2*length(Les),1);




if ifp==-10
figure, plot(Lint_int,nint_int), hold on, plot(Lint_int,f2_int,'r'), pausak
end

%'DTOTO', keyboard

Dtot=[0; cumsum(LtopTJ); cumsum(Ltop)];
[Ds,iso]=sort(Dtot);

dL=diff(Ds);
fiVeri=find(dL>0);


Dfinal=Ds(fiVeri);
Lfinal=diff(Dfinal);

n_intI=interp1(Lint_int,nint_int,Dfinal);
n_extI=interp1(Lint_ext,nint_ext,Dfinal);




f1_I=round(interp1(Lint_intF,f1_int,Dfinal,'nearest'));
f2_I=round(interp1(Lint_intF,f2_int,Dfinal,'nearest'));

f1_Id=(interp1(Lint_intF,f1_int,Dfinal,'next'));
f2_Id=(interp1(Lint_intF,f2_int,Dfinal,'next'));

%de=-1
%n_intIp=interp1(Lint_int,nint_int,1564+de)
%n_extIp=interp1(Lint_ext,nint_ext,1564+de)

%figure, plot(cumsum(L_final),real(n_int))
%hold on, plot(cumsum(L_final),real(n_ext),'r')

LeInter=[[Dfinal(1:end-1)] Dfinal(2:end)];
neInter=[n_intI(2:end) n_intI(2:end) ];
Lplot=reshape(LeInter',2*length(LeInter),1);
np_int=reshape(neInter.',2*length(LeInter),1);

neInter=[n_extI(2:end) n_extI(2:end) ];
np_ext=reshape(neInter.',2*length(LeInter),1);

fInter=[f1_I(1:end-1) f1_I(1:end-1) ];
f1p=reshape(fInter.',2*length(LeInter),1);

fInter=[f2_I(1:end-1) f2_I(1:end-1) ];
f2p=reshape(fInter.',2*length(LeInter),1);

if ifp==-10
figure, plot(Lplot,f1p,Lint_int,f1_int,'o')
hold on, plot(Lplot,f2p,Lint_int,f2_int,'o')

pausak


figure, plot(Lint_ext,real(nint_ext))
hold on, plot(Lint_int,real(nint_int),'r')
pausak
%figure, 
plot(Lplot,real(np_ext),'g.')
hold on, plot(Lplot,real(np_int),'m.')

pausak

plot(Lplot,real(np_ext),'g')
hold on, plot(Lplot,real(np_int),'m')

end

%'fine TJ', keyboard

%PaTJ=radii.TJ;
Lf=length(Lfinal);
puAdd=1+[1:Lf];


puB=[puTJ(end)+1:length(L_i)];

puBm2=[puTJ(end)+1:length(L_i)-2];

%[dv,nv,iauto, shavet, fmlstot, anyf, xm, radii]

% iauto
loc=sum(iauto(pui,:),2);
if length(find(loc>0))>0
 'iauto', keyboard
end
iauto_out=[iauto(1,:); zeros(Lf,size(iauto,2)); iauto(puB,:)];

% shavet
loc=sum(shavet(pui,:),2);
if length(find(loc>0))>0
% 'shavet', keyboard
end

%shaloc_add=zeros(length(iauto_out),1);
%shaloc_add(puAdd,1)=PaTJ.shape;
shavet_out=[shavet(1,:); zeros(Lf,size(loc,2)); shavet(puB,:)];
shavet_out(puAdd,1)=PaTJ.shape;


% xm
loc=sum(xm(pui,:),2);
if length(find(loc>0))>0
 'xm', keyboard
end
xm_o=[xm(1,:); -10*ones(Lf,size(loc,2)); xm(puB,:)];
xm_out=[xm_o xm_o(:,end)];

%xm_add=-10*ones(length(iauto_out),1);



% anyf
loc=sum(anyf(pui,:),2);
if length(find(loc>0))>0
% 'anyf', keyboard
end
anyf_out=[anyf(1,:); 0*ones(Lf,size(loc,2)); anyf(puB,:)];

% radii
sal=size(aitot,2);
if sal>1
% 'raggi da sistemare', keyboard
end

aloc=[aitot(1,:)*0; aitot; aitot(1,:)*0];
loc=sum(aloc(pui,:),2);
if length(find(loc>0))>0
% 'raggi', keyboard
end

% 'raggi', keyboard

R_tj=PaTJ.Ram;
av_du=[aloc(1,:); R_tj*ones(Lf,size(aloc,2)); aloc(puB,:)];
%av_du=[aloc(1,:); R_tj*ones(size(aloc)); aloc(puB,:)];
av_out=av_du(2:end-1,:);



%'avdy', keyboard

% fmlstot DA FIT

[pma,fima]=max(f2_I);
fimir=find(f2_I==pma);
if fimir(end)==length(f2_I)
f1_I(fimir(end-1:end))=0;
f2_I(fimir(end-1:end))=1;
f1_I(fimir-2)=length(fimir);
f2_I(fimir-2)=pma;
else
f1_I(fimir)=length(fimir);
f2_I(fimir)=pma;
end
fmFIT=[flipud(f1_I) flipud(f2_I)];
%fmlstot_out=[fmFIT; fmlstot(puBm2,:)];
fmlstot_out=[fmFIT; fmlstot(puB,:)];

% nv_tot DA FIT
nFIT=[flipud(n_intI) flipud(n_extI)];
if size(n_in,2)==1
 nRest=[n_i(puB,:) zeros(size(puB'))];
else
 nRest=n_in(puB,1:2);
end
nv_out=[[n_i(1) 0]; nFIT(1:end-1,:); nRest];

Lff=flipud(Lfinal);
dv_out=[L_i(1); Lff; L_i(puB)];

ifield_out=[ifield; zeros(length(dv_out)-length(ifield),1)];

avsave= av_out;
nvsave= nv_out;
shsave= shavet_out;
az=aitot_in(1,:)*0;
aitot_inM=[az; aitot_in; az];
n_inM=n_in;
sh_inM=shavet_in;
aitot_inM(puB,1:end-1)=aitot_inM(puB,2:end);
sh_inM(puB,1:end-1)=sh_inM(puB,2:end);
n_inM(puB,1:end-1)=n_inM(puB,2:end);

 lext=size(n_in,2);
% 'fine TJ par', keyboard 

Lfit=Lint_int-Detch;
%Lfit=Lint_int;
 if lext>1
%  av_out=av_out(:,1);
%  nv_out=nv_out(:,1:2);
%  shavet_out=shavet_out(:,1);
  av_out=avsave(:,1);
  nv_out=nvsave(:,1:2);
  shavet_out=shavet_out(:,1);

% 'PRIMA TJ for', keyboard 

  for kt=1:lext-1
%  'kt', keyboard
  n_i=n_inM(:,1+kt);
  aitot=aitot_inM(:,kt);
  shavet=sh_inM(:,kt);
  ntopTJ=flipud(n_i(puiTJ));
  ragTJ=flipud(aitot(puiTJ));
  shTJ=flipud(shavet(puiTJ));
  
  ragInter=[ragTJ ragTJ];
  neInter=[ntopTJ ntopTJ];
  shInter=[shTJ shTJ];
  
  ae_int=reshape(ragInter.',2*length(Les),1);
  ne_int=reshape(neInter.',2*length(Les),1);
  sh_int=reshape(shInter.',2*length(Les),1);


  n_extI=interp1(Lfit,ne_int,Dfinal);
  
  nFIT=flipud(n_extI);
  nRest=n_i(puB,:);
  nv_add=[n_i(1); nFIT(1:end-1,:); nRest];
  
  sh_extI=round(interp1(Lfit,sh_int,Dfinal));
  
  nFIT=flipud(sh_extI);
  nRest=shavet(puB,:);
  sh_add=[shavet(1); nFIT(1:end-1,:); nRest];  
%  'shavet', keyboard
  rag_I=interp1(Lfit,ae_int,Dfinal,'nearest');  
  nFIT=flipud(rag_I);
  nRest=aitot(puB(1:end-1),:);
  av_add=[nFIT(1:end-1,:); nRest];  
% 'rag', keyboard    
% 'rag', keyboard    
% 'rag', keyboard    
  av_out=[av_out av_add];
  nv_out=[nv_out nv_add];
  shavet_out=[shavet_out sh_add];
% 'rag', keyboard     
  
 end

fie=[1 3:length(dv_out)];
fie1=[2:length(av_out)];

%fie=[1:length(dv_out)];
%fie1=[1:length(av_out)];

%  'rag U', keyboard   

dv_out=dv_out(fie);

nv_out=nv_out(fie,:);
shavet_out=shavet_out(fie,:);
iauto_out=iauto_out(fie,:);
fmlstot_out=fmlstot_out(fie,:);
xm_out=xm_out(fie,:);
ifield_out=ifield_out(fie,:);
av_out=av_out(fie1,:);
%dv_out,nv_out,av_out,iauto_out, shavet_out, fmlstot_out, anyf_out, xm_out, radii_out, ifield_out]=..

 

 av_dum=[av_out(1,:); av_out; av_out(end,:)];
 av_dum1=av_dum*0;
 ra_old=av_out;
 nv_old=nv_out;
 sv_old=shavet_out;
 kl=0
 for kl1=1:length(shavet_out)
  kl=kl+1;
  ra=av_dum(kl1,:);
%  pausak
  rn=nv_out(kl1,2:end);
  sn=shavet_out(kl1,:);

  fima=find(ra>0);
  raa=ra(fima);
  if length(fima)>1
  if raa(2)==raa(1);
   fima=fima(1:end-1);
  end  
  end
  
  rs=0;
  if length(fima)>0
   [rs,iso]=sort(ra(fima));
%   if kl>150
%   kl
%  pausak
%  end

   fiso=fima(iso);
   av_dum1(kl,1:length(fiso))=rs;   
   shavet_out(kl,:)=0;   
   nv_out(kl,2:end)=0;   
   shavet_out(kl,[1:length(fiso)])=sn(fiso);   
   nv_out(kl,1+[1:length(fiso)])=rn(fiso);
  end
 end
 
  av_out=av_dum1(2:end-1,:);
 end
% 'fine', keyboard
if exist('av_dum1')
av_du=av_dum1;
end
 
 radii_out.a=av_du;
 radii_out.b=av_du*0;
 radii_out.c=av_du*0;
 for k=1:length(fie)
  radii_out.array{k,1}{1}=0;
 end
 
%if length(find(aitot_in(1,:)>0))>1
 
%nv_out([2 4],2)=nv_out([4 2],2);
%nv_out(3,3)=nv_out(3,2);

%end
%'chiuso NEW', keyboard






if ifp==-10

Ls=cumsum(L_i);
du=[[0; Ls(1:end-1)] Ls];
nu=[real(n_in(:,1)) real(n_in(:,1))];
fu=[fmlstot(:,2) fmlstot(:,2)];
 duri=reshape(du',1,prod(size(du)));
 nuri=reshape(nu',1,prod(size(du)));
 furi=reshape(fu',1,prod(size(du)));

Ls=cumsum(dv_out);
du=[[0; Ls(1:end-1)] Ls];
nu=[real(nv_out(:,1)) real(nv_out(:,1))];
fu=[fmlstot_out(:,2) fmlstot_out(:,2)];
 dur=reshape(du',1,prod(size(du)));
 nur=reshape(nu',1,prod(size(du))); 
 fur=reshape(fu',1,prod(size(du))); 
 
figure, plot(duri,nuri,'o-',dur,nur,'p--')
figure, plot(duri,furi,'o-',dur,fur,'p--')

'chiuso', keyboard
end

