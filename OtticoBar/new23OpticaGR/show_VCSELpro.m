%'in show_VCSEL', keyboard


n_type='real';
coord_num='';

aitotM=aitot;
  if isfield(Ps,'imesa')
    if Ps.imesa==0
    RmeMa=Ps.imesaRag;
    fiRa=find(aitotM>RmeMa);
%    aitotM(fiRa)=0;
   end
 end

Amax=max(max(aitotM))*1.5;

if fish>0
 Ngra=fish-1;
 if isfield(par_grat,'th')==1
  Lg=par_grat.th'/1000;
 else
  Lg=Litot(Ngra);
 end
 perio=par_grat.per;
 dg1=perio*par_grat.DC;
 dg2=perio*(1-par_grat.DC);
 ng1=par_grat.n1.';
 ng2=par_grat.n2.';
 O=ones(size(Lg));
 aG0=[dg1*O (dg1+dg2)*O];
 nG0=[ng1 ng2];

 aG=[dg1/2*O (dg1/2+dg2)*O];
 nG=nG0;
 rlat=sum(aG(1,:));
%'ver;',
%keyboard

  if max(aitot(Ngra,:))>0
   AmaxP=max(aitot(Ngra,:));
  else
   AmaxP=Amax;
  end 
 
 if AmaxP/perio>500
  AmaxP=100*perio;
 end
 
 while aG(1,end)<=AmaxP
  aG=[aG aG(1,end)+aG0];
  nG=[nG nG0];
  rlat=rlat+perio;
 end
nG=[nG nG0(:,2)];
 ncol=size(nG,2); 
 nrig=size(Litot,1)-1+length(Lg); 
  Litotp=zeros(nrig,1);
  aitotp=zeros(nrig,ncol-1);
  nitotp=zeros(nrig,ncol);
 if Ngra==1
  pug=[1:length(Lg)];
  puv=length(Lg)+1:nrig;
  Litotp(pug,:)=Lg;
  aitotp(pug,:)=aG;
  nitotp(pug,:)=nG; 
  Litotp(puv,1)=Litot(2:end);
  aitotp(puv,1:size(aitot,2))=aitot(2:end,:);
  nitotp(puv,1:size(nitot,2))=nitot(2:end,:);   
 else
  fd=fst(:,2);
  Ngrav=Ngra;
  Ngra=sum(fd(1:Ngrav));
  Ngras=Ngra;
  pug=[1:length(Lg)]+Ngra-1;
  puv1=1:Ngra-1;
  puv2=Ngra+length(Lg):nrig;
  puv1o=1:Ngra-1;
  puv2o=(Ngra+1):length(Litot);  
  Litotp(pug,:)=Lg;
  aitotp(pug,:)=aG;
  nitotp(pug,:)=nG; 
  Litotp(puv1,1)=Litot(puv1o);
  aitotp(puv1,1:size(aitot,2))=aitot(puv1o,:);
  nitotp(puv1,1:size(nitot,2))=nitot(puv1o,:); 
  Litotp(puv2,1)=Litot(puv2o);
  aitotp(puv2,1:size(aitot,2))=aitot(puv2o,:);
  nitotp(puv2,1:size(nitot,2))=nitot(puv2o,:);    
 end 
else
 Ngra=0;
 par_grat=0;
 Litotp=Litot;
 aitotp=aitot;
 nitotp=nitot;
 rab=radii.b;
 fi11=find(shavet(:,1)==11);
 if length(fi11)>0
  rab(fi11)=0;
 end
 ratot=[radii.a rab];
 ratotc=ratot(2:end-1,:);
 shatotc=shavet(2:end-1,:);
 fi=find(shatotc==11); 
 if length(fi)>0
  aitotp=[aitotp; zeros(size(aitotp))];
  aitotp(fi,2)=ratotc(fi,2);
  nitotp(fi,3)=nitotp(fi,1);
 end 
 
end

%'qui', keyboard
%fi8=find(shavet(:,1)==-8);
fi8=find(shavet(:,1)==-8);
fi8m=fi8-1;
if length(fi8)==1
%'pro fermo IK', keyboard
 %rav=radii.array{fi8}{12};
 pgra=radii.array{fi8}{13};
 rav=pgra.rDOE;

 nav=nitot(fi8m,1:2);
 ic=0;
 ag=[];
 ng=[];
 if pgra.igra==1
  ras=0;
  Pes=0;
  Pe=pgra.period;
  DC=pgra.DC;
  t1=Pe*DC;
  t2=Pe*(1-DC);
  if DC>0
  while Pes<rav(1)-Pe 
   ic=ic+1;
   if ic==1
    Pes=Pes+t1/2;
   else
    Pes=Pes+t1;
   end
   ag=[ag Pes* [1 1] ];
   ng=[ng nav];
   ic=ic+1;
   Pes=Pes+t2;
   ag=[ag Pes* [1 1] ];
   ng=[ng fliplr(nav)];
  end
  else
   ng=nav;
   ag=rav;
  end
 end
 ic=0;
 while ic<=length(rav)-2
   ic=ic+1;
   rai=rav(ic);
   ag=[ag rai* [1 1] ];
   ng=[ng nav];
   ic=ic+1;
   rai=rav(ic);
   ag=[ag rai* [1 1] ];
   ng=[ng fliplr(nav)];
 end
 %'ap1', keyboard
 ng(end)=nav(end);
 
 if isfield(pgra,'ng')
  ng=pgra.ng;
  ag=pgra.ag;
  ng_p=ng;
  ag_p=ag;
  ng_p0=ng;
  ag_p0=ag;  
  nitotp(fi8m,1:floor(length(ng)/2))=ng(2:2:end);
  aitotp(fi8m,1:floor(length(ng)/2))=ag(2:2:end);  
%  nitotp(fi8m,1:length(ng))=ng;
%  aitotp(fi8m,1:length(ng))=ag;    
  ng_p=ng(2:2:end);
  ag_p=ag(2:2:end);  
%  'passo', keyboard
 end


 
% nitotp(fi8m,1:length(ng))=ng;
% aitotp(fi8m,1:length(ng))=ag;
% ng_p=ng;
% ag_p=ag;
end

%'vedo ret', keyboard
Litotp1=Litotp;
nitotp1=nitotp;
aitotp1=aitotp;

if (ifp==-10 | Ps.ifpstop==1)
  if isfield(Ps,'imesa')
    if Ps.imesa==1
    RmeMa=Ps.imesaRag;
    fiRa=find(aitotp>RmeMa);
%    aitotp(fiRa)=0;
   end
 end
%'prima', keyboard
if ifp~=-4
n_distrG(Litotp,aitotp,real(nitotp),n_type,coord_num,'v',Amax);
%n_distrG(Litotp,aitotp,abs(nitotp),n_type,coord_num,'v',Amax);
 title('VCSEL transverse ref. index (COMPACT)')
 xlabel('radius (um)')
 ylabel('depth (um)')
 caxis auto
 pausak
 end
end 

[shu,Lu,nu,fiQ,fiCav,fiQW,au,ntot]=unpack_fun(dv,nv,xm,radii,fst,ifp,shavet,iauto,anyf);
Lto=Lu/1000;
Lu=Lu(2:end-1)/1000;
shu=shu(2:end-1);
nu=ntot(2:end-1,:);
au=au(2:end-1,:);

 %'qui controllo', keyboard
 ireturno=0;
 if ifp==-10 | Ps.ifpstop==1
  ireturno=0;
 end 
 if ireturno==1
  return
 end
 
% 'VETE{', keyboard
 
if fish>0 
 Ngra=find(shu(:,1)==6 | shu(:,1)==-4 );
 if length(Ngra)==0 & length(find(anyf==1))==1
  Ngra=find(shu(:,1)==4);
 end
% Ngra=fish-1;
% if exist('Ngras')
%  Ngra=Ngras;
% else
%  Ngras=fish;
% end 
 
 if isfield(par_grat,'th')==1
  Lg=par_grat.th'/1000;
 else
  Lg=Litot(Ngra);
 end
 perio=par_grat.per;
 dg1=perio*par_grat.DC;
 dg2=perio*(1-par_grat.DC);
 ng1=par_grat.n1.';
 ng2=par_grat.n2.';
 O=ones(size(Lg));
 aG0=[dg1*O (dg1+dg2)*O];
 nG0=[ng1 ng2];

 aG=[dg1/2*O (dg1/2+dg2)*O];
 nG=nG0;
 rlat=sum(aG(1,:));

  if max(aitot(Ngra,:))>0
   AmaxP=max(aitot(Ngra,:));
  else
   AmaxP=Amax;
  end 
 
 if AmaxP/perio>500
  AmaxP=100*perio;
 end  
  
  
 while aG(1,end)<=AmaxP
  aG=[aG aG(1,end)+aG0];
  nG=[nG nG0];
  rlat=rlat+perio;
 end
nG=[nG nG0(:,2)];
 ncol=size(nG,2); 
 nrig=size(Lu,1)-1+length(Lg); 
  Litotp=zeros(nrig,1);
  aitotp=zeros(nrig,ncol-1);
  nitotp=zeros(nrig,ncol);
%  'nrig', keyboard
 if Ngra==1
  pug=[1:length(Lg)];
  puv=length(Lg)+1:nrig;
  Litotp(pug,:)=Lg;
  aitotp(pug,:)=aG;
  nitotp(pug,:)=nG; 
  Litotp(puv,1)=Lu(2:end);
  aitotp(puv,1:size(au,2))=au(2:end,:);
  nitotp(puv,1:size(nu,2))=nu(2:end,:);   
 else
  pug=[1:length(Lg)]+Ngra-1;
  puv1=1:Ngra-1;
  puv2=Ngra+length(Lg):nrig;
  puv1o=1:Ngra-1;
  puv2o=(Ngra+1):length(Lu);  
  Litotp(pug,:)=Lg;
  aitotp(pug,:)=aG;
  nitotp(pug,:)=nG; 
  Litotp(puv1,1)=Lu(puv1o);
  aitotp(puv1,1:size(aitot,2))=au(puv1o,:);
  nitotp(puv1,1:size(nitot,2))=nu(puv1o,:); 
  Litotp(puv2,1)=Lu(puv2o);
  aitotp(puv2,1:size(aitot,2))=au(puv2o,:);
  nitotp(puv2,1:size(nitot,2))=nu(puv2o,:);    
 end 
else
 Ngra=0;
 par_grat=0;
 Litotp=Lu;
 aitotp=au;
 nitotp=nu;
  rab=radii.b;
  fi11=find(shavet==11);
  if length(fi11)>0
   rab(fi11)=0;
  end
 ratot=[radii.a rab];
 ratotc=ratot(2:end-1,:);
 shatotc=shavet(2:end-1);
 fi=find(shatotc==11); 
 if length(fi)>0
  aitotp=[aitotp; zeros(size(aitotp))];
  aitotp(fi,2)=ratotc(fi,2);
  nitotp(fi,3)=nitotp(fi,1);
 end 
 
end

%fi8=find(shu==-8);
%if length(fi8)==1 & length(fish)==0
% fi8m=fi8;
% nitotp(fi8m,1:length(ng))=ng;
% aitotp(fi8m,1:length(ag))=ag;
%end

fi8=find(shavet(:,1)==-8);
fi8mu=find(shu==-8);
fi8m=fi8-1;

cond8=0;
if length(fi8)==1 
 if dv(fi8)>0
  cond8=1;
 end 
end
if cond8==1
%'pro fermo IK', keyboard
% rav=radii.array{fi8}{12};
 pgra=radii.array{fi8}{13};
 rav=pgra.rDOE;

 nav=nitot(fi8m,1:2);
 nitotp(fi8mu,:)=ng_p(end);
 aitotp(fi8mu,:)=0;
 nitotp(fi8mu,1:length(ng_p))=ng_p;
 aitotp(fi8mu,1:length(ng_p))=ag_p;
 %Amax=max(rav)+3;
end

  if isfield(Ps,'imesa')
    if Ps.imesa==0
    RmeMa=Ps.imesaRag;
    fiRa=find(aitotp>RmeMa);
%    aitotp(fiRa)=0;
   end
 end
 
STR.Litotp=Litotp;
STR.aitotp=aitotp;
STR.nitotp=nitotp;
STR.n_type=n_type;
STR.coord_num=coord_num;
STR.Amax=Amax;
STR.v='v';
%'qui str', keyboard
 
if Ps.ifpstop==1 %& i1D==0
 n_distrGs(STR);
%handl=n_distrGs(STR);
%n_distrG(Litotp,aitotp,abs(nitotp),n_type,coord_num,'v',Amax);

shading flat
 title('VCSEL transverse ref. index ')
 xlabel('radius (um)')
 ylabel('depth (um)')
 colormap jet
 caxis auto
pausak
end

%'fine mappa', keyboard
%'Qui per zoom', keyboard

return

if fish>0
 Ngra=fish-1;
 if isfield(par_grat,'th')==1
  Lg=par_grat.th'/1000;
 else
  Lg=Litot(Ngra);
 end
 perio=par_grat.per;
 dg1=perio*par_grat.DC;
 dg2=perio*(1-par_grat.DC);
 ng1=par_grat.n1.';
 ng2=par_grat.n2.';
 O=ones(size(Lg));
 ag0=[dg1*O (dg1+dg2)*O];
 ng0=[ng1 ng2];

 ag=[dg1/2*O (dg1/2+dg2)*O];
 ng=ng0;
 rlat=sum(ag(1,:));
 if Amax/perio>500
  AmaxP=100*perio;
 else
  if max(aitot(fish,:))>0
   AmaxP=max(aitot(fish,:));
  else
   AmaxP=Amax;
  end 
 end
 while ag(1,end)<=AmaxP
  ag=[ag ag(1,end)+ag0];
  ng=[ng ng0];
  rlat=rlat+perio;
 end
ng=[ng ng0(:,2)];
 ncol=size(ng,2); 
 nrig=size(Litot,1)-1+length(Lg); 
  Litotp=zeros(nrig,1);
  aitotp=zeros(nrig,ncol-1);
  nitotp=zeros(nrig,ncol);
 if Ngra==1
  pug=[1:length(Lg)];
  puv=length(Lg)+1:nrig;
  Litotp(pug,:)=Lg;
  aitotp(pug,:)=ag;
  nitotp(pug,:)=ng; 
  Litotp(puv,1)=Litot(2:end);
  aitotp(puv,1:size(aitot,2))=aitot(2:end,:);
  nitotp(puv,1:size(nitot,2))=nitot(2:end,:);   
 else
  pug=[1:length(Lg)]+Ngra-1;
  puv1=1:Ngra-1;
  puv2=Ngra+length(Lg):nrig;
  puv1o=1:Ngra-1;
  puv2o=(Ngra+1):length(Litot);  
  Litotp(pug,:)=Lg;
  aitotp(pug,:)=ag;
  nitotp(pug,:)=ng; 
  Litotp(puv1,1)=Litot(puv1o);
  aitotp(puv1,1:size(aitot,2))=aitot(puv1o,:);
  nitotp(puv1,1:size(nitot,2))=nitot(puv1o,:); 
  Litotp(puv2,1)=Litot(puv2o);
  aitotp(puv2,1:size(aitot,2))=aitot(puv2o,:);
  nitotp(puv2,1:size(nitot,2))=nitot(puv2o,:);    
 end 
else
 Ngra=0;
 par_grat=0;
 Litotp=Litot;
 aitotp=aitot;
 nitotp=nitot;
 ratot=[radii.a radii.b];
 ratotc=ratot(2:end-1,:);
 shatotc=shavet(2:end-1);
 fi=find(shatotc==11); 
 if length(fi)>0
  aitotp=[aitotp; zeros(size(aitotp))];
  aitotp(fi,2)=ratotc(fi,2);
  nitotp(fi,3)=nitotp(fi,1);
 end 
 
end

'dopp', keyboard
STR.Litotp=Litotp;
STR.aitotp=aitotp;
STR.nitotp=nitotp;
STR.n_type=n_type;
STR.coord_num=coord_num;
STR.Amax=Amax;

n_distrG(Litotp,aitotp,abs(nitotp),n_type,coord_num,'v',Amax);
%n_distrG(Lto,auP,real(nuP),n_type,coord_num,'v',Amax);
shading flat
 title('VCSEL transverse ref. index ')
 xlabel('radius (um)')
 ylabel('depth (um)')
pausak