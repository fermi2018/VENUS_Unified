function ParDD=Stair(ParDD)

mesh_O=ParDD.mesh;
Na_O=ParDD.Na;
Nd_O=ParDD.Nd;
x_O=ParDD.x;
xm_O=ParDD.xm;
IPar_O=ParDD.IPar;
d_O=ParDD.d;
rag_O=ParDD.rag;
rep_O=ParDD.rep;
mat_O=ParDD.mat;

dD=diff(Nd_O,1,2);
dA=diff(Na_O,1,2);
dx=diff(x_O,1,2);
G=dD+dA+dx;

ic=0;
nlongcnt=0;
nlongcntS=0;
iNEGF=[];
for ll=1:length(mesh_O)
iNa=0;
iNd=0;
ix=0;
REp=abs(rep_O(ll,1));
nlongcnt=nlongcnt+REp;

 if G(ll)~=0
  nn=mesh_O(ll);
  REpS=REp*nn;
  di=d_O(ll)/nn;
  if dD(ll)~=0
   iNd=1;
   dy=dD(ll)/nn;
   du=[Nd_O(ll,1):dy:Nd_O(ll,2)]+dy/2;
   du=du(1:end-1);
   Ndi=[du' du'];
  end
  if dA(ll)~=0
   iNa=1;
   dy=dA(ll)/nn;
   du=[Na_O(ll,1):dy:Na_O(ll,2)]+dy/2;
   du=du(1:end-1);   
   Nai=[du' du'];  
  end
  if dx(ll)~=0
   ix=1;
   dy=dx(ll)/nn;
   du=[x_O(ll,1):dy:x_O(ll,2)]+dy/2;
   du=du(1:end-1);
   xi=[du' du'];  
  end
  if nn>1000
   'lllkk'
   keyboard,
  end 
  
  for ici=1:nn
   ic=ic+1;
   if iNa
   Na(ic,:)=Nai(ici,:);
   else 
   Na(ic,:)=Na_O(ll,:);
   end
   if iNd   
    Nd(ic,:)=Ndi(ici,:);
   else
    Nd(ic,:)=Nd_O(ll,:);
   end  
   
   if ix 
    x(ic,:)=xi(ici,:);
   else
    x(ic,:)=x_O(ll,:);
   end   
   d(ic,1)=di;

   xm(ic,:)=xm_O(ll,:);
   IPar(ic,:)=IPar_O(ll,:);
   rag(ic,:)=rag_O(ll,:);
   rep(ic,:)=rep_O(ll,:);
   mat(ic)=mat_O(ll);     
   mesh(ic,1)=1;

  end

  nlongcntS=nlongcntS+REpS;  
 
 else
  ic=ic+1;
  Na(ic,:)=Na_O(ll,:);
  Nd(ic,:)=Nd_O(ll,:);
  x(ic,:)=x_O(ll,:);
  xm(ic,:)=xm_O(ll,:);
  IPar(ic,:)=IPar_O(ll,:);
  d(ic,1)=d_O(ll,1);
  rag(ic,:)=rag_O(ll,:);
  rep(ic,:)=rep_O(ll,:);
  mat(ic)=mat_O(ll);  
  mesh(ic,1)=mesh_O(ll);   
  nlongcntS=nlongcntS+REp;  

 end
 
 if isfield(ParDD,'iNEGF') & sum(nlongcnt==ParDD.iNEGF)
  iNEGF=[iNEGF nlongcntS];
%  ll, keyboard
 end 

end


%keyboard
% dD=diff(Nd,1,2);
% dA=diff(Na,1,2);
% dx=diff(x,1,2);
% G=dD+dA+dx;

 if isfield(ParDD,'iNEGF')
  ParDD.iNEGF=iNEGF;
%  ll, keyboard
 end 
 
rere=find(rep(:,1)>1);
PUo=find(diff(rere)>1)';
PUo=[PUo PUo(end)+1];
for ll=PUo
 PUf=rere(ll);
 fin=find(rep(:,1)==rep(ll,1));
 rep(fin,2)=length(fin);
end 
 
 

ParDD.mesh=mesh;
ParDD.Na=Na;
ParDD.Nd=Nd;
ParDD.x=x;
ParDD.xm=xm;
ParDD.IPar=IPar;
ParDD.d=d;
ParDD.rag=rag;
ParDD.rep=rep;
ParDD.mat=mat;

%'fine', keyboard