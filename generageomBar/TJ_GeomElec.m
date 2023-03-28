function ParDD=TJ_GeomElec(ifp,paTJ,ParDD0)


L_i=ParDD0.d;
%L_i=abs(ParDD0.d);
Netch_TJ=paTJ.Nlay_ethc;

ifgSt=0;
if ifgSt==1
'dentro TJ elec', keyboard
end

% prima cosa da fare è aggiungere una periodicità del DBR per tener conto dello scalino litografico,
% che disallinea i DBR. Lo faccio con un puntatore.

rep=ParDD0.rep;
repsum=sum(abs(rep),2);
fid=find(repsum>1 & repsum~=500);
fi=find(repsum==repsum(fid(1)));
NP=rep(fi(1),1);
fis=find(rep(:,1)==NP);
Lpair=sum(L_i(fis));

fiDBR=find(repsum==repsum(fi(1)));
puAddDBR=[1:fi(end) fiDBR' fi(end)+1:length(repsum)];
repO=rep;
rep=repO(puAddDBR,:);
rep(fiDBR+length(fiDBR),:)=repmat([1 0],length(fiDBR),1);
fiNP=find(rep(:,1)==NP);
rep(fiNP,1)=NP-1;

L_i=L_i(puAddDBR);


%puTJ=ParDD0.iNEGFdd(1:length(ParDD0.iNEGFdd)+Netch_TJ)+length(fiDBR);
%puTJ=find(repsum==500)+length(fiDBR);

puTJ=ParDD0.iNEGFdd(1:end+Netch_TJ)+length(fiDBR);

if L_i(2)<sum(L_i(puTJ))
 'Relief TROPPO spesso! Caso non ancora implementato', keyboard
 'Relief TROPPO spesso! Caso non ancora implementato', keyboard
 'Relief TROPPO spesso! Caso non ancora implementato', keyboard
 'Relief TROPPO spesso! Caso non ancora implementato', keyboard
end


pu_Below=[puTJ(end)+1:length(L_i)];



%'PAET]', keyboard
% ora ricalcolo la periodicità
if ifgSt==1
'periodicità', keyboard
end

repS=rep;
%rep(puTJ,1)=-500;

[rep_outd,d]=SubDD_litho(rep(:,1),L_i,puTJ,pu_Below,ifp);

repe=rep_outd(:,2);
fi=find(abs(repe)<.5);
repe(fi)=1;
fi=find(rep_outd(:,1)==-500);
repe(fi)=-500;

[rep_outd]=SubDD_litho(rep(:,2),L_i,puTJ,pu_Below,ifp);
clear rep
nlayr=rep_outd(:,2);
fi=find(repe==NP-1);
nlayr(fi)=length(fi);
fi=find(abs(nlayr)<.5);
nlayr(fi)=0;


LpairNew=sum(d(fi));
rep(:,1)=repe;
rep(:,2)=nlayr;

if ifgSt==1
	'fine rep TJ', 
	'Inizio MESH', 
	keyboard
end

% ora ricalcolo tutte le variabili in ParDD0

% spessori fettine

mesh=ParDD0.mesh;
mesh=mesh(puAddDBR);
[mesh_outd,d]=SubDD_lithoMesh(mesh,L_i,puTJ,pu_Below,ifp);

meshd=max(mesh_outd,[],2);
fi0=find(mesh_outd(:,2)==0);
meshd(fi0)=mesh_outd(fi0,2);
mesh=meshd;

if ifgSt==1

% frazione molare

	'FINE mesh', 
	'INizio x', keyboard
end

x_in=ParDD0.x(puAddDBR,1);
xTJ=x_in(puTJ(1));

%[xi_out,L_out]=SubDD_litho(x_in,L_i,puTJ,pu_Below,ifp);
%xi_out(find(xi_out(:,2)<0))=-1;
%fi_ext=find(xi_out(:,2)<0);
%fiTJ=find(abs(xi_out(:,1)-xTJ)<1e-10);  %puntatore TJ nuova discretizzazione

%x_in=ParDD0.x(puAddDBR,2);
%[xf_out,L_out]=SubDD_litho(x_in,L_i,puTJ,pu_Below,ifp);
%xf_out(find(xf_out(:,2)<0))=-1;

x_in=ParDD0.x(puAddDBR,1);
%[xi_out,L_out]=SubDD_litho(x_in,L_i,puTJ,pu_Below,ifp);
x_fi=ParDD0.x(puAddDBR,2);
%[xf_out,L_out]=SubDD_litho(x_fi,L_i,puTJ,pu_Below,ifp);

% save sad x_in x_fi L_i puTJ pu_Below ifp

[xi_out,xf_out]=SubDD_lithoGraded(x_in,x_fi,L_i,puTJ,pu_Below,ifp);
%xi_out(find(xi_out(:,2)<0))=-1;
%xf_out(find(xf_out(:,2)<0))=-1;
%shavet=ParDD0.Na(:,1);
%shavet_in=shavet;

x(:,1,:)=xi_out;
x(:,2,:)=xf_out;


if ifgSt==1
	'fine x', keyboard
end
%'fine x', keyboard
% xm

XM=ParDD0.xm(puAddDBR,:);
x_in=XM(:,1);
[xm_out]=SubDD_litho(x_in,L_i,puTJ,pu_Below,ifp);

for kk=2:size(XM,2)
x_in=XM(:,kk);
[xm_d]=SubDD_litho(x_in,L_i,puTJ,pu_Below,ifp);
xm_out=[xm_out xm_d(:,2)] ;
end
xm=xm_out;



if ifgSt==1
	'fine xm', keyboard
end	
% doping Nadb


%x_in=ParDD0.Na(puAddDBR,1);
%[xi_out,L_out]=SubDD_litho(x_in,L_i,puTJ,pu_Below,ifp);
%fi_ext=find(xi_out(:,2)<0);


%x_in=ParDD0.Na(puAddDBR,2);
%[xf_out,L_out]=SubDD_litho(x_in,L_i,puTJ,pu_Below,ifp);

%shavet=ParDD0.Na(:,1);
%shavet_in=shavet;

x_in=ParDD0.Na(puAddDBR,1);
%[xi_out,L_out]=SubDD_litho(x_in,L_i,puTJ,pu_Below,ifp);
x_fi=ParDD0.Na(puAddDBR,2);
%[xf_out,L_out]=SubDD_litho(x_fi,L_i,puTJ,pu_Below,ifp);

[xi_out,xf_out,L_out]=SubDD_lithoGraded(x_in,x_fi,L_i,puTJ,pu_Below,ifp);

Na(:,1,:)=xi_out;
Na(:,2,:)=xf_out;



% IPar


x_in=ParDD0.IPar(puAddDBR,1);
[xi_out,L_out]=SubDD_litho(x_in,L_i,puTJ,pu_Below,ifp);
IPar(:,1)=xi_out(:,1);
x_in=ParDD0.IPar(puAddDBR,2);
[xi_out,L_out]=SubDD_litho(x_in,L_i,puTJ,pu_Below,ifp);
IPar(:,2)=xi_out(:,2);



% doping Nd



x_in=ParDD0.Nd(puAddDBR,1);
%[xi_out,L_out]=SubDD_litho(x_in,L_i,puTJ,pu_Below,ifp);
x_fi=ParDD0.Nd(puAddDBR,2);
%[xf_out,L_out]=SubDD_litho(x_fi,L_i,puTJ,pu_Below,ifp);

[xi_out,xf_out]=SubDD_lithoGraded(x_in,x_fi,L_i,puTJ,pu_Below,ifp);

%Lp=cumsum(L_out);
%figure, plot(Lp,xi_out(:,1),'.',Lp,xi_out1(:,1),'o')
%'RIN fuori graded', keyboard


%'DOPI', keyboard

%ScriptInter

%'DOPo ScriptInter', keyboard

%shavet=ParDD0.Na(:,1);
%shavet_in=shavet;
%	'fine TJ', keyboard

Nd(:,1,:)=xi_out;
Nd(:,2,:)=xf_out;

% sistemo materiali  0, no cella; 1, semiconduttore; -1, line; -2, vuoto; -3, oxide
MM=ParDD0.mat(puAddDBR);

TabellaSemic{1}='AlGaAs';
TabellaSemic{2}='InGaAs';

clear maType
for kk=1:length(MM)
mmi=MM{kk};
 for ii=1:length(mmi)
  if strcmp(mmi(ii),'vacuum')
   maType(kk,ii)=-2;
  elseif strcmp(mmi(ii),'Line')
   maType(kk,ii)=-1;
  elseif strcmp(mmi(ii),'AlOx')
   maType(kk,ii)=-3;
  else
   clear fiMat
   for kmat=1:length(TabellaSemic)
    if strcmp(mmi(ii),TabellaSemic{kmat})
     fiMat=kmat;
	end
   end
   if ~exist('fiMat')
    'errore materiale'
   else
    maType(kk,ii)=fiMat;
   end   
  end
 end
end

%'LATT', keyboard

x_in=maType(:,1);
[xi_out]=SubDD_litho(x_in,L_i,puTJ,pu_Below,ifp);
maTout=xi_out;

for kk=2:size(maType,2)
x_in=maType(:,kk);
[xi_out]=SubDD_litho(x_in,L_i,puTJ,pu_Below,ifp);
maTout=[maTout xi_out(:,2)];
end
maTout=round(maTout);

[R,C]=find(maTout==-1);
if length(R)>0
puDD=R(end):length(maTout);
else
puDD=1:length(maTout);
end


for kk=1:length(maTout)
mmi=maTout(kk,:);
clear Mati
 for ii=1:length(mmi)
  if mmi(ii)==-2
   Mati{ii}='vacuum';
  elseif mmi(ii)==-1
   Mati{ii}='Line';
  elseif mmi(ii)==-3
   Mati{ii}='AlOx';
  else
   fisem=find(mmi>0);
   if length(fisem)>0
    for kks=fisem
     Mati{kks}=TabellaSemic{mmi(kks)};
	end 
%    fisem=find(mmi==0)
%    if length(fisem>0);
%	 for kks=fisem
%      Mati{kks}=TabellaSemic{mmi(kks)};
%	 end
%    end
   end	
  end
 end

 Mater{kk,:}=Mati;
end

if ifgSt==1
	'fine Materiali', keyboard
%	'fine Materiali', keyboard
%	'fine Materiali', keyboard
end
% ' % raggi $$$$$$$$ ', keyboard
rag=ParDD0.rag;
RA0=ParDD0.rag(puAddDBR,:);
x_in=RA0(:,1);
% 'IERI',keyboard
[xi_out]=SubDD_litho(x_in,L_i,puTJ,pu_Below,ifp);
fi_ext=find(xi_out(:,2)<0);
xi_out(fi_ext,2)=xi_out(fi_ext,1);
RaTJ=rag(ParDD0.iNEGFdd(1),1);
Ra=xi_out;
fi=find(Ra==RaTJ);
Ra(fi)=0;

ColTJ=1;
fiSemic=find(maType(:,ColTJ)==1);


RaVero=paTJ.Ram;
Ra(fiSemic:fi(1),1)=RaVero;


for kk=2:size(RA0,2)
x_in=RA0(:,kk);
[xi_out]=SubDD_litho(x_in,L_i,puTJ,pu_Below,ifp);
Ra=[Ra xi_out(:,2)];
end

for kk=1:length(Ra)
 ri=Ra(kk,:);
 rid=[Ra(kk,:) 100];
 dr=diff(rid);
 rid=ri*0;
 fi0=find(dr~=0);
 rid(1:length(fi0))=ri(fi0);
 Ra(kk,:)=rid;
end

Su=sum(Ra);
fima=find(Su~=0);
Ra=Ra(:,fima);

%'fine Raggi', keyboard
if ifgSt==1

	'fine Raggi', keyboard
%	'fine Doping', keyboard
end
%	'fine Raggi', keyboard


ParDD.d=d(puDD);
ParDD.Na=Na(puDD,:,:);
ParDD.Nd=Nd(puDD,:,:);
ParDD.x=x(puDD,:,:);
ParDD.xm=xm(puDD,:);

ParDD.rag=Ra(puDD,:);
ParDD.ragAdd=ParDD0.ragAdd;
ParDD.rep=rep(puDD,:);
ParDD.mesh=mesh(puDD);
ParDD.mat=Mater(puDD);
ParDD.iDD=-ParDD0.iDD;
ParDD.IPar=IPar(puDD,:);
% ParDD.iNEGF_multiple=ParDD0.iNEGF_multiple;
% Needs adjustment for MTJ case
%	'fine Raggi', keyboard
fiNegf=find(repe(puDD)==-500);

for itj=1:ParDD0.nBTJ
    ParDD.iNEGF_multiple{itj}=fiNegf;
end
ParDD.nBTJ=ParDD0.nBTJ;
ParDD.OxVero=ParDD0.OxVero;
ParDD.iNEGFdd=fiNegf;



%	'fine TJ elec', keyboard



