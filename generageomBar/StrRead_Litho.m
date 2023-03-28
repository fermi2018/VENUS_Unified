function [StrDD,StrTT,ParOpt]=StrRead_Litho(ParVet,str,iplot,Str1D)

filename=[str,'.str'];

global flgStop


[ParDD0,ParOpt,ParMore,paTJ]=Lay_ddLitho(filename);


d=ParDD0.d;

x=ParDD0.x;
rag=ParDD0.rag;

DU=ParDD0.IPar;
fize=find(DU==0);
DU(fize)=-1e-5;

IPar=-1./(sort(-fliplr(1./DU),2));
fize=find(IPar==-1e-5);
IPar(fize)=0;

sP=sum(abs(IPar));
fike=find(sP>0);


fiP=find(IPar>0);
puso=sort(IPar(fiP));
puf=[1; diff(puso)];
putot=puso(find(puf==1));
IPar=ParDD0.IPar;
%'Lay_dd  entro ', keyboard
Netch=ParVet(find(ParVet<0));
IPar(ParDD0.iNEGFdd(1:end+Netch),1)=IPar(ParDD0.iNEGFdd(1));

for k=putot'
    fi=find(ParDD0.IPar==k);
    ParDD0.rag(fi)=ParVet(k);
end

%'Prima di TJ_GeomElec ', keyboard

partjdu=squeeze(paTJ.ipar0)';
i_etch=partjdu(find(partjdu(:,2)==-9),1);
if length(ParVet)>=i_etch
 Netch_TJ=ParVet(i_etch);
else
 Netch_TJ=paTJ.Nlay_ethc;
end 
paTJ.Nlay_ethc=Netch_TJ;

i_siz=partjdu(find(partjdu(:,2)==-6),1);
if length(ParVet)>=i_siz
 Siz_TJ=ParVet(i_siz);
 paTJ.Ram=Siz_TJ;
end 

ParOpt.radii.PaTJ=paTJ;

% save LIT
ParDD=TJ_GeomElec(iplot,paTJ,ParDD0);
d=ParDD.d;

x=ParDD.x;
rag=ParDD.rag;

% 'dopo TJ_GeomElec', keyboard

if Str1D.stairs==1
    ParDD_old=ParDD;
    ParDD=Stair(ParDD);
end

StrTT.Bcond=ParMore.Bcond;
if flgStop==1
    'ParDD in StrRead', keyboard
end

StrDD.WMQW=ParOpt.d(find(ParOpt.iauto(:,1)==2))*1e-7;

if isfield(Str1D,'quasi1D') && Str1D.quasi1D==1
    ParDD.mat{1}(1)={'Line'};
    
    for imat=1:length(ParDD.mat)
        if length(ParDD.mat{imat})>1
            ParDD.mat1{imat}=ParDD.mat{imat}(1:2);
            if strcmp(ParDD.mat{imat}(1),'vacuum')
                ParDD.mat1{imat}(1)=ParDD.mat1{imat}(2);
            end
        else
            ParDD.mat1{imat}=ParDD.mat{imat}(1);
        end
    end
    ParDD.mat=ParDD.mat1';
end

if length(ParDD.iDD)==2
    mat=ParDD.mat;
    Me=mat(end);
    Me=Me{1};
    entrato=0;
    for km=1:length(Me)
        fiM=strcmp(Me(km),'Line');
        if fiM==1
            entrato=1;
            Me{km}='Ground';
        end
    end
    if entrato==1
        ME{1}=Me;
        mat(end)=ME;
        ParDD.mat=mat;
    end
end


if isfield(Str1D,'quasi1D')
    if Str1D.quasi1D==1
        rag(:,1)=ParVet(1);
        rag=rag(:,1);
        
        ParDD.xm=ParDD.xm(:,1);
    end
end

if flgStop==1
 'qui rag', keyboard
end

mesh=ParDD.mesh;
Na=ParDD.Na;
Nd=ParDD.Nd;
rep=ParDD.rep;
fi=find(rep(:,1)==-1);
rep(fi,1)=1i;
mireq=ParMore.mireq;


%  	save LITO22 % use fastUnif for a quick debug!
%  	'salvato LITO22',  keyboard
fprintf('StrDD definition\n')
if isfield(ParDD,'iNEGF_multiple')
    StrDD.iNEGF=1;  % assigned inside geom.iNEGF and used in rectmesh
end
StrReadVcselUnifiedLITHO	
 

% spessori z
%T_DD=max(zf);
%Tbuf=350; % spessore buffer (um)

%TBuf_DD=1;
%Tdbr_inf=4.7; % spessore specchio dbr inferiore (um)
%Tdbr_sup=2.721;  % spessore specchio dbr superiore
%Tmetallo=.2; % spessore contatto metallico
%Tcentrale=T_DD-Tdbr_inf-Tmetallo-Tdbr_sup-TBuf_DD;
%Tcav=Tcentrale;



iauto=ParOpt.iauto;
do=ParOpt.d;
repop=ParOpt.flsv;
doMul=do.*abs(repop(:,2));

fiDD=find(iauto(:,2)==-100);
fiCav=find(iauto(:,2)==-4);
if isempty(fiCav)
    pu_Mirup=[];
    pu_Mirdown=[];
    pu_cav=[];
else
    pu_Mirdown=fiCav(2)+1:length(do)-1;
    pu_Mirup=2:fiCav(1)-1;
    
    
    if length(fiDD)>0
        if fiDD(1)<fiCav(1)
            pu_Mirup=fiDD(1):fiCav(1)-1;
        end
        if fiDD(end)>fiCav(2)
            pu_Mirdown=fiCav(2)+1:fiDD(end);
        end
    end
    pu_cav=fiCav(1):fiCav(2);
    
end

Buf=ParMore.Buffer;

TBuf=Buf.th;
T_DD=sum(doMul(2:end))/1000;
Tdbr_sup=sum(doMul(pu_Mirup));
Tcav=sum(doMul(pu_cav));
Tdbr_inf=sum(doMul(pu_Mirdown));

%StrDD.d_dd=d_dd;

StrTT.Tbuf=TBuf;
StrTT.Tbuf_dd=ParMore.Buffer.thdd/1000;
StrTT.Tdbr_inf=Tdbr_inf/1000;
StrTT.Tdbr_sup=Tdbr_sup/1000;
StrTT.Tcav=Tcav/1000;

Mesa=ParMore.mesa;
Dpass=Mesa.Dpass;
StrTT.TrenchPass=0;
if ParMore.r_passiv<0
	StrTT.TrenchPass=1;
	Dpass=abs(ParMore.r_passiv);
end


ro_mesa=R_contact(end);
ro_pass=ro_mesa+Dpass;
ro_max=Buf.size;
StrTT.ro_pass=ro_pass; % questi ro sono tutti radiali (um)
StrTT.ro_met=R_contact(1);
StrTT.ro_mesa=ro_mesa;
StrTT.ro_max=ro_max;
StrTT.Rox=Rox;
% Contact.radius_i Contact.radius_e]);
% if iplot==1
%  'quo termico', keyboard
% end


if isfield(Str1D,'flgBTJ_lithographic') && Str1D.flgBTJ_lithographic==2
nref=ParOpt.n(:,2);
fiox=find(real(nref==1.6));
else
fiox=[];
end

if length(fiox)>0
    T_ox=sum(doMul(2:fiox(end)))/1000;
else
    T_ox=0;
end

StrDD.r_passiv=ParMore.r_passiv;
StrDD.mesh_passiv=ParMore.mesh_passiv;
StrDD.zox=T_ox;
if isempty(pu_cav)
    StrDD.fiPassiv=[];
else
 if length(fiox)>0
    StrDD.fiPassiv=length(Prag)-pu_cav(1);
 else
    StrDD.fiPassiv=0;
 end 
end

if flgStop==1
'fine Str', keyboard
end