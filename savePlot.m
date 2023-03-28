%'ididn', keyboard

%if exist('iDyn')
% 'ididn', keyboard
% if iDyn==1
%  Cf=(1000*mode.ii_dd(end));
%  Cfi=fix(Cf);
%  fi=find(CurDynRef==Cfi);
%  if length(fi)==1
%   CurDynRef(fi)=1000;
%   CurDyn=[CurDyn Cf];
%   dinamica
%  end
% end
%end


%'entro savePlot', keyboard
modePlot.ii_dd=mode.ii_dd;
%                modePlot.len=len;
modePlot.nMaxVet=mode.nMaxVet;
modePlot.pMaxVet=mode.pMaxVet;
modePlot.n3MaxVet=mode.n3MaxVet;
modePlot.p3MaxVet=mode.p3MaxVet;
modePlot.nQW=mode.nQW;
modePlot.NMQW=mesh.NMQW;
%                modePlot.Vmeas=Vmeas;
%                modePlot.Imeas=Imeas;
%                modePlot.Imeas_res=Imeas_res;
%                modePlot.Rmeas=Rmeas;
modePlot.indv=indv;
modePlot.vv0_dd=mode.vv0_dd;
modePlot.Pst_dd=mode.Pst_dd;
modePlot.Psp_dd=mode.Psp_dd;
modePlot.vind=vind;
modePlot.DeltaTmax=mode.DeltaTmax;
modePlot.DeltaTmax_Joule=mode.DeltaTmax_Joule;
modePlot.DeltaTmax_OptAbs=mode.DeltaTmax_OptAbs;
modePlot.DeltaTmax_Ccap=mode.DeltaTmax_Ccap;
modePlot.DeltaTmax_RAD=mode.DeltaTmax_RAD;
modePlot.DeltaTmax_srhAu=mode.DeltaTmax_srhAu;

mode.DeltaTmax_Joule(indv)=DeltaTmax_Joule;
mode.DeltaTmax_srhAu(indv)=DeltaTmax_srhAu;
mode.DeltaTmax_Ccap(indv)=DeltaTmax_Ccap;
mode.DeltaTmax_RAD(indv)=DeltaTmax_RAD;
mode.DeltaTmax_OptAbs(indv)=DeltaTmax_OptAbs;

modePlot.DeltaTmax=mode.DeltaTmax;
modePlot.Gmod=mode.Gmod;
modePlot.Lmod=mode.Lmod;
modePlot.oflg=mode.oflg;
modePlot.FLos=mode.FLos;
modePlot.nmodes=mode.nmodes;
modePlot.lambda=mode.lambda;
modePlot.vv_dd=mode.vv_dd;
modePlot.PTherm=mode.PTherm;
modePlot.NVbias=NVbias;
modePlot.matgain(indv,:)=mode.matgain;
modePlot.fPdif(indv,:)=mode.fPdif;
modePlot.E2(indv,:,:)=mode.E2;
%                'in savePlot', keyboard
modePlot.dn(indv,1:length(mode.DeltaN))=mode.DeltaN';
modePlot.Tqw(indv,:)=mesh.T(mesh.inMQW{1})';
modePlot.Elqw(indv,:)=mode.nQW{end}{2}';
modePlot.Hoqw(indv,:)=mode.pQW{end}{2}';
if isfield(mode,'Fat_cap_e')
    modePlot.Fate(indv,:)=mode.Fat_cap_e';
    modePlot.Fath(indv,:)=mode.Fat_cap_h';
else
    modePlot.Fate(indv,:)=ones(size(mode.nQW{end}{2}'));
    modePlot.Fath(indv,:)=ones(size(mode.nQW{end}{2}'));
end

if isfield(mode,'Fat_cap')
    modePlot.FatCap(indv,:)=mode.Fat_cap';
else
    modePlot.FatCap(indv,:)=ones(size(mode.nQW{end}{2}'));
end
modePlot.y=mesh.ygrid*1e4;
modePlot.x=mesh.xgrid*1e4;
modePlot.yQW=mesh.yMQW{2}*1e4;

Rcdu=reshape(mode.ecb,mesh.nny,mesh.nnx);
Rvdu=reshape(mode.evb,mesh.nny,mesh.nnx);
REfcdu=reshape(mode.EFn,mesh.nny,mesh.nnx);
REfvdu=reshape(mode.EFp,mesh.nny,mesh.nnx);
REldu=reshape(mode.elec,mesh.nny,mesh.nnx);
RHodu=reshape(mode.hole,mesh.nny,mesh.nnx);
REl2du=reshape(mode.N2D,mesh.nny,mesh.nnx);
RHo2du=reshape(mode.H2D,mesh.nny,mesh.nnx);

RRaug=sum(reshape(mode.RAugerQW,mesh.nny,mesh.nnx));
RRsrh=sum(reshape(mode.RSRHQW,mesh.nny,mesh.nnx));
RRlea=sum(reshape(mode.RLeakageQW,mesh.nny,mesh.nnx));
Cn=sum(reshape(mode.Ccapn3D,mesh.nny,mesh.nnx))/NQW*1e-9*mesh.vWMQW{1};
Cp=sum(reshape(mode.Ccapp3D,mesh.nny,mesh.nnx))/NQW*1e-9*mesh.vWMQW{1};
JN_X=reshape(mode.Jn_x,mesh.nny,mesh.nnx);
JN_Y=reshape(mode.Jn_y,mesh.nny,mesh.nnx);
JP_X=reshape(mode.Jp_x,mesh.nny,mesh.nnx);
JP_Y=reshape(mode.Jp_y,mesh.nny,mesh.nnx);

%		'in savePlot', keyboard
%		RRaug=(reshape(mode.RAugerQW,mesh.nny,mesh.nnx));
%		RRaug=RRaug+(reshape(mode.RAuger,mesh.nny,mesh.nnx));
%		RRsrh=(reshape(mode.RSRHQW,mesh.nny,mesh.nnx));
%		RRsrh=RRsrh+(reshape(mode.RSRH,mesh.nny,mesh.nnx));
%		RRrad=(reshape(mode.RradQW,mesh.nny,mesh.nnx));
%		RRrad=RRrad+(reshape(mode.Rrad,mesh.nny,mesh.nnx));
%		Cnd=(reshape(mode.Ccapn3D,mesh.nny,mesh.nnx))/NQW*1e-9*mesh.vWMQW{1};
%		Cpd=(reshape(mode.Ccapp3D,mesh.nny,mesh.nnx))/NQW*1e-9*mesh.vWMQW{1};

RRaug=(reshape(mode.RAuger,mesh.nny,mesh.nnx));
RRsrh=(reshape(mode.RSRH,mesh.nny,mesh.nnx));
RRrad=(reshape(mode.Rrad,mesh.nny,mesh.nnx));
Cnd=(reshape(mode.Ccapn3D,mesh.nny,mesh.nnx))/NQW*1e-9*mesh.vWMQW{1};
Cpd=(reshape(mode.Ccapp3D,mesh.nny,mesh.nnx))/NQW*1e-9*mesh.vWMQW{1};

Rtot=RRaug+RRsrh+RRrad;
Rtn=Rtot-Cnd;
Rtp=Rtot-Cpd;


J_X=JN_X+JP_X;
J_Y=JN_Y+JP_Y;
J_Xn=JN_X;
J_Yn=JN_Y;

% figure, plot(x,RRaug,x,RRsrh,'--',x,RRlea,'.')
%xm=x*1e-2;
% WQW=7.9e-7;
% Iaug=qel*2*pi*trapz(x,x.*RRaug)*1000*WQW;
% Ilea=qel*2*pi*trapz(x,x.*RRlea)*1000*WQW;
% Isrh=qel*2*pi*trapz(x,x.*RRsrh)*1000*WQW;
% Icn=qel*2*pi*trapz(x,x.*Cn)*1000*WQW;
% Icp=qel*2*pi*trapz(x,x.*Cp)*1000*WQW;


ifigura=0;
if ifigura==1
    figure, plot(y,REfcdu(:,[1 end]),'linewidth',2)
    hold on
    ax = gca;
    ax.ColorOrderIndex = 1;
    plot(y,Rcdu(:,[1 end]))
    title('Elettroni')
    
    figure, plot(y,REfvdu(:,[1 end]),'linewidth',2)
    hold on
    ax = gca;
    ax.ColorOrderIndex = 1;
    plot(y,Rvdu(:,[1 end]))
    title('Lacune')
    
    figure,plot(mesh.ygrid*10000,mesh.xmol(1:mesh.nny))
    figure,plot(mesh.ygrid*10000,mode.elec(1:mesh.nny))
    figure,plot(mesh.ygrid*10000,mode.hole(1:mesh.nny))
end

modePlot.Rtn(indv,:,:)=Rtn;
modePlot.Rtp(indv,:,:)=Rtp;

modePlot.JX(indv,:,:)=J_X;
modePlot.JY(indv,:,:)=J_Y;
modePlot.JXn(indv,:,:)=JN_X;
modePlot.JYn(indv,:,:)=JN_Y;
modePlot.JXp(indv,:,:)=JP_X;
modePlot.JYp(indv,:,:)=JP_Y;
modePlot.Cn(indv,:)=Cn;
modePlot.Cp(indv,:)=Cp;
modePlot.Ec(indv,:)=Rcdu(:,1)';
modePlot.Ev(indv,:)=Rvdu(:,1)';
modePlot.EcH(indv,:)=Rcdu(:,fix(end/2))';
modePlot.EvH(indv,:)=Rvdu(:,fix(end/2))';
modePlot.EFc(indv,:)=REfcdu(:,1)';
modePlot.EFv(indv,:)=REfvdu(:,1)';
modePlot.El(indv,:)=REldu(:,1)';
modePlot.Ho(indv,:)=RHodu(:,1)';
%                'stop batbn', keyboard

modePlot.xmol=mesh.xmol(1:mesh.nny);

if indv>1
    modePlot.Ic(indv,:)=curr_n;
    modePlot.Iv(indv,:)=curr_p;
    modePlot.Ind=Ind;
    modePlot.Zl=Zl;
    modePlot.zLim=zLim;
    modePlot.IntCcapN=mode.IntCcapN;
    modePlot.IntCcapP=mode.IntCcapP;
    modePlot.IntRec=mode.IntREc;
    modePlot.IntRecN=mode.IntRecN;
    modePlot.IntRecP=mode.IntRecP;
    modePlot.Fleak=mode.Fleak;
end
% da togliere

modePlot.nQW=mode.nQW;
modePlot.pQW=mode.pQW;
modePlot.irel=irel;

modePlot.fat_gain=mode.fat_gain;
modePlot.Last_Workspace=Last_Workspace;
modePlot.fCondTer=fCondTer;
modePlot.fCondTerZ=fCondTerZ;
modePlot.fatt_dndT=fatt_dndT; %

modePlot.CoeffHilsum= mode.CoeffHilsum; % parte bene
modePlot.Tcap_EXP= mode.Tcap_EXP;
modePlot.ExpH= mode.ExpH;
modePlot.ExpE= mode.ExpE;
modePlot.fat_RAD= mode.fat_RAD;
modePlot.fat_gain= mode.fat_gain;
modePlot.AlTarocco= mode.AlTarocco;
modePlot.T_tauscat= mode.T_tauscat;
modePlot.tauRat= mode.tauRat;
modePlot.Deltalam= mode.Deltalam; % in nm ; % forse 3.5 è l'ideale per lo Standard
modePlot.TAROCCO= TAROCCO;
%
modePlot.idiffusioneQW= mode.idiffusioneQW; % 0 nulla per port qw; 1: solo diffusione;  2 DD; 3 DD con gamma
modePlot.iambipolar= mode.iambipolarQW;
modePlot.FAT_idiffusioneQW_E= mode.FAT_idiffusioneQW_E;
modePlot.FAT_idiffusioneQW_H= mode.FAT_idiffusioneQW_H;
modePlot.tausE= mode.tausE;
modePlot.tausH= mode.tausH;
modePlot.Gamma= sum(mode.Gamma_z);
modePlot.VelmOptions= VelmOptions;
% VelmOptions.ivett=1;
modePlot.Exp_Temp0= mode.Exp_Temp0;

modePlot.fat_ag= mode.fat_ag;     % fattore antiguiding
modePlot.ParVet= ParVet;     % fattore antiguiding


%    modePlot.Contact=Contact;
modePlot.Relief=Relief;
modePlot.Oxide=Oxide;
modePlot.ExpH=ExpH;
modePlot.N_X=N_X;
modePlot.ABSe=mode.ABSe;
modePlot.ABSh=mode.ABSh;
modePlot.ABS_Texp=mode.ABS_Texp;
modePlot.T_tauscat=T_tauscat;
modePlot.Tcap_EXP=Tcap_EXP;
modePlot.tauRat=tauRat;
modePlot.T0=T0;
modePlot.iLUT=iLUT;
modePlot.iStruttura=iStruttura;
modePlot.iTappo=iTappo;
modePlot.FattoreTauE=FattoreTauE;
modePlot.AlTarocco=AlTarocco;
modePlot.CN_Auger=CN_Auger;
modePlot.FatNP_Auger=FatNP_Auger;
modePlot.CTemp_Auger=CTemp_Auger;
modePlot.CTemp_Ion=CTemp_Ion;
modePlot.Fat_regeneration=Fat_regeneration;
modePlot.Auger_broad=Auger_broad;
modePlot.C_Temp=C_Temp;


modePlot.GLUT=mode.GLUT;
modePlot.LUT=LUT;
modePlot.NOMELUT=NOMELUT;
modePlot.structureName=structureName;
modePlot.settings=rad_setting;

modePlot.node=mesh.node;
modePlot.triangle=mesh.triangle;
modePlot.nn=mesh.nn;
modePlot.nt=mesh.nt;
modePlot.zox=StrDD.zox;
iRa=imag(StrDD.raggi);
rox=iRa(iRa>0);
%'qui rox', keyboard
modePlot.rox=StrTT.Rox;
modePlot.Contact_i=StrTT.ro_met;
modePlot.Contact_e=StrTT.ro_mesa;

modePlot.yMQW=mesh.yMQW;
modePlot.vWMQW=mesh.vWMQW;
modePlot.geom=geom;
modePlot.Fasano=mode.Fasano;
if exist('Isize')
    modePlot.Isize=Isize;
end
if exist('Isize')
    modePlot.Isize=Isize;
end

efield_rho = -diff(mode.phi(mesh.inMQW{1}))./diff(mesh.node(1,mesh.inMQW{1}));
modePlot.efield_rho(indv,:)=efield_rho;

if isfield(mode,'JnQW')
    for kp=1:NQW
        modePlot.JnQW(indv,:,kp)=mode.JnQW{kp};
        modePlot.JpQW(indv,:,kp)=mode.JpQW{kp};
    end
end

%figure,plot(mesh.xgrid(1:mesh.nnxQW{1}-1),efield)

%Efie=reshape(mode.efieldy,mesh.nny,mesh.nnx);
%'fine saveplot', keyboard
Epot=reshape(mode.phi,mesh.nny,mesh.nnx);
E0=Epot(:,1);
Ez=-diff(E0)./diff(mesh.ygrid');
modePlot.efield_z(indv,:)=Ez;
modePlot.Phi(indv,:,:)=Epot;
modePlot.elec(indv,:,:)=REldu;
Tedu=reshape(mesh.T,mesh.nny,mesh.nnx);
modePlot.Temp(indv,:,:)=Tedu;
Mob=reshape(mesh.mobn0_n,mesh.nny,mesh.nnx);
modePlot.Mob(indv,:,:)=Mob;
Nc=reshape(mesh.Nc,mesh.nny,mesh.nnx);
modePlot.Nc(indv,:,:)=Nc;
EFn=reshape(mode.EFn,mesh.nny,mesh.nnx);
modePlot.EFn(indv,:,:)=EFn;
ecb=reshape(mode.ecb,mesh.nny,mesh.nnx);
modePlot.ecb(indv,:,:)=ecb;
