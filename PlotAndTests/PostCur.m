iga=1;   %non considera Gamma fermi in corrente
Tar=[1 1 1];  % fattori moltiplicativi Elec, Hol, Phi

p3=reshape(mode.hole,mesh.nny,mesh.nnx);
n3=reshape(mode.elec,mesh.nny,mesh.nnx);
p2=reshape(mode.P2D,mesh.nny,mesh.nnx);
n2=reshape(mode.N2D,mesh.nny,mesh.nnx);
Pd=reshape(mode.dop_am,mesh.nny,mesh.nnx);
Nd=reshape(mode.dop_dp,mesh.nny,mesh.nnx);
Eps=reshape(mesh.epsxx_n,mesh.nny,mesh.nnx);


p30=p3(:,1);
n30=n3(:,1);
p20=p2(:,1);
n20=n2(:,1);
Pd0=Pd(:,1);
Nd0=Nd(:,1);

zcm=mesh.ygrid*1e4;
xcm=mesh.xgrid*1e4;
s_LoadConstants

Ep0=Eps(:,1);
rhot= p20+p30-Pd0-(n20+n30-Nd0);
dz=[0 diff(zcm)]';
cs=qel.*rhot.*dz./Ep0;
fi=find(isnan(cs)==1);
cs(fi)=0;
Ez=cumsum(cs);

modePlot=MODEplot{1};
%mode=modePlot;
'prima Corrrenti', pausak


%[Jn_xR,Jn_yR,Jp_xR,Jp_yR] = f_EvalCurrentDensityPost(geom,mesh,mode,iga,Tar);
%[Jn_x,Jn_y,Jp_x,Jp_y] = f_EvalCurrentDensity(geom,mesh,mode);

[Jn_xR,Jn_yR,Jp_xR,Jp_yR] = f_EvalCurrentDensity(geom,mesh,mode);
J_XNv=reshape(Jn_xR,mesh.nny,mesh.nnx);
J_YNv=reshape(Jn_yR,mesh.nny,mesh.nnx);
J_XPv=reshape(Jp_xR,mesh.nny,mesh.nnx);
J_YPv=reshape(Jp_yR,mesh.nny,mesh.nnx);

load SAP

'dopo Corrrenti', pausak
pu=length(mode.ii_dd);
pu=1;
hm=figure;
set(hm,'pos',[ 356  412  1150  525])
QVpost

fipst=length(mode.ii_dd);
leakagePost