clear all 
close all

clear
clear global
close all
colordef white
dbstop if error
addpath('termico')
addpath('rumenta')
addpath('generageom')
addpath('rumenta/new17Optica')


load fit
Y=T3D;

xf=linspace(xx(1),xx(end),10000);
zf=linspace(zz(1),zz(end),30000);

[X,Z]=meshgrid(xf,zf);  

Yf=interp2(xxc,zzc,Y',X,Z);

Yf=Yf';

dx=diff(xf(1:2))'
dz=diff(zf(1:2))'


dY_ro=diff(Yf,1,1)/dx;
dY_z=diff(Yf,1,2)/dz;

dY_ro2=diff(dY_ro,1,1)/dx;
dY_z2=diff(dY_z,1,2)/dz;


Mis=diag(1./xf(2:end))*dY_ro;

dYdro2=dY_ro2;
dYdz2=dY_z2;

Lap=-B0*(dYdro2(:,1:end-2)+Mis(1:end-1,1:end-2)+dYdz2(1:end-2,:));

zz2=zf(3:end);
rr2=xf(2:end-1);

figure, 
subplot(211)
plot(zz,qtot)
xlim([350 360])
a=ylim;
subplot(212)
plot(zz2,Lap)
xlim([350 360])
a1=ylim;
a1(2)=a(2);
a1(1)=-.2*a(2);
ylim(a1)
subplot(212)

pausak

figure, 
subplot(211)
plot(xx,qtot)
xlim([0 15])
a=ylim;
subplot(212)
plot(rr2,Lap)
a1=ylim;
a1(2)=a(2);
a1(1)=-.2*a(2);
ylim(a1)
xlim([0 15])