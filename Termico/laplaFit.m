clear all 
close all



%load fit1
%load fit2
load fit3

Y=T3D;

%xf=linspace(xx(1)

ro=xx;
dro=diff(xx);
dz=diff(zz);

Iro1=diag(1./dro);
Iro2=diag(1./dro(1:end-1));
Iz1=diag(1./dz);
Iz2=diag(1./dz(1:end-1));

dY_ro=Iro1*diff(Y,1,1);
dY_z=diff(Y,1,2)*Iz1;

dY_ro2=Iro2*diff(dY_ro,1,1);
dY_z2=diff(dY_z,1,2)*Iz2;


Mis=diag(1./ro(2:end))*dY_ro;

dYdro2=dY_ro2;
dYdz2=dY_z2;

Lap=-B0*(dYdro2(:,1:end-2)+Mis(1:end-1,1:end-2)+dYdz2(1:end-2,:));

zz2=zz(3:end);
rr2=xx(2:end-1);

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