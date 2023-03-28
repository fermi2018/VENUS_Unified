clear 
close all

x=linspace(-50,50,1e6);

ybern=bern(x);
ybern1=bern1(x);

addpath('mp')
addpath('C:\Users\albig\Documents\Multiprecision Computing Toolbox 2019\')
mp.Digits(34); % number of significat digits
mp.FollowMatlabNumericFormat(true); % if "true" uses standard MATLAB format
xmp=linspace(mp('-50'),mp('50'),1e6);

ybernmp=mpbern(xmp);

figure,plot(x,ybern,'o')
hold on,plot(x,ybern1,'.')
hold on,plot(xmp,ybernmp,'--')
set(gca,'yscale','log')
legend('bern','bern1','bern mp')

figure,plot(x,abs(1-ybern./ybernmp),'.-')
set(gca,'yscale','log')
title('Rel. err. bern/bernmp')

figure,plot(x,abs(1-ybern1./ybernmp),'.-')
set(gca,'yscale','log')
title('Rel. err. bern1/bernmp')

figure,plot(x,abs(1-ybern./ybern1),'.-')
set(gca,'yscale','log')
title('Rel. err. bern/bern1')

% dbern
ydbern=dbern(x);
ydbern1=dbern1(x);
ydbernmp=mpdbern(xmp);

figure,plot(x,abs(ydbern),'o')
hold on,plot(x,abs(ydbern1),'.')
hold on,plot(x,abs(ydbernmp),'--')
set(gca,'yscale','log')
legend('bern','bern1','bern mp')

figure,plot(x,abs(1-ydbern1./ydbernmp),'.-')
set(gca,'yscale','log')
title('Rel. err. dbern1/dbernmp')

figure,plot(x,abs(1-ydbern./ydbernmp),'.-')
set(gca,'yscale','log')
title('Rel. err. dbern/dbernmp')

figure,plot(x,abs(1-ydbern./ydbern1),'.-')
set(gca,'yscale','log')
title('Rel. err. dbern/dbern1')
