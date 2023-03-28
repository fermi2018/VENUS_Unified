
clear
% close all

addpathVENUS

P=load('out\LW_MarkusN_FINALE_LgDBR_80C_FatPerd=2.mat');
A=load('out\LW_MarkusN_TJ_oxAbove_LgDBR_80C_FatPerd=2');
B=load('out\LW_MarkusN_TJ_oxBelow_LgDBR_80C_FatPerd=2.mat');

T0=P.mode.T0-273;

if double(T0)==20
    load(['MarkusN_4_T',num2str(double(T0)),'.mat'])
elseif double(T0)==50
    load(['MarkusN_4_T',num2str(double(T0)),'.mat'])
elseif double(T0)==80
    load(['MarkusN_4_T',num2str(double(T0)),'.mat'])
elseif double(T0)==110
    load(['MarkusN_4_T',num2str(double(T0)),'.mat'])
end
pa=1:10:length(Imeas);

%% IV
figure(80)
hold on,grid on,box on
plot(P.mode.vv_dd,P.mode.ii_dd*1000,'r.-','LineWidth',2)
plot(A.mode.vv_dd,A.mode.ii_dd*1000,'g.-','LineWidth',2)
plot(B.mode.vv_dd,B.mode.ii_dd*1000,'b.-','LineWidth',2)
plot(Vmeas(pa),Imeas(pa),'r.','LineWidth',2)

set(gca,'FontSize',14,'FontName','Times new roman')
axis([1.4 2.8 0 9])
xlabel('Voltage, V'),ylabel('Current, mA')
legend('pin - Ox. confined','Infinite TJ - Ox. Above','Infinite TJ - Ox. Below','location','northwest')

%% LI
figure(81)
hold on,grid on,box on
plot(P.mode.ii_dd*1000,sum(P.mode.Pst_dd,1),'r.-','LineWidth',2)
plot(A.mode.ii_dd*1000,sum(A.mode.Pst_dd,1),'g.-','LineWidth',2)
plot(B.mode.ii_dd*1000,sum(B.mode.Pst_dd,1),'b.-','LineWidth',2)
plot(Imeas(pa),Lmeas(pa),'r.','LineWidth',2)

plot(P.mode.ii_dd*1000,P.mode.Pst_dd(1,:),'r+','LineWidth',1,'markersize',4)
plot(P.mode.ii_dd*1000,P.mode.Pst_dd(2,:),'rv','LineWidth',1,'markersize',4)

plot(A.mode.ii_dd*1000,A.mode.Pst_dd(1,:),'g+','LineWidth',1,'markersize',4)
plot(A.mode.ii_dd*1000,A.mode.Pst_dd(2,:),'gv','LineWidth',1,'markersize',4)

plot(B.mode.ii_dd*1000,B.mode.Pst_dd(1,:),'b+','LineWidth',1,'markersize',4)
plot(B.mode.ii_dd*1000,B.mode.Pst_dd(2,:),'bv','LineWidth',1,'markersize',4)

set(gca,'FontSize',14,'FontName','Times new roman')
axis([0 9 0 2])
ylabel('Output optical power, mW'),xlabel('Current, mA')
legend('pin - Ox. confined','Infinite TJ - Ox. Above','Infinite TJ - Ox. Below','location','northwest')

%% WPE
figure(82)
hold on,grid on,box on
plot(P.mode.ii_dd*1000.,sum(P.mode.Pst_dd,1)./(P.mode.ii_dd*1000.*P.mode.vv_dd)*100,'r.-','LineWidth',2)
plot(A.mode.ii_dd*1000.,sum(A.mode.Pst_dd,1)./(A.mode.ii_dd*1000.*A.mode.vv_dd)*100,'g.-','LineWidth',2)
plot(B.mode.ii_dd*1000.,sum(B.mode.Pst_dd,1)./(B.mode.ii_dd*1000.*B.mode.vv_dd)*100,'b.-','LineWidth',2)
plot(Imeas(pa),Lmeas(pa)./(Vmeas(pa).*Imeas(pa))*100,'r.','LineWidth',2)

set(gca,'FontSize',14,'FontName','Times new roman')
axis([0 9 0 20])
ylabel('\eta_{WP}, %'),xlabel('Current, mA')
legend('pin - Ox. confined','Infinite TJ - Ox. Above','Infinite TJ - Ox. Below','location','south')

%% Wavelength
figure(83)
hold on,grid on,box on
plot(P.mode.ii_dd*1000,P.mode.lambda(1,:),'r','LineWidth',1.5)
plot(A.mode.ii_dd*1000,A.mode.lambda(1,:),'g','LineWidth',1.5)
plot(B.mode.ii_dd*1000,B.mode.lambda(1,:),'b','LineWidth',1.5)
plot(Cur,LAM,'ro','LineWidth',1,'markersize',4)

thP11=find(P.mode.ii_dd*1e3>3);
thA11=find(A.mode.ii_dd*1e3>3);
thB11=find(B.mode.ii_dd*1e3>3);
plot(P.mode.ii_dd(thP11:end)*1000,P.mode.lambda(2,thP11:end),'r--','LineWidth',1.5)
plot(A.mode.ii_dd(thA11:end)*1000,A.mode.lambda(2,thA11:end),'g--','LineWidth',1.5)
plot(B.mode.ii_dd(thB11:end)*1000,B.mode.lambda(2,thB11:end),'b--','LineWidth',1.5)

set(gca,'FontSize',14,'FontName','Times new roman')
axis([1 9 850 856])

xlabel('Current, mA'),ylabel('Wavelength, nm')
legend('pin - Ox. confined','Infinite TJ - Ox. Above','Infinite TJ - Ox. Below','location','northwest')

