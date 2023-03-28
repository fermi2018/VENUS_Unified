addpathVENUS
colordef white 

clear Imeas Lmeas Rmeas Vmeas T0

T0=293-273;

% if mode.flgBTJ==0
    load(['MarkusN_4_T',num2str(double(T0)),'.mat'])
    pa=1:10:length(Imeas);
% else
%     load VELM3-TJ
%     pa=1:length(Imeas);
% end

vind=[];
for indplot=1:length(VELMInfo)
    vind=[vind,VELMInfo(indplot).indVoltage];
end

figure(80)
hold on
grid on
box on
plot(Vmeas(pa),Imeas(pa),'ko','markersize',4,'linewidth',1)

plot(mode.vv_dd,mode.ii_dd*1000,'.-','LineWidth',2,'markersize',8)
plot(mode.vv_dd(vind),mode.ii_dd(vind)*1000,'bo','markersize',4)
axis([1.4 mode.vv_dd(end)+.2 0 1000*mode.ii_dd(end)+1])
xlabel('Voltage, V')
ylabel('Current, mA')

figure(81)
hold on
grid on
box on
plot(Imeas(pa),Lmeas(pa),'ko','markersize',4,'linewidth',1)
PPst=sum(mode.Pst_dd,1)+mode.Psp_dd;

plot(mode.ii_dd*1000,PPst,'.-','LineWidth',2,'markersize',8)
plot(mode.ii_dd(vind)*1000,PPst(vind),'bo','markersize',4)

axis([0 1000*mode.ii_dd(end)+1 0 max(sum(mode.Pst_dd,1))*1.1+.1])
xlabel('Current, mA')
ylabel('Optical power, mW')

% figure(80)
% axis on
% grid on
% hold on
% box on
% plot(Vmeas(pa),Imeas(pa),'ko','markersize',4,'linewidth',1)
% ylabel('Current, mA')
% xlabel('Voltage, V')
% 
% axis([1.2 3.2 0 mode.ii_dd(end)*1e3])
% 
% figure(81)
% axis on
% grid on
% hold on
% box on
% plot(Imeas(pa),Lmeas(pa),'ko','markersize',4,'linewidth',1)
% xlabel('Current, mA')
% ylabel('Optical power, mW')
% xlim([0 mode.ii_dd(end)*1e3])
% 
% 
% figure(80),hold on,plot(mode.vv_dd,mode.ii_dd*1e3)
% figure(81),hold on,plot(mode.ii_dd*1e3,mode.Pst_dd)