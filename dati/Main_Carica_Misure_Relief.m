
clear
close all
clc

load Meas_ReliefVCSEL
% load Meas_StandardVCSEL_fbu3

% Resistence=diff(mode.vv0_dd)./diff(mode.ii_dd);
% Resistence=[Resistence(1),Resistence];
% figure
% hold on
% grid on
% plot(Imeasinterp,Resmeasinterp)
% xlabel('Current, mA')
% ylabel('Resistance, \Omega')
% % xlim([2,mode.ii_dd(end)*1000])
% ylim([80,120])
% legend('Simulation','Measurements','Location','Best')

ResFake=0;

figure(4873)
set(gcf,'Position',[610 470 1201 483])
subplot(1,2,1)
hold on
grid on
plot(Vmeas,Imeas)
xlabel('Voltage, V')
xlabel('Current, mA')
legend('Simulation','Measurements','Location','Best')
subplot(1,2,2)
hold on
grid on
plot(Imeas,Lmeas)
xlabel('Voltage, V')
xlabel('Current, mA')
legend('Simulation','Measurements','Location','Best')

load('Meas_StandardVCSEL_ebr3.mat')
figure(4873)
set(gcf,'Position',[610 470 1201 483])
subplot(1,2,1)
hold on
grid on
plot(Vmeas,Imeas)
xlabel('Voltage, V')
xlabel('Current, mA')
legend('Simulation','Measurements','Location','Best')
subplot(1,2,2)
hold on
grid on
plot(Imeas,Lmeas)
xlabel('Voltage, V')
xlabel('Current, mA')
legend('Simulation','Measurements','Location','Best')