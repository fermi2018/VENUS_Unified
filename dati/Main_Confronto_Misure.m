
clear
close all
clc

load Meas_ReliefVCSEL
LAM_rel=LAM;
C_rel=Cur;
V_rel=Vmeas;
I_rel=Imeas;
L_rel=Lmeas;
R_rel=diff(Vmeas)./diff(Imeas/1000);
R_rel=[R_rel(1) R_rel];

load Meas_StandardVCSEL_ebr3
LAM_st=LAM;
C_st=Cur;
V_st=Vmeas;
I_st=Imeas;
L_st=Lmeas;
R_st=diff(Vmeas)./diff(Imeas/1000);
R_st=[R_st(1) R_st];

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

set(gcf,'Position',[265          55        1201         929])
subplot(2,2,1)
hold on
grid on
plot(V_rel,I_rel,V_st,I_st,'linewidth',2)
ylabel('Voltage, V')
xlabel('Current, mA')
legend('Relief','Standard','Location','Best')
subplot(2,2,2)
hold on
grid on
plot(I_rel,L_rel,I_st,L_st,'linewidth',2)
ylabel('Voltage, V')
xlabel('Current, mA')
%legend('Simulation','Measurements','Location','Best')
subplot(2,2,3)
hold on
grid on
plot(C_rel,LAM_rel,'o',C_st,LAM_st,'o','linewidth',2)
ylabel('Voltage, V')
xlabel('Current, mA')
subplot(2,2,4)
hold on
grid on
plot(I_rel,R_rel,'.',I_st,R_st,'.','linewidth',2)
ylabel('Voltage, V')
xlabel('Current, mA')
ylim([60 150])
%legend('Simulation','Measurements','Location','Best')
