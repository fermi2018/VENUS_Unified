% close all

if isfield(mode,'AM')

G_ss=mode.Gss;
C_ss=mode.Css;

% R_ss=1./G_ss;
% Z_ss=zeros(size(R_ss));
Rs =46; 
Cp=130e-15;
% Rs =20; 
% C_p=180e-15;

omega = mode.fvet*2*pi;

% for iF = 1:length(mode.fvet)
%     Z_ss(:,iF) = 1./(1./(1./(G_ss(:,iF)+1i*omega(iF)*C_ss(:,iF))+Rs)+1i*omega(iF)*Cp);
% end

Z0=50; 
% S11=(Z_ss-Z0)./(Z_ss+Z0);

AM=mode.AM;
fGH=mode.fvet/1e9;

figure
set(gcf,'Position',[  379         522        1307         420])
subplot(1,3,1)
semilogx(fGH,G_ss,'.-')
grid on
hold on
xlabel('Frequency, GHz')
ylabel('Differential conductance, S')
% title(['Current: ',num2str(mode.CurDyn,' %0.1f '),' mA'])
% 
subplot(1,3,2)
semilogx(fGH,C_ss*1e15,'.-')
grid on
hold on
xlabel('Frequency, GHz')
ylabel('Differential capacitance, fF')
% title(['Current: ',num2str(mode.CurDyn,' %0.1f '),' mA'])
% 

subplot(1,3,3)
semilogx(fGH,10*log10(abs(AM)),'.-')
grid on
hold on
xlabel('Frequency, GHz')
ylabel('AM response, normalized')
set(gca,'xscale','linear'),xlim([0 20])
legend(num2str(mode.CurDyn',' %0.1f '),'location','best')

 %%
% figure(1),subplot(2,3,4),hold on
% plot(mode.CurDyn,real(Z_ss(:,1)),'k*')

% figure
% subplot(1,3,1)
% % plot(fGH,10*log10(abs(S11/S11(1)).^2),'.-')
% plot(fGH,real(S11),'.-')
% grid on,box on,hold on
% % xlabel('frequency, GHz'),ylabel('|S_{11}|^2')
% xlabel('frequency, GHz'),ylabel('Re\{S_{11}\}')
% xlim([0 20])
% subplot(1,3,2)
% plot(fGH,imag(S11),'.-')
% grid on,box on,hold on
% xlabel('frequency, GHz'),ylabel('Im\{S_{11}\}')
% xlim([0 20])
% 
% for iI = 1:length(mode.CurDyn)
%     subplot(1,3,3)
%     P=polar(angle(S11(iI,:)),abs(S11(iI,:)));
%     title('S_{11}')
%     set(P,'linewidth',2)
%     hold on
% end

% figure
% subplot(1,2,1)
% plot(fGH,real(Z_ss),'.-')
% grid on
% xlim([0 20])
% xlabel('frequency, GHz'),ylabel('Re\{Z\}, \Omega')
% subplot(1,2,2)
% plot(fGH,imag(Z_ss),'.-')
% grid on
% xlim([0 20])
% xlabel('frequency, GHz'),ylabel('Im\{Z\}, \Omega')

% subplot(1,3,3)
% grid on,box on,hold on
% plot(fGH,abs(Z_ss),'.-')
% xlim([0 20])
% xlabel('frequency, GHz'),ylabel('|Z|, \Omega')
% legend(num2str(mode.CurDyn',' %0.1f '),'location','best')

% figure
% grid on, hold on, box on 
% plot(fGH,abs(Z_ss),'.-','linewidth',2)
% xlim([0 20])
% xlabel('frequency, GHz'),ylabel('|Z|, \Omega')
% legend(num2str(mode.CurDyn',' %0.1f '),'location','best')

% subplot(2,2,4)
% grid on,box on,hold on
% plot(fGH,10*log10(abs(S21/S21(1)).^2),'.-')
% xlabel('frequency, GHz'),ylabel('|S_{21}|^2')

end