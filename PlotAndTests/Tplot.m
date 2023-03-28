function []=Tplot(xx,zz,fiat,ro_mesa,Tbuf,Temp)

figure(38992)
subplot(1,3,1)
grid on
hold on
box on
plot(xx,Temp(:,fiat),'linewidth',2)
ylabel('Temp rise (K)')
xlabel('\rho (um)')
a=axis;
a(2)=2*ro_mesa;
a(4)=a(4)*1.2;
axis(a)
subplot(1,3,2)
grid on
hold on
box on
plot(zz,Temp(1,:),'linewidth',2)
a=axis;
xlabel('z (um)')

a(1)=Tbuf*.95;
a(2)=zz(end);
a(4)=a(4)*1.2;
axis(a) 
subplot(1,3,3)
grid on
hold on
box on
plot(zz,Temp(1,:),'linewidth',2)
xlabel('z (um)')
% pausak
drawnow