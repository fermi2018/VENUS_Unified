h=figure,
set(h,'pos',[247         414        1427         525])

Tv=20:20:80
for T=Tv

 eval(['load Markus4_T',num2str(T)])
 subplot(131)
hold on
 plot(Vmeas,Imeas,'linewidth',2)
ylabel('Current')
xlabel(' Voltage')
 axis([1.4 2.6 0 12])
  subplot(132)
hold on
 plot(Imeas,Lmeas,'linewidth',2)
xlabel('Current')
ylabel(' Power')
xlim([0 12])
  subplot(133)

  Rmeas=1000*diff(Vmeas)./diff(Imeas);
 %plot(Imeas(1:length(Rmeas)),log10(Rmeas),'.','linewidth',2)
 semilogy(Imeas(1:length(Rmeas)),(Rmeas),'.','linewidth',2)
hold on
 xlabel('Current')
ylabel(' Res')
 axis([0 4 50 1000])
end 

legend(num2str(Tv'),'location','best')

%axis([0 4 0 2])