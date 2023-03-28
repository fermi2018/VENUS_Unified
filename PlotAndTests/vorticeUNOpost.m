 figure(2002)
 set(2002,'pos',[610    30   526   923])
 subplot(411)

%JN_Y=mode.JYn;

 XL=[-30 30];
 XL=[-100 100];
 Vcori=fipst;
 plot(yqw,(JN_Y(:,1)),'linewidth',2)
  hold on
grid on 
   ax = gca;
   ax.ColorOrderIndex = 1;
%  plot(yqw,squeeze(JP_Y(Vcori,:,1)),'--','linewidth',2)  
  xlim(XL)
set(gca,'XTick',[])
   ylabel(' Jz, mA/cm2')

 subplot(412)
 

 plot(yqw,modePlot.El(Vcori,:),'linewidth',2)
grid on
   hold on
    ax = gca;
   ax.ColorOrderIndex = 1;
% semilogy(yqw,mode.Ho(Vcori,:),'--','linewidth',2)
 %hold on
 %ylim([1e14 1e19])
 %xlim(y(end)*[.98 1])
   xlim(XL)
 ylabel(' N-P  (1/cm^3)')
set(gca,'XTick',[])

   
 subplot(413)
   pp=plot(yqw(1:end-1),(modePlot.efield_z(Vcori,:)),'linewidth',2)
   hold on
   plot(yqw,Ez,'r','linewidth',2)
   grid
 
 ylabel(' E_z (r=0), volt/cm') 
 legend('From DD','Integrale rho','location','best')
  xlim(XL) 
set(gca,'XTick',[]) 
  subplot(414)
    plot(yqw,(modePlot.Phi(Vcori,:,1)),'linewidth',2)
    grid
     xlabel(' z, nm')
 ylabel(' Phi (r=0), volt') 
 
  xlim(XL)
 xlabel(' z around QW ,nm ')
