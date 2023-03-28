 figure(1002)
 set(1002,'pos',[610    30   526   923])
 subplot(311)
 
 plot(yqw,squeeze(JN_Y(Vcor,:,1)),'linewidth',2)
  hold on
   ax = gca;
   ax.ColorOrderIndex = 1;
  plot(yqw,squeeze(JP_Y(Vcor,:,1)),'--','linewidth',2)  
  xlim([-550 150])

   ylabel(' J_z, mA/cm^2')
grid
 subplot(312)
 

 semilogy(yqw,mode.El(Vcor,:),'linewidth',2)
   hold on
    ax = gca;
   ax.ColorOrderIndex = 1;
 semilogy(yqw,mode.Ho(Vcor,:),'--','linewidth',2)
 %hold on
 ylim([1e14 1e19])
 %xlim(y(end)*[.98 1])
   xlim([-550 150])
 ylabel(' N-P  (1/cm^3)')

 grid

   
 subplot(313)
   plot(yqw(1:end-1),cumsum(mode.efield_z(Vcor,:)),'linewidth',2)
   grid
    xlabel(' z, nm')
 ylabel(' \phi (\rho=0), V') 
  xlim([-550 150])
 xlabel(' z around QW, nm ')
