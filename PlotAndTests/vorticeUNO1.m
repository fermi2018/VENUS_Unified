

 figure(F1003)
 set(F1003,'pos',[810    30   526   923])
 subplot(411)

JN_Y=mode.JYn;


 Vcori=Vcor([1 end]);
 Vcori=Vcor;
 plot(yqw,squeeze(JN_Y(Vcori,:,1)),'linewidth',2)
  hold on
grid on 
   %ax = gca;
   %ax.ColorOrderIndex = 1;
%  plot(yqw,squeeze(JP_Y(Vcori,:,1)),'--','linewidth',2)  
  xlim(XL)

   ylabel(' J_z, mA/cm^2')

 subplot(412)
 

 plot(yqw,mode.El(Vcori,:),'linewidth',2)
grid on
   hold on
   % ax = gca;
   %ax.ColorOrderIndex = 1;
% semilogy(yqw,mode.Ho(Vcori,:),'--','linewidth',2)
 %hold on
 %ylim([1e14 1e19])
 %xlim(y(end)*[.98 1])
   xlim(XL)
 ylabel(' N-P  (1/cm^3)')


   
 subplot(413)
   plot(yqw(1:end-1),(mode.efield_z(Vcori,:)),'linewidth',2)
   hold on
   %   ax = gca;
   %ax.ColorOrderIndex = 1;
   grid
    xlabel(' z, nm')
 ylabel(' E_z (\rho=0), V/cm') 
  xlim(XL) 
 
  subplot(414)
    plot(yqw,(mode.Phi(Vcori,:,1)),'linewidth',2)
    hold on
   %    ax = gca;
   %ax.ColorOrderIndex = 1;
    grid
     xlabel(' z, nm')
 ylabel(' \phi (\rho=0), V') 
 
  xlim(XL)
 xlabel(' z around QW, nm ')
