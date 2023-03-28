load(['C:\Users\albig\Politecnico\Dottorato\3b_VENUS\03_20200429_D1ANA_SS_quasi1D\IntermediateResults\D1ANA',num2str(mode.vv_dd(indv)),'V.mat'])
    
y=modePlot.y;
lim = [114.2 116.2];

figure(2)
set(2,'pos',[141         136        1275         829 ])
subplot(223)
 plot(y,modePlot.Ec(end,:),'.','linewidth',2)
 hold on, grid on, box on
 plot(y,modePlot.Ev(end,:),'.','linewidth',2)
 plot(exportMesh.node*1e4,exportMode.Ec,'.','linewidth',2)
 plot(exportMesh.node*1e4,exportMode.Ev,'.','linewidth',2)

 legend('E_c VENUS1', 'E_v', 'E_c D1ANA', 'E_v','location','northwest') 
 
 plot(y,modePlot.EFc(end,:),'k.',y,modePlot.EFv(end,:),'k--')
 plot(exportMesh.node*1e4,exportMode.EFn,'g.',exportMesh.node*1e4,exportMode.EFp,'g--')
 xlim([lim(1),lim(2)])
 title(' Energy band diagram - Quasi 1D')
 xlabel(' z around QW (\mum) ')


subplot(222)
hold on, grid on, box on
   xlabel(' z, um')
 ylabel(' Doping level, cm^{-3}')
 ylim([2e17 2e19])
 xlim([lim(1),lim(2)])
plot(1e4*mesh.node(2,1:mesh.nny),mesh.dop_d(1:mesh.nny),'o')
plot(1e4*mesh.node(2,1:mesh.nny),mesh.dop_a(1:mesh.nny),'o')

plot(exportMesh.node*1e4,exportMesh.dop_d,'o','markersize',3)
plot(exportMesh.node*1e4,exportMesh.dop_a,'o','markersize',3)
set(gca,'yscale','log')
  
  legend('N_D VENUS1','N_A','N_D D1ANA','N_D','location','southeast')
  xlim([lim(1),lim(2)])
 
 subplot(224)
 hold on, grid on, box on
plot(y,modePlot.El(end,:),'o','markersize',3)
plot(y,modePlot.Ho(end,:),'o','markersize',3)
plot(exportMesh.node*1e4,exportMode.elec,'o','markersize',3)
plot(exportMesh.node*1e4,exportMode.hole,'o','markersize',3)

legend('n, VENUS1','p','n, D1ANA','p')
ylim([1e17 1e20])
 xlim([lim(1),lim(2)])
 ylabel(' N-P  (1/cm)')
 xlabel(' z around QW (um) ')
 
 set(gca,'yscale','log')
 
  subplot(221)
 hold on, grid on, box on
 plot(y,modePlot.Phi(end,:,end),'.')
 
 plot(exportMesh.node*1e4,exportMode.phi,'.')
   legend('VENUS1','D1ANA')
 xlim([lim(1),lim(2)])
 ylabel(' Electrostatic Potential (eV)')
 xlabel(' z (um) ')

  figure
 hold on, grid on, box on
 plot(exportMesh.node*1e4,exportMesh.xmol,'.')
 plot(1e4*mesh.node(2,1:mesh.nny),mesh.xmol(1:mesh.nny),'o')
 legend('D1ANA','VENUS1')
 xlim([110 120])
 
 keyboard