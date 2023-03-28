mode.rAperture=ParVet(1);   % needed to define the TJ radius

[~,iRagTJ]=min(abs(mesh.xgrid*1e4-mode.rAperture));

dopD_rag=reshape(mesh.dop_d,mesh.nny,mesh.nnx);
dopA_rag=reshape(mesh.dop_a,mesh.nny,mesh.nnx);
xmol_rag=reshape(mesh.xmol,mesh.nny,mesh.nnx);

figure(1111)
set(gcf,'Position',[144  179 1107  787])
subplot(1,2,1)
axis on
grid on
hold on
box on
plot(1e7*mesh.node(2,1:mesh.nny),xmol_rag(:,1),'b','LineWidth',2)
plot(1e7*mesh.node(2,1:mesh.nny),xmol_rag(:,iRagTJ+1),'k--','LineWidth',2)
xlabel('z (nm)')
ylabel('Molar fraction')
subplot(1,2,2)
axis on
grid on
hold on
box on
plot(1e7*mesh.node(2,1:mesh.nny),dopD_rag(:,1),'b','LineWidth',2)
plot(1e7*mesh.node(2,1:mesh.nny),dopA_rag(:,1),'r','LineWidth',2)

plot(1e7*mesh.node(2,1:mesh.nny),dopD_rag(:,iRagTJ+1),'k--','LineWidth',2)
plot(1e7*mesh.node(2,1:mesh.nny),dopA_rag(:,iRagTJ+1),'k-.','LineWidth',2)

xlabel('z (nm)')
ylabel('Doping profile cm^{-3}')
legend('N_D','N_A','location','Best')
set(gca,'yscale','log')
drawnow
