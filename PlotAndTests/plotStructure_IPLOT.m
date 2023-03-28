node='off'; triangle='off'; color='on'; vpath='off'; arrows='off';
sd=1:geom.nd; scale=[]; cmap=[]; grido='on'; cbar='off';
%sd=1:geom.nd; scale=[]; cmap=[]; grido='off'; cbar='off';
figure
plot_tri_meshPLOT(geom,mesh,[],sd,color,triangle,node,grido,cbar,cmap,scale,vpath,arrows)
title([' Grigliato: Punti X Y   ', num2str(length(mesh.xgrid)),'  x  ',num2str(length(mesh.ygrid))])
axis normal
drawnow
pausak
node='off'; triangle='off'; color='on'; vpath='off'; arrows='off';
%sd=1:geom.nd; scale=[]; cmap=[]; grido='on'; cbar='off';
sd=1:geom.nd; scale=[]; cmap=[]; grido='off'; cbar='off';
figure
plot_tri_meshPLOT(geom,mesh,[],sd,color,triangle,node,grido,cbar,cmap,scale,vpath,arrows)
title([' Grigliato: Punti XxY   ', num2str(length(mesh.brad))])

drawnow
pausak


figure(1001)
set(gcf,'Position',[144  179 1107  787])
subplot(2,2,1)
axis on
grid on
hold on
box on
plot(1e7*mesh.node(2,1:mesh.nny),mesh.xmol(1:mesh.nny),'b','LineWidth',2)
xlabel('z (nm)')
ylabel('Molar fraction')
subplot(2,2,2)
axis on
grid on
hold on
box on
plot(1e7*mesh.node(2,1:mesh.nny),mesh.dop_d(1:mesh.nny),'b','LineWidth',2)
plot(1e7*mesh.node(2,1:mesh.nny),mesh.dop_a(1:mesh.nny),'r','LineWidth',2)
xlabel('z (nm)')
ylabel('Doping profile cm^{-3}')
legend('N_D','N_A','location','Best')
set(gca,'yscale','log')
drawnow
%
disp([' ']);
disp(['Verifying molar fraction and doping properties']);
disp([' ']);
subplot(2,2,3)
axis on
grid on
hold on
box on
% Eg=mesh.Eg;
% plot(1e7*mesh.ygrid,(Eg(1:mesh.nny)),'b','LineWidth',2)
%Ec=-mesh.affinity;
%Ec = mesh.Eg/2 - mesh.T*kB/qel./2.*log(mesh.Nv./mesh.Nc) + mesh.phi_r;
Ec = mesh.Eg/2  + mesh.phi_r;
Ev=Ec-mesh.Eg;
plot(1e7*mesh.ygrid,(Ec(1:mesh.nny)),'b','LineWidth',2)
plot(1e7*mesh.ygrid,(Ev(1:mesh.nny)),'r','LineWidth',2)
xlabel('z (nm)')
ylabel('Cold Band diagram, eV')
legend('E_c','E_v','location','Best')
%set(gca,'yscale','log')
drawnow
%
if isfield(mode,'flgBTJ_lithographic')==1 && mode.flgBTJ_lithographic==2
    indLITHO=sum(geom.div_x)-geom.div_x(end)-geom.div_x(end-1)+1;
    meshLITHO=(indLITHO*mesh.nny+1):(indLITHO+1)*mesh.nny;
    
    figure(1002)
    set(gcf,'Position',[144  179 1107  787])
    subplot(2,2,1)
    axis on
    grid on
    hold on
    box on
    plot(1e7*mesh.node(2,1:mesh.nny),mesh.xmol(1:mesh.nny),'b','LineWidth',1.5)
    plot(1e7*mesh.node(2,meshLITHO),mesh.xmol(meshLITHO),'b--','LineWidth',1.5)
    XX=xlim;
    xlim(XX(end)*[.9 1])
    xlabel('z (nm)')
    ylabel('Molar fraction')
    subplot(2,2,2)
    axis on
    grid on
    hold on
    box on
    plot(1e7*mesh.node(2,1:mesh.nny),mesh.dop_d(1:mesh.nny),'b','LineWidth',1.5)
    plot(1e7*mesh.node(2,1:mesh.nny),mesh.dop_a(1:mesh.nny),'r','LineWidth',1.5)
    plot(1e7*mesh.node(2,meshLITHO),mesh.dop_d(meshLITHO),'b--','LineWidth',1.5)
    plot(1e7*mesh.node(2,meshLITHO),mesh.dop_a(meshLITHO),'r--','LineWidth',1.5)
    xlim(XX(end)*[.9 1])
    xlabel('z (nm)')
    ylabel('Doping profile cm^{-3}')
    legend('N_D','N_A','location','Best')
    set(gca,'yscale','log')
    drawnow
    %
    disp([' ']);
    disp(['Verifying molar fraction and doping properties of the two columns']);
    disp([' ']);
    subplot(2,2,3)
    axis on
    grid on
    hold on
    box on
    Ec = mesh.Eg/2  + mesh.phi_r;
    Ev=Ec-mesh.Eg;
    plot(1e7*mesh.ygrid,(Ec(1:mesh.nny)),'b','LineWidth',1.5)
    plot(1e7*mesh.ygrid,(Ev(1:mesh.nny)),'r','LineWidth',1.5)
    plot(1e7*mesh.ygrid,(Ec(meshLITHO)),'b--','LineWidth',1.5)
    plot(1e7*mesh.ygrid,(Ev(meshLITHO)),'r--','LineWidth',1.5)
    xlim(XX(end)*[.9 1])
    xlabel('z (nm)')
    ylabel('Cold Band diagram, eV')
    legend('E_c','E_v','location','Best')
end
