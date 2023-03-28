% modeD=mode;
% meshD=mesh;
% save d1ana20 modeD meshD

figure(80),plot(mode.vv_dd,mode.ii_dd*1e3)

load d1ana20
figure(80),hold on,plot(modeD.vv0_dd,modeD.ii_dd*1e3,'--')
if modeD.Optical==1
    figure(81),plot(mode.ii_dd*1e3,mode.Pst_dd)
    figure(81),hold on,plot(modeD.ii_dd*1e3,modeD.Pst_dd,'--')
end

% puntatore=mesh.puntatore;
puntatore=mesh.puntatore+mesh.nny;

% Doping
figure
subplot(2,1,1)
hold on,grid on,box on
plot(mesh.node(2,puntatore)*1e4,mesh.dop(puntatore),'.-')
plot(meshD.node*1e4,meshD.dop,'o--','markersize',4)
xlabel('z, \mum'),ylabel('doping')
% set(gca,'yscale','log')
subplot(2,1,2)
hold on,grid on,box on
plot(mesh.dop(puntatore)-meshD.dop,'.-')

% Molar fraction
figure
subplot(2,1,1)
hold on,grid on,box on
plot(mesh.node(2,puntatore)*1e4,mesh.xmol(puntatore),'.-')
plot(meshD.node*1e4,meshD.xmol,'o--','markersize',4)
xlabel('z, \mum'),ylabel('xmol')

subplot(2,1,2)
hold on,grid on,box on
plot(meshD.node*1e4,mesh.xmol(puntatore)-meshD.xmol,'.-')

% Mesh
figure
subplot(2,1,1)
hold on,grid on,box on
plot(mesh.node(2,puntatore)*1e4,'.-')
plot(meshD.node*1e4,'o--','markersize',4)

subplot(2,1,2)
hold on,grid on,box on
plot(mesh.node(2,puntatore)*1e4-meshD.node*1e4,'.-')

figure,
plot(meshD.node,modeD.Ec,'.-')
hold on,plot(meshD.node,modeD.Ev,'.-')

plot(mesh.node(2,puntatore),mode.ecb(puntatore),'.-')
plot(mesh.node(2,puntatore),mode.evb(puntatore),'.-')