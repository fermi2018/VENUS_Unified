dop_d=reshape(mesh.dop_d,mesh.nny,mesh.nnx);
dop_a=reshape(mesh.dop_a,mesh.nny,mesh.nnx);
xmol=reshape(mesh.xmol,mesh.nny,mesh.nnx);

y=mesh.ygrid*1e4;
x=mesh.xgrid*1e4;

keyboard

figure
hold on,grid on
plot(y,dop_d(:,1),'.-')
plot(y,dop_a(:,1),'.-')
% plot(y,dop_d(:,1),'ko--')
% plot(y,dop_a(:,1),'ko--')
chold
if isfield(mesh,'yBTJ')
    plot(y(ITJ),dop_d(ITJ,1),'o')
    plot(y(ITJ),dop_a(ITJ,1),'o')
end
set(gca,'yscale','log')
xlim([StrTT.Tbuf+4.5 StrTT.Tbuf+6.5])

figure
hold on,grid on
plot(y,xmol(:,1),'.-')
plot(y,xmol(:,1),'ko--')
chold
if isfield(mesh,'yBTJ')
    plot(y(ITJ),dop_d(ITJ,1),'o')
    plot(y(ITJ),dop_a(ITJ,1),'o')
end
xlim([StrTT.Tbuf+4.5 StrTT.Tbuf+6.5])

% figure
% chold,grid on
% plot(y,dop_d(:,1)/max(dop_d(:,1)),'o')
% plot(y,dop_a(:,1)/max(dop_a(:,1)),'o')
% plot(y,xmol(:,1)/max(xmol(:,1)),'ko')
