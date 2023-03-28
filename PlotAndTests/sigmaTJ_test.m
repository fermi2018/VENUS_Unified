s_LoadConstants
mobn_n=mesh.mobn0_n;
mobp_n=mesh.mobp0_n;

sigman=qel.*(mode.elec.*mobn_n);
sigmap=qel.*(mode.hole.*mobp_n);
sigma=reshape(mode.sigma,mesh.nny,mesh.nnx);

sigman=reshape(sigman,mesh.nny,mesh.nnx);
sigmap=reshape(sigmap,mesh.nny,mesh.nnx);

ITJ_finder

sigmaTJ=zeros(size(sigma));
sigmaTJ(ITJ)=mode.sigmaTJ;

y=mesh.ygrid*1e4;
x=mesh.xgrid*1e4;


[~,iRagTJ]=min(abs(mesh.xgrid*1e4-mode.rAperture));

figure
hold on,box on
plot(y,1./sigma(:,[1 iRagTJ+1]),'.-')
plot(y(ITJ),1./sigma(ITJ,[1 iRagTJ+1]),'ko')
plot(y,1./sigmaTJ(:,1),'g*')
xlim([mesh.yMQW{1}*1e4-1 mesh.yMQW{1}*1e4+1])

set(gca,'yscale','log')

figure,hold on,box on
plot(y,sigma(:,[1 iRagTJ+1]),'.-')
plot(y(ITJ),sigma(ITJ,[1 iRagTJ+1]),'ko')
plot(y,sigmaTJ(:,1),'g*')
xlim([mesh.yMQW{1}*1e4-1 mesh.yMQW{1}*1e4+1])
set(gca,'yscale','log')


figure,hold on,box on
plot(y,sigman(:,1),'.-')
plot(y,sigmap(:,1),'.-')
plot(y(ITJ),sigman(ITJ,1),'ko')
plot(y(ITJ),sigmap(ITJ,1),'k*')
xlim([mesh.yMQW{1}*1e4-1 mesh.yMQW{1}*1e4+1])
legend('\sigma_n','\sigma_p')
set(gca,'yscale','log')