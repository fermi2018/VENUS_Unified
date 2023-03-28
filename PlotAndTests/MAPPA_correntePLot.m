
% quiver plot
mesh.X=reshape(mesh.node(1,:),mesh.nny,mesh.nnx);
x=mesh.X*1e4;
mesh.Y=reshape(mesh.node(2,:),mesh.nny,mesh.nnx);
y=mesh.Y*1e4;

step_x=1;
step_y=1;
sX=1:step_x:size(mesh.X,2);
sY=30:step_y:size(mesh.X,1);
X0=mesh.X(1,:)*1e4;
Y0=mesh.Y(:,1)*1e4;
dx=diff(X0);
sogx=.3;
su=0;
sX=[];
for k=1:length(dx)
    su=su+dx(k);
    if su>sogx
        su=0;
        sX=[sX k];
    end
end
figure, plot(X0,'.')
hold on
plot(sX,X0(sX),'ro'), pausak

dx=diff(Y0);
sogx=.3;
su=0;
sY=[];
for k=1:length(dx)
    su=su+dx(k);
    if su>sogx
        su=0;
        sY=[sY k];
    end
end
figure, plot(Y0,'.')
hold on
plot(sY,Y0(sY),'ro'), pausak


SX=mesh.X(sY,sX)*1e4;
SY=mesh.Y(sY,sX)*1e4;


JN_X=reshape(mode.Jn_x,mesh.nny,mesh.nnx);
JN_Y=reshape(mode.Jn_y,mesh.nny,mesh.nnx);
JP_X=reshape(mode.Jp_x,mesh.nny,mesh.nnx);
JP_Y=reshape(mode.Jp_y,mesh.nny,mesh.nnx);


J_X=JN_X+JP_X;
J_Y=JN_Y+JP_Y;
% 
modep=MODEplot{1};
% 
% % Current densities
% node='off'; triangle='off'; color='on'; vpath='off'; arrows='off';
% sd=1:geom.nd; scale=[]; cmap=[]; grido='off'; cbar='off';
% 
% figure
% plot_tri_mesh(geom,mesh,[],sd,color,triangle,node,grido,cbar,cmap,scale,vpath,arrows)
% hold on
% yMQW=modep.yMQW{end}*1e4;
% Xqw=[0 10 10 0];
% thQ=((modep.yMQW{1}-modep.yMQW{end})+modep.vWMQW{1})*1e4;
% Yqw=[yMQW yMQW  yMQW+thQ yMQW+thQ];
% hhco=patch(Xqw,Yqw,'c'); set(hhco,'EdgeColor','none');
% quiver(SX,SY,J_X(sY,sX),J_Y(sY,sX),.1)
% quiver(mesh.xgrid*1e4,mesh.ygrid*1e4,(J_XN),(J_YN),.2)
% quiver(mesh.xgrid*1e4,mesh.ygrid*1e4,(J_XP),(J_YP),.2)

% axis([0 10 110 120])
% 
% xlabel('z, \mum')
% ylabel('\rho, \mum')
% axis normal
% 
% % QW 
% Cur_2D=1;
% xcm=x*1e-4; xQW=x(1:mesh.nnxQW{1}); WQW=mesh.vWMQW{1};
% xDiv=(xcm+diff(xcm(1:2))/2)*2*pi;
% jnQ_x=sum((modep.JnQW),2)*1000;
% jpQ_x=sum((modep.JpQW),2)*1000;
% jQn=SEGNO_elettroni*[jnQ_x' zeros(1,length(xcm)-length(xQW))]./xDiv;
% jQp=[jpQ_x' zeros(1,length(xcm)-length(xQW))]./xDiv;
% 
% if Cur_2D==1
%     yQWmedioE=modep.yMQW{3}*ones(size(sX0))*1e4;
%     yQWmedioH=modep.yMQW{1}*ones(size(sX0))*1e4;
%     quiver(x(sX0),yQWmedioE,jQn(sX0)*WQW,zeros(size(sX0)),.1,'b')
%     quiver(x(sX0),yQWmedioH,jQp(sX0)*WQW,zeros(size(sX0)),.1,'r')
% end

iV=input('Which current point (index)?\n');
pcor=iV;

figure,contourf(x,y,(J_X)),title('J_\rho'),shading interp,xlabel('\rho'),ylabel('z')%,set(gca,'zscale','log')
figure,contourf(x,y,(J_Y)),title('J_z'),shading interp,xlabel('\rho'),ylabel('z')%,set(gca,'zscale','log')

figure,[hC,hC]=contourf(mesh.X,mesh.Y,(J_X),-1e5:1e3:1e5),title('J_\rho')
set(hC,'linestyle','none')
caxis([-5e3,5e3])

figure,[hC,hC]=contourf(mesh.X,mesh.Y,(J_Y),-1e5:1e3:1e5),title('J_\rho')
set(hC,'linestyle','none')
caxis([-5e3,5e3])

x=mesh.xgrid;
for indy=1:mesh.nny
    
    jn_y=(JN_Y(indy,:)+JP_Y(indy,:));
    curr(indy)=trapz(x,2.*pi.*x.*jn_y);
    
    
end
figure,plot(mode.vv0_dd,1000*mode.ii_dd),xlabel('V'),title('current vs V')

figure,plot(1e7*mesh.ygrid,curr*1000),xlabel('z'),title('current vs z')

figure,plot(1e7*mesh.node(2,:),mode.HeatJoule_y),xlabel('z'),title('Joule heat z')

figure,plot(1e7*mesh.node(2,:),mode.HeatJoule_x),xlabel('z'),title('Joule heat \rho')

figure,plot(1e7*mesh.node(2,:),mode.HeatRec_nr),xlabel('z'),title('NR recomb heat')

el=mode.elec; el=el/max(el).*5;
ho=mode.hole; ho=ho/max(ho).*5;
% el=mesh.mobn_n.*mode.elec; el=el/max(el).*5;
% ho=mesh.mobp_n.*mode.hole; ho=ho/max(ho).*5;

figure,plot(1e7*mesh.node(2,:),1./mode.sigma,1e7*mesh.node(2,:),el,'r--',1e7*mesh.node(2,:),ho,'g--'),xlabel('z'),title('\rho = 1/\sigma')

ind=1:mesh.nny;
ind=1:mesh.nn;
figure,hold on,
plot(1e7*mesh.node(2,ind),mode.elec(ind),'b',1e7*mesh.node(2,ind),mode.hole(ind),'r','Linewidth',2),
% plot(1e7*mesh.node(2,ind),mode.N2D(ind),'k',1e7*mesh.node(2,ind),mode.P2D(ind),'m','Linewidth',2),
plot(1e7*mesh.node(2,ind),mode.elec(ind)+mode.N2D(ind),'k',1e7*mesh.node(2,ind),mode.hole(ind)+mode.P2D(ind),'m','Linewidth',2),
set(gca,'yscale','log')

figure,plot((mesh.node(2,end)-mesh.node(2,:))*1e7,mode.sigma)
SIGMA=reshape(mode.sigma,mesh.nny,mesh.nnx);


figure,[hC,hC]=contourf(mesh.X,mesh.Y,log10(1./SIGMA),-5:0.1:1),title('\sigma')
set(hC,'linestyle','none')
caxis([-3,0.5])


JOULE=reshape(mode.HeatJoule,mesh.nny,mesh.nnx);
figure,[hC,hC]=contourf(mesh.X,mesh.Y,log10(JOULE),-16:0.1:-3),title('Joule')
set(hC,'linestyle','none')
caxis([-9,-3])


NONRAD=reshape(mode.HeatRec_nr,mesh.nny,mesh.nnx);
figure,[hC,hC]=contourf(mesh.X,mesh.Y,log10(abs(NONRAD)),-30:0.1:-3),title('Non-rad sum')
set(hC,'linestyle','none')
caxis([-9,-3])

NONRAD_QW=reshape(mode.HeatRec_nr_QW,mesh.nny,mesh.nnx);
figure,[hC,hC]=contourf(mesh.X,mesh.Y,log10(abs(NONRAD_QW)),-30:0.1:-3),title('Non-rad QW')
set(hC,'linestyle','none')
caxis([-9,-3])

NONRAD_bulk=reshape(mode.HeatRec_nr_bulk,mesh.nny,mesh.nnx);
figure,[hC,hC]=contourf(mesh.X,mesh.Y,log10(abs(NONRAD_bulk)),-30:0.1:-3),title('Non-rad bulk')
set(hC,'linestyle','none')
caxis([-9,-3])



% save('Workspace_2_5_V.mat')

% vlam=[]; vind=[];
% for indplot=1:length(VELMInfo)
%     vind=[vind,VELMInfo(indplot).indVoltage];
%     vlam=[vlam,VELMInfo(indplot).vlambda];
% end
%
%
% voltind=1:10:136;
% for indplot=1:length(voltind)
%
%     figure(289)
%     hold on
%     grid on
%     box on
%     ind=voltind(indplot);
%     np=(mode.nQW{ind}+mode.pQW{ind})./2./mesh.WQW;
%     plot(1e7*mesh.xgrid,np)
%
% end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Saving output data
% if((isfield(mode,'symmetry'))&&(strcmp(mode.symmetry,'Cylindrical-X')||strcmp(mode.symmetry,'Cylindrical-Y')))
%     save(['output_' structureName '_cylindrical.mat'],'mesh','mode','geom')
% else
%     save(['output_' structureName '.mat'],'mesh','mode','geom')
% end


return

figure(2)
set(gcf,'Position',[425 164 814 646])
axis on
grid on
hold on
box on
plot(mesh.node(2,:)*1e7,mode.ecb,'b','LineWidth',2)
plot(mesh.node(2,:)*1e7,mode.evb,'b','LineWidth',2)
plot(mesh.node(2,:)*1e7,mode.EFn,'k-.','LineWidth',2)
plot(mesh.node(2,:)*1e7,mode.EFp,'k-.','LineWidth',2)
xlabel('z (nm)')
ylabel('Band diagram (eV)')
set(gca,'FontSize',14)

figure(3)
set(gcf,'Position',[369 320 1031 420])
subplot(1,2,1)
axis on
grid on
hold on
box on
plot(mesh.node(2,:)*1e7,mode.elec,'LineWidth',2)
if(mode.oflg)
    plot(mesh.node(2,:)*1e7,mode.elec+mode.N2D,'k-','LineWidth',2)
end
xlabel('z (nm)')
ylabel('n(z) (1/cm^3)')
set(gca,'FontSize',14)
set(gca,'YScale','log')
subplot(1,2,2)
axis on
grid on
hold on
box on
plot(mesh.node(2,:)*1e7,mode.hole,'LineWidth',2)
if(mode.oflg)
    plot(mesh.node(2,:)*1e7,mode.hole+mode.P2D,'k-','LineWidth',2)
end
xlabel('z (nm)')
ylabel('p(z) (1/cm^3)')
set(gca,'FontSize',14)
set(gca,'YScale','log')

