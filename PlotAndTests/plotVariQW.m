% QW plots

figure,plot(mesh.xgrid,mesh.T(mesh.inMQW{1}),mesh.xgrid,mesh.T(mesh.inMQW{3}))
title('Temp in QWs')
pausak
figure,plot(mesh.xgrid,mode.nQW{end}{1},mesh.xgrid,mode.nQW{end}{2},mesh.xgrid,mode.nQW{end}{3})
title('Electrons in QWs')
pausak
figure,plot(mesh.xgrid,mode.pQW{end}{1},mesh.xgrid,mode.pQW{end}{2},mesh.xgrid,mode.pQW{end}{3})
title('Holes in QWs')
pausak
Rc=reshape(mode.ecb,mesh.nny,mesh.nnx);
Rv=reshape(mode.evb,mesh.nny,mesh.nnx);
y=mesh.ygrid*1e4;
x=mesh.xgrid*1e4;
figure, plot(y,Rv(:,[1 end])), hold on, plot(y,Rc(:,[1 end]))

y=mesh.ygrid*1e4;
El2=reshape(mode.N2D,mesh.nny,mesh.nnx);
Ho2=reshape(mode.P2D,mesh.nny,mesh.nnx);
figure, semilogy(y,El2(:,1),y,Ho2(:,1))

El=reshape(mode.elec,mesh.nny,mesh.nnx);
Ho=reshape(mode.hole,mesh.nny,mesh.nnx);
figure, semilogy(y,El2(:,1)+El(:,1),y,Ho2(:,1)+Ho(:,1))

m=sum(geom.div_x(1:3))

[du,fim]=min(abs(y-mesh.yMQW{2}*1e4));
figure, plot(x,El(fim,:),x,Ho(fim,:)), pausak
 figure,plot(m*mesh.node(2,mesh.ICAV),'ro')
 zmin=min(mesh.node(2,mesh.ICAV));
 zmax=max(mesh.node(2,mesh.ICAV));
 fi=find(mesh.ygrid>=zmin & mesh.ygrid<=zmax );
 figure, plot(y), hold on, plot(fi,y(fi),'r.') 


% tag punti cavità
nnxCAV = mesh.nnxQW{1};
nnyCAV = length(mesh.ICAV)/nnxCAV;
ICAV = reshape(mesh.ICAV,nnyCAV,nnxCAV);
figure,plot(mesh.node(1,mesh.inMQW{1}),mode.elec(ICAV)) % lungo rho
figure,plot(mesh.node(1,mesh.inMQW{1}),mode.hole(ICAV)) % lungo rho
figure,plot(mesh.node(2,mesh.inMQW{1}),mode.elec(ICAV),'.') % lungo z


ficav=fi;
figure, plot(x,El(ficav,:)), pausak
figure, plot(x,Ho(ficav,:)), pausak


JN_X=reshape(mode.Jn_x,mesh.nny,mesh.nnx);
JN_Y=reshape(mode.Jn_y,mesh.nny,mesh.nnx);
JP_X=reshape(mode.Jp_x,mesh.nny,mesh.nnx);
JP_Y=reshape(mode.Jp_y,mesh.nny,mesh.nnx);

NQW=mesh.NMQW;

RRaug=(reshape(mode.RAuger,mesh.nny,mesh.nnx));
RRsrh=(reshape(mode.RSRH,mesh.nny,mesh.nnx));
RRrad=(reshape(mode.Rrad,mesh.nny,mesh.nnx));
Cnd=(reshape(mode.Ccapn3D,mesh.nny,mesh.nnx))/NQW*1e-9*mesh.vWMQW{1};
Cpd=(reshape(mode.Ccapp3D,mesh.nny,mesh.nnx))/NQW*1e-9*mesh.vWMQW{1};

Rtot=RRaug+RRsrh+RRrad;
Rtn=Rtot-Cnd;
Rtp=Rtot-Cpd;

J_X=JN_X+JP_X;
J_Y=JN_Y+JP_Y;
J_Xn=JN_X;
J_Yn=JN_Y;

leQW=length(mode.nQW{end}{2});
xQW=x(1:leQW);
    
figure, plot(x,-J_Xn(fim,:),'linewidth',2),
hold on
plot(xQW,-mode.JnQW{2},'linewidth',2), 
xlabel('radial coord.')
ylabel('Electron flux')
legend(' 3D ',' 2D')
pausak
%figure, plot(y,Rc-Rv)
dX=100e-7;
xlim([mesh.yMQW{3}-dX,mesh.yMQW{1}+dX]*1e4)
title(' Bands in Cavity')

figure,hold on,plot(mesh.xgrid,VELMInfo(end).E2), pausak

% plots doping e Free Carrier
% figure,plot(mesh.ygrid,mode.dop_dp(1:mesh.nny),mesh.ygrid,mesh.dop_d(1:mesh.nny),mesh.ygrid,mode.elec(1:mesh.nny)+10,'k'),set(gca,'yscale','log')
 %figure,plot(mesh.ygrid,mode.dop_dp(1:mesh.nny),mesh.ygrid,mesh.dop_d(1:mesh.nny),mesh.ygrid,mode.elec(1:mesh.nny)+10,'k'),legend('ND+','ND','elec'),set(gca,'yscale','log')
% figure,plot(mesh.ygrid,mode.dop_dp(1:mesh.nny),mesh.ygrid,mesh.dop_d(1:mesh.nny),mesh.ygrid,mode.elec(1:mesh.nny)+10,'k'),set(gca,'yscale','log'),ylim([1e17,1e20])
 
 figure,plot(mesh.ygrid,mode.dop_dp(1:mesh.nny),mesh.ygrid,mesh.dop_d(1:mesh.nny),mesh.ygrid,mode.elec(1:mesh.nny)+10,'k'),
 legend('ND+','ND','elec','location','best'),set(gca,'yscale','log')
 figure,plot(mesh.ygrid,mode.dop_am(1:mesh.nny),mesh.ygrid,mesh.dop_a(1:mesh.nny),mesh.ygrid,mode.hole(1:mesh.nny)+10,'k'),
 legend('AM+','AM','hole','location','best'),set(gca,'yscale','log')
 
  figure,ind=25;plot(mesh.ygrid,mode.ecb(mesh.nny*ind+[1:mesh.nny]),'b',mesh.ygrid,mode.evb(mesh.nny*ind+[1:mesh.nny]),'r',mesh.ygrid,mode.EFn(mesh.nny*ind+[1:mesh.nny]),'k--',mesh.ygrid,mode.EFp(mesh.nny*ind+[1:mesh.nny]),'k-.','LineWidth',2)


% potenziale

PHI=reshape(mode.phi,mesh.nny,mesh.nnx);
 figure, contourf(mesh.X*1e4,mesh.Y*1e4,PHI)