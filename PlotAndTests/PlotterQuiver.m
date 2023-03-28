

if mode.Elementi==1 
[Jn_x,Jn_y,Jp_x,Jp_y] = f_EvalCurrentDensityElementi(geom,mesh,mode);
else
[Jn_x,Jn_y,Jp_x,Jp_y] = f_EvalCurrentDensity(geom,mesh,mode);
end



% ('Total_25Gradi') : mode.Tflg=0 , mode.Elementi=0 % vortice ?
% 'Total_25Gradi_Elementi_termico' : mode.Tflg=1 , mode.Elementi=1 % vortice ? 
% '': mode.Tflg=1 , mode.Elementi=0 % vortice ?

Jnx=reshape(modePlot.JXn(iV,:,:),mesh.nny,mesh.nnx) ; 
Jny=reshape(modePlot.JYn(iV,:,:),mesh.nny,mesh.nnx) ; 
Jpx=reshape(modePlot.JXp(iV,:,:),mesh.nny,mesh.nnx); 
Jpy=reshape(modePlot.JYp(iV,:,:),mesh.nny,mesh.nnx) ; 


for indRow=1:mesh.nny
    
Jnytot(indRow)=trapz(mesh.xgrid,Jny(indRow,:).*mesh.xgrid*2*pi) ; 
Jpytot(indRow)=trapz(mesh.xgrid,Jpy(indRow,:).*mesh.xgrid*2*pi) ;

end

 Jtot=-(Jpytot)-(Jnytot) ; 
 figure
 plot(mesh.ygrid*1e4,-(Jpytot))
 hold on, semilogy(mesh.ygrid*1e4,(Jtot),'r o','Markersize',2) ; 
 hold on, semilogy(mesh.ygrid*1e4,-(Jnytot))
 
 set(gca,'yscale','lin','fontsize',12) ; %<- la scala logaritmica fa perdere
 

%il dettaglio!
% Jx=Jnx+Jpx ;
% Jy=Jny+Jpy ; 
% Jnx=reshape(Jn_x,mesh.nny,mesh.nnx)./sqrt(Jnx.^2+Jny.^2) ; 
% Jny=reshape(Jn_y,mesh.nny,mesh.nnx)./sqrt(Jnx.^2.+Jny.^2) ; 
% 
% Jpx=reshape(Jp_x,mesh.nny,mesh.nnx)./sqrt(Jpx.^2+Jpy.^2) ; 
% Jpy=reshape(Jp_y,mesh.nny,mesh.nnx)./sqrt(Jpx.^2+Jpy.^2) ; 


% RigheSamplingFactor = 2 ;
% ColoSamplingFactor=2 ; 
% JNX=zeros(mesh.nny/ColoSamplingFactor,mesh.nnx/RigheSamplingFactor) ;
% %dirada le righe
% for indRow=1:RigheSamplingFactor:nnx
%     JNX
%     
% end
%Interpolation
xstep=0.25 ; %micron
ystep=0.25; %micron
Xmax=10 ;
Ymax=160 ; 
Xq=(0:xstep:Xmax) ; 
Yq=(0:ystep:Ymax) ; 
[XQ,YQ]=meshgrid(Xq,Yq); 
[X,Y]=meshgrid(mesh.xgrid*1e4,mesh.ygrid*1e4); 
Jx=interp2(X,Y,Jnx,XQ,YQ) ; 
Jy=interp2(X,Y,Jny,XQ,YQ) ; % Jy+

Jx=Jx+interp2(X,Y,Jpx,XQ,YQ) ; 
Jy=Jy+interp2(X,Y,Jpy,XQ,YQ)  ;
% 
%  Jy=Jy./sqrt(Jy.^2+Jx.^2) ; 
%  Jx=Jx./sqrt(Jy.^2+Jx.^2) ; 
 
Xplot=XQ(1,:) ;
Yplot=YQ(:,1) ; 

% for i=1:mesh.nny
% X=[X mesh.xgrid ] ; 
% 
% end
% 
% for i=1:mesh.nnx
%    Y=[Y mesh.ygrid] ; 
% end





% SamplingFactor=8 ;
% [X,Y]=meshgrid(mesh.xgrid*1e4,mesh.ygrid*1e4) ; 
node='off'; triangle='off'; color='on'; vpath='off'; arrows='off';
sd=1:geom.nd; scale=[]; cmap=[]; grido='off'; cbar='off';
figure
plot_tri_mesh(geom,mesh,[],sd,color,triangle,node,grido,cbar,cmap,scale,vpath,arrows)
hold on
quiver(Xplot,Yplot,Jx,Jy,0.6,'b')

% quiver(mesh.xgrid(1:SamplingFactor:end)*1e4,mesh.ygrid(1:SamplingFactor:end)*1e4,Jnx(1:SamplingFactor:end,1:SamplingFactor:end),Jny(1:SamplingFactor:end,1:SamplingFactor:end),'b')
hold on
% quiver(mesh.xgrid(1:SamplingFactor:end)*1e4,mesh.ygrid(1:SamplingFactor:end)*1e4,Jpx(1:SamplingFactor:end,1:SamplingFactor:end),Jpy(1:SamplingFactor:end,1:SamplingFactor:end),'r')
 xlim([0 10])
 ylim([112 118])
 
 figure,surf(XQ,YQ,Jx)
view(2)
shading interp
xlim([0 10])
ylim([112 118])

 figure,surf(XQ,YQ,Jy)
view(2)
shading interp
xlim([0 10])
ylim([112 118])


return 
 
 iiTJ=168:188 ;
 Ncol=0 ; 
figure,plot(mesh.ygrid*1e4,Jn_y(Ncol*mesh.nny+1:(Ncol+1)*mesh.nny))
hold on, plot(mesh.ygrid(iiTJ)*1e4,mesh.nny+Jn_y(Ncol*mesh.nny+1+iiTJ),' r o ')


Ncol=10 ; 
rvet=abs(mode.rvet) ;
figure,plot(mesh.ygrid*1e4,mesh.xmol(Ncol*mesh.nny+1:(Ncol+1)*mesh.nny))
xlim([114 116])
set(gca,'yscale','log')