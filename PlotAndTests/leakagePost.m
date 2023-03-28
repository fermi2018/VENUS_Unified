%JN_X=mode.JXn;
%JN_Y=(mode.JYn);
%JP_X=mode.JXp;
%JP_Y=(mode.JYp);


Vcor=fipst;

z=y;
ycm=z*1e-4;
xcm=x*1e-4;

ipcor=0;
hc=figure;
%set(hc,'pos',[80         421        1800         525 ])

set(hc,'pos',[141         136        1275         829 ])
subplot(221)
xqw=1000*(modePlot.y-modePlot.yQW);
 plot(xqw,modePlot.Ec(fipst,:),'linewidth',2)
 %plot(1000*(modePlot.y-modePlot.yQW),modePlot.Ec(fipst,:)-modePlot.Ev(fipst,:),'linewidth',2)
  ax = gca;
  ax.ColorOrderIndex = 1;
 hold on, plot(xqw,modePlot.Ev(fipst,:)+1.8,'--','linewidth',2)
 
 % figure
 % xqw=(modePlot.y)-350;
 % plot(xqw,modePlot.Ec(fipst,:)+1.8,'linewidth',1)
 %  ax = gca;
 %  ax.ColorOrderIndex = 1;
 % hold on, plot(xqw,modePlot.Ev(fipst,:)+1.8,'r','linewidth',1)
 % xlim([-1 9])
 % plot(xqw,modePlot.EFc(fipst,:)+1.8,'b','linewidth',2)
 % plot(xqw,modePlot.EFv(fipst,:)+1.8,'r','linewidth',2)
 
 %figure, plot(y,Rc-Rv)
 dX=80;
 xlim([-dX,dX])
 title(' Ec (cont), Ev+1.8 (dashed) in Cavity')
 xlabel(' z around QW (nm) ')
grid

%keyboard

for pcor=Vcor
 ipcor=ipcor+1; 
for indy=1:length(z)
    
    jn_y=(JN_Y(indy,:));
    jp_y=(JP_Y(indy,:));
%    figure, plot(x,jn_y,x,jp_y), pausak
    curr_n(indy)=1000*trapz(xcm,2.*pi.*xcm.*jn_y);
    curr_p(indy)=1000*trapz(xcm,2.*pi.*xcm.*jp_y);
    
end

zcm=ycm;

if Cur_2D==1
    jnQ_x=sum(squeeze(mode.JnQW(pcor,:,:)),2)*WQW*1000;
    jpQ_x=sum(squeeze(mode.JpQW(pcor,:,:)),2)*WQW*1000;
    jQn=[jnQ_x' zeros(1,length(xcm)-length(xQW))];
    jQp=[jpQ_x' zeros(1,length(xcm)-length(xQW))];
end    
zcm=ycm;
clear curx_n curx_p
for indx=1:length(x)
    
    jn_x=(JN_X(:,indx))';
    jp_x=(JP_X(:,indx))';
%    figure, plot(zcm,jn_y,zcm,jp_y), pausak
    curx_n(indx)=1000*trapz(zcm,jn_x);
    curx_p(indx)=1000*trapz(zcm,jp_x);
    
end
    


inoplot=1
if inoplot==0
 xcD=xcm*1e4;
 figure(199), 
 set(99,'pos',[ 358         419        1500         525])
 subplot(131)
 plot(xcD,2*pi*xcm.*curx_n,xcD,2*pi*xcm.*curx_p),
 if Cur_2D==1
 hold on
 ax = gca;
 ax.ColorOrderIndex = 1;
 plot(xcD,jQn,'--',xcD,jQp,'--','linewidth',2)
 legend('3D elec','3D hole','2D elec','2D hole')
 else
 legend('3D elec','3D hole') 
 end
 
 xlabel(' \rho, um')
 ylabel(' Current \rho, mA')
 subplot(132)
 plot(xQW(1:end-1),mode.efield_rho)
  xlabel(' \rho, um')
 ylabel(' E \rho, V/cm')
 
  subplot(133)
  yqw=y-modePlot.yQW;
  plot(yqw(1:end-1),mode.efield_z)
   xlabel(' \rho, um')
 ylabel(' E \rho, V/cm')
 pausak
 
 JTcil=2*pi*xcm.*curx_n+2*pi*xcm.*curx_p+jQn+jQp;
 'Cot'
 keyboard
 
 figure(188), plot(xcm,JTcil), 
 title(' Total Current normal to cilinder vs. cil. radius')
 
 pausak
 Jc=trapz(xcm,JTcil)
 figure(188), plot(xcm,JTcil), 
 pausak

end 


 ze=y(end);
 yo=ze-modePlot.zox;
 pcor=length(modePlot.ii_dd);
subplot(222)
yqw=y-modePlot.yQW;
 semilogy(yqw,modePlot.El(pcor,:),yqw,modePlot.Ho(pcor,:))
 %hold on
 ylim([1e14 1e19])
 %xlim(y(end)*[.98 1])
 xlim([-1 1])
 ylabel(' N-P  (1/cm)')
 xlabel(' z around QW (um) ')
 grid



NQW=modePlot.NMQW;
yMQW=modePlot.yMQW;
vWMQW=modePlot.vWMQW;

Zl=[];
Ind=[];
%dX=50e-7;
for kn=1:NQW
[val,ind]=min(abs(ycm-(yMQW{kn}-vWMQW{kn}/2)));
Zl=[Zl xqw(ind)];
Ind=[Ind ind];
[val,ind]=min(abs(ycm-(yMQW{kn}+vWMQW{kn}/2)));
ind=ind+1;
Zl=[Zl xqw(ind)];
Ind=[Ind ind];
%zIn=(yMQW{3}-dX)*1e4;
%zFi=(yMQW{1}+dX)*1e4;
%zLim=[zIn zFi];
end
fP=abs((curr_n)+(curr_p));
figure(hc)
subplot(224)
%semilogy(z,curr_n*1000,z,curr_p*1000), hold on, 
%semilogy(Zl,curr_n(Ind)*1000,'o',Zl,curr_p(Ind)*1000,'o'),
%semilogy(z,fP,'g','linewidth',2), 
semilogy(xqw,abs(curr_n),xqw,abs(curr_p)), hold on, 
semilogy(Zl,abs(curr_n(Ind)),'o',Zl,abs(curr_p(Ind)),'o'),
semilogy(xqw,fP,'g','linewidth',2), 
xlim([-dX,dX])
hold off
 title(['Current = ',num2str(1000*modePlot.ii_dd(end),2), ' (mA)'])
grid
%xlim(zLim)
 ylabel(' N-P currents (mA)')
 xlabel(' z around QW (nm) ')
 
 leQW=length(modePlot.nQW{1}{1});
 xQW=X0(1:leQW);


     subplot(223)

%     nQW=0; pQW=0;
%     for indQW=1:mesh.NMQW
%         nQW=nQW+modePlot.nQW{end}{indQW};
%         pQW=pQW+modePlot.pQW{end}{indQW};
%     end

    plot(xQW,modePlot.nQW{pcor}{2},'b.-',xQW,modePlot.pQW{pcor}{2},'r.-')
    %hold on
    %plot(x,modePlot.Cn(end,:),'c.',x,modePlot.Cn(end,:),'m.')
%    K>> figure, plot(x,mode.Cp*8e-7*1e-9/9)
%K>> figure, plot(x,mode.Cn*8e-7*1e-9/9)

    xlim([xQW(1),xQW(end)])
    xlabel('\rho, nm')
    ylabel('N-P 2D density 1e12/cm^2)')
    hold off
    grid
pausak
end

xo=1000*(modePlot.y-yo);
figure(201)
 plot(xo,modePlot.EcH(fipst,:),'linewidth',2)

  ax = gca;
  ax.ColorOrderIndex = 1;
 hold on, plot(xo,modePlot.EvH(fipst,:)+1.8,'--','linewidth',2)
 

 dX=180;
 xlim([-dX,dX])
 title(' Ec (cont), Ev+1.8 (dashed) below Oxide')
 xlabel(' z around Oxide (nm) ')

figure(200)
 set(200,'pos',[ 188          84        1100         893]) 
 subplot(221)
  yqw=1000*(y-modePlot.yQW);
  [du,iz]=min(abs(yqw));
  p=diag(modePlot.EFc(Vcor,iz))*ones(size(modePlot.EFc(Vcor,:)));
 plot(yqw,modePlot.EFc(Vcor,:)-p,'-','linewidth',2)
 grid
 hold on
  ax = gca;
  ax.ColorOrderIndex = 1;
% bande
%    p=diag(mode.EFv(Vcor,iz))*ones(size(mode.EFv(Vcor,:)));
  plot(yqw,modePlot.Ec(Vcor,:)-p,'.','linewidth',2)
% valenza
%    p=diag(mode.EFv(Vcor,iz))*ones(size(mode.EFv(Vcor,:)));
%  plot(yqw,mode.EFv(Vcor,:)-p,'--','linewidth',2)
 xlabel(' z, nm')
 ylabel(' Cond. Band diagram related to zQW, eV') 
 title('Dots Ec, Cont. Ef, both at (r=0), eV') 
 xlim([-150 150])

 
 subplot(222)
  yqw=1000*(y-modePlot.yQW);
  [du,iz]=min(abs(yqw));
  p=diag(modePlot.EFv(Vcor,iz))*ones(size(modePlot.EFv(Vcor,:)));
 plot(yqw,modePlot.EFv(Vcor,:)-p,'-','linewidth',2)
 grid
 hold on
  ax = gca;
  ax.ColorOrderIndex = 1;
% bande
%    p=diag(mode.EFv(Vcor,iz))*ones(size(mode.EFv(Vcor,:)));
  plot(yqw,modePlot.Ev(Vcor,:)-p,'.','linewidth',2)
% valenza
%    p=diag(mode.EFv(Vcor,iz))*ones(size(mode.EFv(Vcor,:)));
%  plot(yqw,mode.EFv(Vcor,:)-p,'--','linewidth',2)
 xlabel(' z, nm')
 ylabel(' Val. Band diagram related to zQW, eV') 
 title('Dots Ev, Cont. Efv, both at (r=0), eV') 
 xlim([-150 150])

subplot(223)
 semilogy(yqw,modePlot.El(Vcor,:),'linewidth',2)
 %hold on
% ylim([1e16 2e18])
 ylim([1e17 1e19])
 %xlim(y(end)*[.98 1])
 xlim([-150 150])
 ylabel(' N  (1/cm)')
 xlabel(' z around QW (um) ')
 grid

subplot(224)
 semilogy(yqw,modePlot.Ho(Vcor,:),'linewidth',2)
 %hold on
 ylim([1e17 1e19])
 %xlim(y(end)*[.98 1])
 xlim([-150 150])
 ylabel(' P  (1/cm)')
 xlabel(' z around QW (um) ')
 grid

pausak

figure(201)
 set(201,'pos',[ 188         70        1100         381]) 
 subplot(121)
 plot(xQW(1:end-1),modePlot.efield_rho(Vcor,:),'linewidth',2)
  xlabel(' \rho, um')
 ylabel(' E\rho (z=QW), V/cm')
 %title(['Correnti = ',num2str(CORRENTI)])
grid
  subplot(122)
  %yqw=y-modePlot.yQW;
  plot(yqw(1:end-1),modePlot.efield_z(Vcor,:),'linewidth',2)
  grid
   xlabel(' z, nm')
 ylabel(' Ez (r=0), V/cm') 
 xlim([-150 150])

pausak
 figure(199), 
 set(199,'pos',[  118         171        1258         431])
 subplot(121)
 xcD=xcm*1e4;
 plot(xcD,2*pi*xcm.*curx_n,xcD,2*pi*xcm.*curx_p), 
 if Cur_2D==1
 hold on
 ax = gca;
 ax.ColorOrderIndex = 1;
 plot(xcD,jQn,'--',xcD,jQp,'--','linewidth',2)
 legend('3D elec','3D hole','2D elec','2D hole')
 else
  legend('3D elec','3D hole')
 end
 title(' @ last current')
 xlabel(' \rho, um')
 ylabel(' Integral over z of Current \rho, mA')
 grid

 subplot(122)
 
 plot(yqw,(JN_Y(:,1)),'linewidth',2)
  hold on
   ax = gca;
   ax.ColorOrderIndex = 1;
  plot(yqw,(JP_Y(:,1)),'--','linewidth',2)  
  xlim([-550 150])
   xlabel(' zQW, um')
   ylabel(' Jz, mA/cm2')
grid
pausak

vorticeUNOpost
pausak