 Cn=mode.Cn;
 Cn(:,1)=Cn(:,2);
 Cp=mode.Cp;
 Cp(:,1)=Cp(:,2);
 %puc=[8 19];
 %puc=[4 8];
 puc=pu;
 hm=figure;
 set(hm,'pos',[349         398        1525         525 ])
subplot(121)
 plot(x,Cn(puc,:),'o')
         ax = gca;
       ax.ColorOrderIndex = 1
  hold on, plot(x,Cp(puc,:),'.')
  legend(num2str(Cpl(puc)',2))
  xlabel('\rho, \mum')
  ylabel('2D Injection, (cm^2 ns)^{-1}')
    
    pausak
    hold off

J_X=mode.JX;
J_Y=mode.JY;
y=modePlot.y;
[X,Y]=meshgrid(x,y);


 X0=x;
 Y0=y;
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
 sX0=sX;
  ifiver=0;
  
  if ifiver==1
   figure, plot(X0,'.')
   hold on
   plot(sX,X0(sX),'ro'), pausak
  end 
 
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
  sY0=sY;
  
   X0=x;
   Y0=y;
   dx=diff(X0);
   sogx=.1;
   su=0;
   sX=[];
   for k=1:length(dx)
    su=su+dx(k);
    if su>sogx
     su=0;
     sX=[sX k];
    end
   end 
   sX1=sX;
    ifiver=0;
    
    if ifiver==1
     figure, plot(X0,'.')
     hold on
     plot(sX,X0(sX),'ro'), pausak
    end 
   
    dx=diff(Y0);
    sogx=.04;
    su=0;
    sY=[];
    for k=1:length(dx)
     su=su+dx(k);
     if su>sogx
      su=0;
      sY=[sY k];
     end
    end 
  sY1=sY;
  
  
  if ifiver==1
   figure, plot(Y0,'.')
   hold on
   plot(sY,Y0(sY),'ro'), pausak
  end 
 
 
 SX=X(sY0,sX0);
 SY=Y(sY0,sX0);
 SX1=X(sY1,sX1);
 SY1=Y(sY1,sX1);
 
 geom=mode.geom;
 mesh.node=mode.node;
 mesh.triangle=mode.triangle;
 mesh.nn=mode.nn;
 mesh.nt=mode.nt;
 
 xMax=10;
  node='off'; triangle='off'; color='on'; vpath='off'; arrows='off';
 sd=1:geom.nd; scale=0.000002; cmap=[]; grido='off'; cbar='off';
 
 iar=0
 if iar==1
 Jt_x=pdeintrp(node,triangle(1:4,:),J_X.'); % T on triangles
 Jt_y=pdeintrp(node,triangle(1:4,:),J_Y.'); % T on triangles
 uvet=[Jt_x;Jt_y];
else
 uvet=[];
end  
figure(hm)
subplot(122)

 pcor=8;
%for pcor=[8 15 25]
Vcor=pu;
ipcor=0;
for pcor=Vcor
 ipcor=ipcor+1; 
subplot(121)

 plot(x,Cn(pcor,:),'o')
         ax = gca;
       ax.ColorOrderIndex = 1
  hold on, plot(x,Cp(pcor,:),'.')
  legend(num2str(Cpl(pcor)',2))
  xlabel('\rho, \mum')
  ylabel('2D Injection, (cm^2 ns)^{-1}') 
  hold off
subplot(122)
 PlotStru
 quiver(SX,SY,squeeze(J_X(pcor,sY0,sX0)),squeeze(J_Y(pcor,sY0,sX0)),.2),
 hold on

 title(['Current = ',num2str(Cpl(Vcor(1:ipcor)),2),' (mA)'])
 ylim([352 359]), xlim([0 xMax]),  hold off, pausak

 
 iZoom=input(' iZoom ');
 if length(iZoom)==0
  iZoom=0;
 else
  iZoom=1;
 end
 if iZoom==1
 figure
 plot_tri_meshPLOT(geom,mesh,uvet,sd,color,triangle,node,grido,cbar,cmap,scale,vpath,arrows)
  hold on
  ze=y(end);
  yo=ze-mode.zox;
  xo=mode.rox;
  thvis=.1;
  xc=mode.Contact_i;
  yc=ze;
  wc=mode.Contact_e-xc;
  wo=mode.Contact_e-xo;
  RecCont=[xc yc wc thvis];
  Xc=[xc xc+wc xc+wc xc];
  Yc=[yc yc  yc+thvis yc+thvis];
  RecOx=[xo yo wo thvis];
  Xo=[xo xo+wo xo+wo xo]; 
  Yo=[yo yo  yo+thvis yo+thvis];
  hhco=patch(Xc,Yc,'y'); set(hhco,'EdgeColor','none');
  hhco=patch(Xqw,Yqw,'c'); set(hhco,'EdgeColor','none'); 
  ylim(yo+[-.5 .5]), xlim(xo+[-2 2])
  quiver(SX1,SY1,squeeze(J_X(pcor,sY1,sX1)),squeeze(J_Y(pcor,sY1,sX1)),.15),
  axis normal
  pausak
 end
end 


