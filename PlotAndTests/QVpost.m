in=1;
Cur_2D=0
SEGNO_elettroni=-1;
%[Jn_x,Jn_y,Jp_x,Jp_y]
J_XN=reshape(Jn_xR,mesh.nny,mesh.nnx);
J_YN=reshape(Jn_yR,mesh.nny,mesh.nnx);
J_XP=reshape(Jp_xR,mesh.nny,mesh.nnx);
J_YP=reshape(Jp_yR,mesh.nny,mesh.nnx);

JN_X=J_XN;
JN_Y=J_YN;
JP_X=J_XP;
JP_Y=J_YP;



%J_X=mode.JX;l
%J_Y=mode.JY;
y=modePlot.y;
x=modePlot.x;



SZ=.02

 X0=x;
 Y0=y;
[X,Y]=meshgrid(X0,Y0);
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
   sogx=SZ;
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
    sogx=SZ;
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
 
 if in==0
 geom=mode.geom;
 mesh.node=mode.node;
 mesh.triangle=mode.triangle;
 mesh.nn=mode.nn;
 mesh.nt=mode.nt;
 end
 
 xMax=10;
  node='off'; triangle='off'; color='on'; vpath='off'; arrows='off';
 sd=1:geom.nd; scale=0.000002; cmap=[]; grido='off'; cbar='off';
 
 iar=0
 if iar==1
 Jt_x=pdeintrp(node,triangle(1:4,:),J_XN.'); % T on triangles
 Jt_y=pdeintrp(node,triangle(1:4,:),J_YN.'); % T on triangles
 uvetN=[Jt_x;Jt_y];
 Jt_x=pdeintrp(node,triangle(1:4,:),J_XP.'); % T on triangles
 Jt_y=pdeintrp(node,triangle(1:4,:),J_YP.'); % T on triangles
 uvetP=[Jt_x;Jt_y]; 
else
 uvet=[];
end  


 pcor=8;
%for pcor=[8 15 25]
Vcor=pu;
ipcor=0;
for pcor=Vcor
 ipcor=ipcor+1; 
figure(hm)
subplot(121)

 plot(xcm*1e4,Cn(end,:),'o-',xcm*1e4,Cp(end,:),'.-')
  legend('Electrons','Holes')
  xlabel('\rho, \mum')
  ylabel('2D Injection, (cm^2 ns)^{-1}') 
  hold off
subplot(122)
 PlotStru
 quiver(SX,SY,SEGNO_elettroni*(J_XN(sY0,sX0)),SEGNO_elettroni*(J_YN(sY0,sX0)),.2),
 ylim([yo-2 yc]), xlim([0 xMax]),  
 pausak
 quiver(SX,SY,(J_XP(sY0,sX0)),(J_YP(sY0,sX0)),.2),
 hold on

 title(['Particle Currents @ ',num2str(1000*modePlot.ii_dd(end),2),', mA'])
 

   if Cur_2D==1
    xcm=x*1e-4;
    xDiv=(xcm+diff(xcm(1:2))/2)*2*pi;
    jnQ_x=sum(squeeze(mode.JnQW(pcor,:,:)),2)*1000;
    jpQ_x=sum(squeeze(mode.JpQW(pcor,:,:)),2)*1000;
    jQn=SEGNO_elettroni*[jnQ_x' zeros(1,length(xcm)-length(xQW))]./xDiv;
    jQp=[jpQ_x' zeros(1,length(xcm)-length(xQW))]./xDiv;
    

    yQWmedioE=mode.yMQW{3}*ones(size(sX0))*1e4;
    yQWmedioH=mode.yMQW{1}*ones(size(sX0))*1e4;
    quiver(x(sX0),yQWmedioE,jQn(sX0)*WQW,zeros(size(sX0)),.3,'g')
    quiver(x(sX0),yQWmedioH,jQp(sX0)*WQW,zeros(size(sX0)),.3,'m')
   end
   
    hold off, pausak
 
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
  ylim(yo+[-.8 1]), xlim(xo+[-2 2])
  quiver(SX1,SY1,SEGNO_elettroni*(J_XN(sY1,sX1)),SEGNO_elettroni*(J_YN(sY1,sX1)),.15),
  quiver(SX1,SY1,(J_XP(sY1,sX1)),(J_YP(sY1,sX1)),.15),
  axis normal
   if Cur_2D==1
    yQWmedioE=mode.yMQW{3}*ones(size(sX1))*1e4;
    yQWmedioH=mode.yMQW{1}*ones(size(sX1))*1e4;
    quiver(x(sX1),yQWmedioE,jQn(sX1)*WQW,zeros(size(sX1)),.3,'g')
    quiver(x(sX1),yQWmedioH,jQp(sX1)*WQW,zeros(size(sX1)),.3,'m')
   end 
    hold off, pausak
  pausak
 end
end 


