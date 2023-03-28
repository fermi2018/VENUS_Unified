plot_tri_meshPLOT(geom,mesh,uvet,sd,color,triangle,node,grido,cbar,cmap,scale,vpath,arrows)
 hold on
 ze=y(end);
 yo=ze-modePlot.zox;
 xo=modePlot.rox;
 thvis=.1;
 xc=modePlot.Contact_i;
 yc=ze;
 wc=modePlot.Contact_e-xc;
 wo=modePlot.Contact_e-xo;
 RecCont=[xc yc wc thvis];
 Xc=[xc xc+wc xc+wc xc];
 Yc=[yc yc  yc+thvis yc+thvis];
 RecOx=[xo yo wo thvis];
 Xo=[xo xo+wo xo+wo xo]; 
 Yo=[yo yo  yo+thvis yo+thvis];
 hhco=patch(Xc,Yc,'y'); set(hhco,'EdgeColor','none');
 %hhco=patch(Xo,Yo,'c'); set(hhco,'EdgeColor','none');
   yMQW=modePlot.yMQW{end}*1e4;
   Xqw=[0 xMax xMax 0 ]; 
   thQ=((modePlot.yMQW{1}-modePlot.yMQW{end})+modePlot.vWMQW{1})*1e4;
   Yqw=[yMQW yMQW  yMQW+thQ yMQW+thQ];
  hhco=patch(Xqw,Yqw,'c'); set(hhco,'EdgeColor','none'); 