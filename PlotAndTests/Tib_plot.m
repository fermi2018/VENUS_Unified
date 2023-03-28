%portatori
 figure, semilogy(mesh.ygrid,mode.elec(1:mesh.nny),mesh.ygrid,mode.hole(1:mesh.nny))
 ylim([1e16 1e19])
 xlim(mesh.ygrid(end)*[.98 1])

meshPlot.node=mesh.node;
meshPlot.triangle=mesh.triangle;
meshPlot.nn=mesh.nn;
meshPlot.nt=mesh.nt;

node='off'; triangle='off'; color='on'; vpath='off'; arrows='off';
sd=1:geom.nd; scale=0.000002; cmap=[]; grido='off'; cbar='off';


Jt_x=pdeintrp(mesh.node,mesh.triangle(1:4,:),(mode.Jn_x+mode.Jp_x).'); % T on triangles
Jt_y=pdeintrp(mesh.node,mesh.triangle(1:4,:),(mode.Jn_y+mode.Jp_y).'); % T on triangles
uvet=[Jt_x;Jt_y];
figure,plot_tri_mesh(geom,mesh,uvet,sd,color,triangle,node,grido,cbar,cmap,scale,vpath,arrows)