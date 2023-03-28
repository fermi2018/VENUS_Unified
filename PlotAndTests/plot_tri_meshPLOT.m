function plot_tri_meshPLOT(geom,mesh,uvet,sd,color,triangle,node,grid,cbar,cmap,scale,vpath,arrows)

if(nargin==3), 
node='off'; triangle='off'; color='off'; vpath='off'; arrows='off';
sd=1:geom.nd; scale=[]; cmap='jet'; grid='off'; cbar='on';
end   
    
cla reset, set(gca,'NextPlot','add')

if(strcmp(grid,'on')), str='k'; else, str='none'; end
if(strcmp(color,'on')), color_patch=1; else, color_patch=0; end

% mesh plot 
if(not(isempty(mesh)))
x=mesh.node(1,:)*1e4; y=mesh.node(2,:)*1e4; 
%
in1=mesh.triangle(1,:); in2=mesh.triangle(2,:); in3=mesh.triangle(3,:);  
region=mesh.triangle(4,:); 
x1=x(in1); x2=x(in2); x3=x(in3);  
y1=y(in1); y2=y(in2); y3=y(in3); 
X=[x1; x2; x3]; Y=[y1; y2; y3];
nt=mesh.nt; nn=mesh.nn;
if(strcmp(color,'on'))
for n=1:geom.nd;
it=find(region==n);  
patch([x1(it); x2(it); x3(it)],[y1(it); y2(it); y3(it)], ...
   geom.color(n,:),'EdgeColor',str,'FaceAlpha',color_patch);
end, end, end

% matrix data
if(not(isempty(mesh)))
it=ismember(mesh.triangle(4,:),sd);
in1=mesh.triangle(1,it); in2=mesh.triangle(2,it); in3=mesh.triangle(3,it); 
x1=x(in1); x2=x(in2); x3=x(in3); y1=y(in1); y2=y(in2); y3=y(in3); 
X=[x1; x2; x3]; Y=[y1; y2; y3]; xc=(x1+x2+x3)/3; yc=(y1+y2+y3)/3;
end

% rectangles
if(strcmp(triangle,'on')), text(xc,yc,num2str([1:nt]'),'FontSize',8,'color','r'), end

% nodes
if(strcmp(node,'on')), text(x,y,num2str([1:nn]'),'FontSize',8), end 

% vpath
if(strcmp(vpath,'on'))
[vpath]=find_path(geom,mesh);
x=mesh.node(1,vpath); y=mesh.node(2,vpath);
plot(x,y,'m','linewidth',2)
end

% mark contact points
%ii=find(mesh.contact); plot(mesh.node(1,ii),mesh.node(2,ii),'k.')

% plot solution 
if(not(isempty(uvet))),   
%    
if(size(uvet)==[1 nn]),     
C=([uvet(in1); uvet(in2); uvet(in3)]);
patch(X,Y,C,'EdgeColor',str)
end
%
if(size(uvet)==[1 nt]), 
C=([uvet(it); uvet(it); uvet(it)]);    
patch(X,Y,C,'EdgeColor',str)
end
%
if(size(uvet)==[2 nt]),
it=ismember(mesh.triangle(4,:),sd);
uvet=uvet(:,it); uvet2=sqrt(abs(uvet(1,:)).^2+abs(uvet(2,:)).^2);   
C=[uvet2; uvet2; uvet2];
if(strcmp(color,'on')), patch(X,Y,C,'EdgeColor',str), end
if(strcmp(arrows,'on')),
% plot only arrows with nonzero length 
absuvet2=sqrt(real(uvet(1,:)).^2+real(uvet(2,:)).^2);    
ii=(absuvet2>max(absuvet2)/20); 
xc=xc(ii); yc=yc(ii); 
uvet=uvet(:,ii);
quiver(xc,yc,scale*real(uvet(1,:)),scale*real(uvet(2,:)),0,'b')
quiver(xc,yc,scale*imag(uvet(1,:)),scale*imag(uvet(2,:)),0,'r')
end, end
%
end
%
% geometry
% for bs=1:size(geom.dgm,2)
% x=geom.dgm(2:3,bs)*1e4; y=geom.dgm(4:5,bs)*1e4; hl=plot(x,y);
% if(ismember(bs,geom.bspec)), set(hl,'Color',[1 0 0],'linewidth',2)
% elseif(ismember(bs,geom.bspmc)), set(hl,'Color',[0 0 1],'linewidth',2), 
% else, set(hl,'Color',[0 0 0],'linewidth',0.5), end
% end     
%ii=find(mesh.edge(1,:));
%i1=mesh.edge(1,ii); i2=mesh.edge(2,ii);
%x1=x(i1); x2=x(i2); y1=y(i1); y2=y(i2);
%for n=1:length(x1);
%line([x1(n) x2(n)],[y1(n) y2(n)],'color',[0 0 0])
%end
%
if(strcmp(cbar,'on')), colorbar('vert'), end
if(not(isempty(cmap))), colormap(cmap), end
%
set(gca,'NextPlot','replace','FontSize',14,'FontName','Arial','box','on')
xlabel('\rho, \mum')
ylabel('z, \mum')
%    
% set(gca,'XLabel',text('String','x, \mum'), ...
%         'YLabel',text('String','y, \mucm'))
%     
axis equal, 
zoom on