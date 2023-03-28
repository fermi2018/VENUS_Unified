JN_X=mode.JXn;
JN_Y=abs(mode.JYn);
JP_X=mode.JXp;
JP_Y=abs(mode.JYp);




z=y;
ycm=z*1e-4;
xcm=x*1e-4;

ipcor=0;
hc=figure;
set(hc,'pos',[349         398        1525         525 ])
subplot(121)
xqw=1000*(modePlot.y-modePlot.yQW);
 plot(xqw,modePlot.Ec(fipst,:),'linewidth',2)
 %plot(1000*(modePlot.y-modePlot.yQW),modePlot.Ec(fipst,:)-modePlot.Ev(fipst,:),'linewidth',2)
  ax = gca;
  ax.ColorOrderIndex = 1;
 hold on, plot(xqw,modePlot.Ev(fipst,:)+1.8,'--','linewidth',2)
 
 %figure, plot(y,Rc-Rv)
 dX=80;
 xlim([-dX,dX])
 title(' Ec (cont), Ev+1.8 (dashed) in Cavity')
 xlabel(' z around QW (nm) ')
grid

for pcor=Vcor
 ipcor=ipcor+1; 
for indy=1:length(z)
    
    jn_y=squeeze(JN_Y(pcor,indy,:))';
    jp_y=squeeze(JP_Y(pcor,indy,:))';
%    figure, plot(x,jn_y,x,jp_y), pausak
    curr_n(indy)=trapz(xcm,2.*pi.*xcm.*jn_y);
    curr_p(indy)=trapz(xcm,2.*pi.*xcm.*jp_y);
    
end



iplolea=1;

NQW=mode.NMQW;
yMQW=mode.yMQW;
vWMQW=mode.vWMQW;

Zl=[];
Ind=[];
dX=50e-7;
for kn=1:NQW
[val,ind]=min(abs(ycm-(yMQW{kn}-vWMQW{kn}/2)));
Zl=[Zl z(ind)];
Ind=[Ind ind];
[val,ind]=min(abs(ycm-(yMQW{kn}+vWMQW{kn}/2)));
ind=ind+1;
Zl=[Zl z(ind)];
Ind=[Ind ind];
zIn=(yMQW{3}-dX)*1e4;
zFi=(yMQW{1}+dX)*1e4;
zLim=[zIn zFi];
end
fP=((curr_n)+(curr_p))*1000;
figure(hc)
subplot(122)
%semilogy(z,curr_n*1000,z,curr_p*1000), hold on, 
%semilogy(Zl,curr_n(Ind)*1000,'o',Zl,curr_p(Ind)*1000,'o'),
%semilogy(z,fP,'g','linewidth',2), 
semilogy(xqw,curr_n*1000,xqw,curr_p*1000), hold on, 
semilogy(Zl,curr_n(Ind)*1000,'o',Zl,curr_p(Ind)*1000,'o'),
semilogy(z,fP,'g','linewidth',2), 
xlim(zLim)
pausak
grid

end

