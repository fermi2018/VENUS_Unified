XI=20;
XS=XI*.9;
n=pi/2;
lg=XI/6;
xP=0.4;

Dd=25;
Da=25;

iplot=1;

x0=2*XI*xP-XS;
xd0=(x0+XS)/(2*XI);
yd0=(n+atan(x0))/(n*2);
G0=exp(-((x0)/lg).^2);
Yd0=3+(1+xd0/8).*(yd0*30+Dd*G0);

Ya0=25+Da*xd0;


if iplot==1

x=-XI:.1:XI*1.2;
xd=(x+XS)/(2*XI);
yd=(n+atan(x))/(n*2);
G=exp(-((x)/lg).^2);
Yd=3+(1+xd/8).*(yd*30+G*Dd);
Ya=25+Da*xd;




h=figure, 
set(h,'pos',[1034         212         560         631])

plot(xd,Yd,'linewidth',2)
hold on, plot(xd0,Yd0,'bo')
plot(xd,Ya,'r','linewidth',2)
hold on, plot(xd0,Ya0,'ro')
grid
xlabel('x molar fraction')
ylabel('Act. En. (meV)')
xlim([0 1])
ylim([0 120])

end