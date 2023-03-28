close all
irea=0;
if irea==1
 d=xlsread('mob.xlsx','N','B2:D18');
 x=d(:,1);
 n=d(:,2);
 m=d(:,3);
 save ron x n m
 dh=xlsread('mob.xlsx','P','B2:D5');
 xh=dh(:,1);
 nh=dh(:,2);
 mh=dh(:,3);
 save roh xh nh mh 
else
 load ron
 load roh
end

xmol=x/100;


%xmol=linspace(0,1,101);
indReg2=find(xmol>0.425);

        mobnint = 8000-22000.*xmol+10000.*xmol.^2;
        mobnint(indReg2)= -255+1160*xmol(indReg2)-720*xmol(indReg2).^2;
        mobnint = 8000-24000.*xmol+13000.*xmol.^2;
        mobnint(indReg2)= 447*(xmol(indReg2)-.45).^2+148;
        mobnint(indReg2)= 1200*(xmol(indReg2)-.45).^2+148;
%        mobnint4= -255+1160*xmol-720*xmol.^2;
        %
        % Temperature dependence from 2005 Adachi, p. 325
%        mobnint=mobnint.*(300./T).^mesh.ExpE;
        %

%keyboard


N0_H=.1;
N0=n;
 
 mun=mobnint./(1+(N0/N0_H).^.35);


h=figure;
set(h,'pos',[ 159         471        1294         477])
subplot(131)
semilogy(x,m,'ro',x,mun,'gs'), 
grid
xlabel('molar fraction %')
ylabel('Mob. elec')
subplot(132)
loglog(n,m,'ro',n,mun,'gs'),
grid
xlabel('Free carries (1e18)'), 

  meInt=m.*(1+(N0/N0_H).^.35);
  figure, plot(x,mobnint,'go',x,meInt,'ro'), pausak

[xs,is]=sort(x);
xsd=[-1; xs];
fi=find(diff(xsd)>0);
xu=xs(fi);
for k=1:length(xu)
 xi=xu(k);
 fi=find(x==xi);
 mum(k)=mean(meInt(fi));
 muma(k)=max(meInt(fi));
end

  figure, plot(x,mobnint,'go',xu,mum,'ro',xu,muma,'rd',x,meInt,'r.'), pausak


N0_H=1;
N0_H=.2;
N0=nh;
x=xh/100;
xmol=x;

%xmol=linspace(0,1,101);
        % Hole low-field mobility
         mobpint1=400-700*xmol+450*xmol.^2;  %Calciati?  Sentaurus
         mobpint2=370-970*xmol+740*xmol.^2;  % Joffe
         mobpint3=400-775*xmol+535*xmol.^2;  % Joffe         
         mobpint=mobpint3;
 
 muh=mobpint./(1+(N0/N0_H).^.35);
 
 mhInt=mh.*(1+(N0/N0_H).^.35);
 %co=polyfit(x,mhInt,2);
 %va=polyval(co,x);
  figure, plot(xmol,mobpint,x,mhInt,'ro'), keyboard
 %figure

figure(h) 
subplot(133)
 semilogy(x,mh,'ro',x,muh,'gs'), 
grid
 xlabel('molar fraction %')
ylabel('Mob. holes')


  
  
