function [xi_out,xf_out,Lfinal]=SubDD_lithoGraded(x_in,x_fi,L_in,puTJ,pu_Below,ifp)

ifi=0;
if ifp==-10
ifi=1;
end

L_i=(L_in);
%L_i=abs(L_in);

Detch=sum(L_i(puTJ));

if ifi==1
Cs=cumsum(L_i);
figure, plot(Cs,x_in,'.-')


hold on

D=sum(L_i(1:puTJ(1)-1));
%plot(cumsum(L_i(puTJ))+D,x_in(puTJ),'r','linewidth',2)
plot(Cs(puTJ),x_in(puTJ),'or','linewidth',2)
pausak
end

pui=[1:puTJ(1)-1];
puiTJ=[1:puTJ(end)];


xt=flipud([x_fi(puiTJ) x_in(puiTJ)]);
Lt=flipud([1e-10*ones(size(puiTJ')) L_i(puiTJ)]);
Lt(1)=0;

DD=cumsum([0;L_i])-Detch;
L_is=diff(DD);

xv=reshape(xt',prod(size(xt)),1);
Lv=reshape(Lt',prod(size(xt)),1);

Csv=cumsum(Lv);

Le=[Detch; L_i(pui)];


Ltop=flipud(Le);

puiTJ=[1:puTJ(end)];
LtopTJ=flipud(L_i(puiTJ));
ntopTJ=flipud(x_in(puiTJ));
Les=cumsum(LtopTJ);



Dtot=[0; cumsum(LtopTJ); cumsum(Ltop)];
[Ds,iso]=sort(Dtot);

dL=diff([-10; Ds]);
fiVeri=find(dL>0);


Dfinal=Ds(fiVeri);
Lfinal=diff(Dfinal);


if ifp==-10
'VEDI FIT', keyboard
end
Di=[Dfinal(1:end-1); Dfinal(end)+1e-3];
Df=[Dfinal(2)-1e-2;Dfinal(2:end)];
XIi=interp1(Csv,xv,Di);
XFi=interp1(Csv,xv,Df);
XIe=interp1(Csv-Detch,xv,Di);
%'LLLL', keyboard
XIe(find(isnan(XIe)==1))=0;
XFe=interp1(Csv-Detch,xv,Df);
XFe(find(isnan(XFe)==1))=0;
%XI=[XI; XF(end-1)];
%XF=[XI(1); XF];
xi_outd=[XIi XIe];
xf_outd=[XFi XFe];

xi_out=[flipud(xi_outd(2:end,:)); repmat(x_in(pu_Below,:),1,2)];
xf_out=[flipud(xf_outd(2:end,:)); repmat(x_fi(pu_Below,:),1,2)];
%'VEDI FIT', keyboard

ifig=0;

if ifig==0
return
end
Dm=Dfinal;
figure, plot(Dm,XIi,'.',Dm,XIe,'+')
%hold on, plot(Dm,XIe,'o',Dm,XFe,'s')
hold on, plot(Csv,xv)


keyboard
'qui; GRADED', keyboard

figure, plot(Dm,XIi,'.',Dm-Detch,XIe,'o')



%'qui; GRADED', keyboard

