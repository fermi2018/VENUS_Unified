
'Inter', keyboard
ifp=-10;
Cs=cumsum(L_i);

x_in=ParDD0.Nd(puAddDBR,1);
x_fi=ParDD0.Nd(puAddDBR,2);

xt=[x_in x_fi];
Lt=[L_i 1e-10*ones(size(L_i))];

xv=reshape(xt',prod(size(xt)),1);
Lv=reshape(Lt',prod(size(xt)),1);

Cs=cumsum(L_i);
Csv=cumsum(Lv);

xfitto=linspace(Csv(1),Csv(end),10000);
XI=interp1(Csv,xv,xfitto,'linear');

figure, plot(Cs,x_in,'go',Cs,x_fi,'ro'), hold on, plot(xfitto,XI,'k.')
%plot(Csv,xv,'ro')
keyboard

[xi_out,xf_out,L_out]=SubDD_lithoGraded(x_in,x_fi,L_i,puTJ,pu_Below,ifp);
'RIN fuori graded',

keyboard
fi_ext=find(xi_out(:,2)<0);


x_in=ParDD0.Nd(puAddDBR,2);
[xf_out,L_out]=SubDD_litho(x_in,L_i,puTJ,pu_Below,ifp);
