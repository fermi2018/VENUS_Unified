T=linspace(0,300,1001);
T300=300;
Tp=[300 350 400];

Et0=-1.3;
for k=1:length(Tp)
 Et=Et0+1e-3*T;
 
 
 y(k,:)=((Tp(k)+T)/T300).^Et;
 y0(k,:)=((Tp(k)+T)/T300).^Et0;
end

figure, plot(T,y), hold on, 
ResetColor
plot(T,y0,'--')