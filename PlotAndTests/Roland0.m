close all



xmol=linspace(0,1,101);
indReg2=find(xmol>0.425);

        mobnint1 = 8000-22000.*xmol+10000.*xmol.^2;
        mobnint1(indReg2)= -255+1160*xmol(indReg2)-720*xmol(indReg2).^2;
        mobnint = 8000-24000.*xmol+13000.*xmol.^2;
        mobnint(indReg2)= 1200*(xmol(indReg2)-.45).^2+148;
%        mobnint4= -255+1160*xmol-720*xmol.^2;
        %
        % Temperature dependence from 2005 Adachi, p. 325
%        mobnint=mobnint.*(300./T).^mesh.ExpE;
        %

%keyboard



  figure, plot(xmol,mobnint,xmol,mobnint1), 
  xlabel('molar fraction')
  ylabel('Elec. Mob')
  legend('New','Old')
  pausak




        % Hole low-field mobility
         mobpint1=400-700*xmol+450*xmol.^2;  %Calciati?  Sentaurus
         mobpint2=370-970*xmol+740*xmol.^2;  % Joffe
         mobpint3=400-775*xmol+535*xmol.^2;  % Roland       
         mobpint=mobpint3;
 

 %co=polyfit(x,mhInt,2);
 %va=polyval(co,x);
  figure, plot(xmol,mobpint,xmol,mobpint2),
  grid
    xlabel('molar fraction')
    ylabel('Hole Mob')
  legend('New','Old')
 %figure

