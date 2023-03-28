s_LoadConstants


nv=[.1 1 10]*1e18;
T=300;
lam=850e-9;  %in cm
N0=1e18;


mesh.ExpE=1;

xmol=0:.02:1;

nref=sqrt(13-4.*xmol);

indReg2=find(xmol>0.45);

        mobnint = 8000-22000.*xmol+10000.*xmol.^2;
        mobnint(indReg2)= -255+1160*xmol(indReg2)-720*xmol(indReg2).^2;
        %
        % Temperature dependence from 2005 Adachi, p. 325
        mobnint=mobnint.*(300./T).^mesh.ExpE;
        macro.mobnint=mobnint;
        %
        % Hole low-field mobility
         mobpint1=400-700*xmol+450*xmol.^2;  %Calciati?  Sentaurus
         mobpint2=370-970*xmol+740*xmol.^2;  % Joffe
         mobpint=(mobpint1+mobpint2)/2;
         mobpint=mobpint2;
         
	Eg0_G=1.519+1.155.*xmol+0.370.*xmol.^2;
        alpha_G=5.41e-4; beta_G=204;
        Eg_G=Eg0_G-(alpha_G.*T.^2)./(beta_G+T);
        %
        Eg0_X=1.981+0.124.*xmol+0.144.*xmol.^2;
        alpha_X=4.6e-4; beta_X=204;
        Eg_X=Eg0_X-(alpha_X.*T.^2)./(beta_X+T);
        %
        Eg0_L=1.815+0.69.*xmol; % mind the error in Ioffe!
        alpha_L=6.05e-4; beta_L=204;
        Eg_L=Eg0_L-(alpha_L.*T.^2)./(beta_L+T);
        %
        Eg = min(Eg_G,Eg_X);
        %
        % saving bandgaps
        macro.Eg_G = Eg_G;
        macro.Eg_X = Eg_X;
        macro.Eg_L = Eg_L;
        macro.Eg = Eg;
        %
        %-- Effective DoS masses (1985Adachi_JAP, Table II)
        meffn_G = 0.067+0.083.*xmol;
        meffn_X = 0.850-0.140.*xmol;
        meffn_L = 0.560+0.100.*xmol;
        %
        % DoS
        Nc_G = 2*(2*pi*kB*T*m0).^(3/2)/(h.^3).*1e-6.*meffn_G.^(3/2);
        Nc_X = 2*(2*pi*kB*T*m0).^(3/2)/(h.^3).*1e-6.*meffn_X.^(3/2);
        Nc_L = 2*(2*pi*kB*T*m0).^(3/2)/(h.^3).*1e-6.*meffn_L.^(3/2);
        %
        % Equivalent DoS
        Vt=kB.*T./qel;
        DeltaE_XG = Eg_X-Eg_G;
        DeltaE_LG = Eg_L-Eg_G;
        Nc_1 = Nc_G.*(1+Nc_X./Nc_G.*exp(-DeltaE_XG./Vt)+Nc_L./Nc_G.*exp(-DeltaE_LG./Vt));
        DeltaE_GX = Eg_G-Eg_X;
        DeltaE_LX = Eg_L-Eg_X;
        Nc_2 = Nc_X.*(1+Nc_G./Nc_X.*exp(-DeltaE_GX./Vt)+Nc_L./Nc_X.*exp(-DeltaE_LX./Vt));
        %
        indX=find(Eg_G>Eg_X); % index where band X is dominant
        Nc = Nc_1;
        Nc(indX) = Nc_2(indX);
        meffn = (Nc./(2.*(2.*pi.*kB.*T.*m0).^(3/2)/(h.^3).*1e-6)).^(2/3);         
        meffp=0.55+0.26*xmol; % effective hole mass
        Nv=2*(2.*pi.*kB.*T.*m0).^(3/2)/(h.^3)*1e-6.*meffp.^(3/2); 
        meffp = (Nv./(2.*(2.*pi.*kB.*T.*m0).^(3/2)/(h.^3).*1e-6)).^(2/3);        
        
 figure, plot(xmol,meffn,xmol,meffp)   
 xlabel(' molar fraction')
 ylabel(' Eff. mass')
 , pausak

 
 N0_H=1e17;
 
 mun=mobnint./(1+(N0/N0_H).^.35);
 mup=mobpint./(1+(N0/(2*N0_H)).^.35);
 
  figure, semilogy(xmol,mun,xmol,mup)       
   xlabel(' molar fraction')
  ylabel(' Eff. mob.')
 , pausak
 
 Clight=3e8;
 
 
 
 C_Drude=qel^3*lam^2./(4*pi^2*nref*eps0*Clight^3*m0^2);
 
 keyboard
 
 for k=1:length(nv)
 
  N0=nv(k);
  mun=1e-4*mobnint./(1+(N0/N0_H).^.35);
  mup=1e-4*mobpint./(1+(N0/(2*N0_H)).^.35);
  N0=N0*1e6;
  Al_DrudeN(k,:)=C_Drude*N0./(meffn.^2.*mun);
  Al_DrudeP(k,:)=C_Drude*N0./(meffp.^2.*mup);

end  

figure, 
semilogy(xmol,Al_DrudeP/100), hold on

ax = gca;
ax.ColorOrderIndex = 1;
semilogy(xmol,Al_DrudeN/100,'--'),
