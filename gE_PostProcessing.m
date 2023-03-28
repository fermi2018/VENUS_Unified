
%clear all
%clc
%%close all 

%N=load('OX_Nodi_ThFake.mat') ;
s_LoadConstants
LutName='C:\Users\MDale\Downloads\VenusTJ_lithoJan23_1002\VenusTJ_lithoJan23\dati\LUT4D_Jun_Markus_nMark_40_Der.mat' ;
LUT=load(LutName) ;
nindex=3.6000;
vph=Clight/N.mode.nindexQW ;
nnQW=1:51 ;


for indQW=1:3
    
    inQW = N.mesh.inMQW{indQW} ;
    xQW = N.mesh.node(1,inQW) ;
    LeQW = xQW(2:end) - xQW(1:end-1); % edge length, cm
    Lp = zeros(1,51); % box length, cm
    xcQW = (xQW(2:end) + xQW(1:end-1))/2; % edge centers, cm
    iiQW1=1:50 ;
    iiQW2=2:51 ;
    Lp(1:(51-1)) = LeQW/2; Lp(2:51) = Lp(2:51) + LeQW/2;
    Lp1=Lp(iiQW1)/2; Lp1(1)=2*Lp1(1);
    Lp2=Lp(iiQW2)/2; Lp2(end)=2*Lp2(end);
    Lp1=2.*pi.*Lp1.*xcQW;
    Lp2=2.*pi.*Lp2.*xcQW;
 
    for iV=10:40
        
        nQWN=N.modePlot.nQW{iV}{indQW} ;
        pQWN=N.modePlot.pQW{iV}{indQW} ;
        Temp=reshape(N.modePlot.Temp(iV,:,:),[],1) ;
        TQWN = Temp(inQW) ;
        E1=squeeze(N.modePlot.E2(iV,1,nnQW))' ;
        E2=squeeze(N.modePlot.E2(iV,2,nnQW))' ;
        E3=squeeze(N.modePlot.E2(iV,3,nnQW))' ;
        
        
        n=nQWN ;
        p=pQWN ;
        T=TQWN ;
        
        for iLambda=1:3
            lambdaN=N.modePlot.lambda(iLambda,iV) ;
            lambda=ones(1,length(nQWN))*lambdaN*1e-3 ;
            g(iLambda,:)=EvalGains_interp(LUT,n,p,lambda+N.mode.Deltalam*1e-3,T');  % gain Computation
           
            
            
        end
        
           % NON-LINEAR GAIN COMPRESSION
            UU = (N.modePlot.Pst_dd(:,iV)'.*N.modePlot.fPdif(iV,:))*squeeze(N.modePlot.E2(iV,:,:));
            uudu=1./(1+N.mode.epsNLg*UU);
            nlG=uudu(1:51);
        
        
        %g=g*vph ;
        
                gE1Mat(iV,1,:)=reshape(g(1,:).*E1,1,[]) ;
                gE1=reshape(g(1,:).*E1,1,[]).*nlG ;
                gE2Mat(iV,1,:)=reshape(g(2,:).*E2,1,[]) ;
                gE2=reshape(g(2,:).*E2,1,[]).*nlG ;
                gE3Mat(iV,1,:)=reshape(g(3,:).*E3,1,[]) ;
                gE3=reshape(g(3,:).*E3,1,[]).*nlG ;
        %
        %         Gmod1(iV,indQW)=trapz(xQW,xQW.*gE1*2*pi)*N.mode.Gamma_z(1,indQW) ;
        %         Gmod2(iV,indQW)=trapz(xQW,xQW.*gE2*2*pi)*N.mode.Gamma_z(2,indQW) ;
        %         Gmod3(iV,indQW)=trapz(xQW,xQW.*gE3*2*pi)*N.mode.Gamma_z(3,indQW) ;
        
        
        Gmod1(iV,indQW) = sum([Lp1.*gE1(iiQW1) Lp2.*gE1(iiQW2)]);
        Gmod1(iV,indQW)=N.mode.Gamma_z(1,indQW)*Gmod1(iV,indQW) ;
        Gmod2(iV,indQW) = sum([Lp1.*gE2(iiQW1) Lp2.*gE2(iiQW2)]);
        Gmod2(iV,indQW)=N.mode.Gamma_z(2,indQW)*Gmod2(iV,indQW) ;
        Gmod3(iV,indQW) = sum([Lp1.*gE3(iiQW1) Lp2.*gE3(iiQW2)]);
        Gmod3(iV,indQW)=N.mode.Gamma_z(3,indQW)*Gmod3(iV,indQW) ;
        
        
    end


end

figure

hold on, plot(N.mode.ii_dd(10:40)*1e3,sum(Gmod1(10:40,:),2),'r')
hold on, plot(N.modePlot.ii_dd*1e3,N.modePlot.Gmod(1,:),'b')
hold on, plot(N.mode.ii_dd(10:40)*1e3,sum(Gmod2(10:40,:),2),'r')
hold on, plot(N.modePlot.ii_dd*1e3,N.modePlot.Gmod(2,:),'b')
hold on, plot(N.mode.ii_dd(10:40)*1e3,sum(Gmod3(10:40,:),2),'r')
hold on, plot(N.modePlot.ii_dd*1e3,N.modePlot.Gmod(3,:),'b')
grid on
xlabel('Current')
legend('Post Processing','Saved')


vph=Clight/N.mode.nindexQW ;
Gamma_zMean=  0.0101  ; 

for iV=10:40
   
        E1=squeeze(N.modePlot.E2(iV,1,1:51))' ;
        E2=squeeze(N.modePlot.E2(iV,2,1:51))' ;
        E3=squeeze(N.modePlot.E2(iV,3,1:51))' ;
        
        gE1=N.modePlot.matgain(iV,:).*E1.*Gamma_zMean ;
        gE2=N.modePlot.matgain(iV,:).*E2.*Gamma_zMean ;
        gE3=N.modePlot.matgain(iV,:).*E3.*Gamma_zMean ;
        
        GM1(iV)=sum([Lp1.*gE1(iiQW1) Lp2.*gE1(iiQW2)]) ;
        GM2(iV)=sum([Lp1.*gE2(iiQW1) Lp2.*gE2(iiQW2)]) ;
        GM3(iV)=sum([Lp1.*gE3(iiQW1) Lp2.*gE3(iiQW2)]) ;
    
end

figure, plot(N.mode.ii_dd(10:39)*1e3,GM1(10:39)*3,'r')
hold on, plot(N.modePlot.ii_dd*1e3,N.modePlot.Gmod(1,:),'b')
hold on, plot(N.mode.ii_dd(10:39)*1e3,GM2(10:39)*3,'r')
hold on, plot(N.modePlot.ii_dd*1e3,N.modePlot.Gmod(2,:),'b')
hold on, plot(N.mode.ii_dd(10:39)*1e3,GM3(10:39)*3,'r')
hold on, plot(N.modePlot.ii_dd*1e3,N.modePlot.Gmod(3,:),'b')
grid on
xlabel('Current')
legend('Post Processing','Saved')

s_LoadConstants
vph=Clight/N.mode.nindexQW ;
Gmod1=zeros(1,40) ; 
Gmod2=zeros(1,40) ; 
Gmod3=zeros(1,40) ; 

for iV=1:40
        E1=squeeze(N.modePlot.E2(iV,1,1:51))' ;
        E2=squeeze(N.modePlot.E2(iV,2,1:51))' ;
        E3=squeeze(N.modePlot.E2(iV,3,1:51))' ;
        
    for indQW=1:3
        
        
        inQW = N.mesh.inMQW{indQW} ;
        xQW = N.mesh.node(1,inQW) ;
        LeQW = xQW(2:end) - xQW(1:end-1); % edge length, cm
        Lp = zeros(1,51); % box length, cm
        xcQW = (xQW(2:end) + xQW(1:end-1))/2; % edge centers, cm
        iiQW1=1:50 ;
        iiQW2=2:51 ;
        Lp(1:(51-1)) = LeQW/2; Lp(2:51) = Lp(2:51) + LeQW/2;
        Lp1=Lp(iiQW1)/2; Lp1(1)=2*Lp1(1);
        Lp2=Lp(iiQW2)/2; Lp2(end)=2*Lp2(end);
        Lp1=2.*pi.*Lp1.*xcQW;
        Lp2=2.*pi.*Lp2.*xcQW;
        
        
%         UU = (N.modePlot.Pst_dd(:,iV)'.*N.modePlot.fPdif(iV,:))*squeeze(N.modePlot.E2(iV,:,:));
%         uudu=1./(1+N.mode.epsNLg*UU);
%         nlG=uudu(1:51);
%         nlG=ones(1,51) ; 
        gE1=squeeze(N.mode.gE(iV,indQW,1,:))'.*N.mode.Gamma_z(1,indQW)/vph ; 
        gE2=squeeze(N.mode.gE(iV,indQW,2,:))'.*N.mode.Gamma_z(2,indQW)/vph ; 
        gE3=squeeze(N.mode.gE(iV,indQW,3,:))'.*N.mode.Gamma_z(3,indQW)/vph ; 
        Gmod1(iV)=Gmod1(iV)+sum([gE1(iiQW1).*Lp1 gE1(iiQW2).*Lp2]) ; 
        Gmod2(iV)=Gmod2(iV)+sum([gE2(iiQW1).*Lp1 gE2(iiQW2).*Lp2]) ; 
        Gmod3(iV)=Gmod3(iV)+sum([gE3(iiQW1).*Lp1 gE3(iiQW2).*Lp2]) ; 
        % NON-LINEAR GAIN COMPRESSION
       
 
    end
    
end

figure, plot(N.mode.ii_dd(10:39)*1e3,Gmod1(10:39),'r')
hold on, plot(N.modePlot.ii_dd*1e3,N.modePlot.Gmod(1,:),'b')
figure, plot(N.mode.ii_dd(10:39)*1e3,Gmod2(10:39),'r')
hold on, plot(N.modePlot.ii_dd*1e3,N.modePlot.Gmod(2,:),'b')
figure, plot(N.mode.ii_dd(10:39)*1e3,Gmod3(10:39),'r')
hold on, plot(N.modePlot.ii_dd*1e3,N.modePlot.Gmod(3,:),'b')
grid on
xlabel('Current')
legend('Post Processing','Saved')

hold on ,plot(N.modePlot.ii_dd*1e3,N.modePlot.Lmod,'black o')


%% plot g-E 

gE1=zeros(40,51) ; 
gE2=zeros(40,51) ; 
gE3=zeros(40,51) ; 

for iV=1:40
    
        E1=squeeze(N.modePlot.E2(iV,1,1:51))' ;
        E2=squeeze(N.modePlot.E2(iV,2,1:51))' ;
        E3=squeeze(N.modePlot.E2(iV,3,1:51))' ;
        gE1(iV,:)=N.modePlot.matgain(iV,:).*E1.*Gamma_zMean ;
        gE2(iV,:)=N.modePlot.matgain(iV,:).*E2.*Gamma_zMean ;
        gE3(iV,:)=N.modePlot.matgain(iV,:).*E3.*Gamma_zMean ;
        g(iV,:) = N.modePlot.matgain(iV,:)  ;
    
end

I=N.mode.ii_dd*1e3 ; 
Elementi=N.mode.Elementi ; 
V=N.mode.vv_dd ;
P=N.mode.Pst_dd  ; 

save('gE_E','I','P','V','gE1','gE2','gE3','xQW','g') ; 


return
clear all 
load('gE_N') ; 
figure,subplot(2,1,1),plot(xQW*1e4,gE1(40,:),'b')
title('Red: Elements, Blue: nodes')
hold on,plot(xQW*1e4,gE2(15,:),'b')
hold on,plot(xQW*1e4,gE3(15,:),'b')
xlabel('\rho, \mum')
ylabel('g-E prodcut')
grid on 
xlim([0 3.5])
subplot(2,1,2),plot(xQW*1e4,g(40,:),'b')
xlabel('\rho, \mum')
ylabel('gain profile')
grid on
xlim([0 5])
EE=load('gE_E') ; 
subplot(2,1,1),hold on, plot(xQW*1e4,EE.gE1(40,:),'r')
hold on,plot(xQW*1e4,EE.gE2(15,:),'r')
hold on,plot(xQW*1e4,EE.gE3(15,:),'r')
xlabel('\rho, \mum')
ylabel('g-E prodcut')
grid on 
xlim([0 3.5])
subplot(2,1,2),hold on, plot(xQW*1e4,EE.g(40,:),'r')
xlabel('\rho, \mum')
ylabel('gain profile')
grid on
xlim([0 5])
