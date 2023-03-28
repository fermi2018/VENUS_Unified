iPloTer=0;



if iPloTer==1
zp=[z(1) z z(end)];
figure, 
subplot(121)
plot(zp,Tdd,'.')
xlim([zp(end)-10 zp(end)])
subplot(122)
plot(rho,Tdd,'.')
end
% calcolo contributi
qqsave=qq;

clear T3D2

   iRAD=0;   
   qq=(Joule);
   Calcola_TEMP
   T3D2{1}=Ru*CoL*S;   
   Q{1}=qtot;
   
   qq=Rec_srhAu; 
   Calcola_TEMP
   T3D2{2}=Ru*CoL*S;   
   Q{2}=qtot;   
   
   qq=Rec_Cap; 
   Calcola_TEMP
   T3D2{3}=Ru*CoL*S;  
   Q{3}=qtot;   

   qq=Rec_RAD; 
   iRAD=iRAD_spalmato;   
   %iRAD=0;
   Calcola_TEMP
   T3D2{4}=Ru*CoL*S;   
   Q{4}=qtot;   

   iRAD=0;
   qq=OptAbs; 
   Calcola_TEMP
   T3D2{5}=Ru*CoL*S;   
   Q{5}=qtot;   

   qq=Tho; 
   Calcola_TEMP
   T3D2{6}=Ru*CoL*S;      
   Q{6}=qtot;   


if iPloTer==1
	Leg{1}='Joule';
	Leg{2}='SHR + Auger';
	Leg{3}='Ccap';
	Leg{4}='Rad';
	Leg{5}='Opt Abs';
	Leg{6}='Thomson';
   
   ff=figure;
   set(ff,'pos',[  151    10   823   974])
   for kp=1:6
   subplot(3,2,kp)
   plot(xxc',T3D2{kp}), 
   xlim([0 15])
   title(Leg{kp})
   end
   pausak
   
   ff=figure;
   set(ff,'pos',[362    10   860   974])
   for kp=1:6
   subplot(3,2,kp)
   plot(zzc',T3D2{kp}), 
   xlim([zp(end)-10 zp(end)])
   title(Leg{kp})
   end   
   pausak
   
   ff=figure;
   set(ff,'pos',[  551    10   823   974])
   for kp=1:6
   subplot(3,2,kp)
   plot(xxc',Q{kp}), 
   xlim([0 15])
   title(Leg{kp})
   end
   pausak
   
   ff=figure;
   set(ff,'pos',[662    10   860   974])
   for kp=1:6
   subplot(3,2,kp)
   plot(zzc,Q{kp}'), 
   xlim([zp(end)-10 zp(end)])
   title(Leg{kp})
   end   
   pausak   
end   


T_Contributi=T3D2; 
clear T3D2
qq=qqsave;