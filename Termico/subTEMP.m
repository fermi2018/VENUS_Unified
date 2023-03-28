  icnu=0;
  for nua=nu_vett
   icnu=icnu+1;
   fTr=fTrr{icnu};
   fTf=fTrf{icnu};
   fT=fM{icnu};
   Mulr0=fTr*diag(xdx0);
   Mulf0=fTf*diag(xdx0);
   Mulx0=fT*diag(xdx0);
   Mulz0=f*diag(dz);
   Mulpz0=fp*diag(dz);
%   Mulx=repmat(Mulx0,L,1);
%   Mulz=repmat(Mulz0,Js(icnu),1);
   pJL=L*Js(icnu);
   ic=0;
   clear MLi
   for kx=1:Js(icnu)
    for kz=1:L
     ic=ic+1;    
     MLtr=Mulz0*diag(f(kz,:))*Fcondt*diag(fTr(kx,:))*Mulr0';
     MLtf=Mulz0*diag(f(kz,:))*Fcondt*diag(fTf(kx,:))*Mulf0';
     MLz=Mulpz0*diag(fp(kz,:))*Fcondz*diag(fT(kx,:))*Mulx0';
     MLi(ic,:)=reshape(MLtr+MLtf+MLz,pJL,1);     
    end
   end 
     MLm{icnu}=MLi;  
     
  %'quib', keyboard  
  end  %nua
    

 J=Jsu;
 pJL=J*L; 
 
  icnu=0;
 for nua=nu_vett
   icnu=icnu+1;
    ML=MLm{icnu};
    iML=inv(ML);
%   'qui inv', keyboard    
    if nua==0
     iMLt{icnu}=iML/(2*pi);
    else
     iMLt{icnu}=iML/pi;
    end

   end %azimut
   
   
   cijL=[];
   for inu=1:length(nu_vett)
    Ji=Js(inu);
    Jsh=sum(Js(1:inu-1));
    QL=reshape(P((1:Ji)+Jsh,:)',Ji*L,1);
    cijL=[cijL; iMLt{inu}*QL];
   end
   CoL=fattore_correttivo*reshape(cijL,L,J)';

   T3D1=Ru*CoL*S;   

