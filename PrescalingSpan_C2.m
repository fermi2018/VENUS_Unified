        C2vet=logspace(0,22,150); 
        C1=1e10; 
        
        clear condiz_prima
        clear condiz_dopo
        disp('Reconditioning: setting C3 ')
        
        
for indC=1:length(C2vet)
    
        Ctry=C2vet(indC); 
        C2=1; 
  
        
        SMdiag=sparse(1,neq); 
        SMdiag(1:nn)=1; 
        SMdiag(nn+1:pp)=Ctry; 
        SMdiag(pp+1:vv)=Ctry*C2; 
        SMdiag(vv+1:ww)=Ctry; 
        SMdiag(ww+1:ss)=Ctry; 
        SMdiag(ss+1)=1; 
        SMdiag(rr+1)=1; 
        SMdiag(qq+1)=1; 
        
        
        SM=sparse(1:neq,1:neq,SMdiag,neq,neq); 
        SMinv=sparse(1:neq,1:neq,1./SMdiag,neq,neq);
        
        
        
        
        Jmat_try=SMinv*Jmat*SM;  
        condiz_prima(indC)=condest(Jmat_try); 
        
         [R,C]=dgsequ(Jmat_try); 
         condiz_dopo(indC)=condest(R*Jmat_try*C); 
         
%         [P1,R1,C1]=equilibrate(Jmat_try); 
%         condiz_dopo_equilibrate(indC)=condest(R1*P1*Jmat_try*C1); 
     
end

disp('done')
strsave=[num2str(indv),'_2TJ_C2']; 
save(strsave,'condiz_dopo','condiz_prima','C2vet','C1'); 
