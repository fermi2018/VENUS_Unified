function [uvet]=clipping(uvet,nn,pp,CarrierNorm)
Nmin=1e2; 
Pmin=1e2 ; 
phi=uvet(1:nn) ; 
n=uvet(nn+1:pp) ; 
p=uvet(pp+1:3*nn)  ;




indN=find(n<Nmin) ; 
indP=find(p<Pmin) ; 

n(indN)=Nmin/CarrierNorm ; 
p(indP)=Pmin/CarrierNorm ; 
resto=uvet(3*nn+1:end) ;

if indN 
   disp('clipped n!')  ; 
    
end
if indP
   disp('clipped p!') 
end

uvet=[phi  n  p resto] ; 





end