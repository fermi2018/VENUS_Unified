function VtMean =  f_ComputeMean(Vt,rho,z)

nz=size(z,2) ;
nrho=size(rho,2) ; 
Vt_z=zeros(1,nz) ;

RhoMax=rho(end) ; 
Zmax=max(z) ; 

VtMat=reshape(Vt,nz,nrho) ; 

for i=1:nz
Vt_z(i)=trapz(rho,VtMat(i,:)*2*pi.*rho)  ; %integration along rho
end

VtMean=trapz(z,Vt_z)/(pi*RhoMax^2*Zmax) ;

end