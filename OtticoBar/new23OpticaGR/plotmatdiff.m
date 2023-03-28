Nkk=31
DIFFERE=0;
saKa=[];
for k=1:Nkk,
 if sum(sum(TD(:,:,k)))>1e-15
  DIFFERE=DIFFERE+1;
  saKa(DIFFERE)=k;
  map(TD(:,:,k)), title(num2str(k)), pausak,
 end
end 
'numero differenze'
DIFFERE
pausak
saKa
close all