
for kr=1:nrig
 R=real(Rag(kr,:));
 fiR=find(R>0);
 R=R(fiR);
  for kc=1:length(ra_dd)-1
   fiC=find(R==ra_dd(kc));
   if length(fiC)>0
    material{kr,kc+1}=materialO{kr,fiC+1};
    fiCp=fiC;
    fiC0=kc+1;
   end 
  end
  if length(fiR)==1
   for kc=fiC0+1:ncol
    material{kr,kc}=material{kr,fiC0};
   end  
  end 
end 