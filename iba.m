load iba
istop=iba;
istopi=input(' input istop = ? ');
%istopi=1;
if length(istopi)==1

 iba=istopi;
 save iba iba 
else
 'toggle iba'
 if istop==0
  iba=1
 else
  iba=0
 end 
  save iba iba  
end