load stopDD
%istop
istopi=input(' input istop = ? ');
%istopi=1;
if length(istopi)==1

 istop=istopi;
 save stopDD istop    -v4
else
 'toggle istop'
 if istop==0
  istop=1
 else
  istop=0
 end 
  save stopDD istop    -v4
end