load stop
istop
istopsa=istop;
istop=input(' input istop = ? ');
if length(istop)==1
  save stop istop -v4
else
 if istopsa==1
  istop=0
 else
  istop=1
 end
  save stop istop -v4
end
