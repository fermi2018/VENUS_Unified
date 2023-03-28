load sa
FM=fmlsi(ficri,2)';
firep=find(FM>1);

dfirep=diff(firep);
fdif=find(dfirep>1)+1;

 ' sono dentro  gam_crit scatt', keyboard
 
 puesp=ficri';
if length(firep)>0
 firepid=1:fdif(1)-1;
 firepi=firep(firepid);
 fire=puesp(firepi); 
 re=FM(fire);
 puesp=[puesp(1:firepi(1)-1) repmat(fire,1,re(1)) puesp(firepi(end)+1:end)];
 FM=[FM(1:firepi(1)-1) repmat(ones(size(fire)),1,re(1)) FM(firepi(end)+1:end)]; 
 firep=find(FM>1);
 'qui', keyboard
 while length(firep)>0
  dfirep=diff(firep);
  fdif=find(dfirep>1)+1;
  if length(fdif)>0
   firepid=1:fdif(1)-1;
  else
    firepid=1:length(firep);
  end
  firepi=firep(firepid);
  re=FM(firepi);
  fire=puesp(firepi);
  puesp=[puesp(1:firepi(1)-1) repmat(fire,1,re(1)) puesp(firepi(end)+1:end)];
  FM=[FM(1:firepi(1)-1) repmat(ones(size(fire)),1,re(1)) FM(firepi(end)+1:end)];
   firep=find(FM>1);
 end 
end
