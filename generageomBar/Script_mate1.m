if flgStop==1
    'entro mate 1', keyboard
end

for kr=1:nrig
 R=real(Rag(kr,:));
 MM=materialO(kr,:);
  Rv=[];
  rox=0;
 for kl=1:length(MM)
  fiLine=strcmp(MM(kl),'Line');
  if fiLine==1
   klS=kl;
   if length(R)>1
    Rv=R(kl+[-1 0]);
   else 
    Rv=R;
   end
  end
  fiLine=strcmp(MM(kl),'Ground');
  if fiLine==1
   klS=kl;
   Rv=R(kl+[-1 0]);
  end  
  fiLine=strcmp(MM(kl),'AlGaAs');
  if fiLine==1
   if kl>1
%    Rv=R(kl+[-1 0]);
    if length(R)>1
     Rv=R(kl+[-1 0]);
    else 
     Rv=R;
    end    
    fiV=find(Rv>0);
    Rv=Rv(fiV);
%    'alta', keyboard
   end 
  end  
  fiLine=strcmp(MM(kl),'AlOx');
  if fiLine==1
   klS=kl;
   Rv=R(kl-1);
   rox=Rox;
%  'qui ox', keyboard
  end  
 end 
% Rv, 
 if length(Rv)==2
  fi1=find(Rv(1)==ra_dd)+1;
    if rox>0
%     'rox', keyboard    
     fi1=fi1+1;
    end
  fi2=find(Rv(2)==ra_dd); 
  fiC=find(R==Rv(2));
  for kc=fi1:fi2
   if length(fiC)>0
    material{kr,kc}=materialO{kr,fiC};
   end 
  end
  elseif length(Rv)==1
    fi1=find(Rv(1)==ra_dd)+1;
    fi2=ncol;
    fiC=find(R==Rv(1))+1;
%    fiC,  keyboard
   for kc=fi1:fi2
    material{kr,kc}=materialO{kr,fiC};
   end  
  end 
end

if flgStop==1
'fine mate1', keyboard
end