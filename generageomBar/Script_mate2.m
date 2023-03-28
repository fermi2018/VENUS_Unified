%'entro mate 2', keyboard
materialO1=[repmat(materialO(:,1),1,length(raoADD)-1) materialO];
% materialO1=materialO1(2:end,:);
% nrig=nrig-1;
LMa=size(material,2);
material=materialO1(:,1:LMa);

%for kr=1:nrig
for kr=1:length(Rag)
 R=real(Rag(kr,:));
 MM=materialO1(kr,:);
% keyboard
   clear kme
  for kl=1:length(MM)
   Cel=MM{kl};

   if length(Cel)==0
    kme=kl;
	break
   end

  end
   if exist('kme')
%   'qui', keyboard
    if kme<=LMa
	 for kkk=kme:LMa
	  material(kr,kkk)=materialO1(kme-1);
	 end
	end
   end
end
material(end-1,end)=material(end-1,end-1);
 if flgStop==1
'fine mate2', keyboard
 end