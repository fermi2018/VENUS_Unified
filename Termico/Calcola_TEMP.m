% DT=flipud(DeltaT);
 fiN=find(isnan(qq)==1);
 qq(fiN)=0;
 qtotV=interp2(roM,zM,qq.',xTdd,zTdd','linear',0).';
 
 % ' qui inter', keyboard
 qtot=zeros(length(xx),length(zz));
 qtot(fir,fi_dd)=qtotV;



if iRAD==1
 qtot=Rspalmato; 
end
  Pqs=qtot*(diag(dz)*S.');
  P=Rdx*Pqs;

   cijL=[];
   for inu=1:length(nu_vett)
    Ji=Js(inu);
    Jsh=sum(Js(1:inu-1));
    QL=reshape(P((1:Ji)+Jsh,:)',Ji*L,1);
    cijL=[cijL; iMLt{inu}*QL];
   end
   CoL=fattore_correttivo*reshape(cijL,L,J)';