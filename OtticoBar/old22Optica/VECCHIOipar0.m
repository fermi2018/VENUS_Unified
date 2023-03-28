 PaTJ=radii.TJ;
 puTJdu=PaTJ.puTJ0;

 if isfield(PaTJ,'ipar0')
  iPar0=PaTJ.ipar0;
  for kpai=1:size(iPar0,3)
   pai=iPar0(:,:,kpai);
   if pai(2)==-6
    PaTJ.Ram=par_in{pai(1)};
   elseif pai(2)==-9
    thTJ=par_in{pai(1)};
    iNO=0;
    if length(puTJdu)<=thTJ | imag(thTJ)~=0
     'NO TJ', 

     if ifp==-10
      keyboard
     end
     iNO=1;
    end
	
    if thTJ>=0 
     puTJ=puTJdu(1+thTJ:length(puTJdu));
    else
     puTJ=puTJdu(1:length(puTJdu)+thTJ);
    end 
%	'puTJ', keyboard
   end
  end

  ipar1=ipar(:,1,1);
  if max(ipar1)>length(par_in)
   puREL=(find(ipar1>length(par_in)));
   thREL=sum(dv(puREL));
   thTunnel=sum(dv(puTJ));
   if thTunnel>thREL
    'problema da risolvere', keyboard
   else
    Trelie(1)=thREL-thTunnel;   
    Trelie(2)=thTunnel;   
	dv(puREL)=Trelie;
	
   end
   
  end

  radii.TJ=PaTJ;
  PaTJ.puTJ=puTJ;
 end