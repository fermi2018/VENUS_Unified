
if isfield(radii,'TJ')

 if ifp==-10
 'prima TJ', keyboard
 end

 ifTJ=1;
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
 
 
 'sono qui TJ;', keyboard
 Dop_in=dglo.Dop;
 dv_in=dv;
 nv_in=nv0;
 iauto_in=iauto;
 aitot=[ar.t; ar.a; ar.b];
 aitot_in=aitot;
 shavet_in=shavet;
 fst_in=fst;
 anyf_in=anyf;
 xm_in=xm;
 radii_in=radii;
 ifield_in=ifield;
 ifp4=-4;
 ifp4=-ifp;
 %'salva deb', keyboard
 %[dv,nv,aitot,iauto, shavet, fst, anyf, xm, Dop, radii, ifield]=...
 % TJ_funNEWreliefDop(ifp4,PaTJ,dv_in,nv_in,aitot_in,iauto_in, shavet_in, fst_in, anyf_in, xm_in, Dop_in, radii_in, ifield_in);
 [dv,nv,aitot,iauto, shavet, fst, anyf, xm, radii, ifield]=...
  TJ_funNEWrelief(ifp4,PaTJ,dv_in,nv_in,aitot_in,iauto_in, shavet_in, fst_in, anyf_in, xm_in, radii_in, ifield_in);

% TJ_funOLD(ifp4,PaTJ,dv_in,nv_in,aitot_in,iauto_in, shavet_in, fst_in, anyf_in, xm_in, radii_in, ifield_in);


dv_c=dv;
nv_c=nv;
aitot_c=aitot;
iauto_c=iauto;
shavet_c=shavet;
fst_c=fst;
xm_c=xm;
radii_c=radii;
save sa  dv_c nv_c aitot_c iauto_c  shavet_c  fst_c  xm_c  radii_c 







%'dopo', keyboard, keyboard
 dv3D=dv;
 nv0=nv;
  reassign
  any.t=anyf(put,:);
  any.b=anyf(pub,:);
  any.a=anyf(pua,:);
  

% fst=[[0 1]; fmlstot; [0 1]];
 fmlstot=abs(fst(2:end-1,:));
% close(2:end)

%' doppo TJ   aaa ', keyboard 

Litot=dv(2:end-1);
nitot=nv(2:end-1,:);

end 