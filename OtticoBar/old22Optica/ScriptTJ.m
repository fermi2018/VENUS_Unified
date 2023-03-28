
if isfield(radii,'TJ')
 if ifp==-10
 'prima TJ', keyboard
 end

 ifTJ=1;
 if isfield(Pf,'PaTJ')
 PaTJ=Pf.PaTJ;
 PaTJ.puTJ0=radii.TJ.puTJ0;
 thTJ=PaTJ.Nlay_ethc;
 puTJdu=PaTJ.puTJ0;
 puTJ=puTJdu(1:length(puTJdu)+thTJ);
 PaTJ.puTJ=puTJ;
 else
 VECCHIOipar0
 end
 
%  save TJ
 %'sono qui TJ;', keyboard
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
 %ifp4=-ifp;
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
%'dopo', keyboard, keyboard
%'dopo', keyboard, keyboard
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

if Ps.ifpstop==1 | ifp==-10
n_type='real';
 show_VCSELproTYPE
 caxis([0 4])
 
 pausak
n_type='imag'; 
 show_VCSELproTYPE
 title('VCSEL Loss Profile')
 caxis([0 50])
 pausak
end

 nitn=nr.t.'; 
 aitn=aitot(put-1,:).'; 
% aitn=ar.t.'; 
 Litn=dv(put)/1000; 
 fmlsp=abs(frp.t);
 fmlst=abs(frp.t);
  fi1=find(fmlsp(:,1)==0);
  fmlst(fi1,1)=1;
  fi2=find(fmlsp(:,2)==1);
  fmlst(fi1,2)=0; 
 nib=nr.b.'; 
 niat=nr.a.'; 
 nitot=[nitn.'; niat.'; nib.']; 
 nvsave=nv;
 nv3D=nv;

end 