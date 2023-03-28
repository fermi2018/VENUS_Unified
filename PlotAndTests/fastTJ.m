 clear all
 close all
 
 load TJ
 

 
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
 %ifp4=-10; 
 'salva deb OLD', keyboard
 [dv,nv,aitot,iauto, shavet, fst, anyf, xm, radii, ifield]=...
  TJ_funNEWrelief(ifp4,PaTJ,dv_in,nv_in,aitot_in,iauto_in, shavet_in, fst_in, anyf_in, xm_in, radii_in, ifield_in);
  
  Litot=dv(2:end-1);
  nitot=nv(2:end-1,:);
  
  
  n_type='real';
 show_VCSELproTYPE
 caxis([0 4])
 pausak
n_type='imag'; 
 show_VCSELproTYPE
 caxis([0 50])
 pausak