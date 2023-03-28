figure,plot(sum(mode.Pst_dd,1),'o-')
iV=input('Which iV?\n');

StrTT=MODEplot{1}.StrTT;
Tprec=squeeze(MODEplot{1}.Tprec(iV,:,:));

modeFAST=f_EvalHeatingTerms(geom,mesh,mode);

IPLOT=2;

[DeltaT,Tprec,PTherm,T_Contributi]=f_ThermicFunMOD(Tprec,mesh,mode,StrTT,IPLOT);