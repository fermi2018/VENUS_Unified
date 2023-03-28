function [mode]=velmFitting(indVELM,fitStart,fitDegree,VELMInfo,mode,v0_dd)          

Clight=2.99792458e10;
vph=Clight./mode.nindexQW;

if mode.IdriveON==0
    LVELM=[VELMInfo.Lm];
    Lmv=LVELM(indVELM+[-fitStart:0]);
    lamVELM=[VELMInfo.vlambda];
    lambdav=lamVELM(indVELM+[-fitStart:0]);
    fPdVELM=[VELMInfo.fPdif];
    fPdifv=fPdVELM(indVELM+[-fitStart:0]);
    
    ind=[VELMInfo.indVoltage];
    xVolt=mode.vv0_dd(ind(indVELM+[-fitStart:0]));
    coeff=polyfit(xVolt,Lmv,fitDegree);
    mode.Lm=polyval(coeff,v0_dd);
    coeff=polyfit(xVolt,lambdav,fitDegree);
    mode.vlambda=polyval(coeff,v0_dd);
    coeff=polyfit(xVolt,fPdifv,fitDegree);
    
    mode.fPdif=polyval(coeff,v0_dd)/vph/mode.Cpot*mode.fPdifScaling;
    %             'Interp',keyboard
else
    i0_dd=v0_dd;
    
    LVELM=[VELMInfo.Lm];
    Lmv=LVELM(indVELM+[-fitStart:0]);
    lamVELM=[VELMInfo.vlambda];
    lambdav=lamVELM(indVELM+[-fitStart:0]);
    fPdVELM=[VELMInfo.fPdif];
    fPdifv=fPdVELM(indVELM+[-fitStart:0]);
    
    ind=[VELMInfo.indVoltage];
    xCurr=mode.i0_dd(ind(indVELM+[-fitStart:0]));
    coeff=polyfit(xCurr,Lmv,fitDegree);
    mode.Lm=polyval(coeff,i0_dd);
    coeff=polyfit(xCurr,lambdav,fitDegree);
    mode.vlambda=polyval(coeff,i0_dd);
    coeff=polyfit(xCurr,fPdifv,fitDegree);

    mode.fPdif=polyval(coeff,i0_dd)/vph/mode.Cpot*mode.fPdifScaling;
    
      %             'Interp',keyboard
end
