

function [g,rsp,dgdn,Rsp,dRspdn,drspdn]=EvalGains_interp(LUT,vportE,vportH,lambda_i,T)

%%
% vportH : holes density expressed in 1/cm^2
% vportE : electron density expressed in 1/cm^2
% T : temperature in Kelvin
% lambda : wavelength expressed in micron

%%
g=interpn(LUT.por_E,LUT.por_E,LUT.lav,LUT.Tv,LUT.G,vportH,vportE,lambda_i,T,'linear');
rsp=interpn(LUT.por_H,LUT.por_H,LUT.lav,LUT.Tv,LUT.Es,vportH,vportE,lambda_i,T,'linear');
dgdn=interpn(LUT.por_E,LUT.por_H,LUT.lav,LUT.Tv,LUT.dGdE,vportH,vportE,lambda_i,T,'linear');
dgdp=interpn(LUT.por_E,LUT.por_H,LUT.lav,LUT.Tv,LUT.dGdH,vportH,vportE,lambda_i,T,'linear');
dgdn=dgdn+dgdp ;

%% spontaneous emission in mode equation:
drspdn=interpn(LUT.por_E,LUT.por_H,LUT.lav,LUT.Tv,LUT.dEdE,vportH,vportE,lambda_i,T,'linear');
drspdp=interpn(LUT.por_E,LUT.por_H,LUT.lav,LUT.Tv,LUT.dEdH,vportH,vportE,lambda_i,T,'linear');
drspdn=drspdn+drspdp ; 

%% spontaneous recombinations :
Rsp=interpn(LUT.porE_Rsp,LUT.porH_Rsp,LUT.Tv,LUT.Ric,vportH,vportE,T,'linear');
dRspdn=interpn(LUT.porE_Rsp,LUT.porH_Rsp,LUT.Tv,LUT.dRicdE,vportH,vportE,T,'linear');
dRspdp=interpn(LUT.porE_Rsp,LUT.porH_Rsp,LUT.Tv,LUT.dRicdH,vportH,vportE,T,'linear');
dRspdn=dRspdn+dRspdp ;

%%
if sum(isnan([g dgdn Rsp dRspdn rsp drspdn]))
   disp('LUT Out of range') 
   keyboard
end


end






