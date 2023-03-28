clear

%load inputTHERM
load TERMOD
addpath('Termico')

[DeltaTold,Tprec,PTherm,T_Contributi,K,TotalHeat]=f_ThermicFun(Tprec,mesh,mode,StrTT,IPLOT);