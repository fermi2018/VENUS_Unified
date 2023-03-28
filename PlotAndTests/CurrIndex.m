function iCurr = CurrIndex(Ivect,Idd)
% Extract the current index corresponding to Ivect, from bias current Idd!

CorLav=Ivect;

pu=[];
for kC=1:length(CorLav)
    corrente=Idd;
    [~,puk]=min(abs(corrente-CorLav(kC)));
    pu=[pu puk];
end

iCurr=pu;
