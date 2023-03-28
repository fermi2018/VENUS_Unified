% % select puC from the LI (new)
% figure,plot(sum(mode.Pst_dd,1),'o-')
% puC=[20 25 30];
% modePlot=MODEplot{1};
% mode=MODEplot{1};
%

I=mode.ii_dd*1000;
%'pa', keyboard
V=mode.vv_dd;
P=sum(mode.Pst_dd,1)+mode.Psp_dd;


Rec=mode.IntRec;
Ccap=-sum(mode.IntCcapN,2);
Rnor=diag(1./Ccap)*Rec;


T0=mode.T0;
WQW=mode.vWMQW{1};
leQW=length(modePlot.nQW{end}{2});
%    x=modePlot.x;
x=modePlot.x(1:leQW);

%         x=mode.x*1e-4;
qel=1.6e-19;
NumQW=length(mode.vWMQW);
FatNP_Auger= mode.FatNP_Auger;
CN_Auger= mode.CN_Auger;
CTemp_Auger=mode.CTemp_Auger;
for kv=2:length(modePlot.ii_dd)
    n2D=modePlot.Elqw(kv,:);
    p2D=modePlot.Hoqw(kv,:);
    TQW=modePlot.Tqw(kv,:);
    
    Cnnp =CN_Auger*exp(CTemp_Auger.*(TQW-T0)./300)*1e-30/WQW.^2;
    %	Cnnp =CN_Auger*abs(1+CTemp_Auger*(TQW-T0)./300)*1e-30/(WQW.^2);
    Cppn =FatNP_Auger.*Cnnp;
    % xa=TQW-300;
    % fT=(xa/20).^4.*exp(-xa/xad)*1e-30/(WQW.^2);
    % Cnnp =CN_Auger*fT
    % Cppn =FatNP_Auger.*CN_Auger.*fT;
    np2D=n2D.*p2D;
    RAuger = NumQW*qel*(Cnnp.*n2D+Cppn.*p2D).*np2D;
    RAuger_n = NumQW*qel*(Cnnp.*n2D).*np2D;
    RAuger_p = NumQW*qel*(Cppn.*p2D).*np2D;
    Nlea(kv)=2*pi*trapz(x,x.*RAuger)/Ccap(kv);
    Nlea_n(kv)=2*pi*trapz(x,x.*RAuger_n)/Ccap(kv);
    Nlea_p(kv)=2*pi*trapz(x,x.*RAuger_p)/Ccap(kv);
    % 'ver;,', keyboard
end
% 'ver;,', keyboard
fi=find(Nlea>1);
Nlea(fi)=1;
Cd=mode.ii_dd*1000;

