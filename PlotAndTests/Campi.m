%load Campi
%puCs=puC;
%save puC puC
close all

% Old visualization...
% load puC
%
% puC=puC([3 6 9]);

% select puC from the LI (new)
figure,plot(sum(mode.Pst_dd,1),'o-')
% puC=[20 25 30];
puC=pu;
modePlot=MODEplot{1};
mode=MODEplot{1};
%
I=mode.ii_dd*1e3;
LW=[1 1.8 2.5];     % linewidth vector


x=modePlot.x;
Te=mode.Temp;
Tqw=squeeze(Te(puC,150,:));
figure(30),set(gcf,'pos',[471   369   945   540])
subplot(121)
for kcur=1:length(puC)
    plot(x,Tqw(kcur,:),'r','linewidth',LW(kcur))
    hold on
    
    lgd{kcur}=['I=',num2str(I(puC(kcur))),' mA'];
end
legend(lgd)
xlabel('\rho, \mum')
ylabel('Active region temp., K')

xlim([0 20])
xfi=mode.x;
Ca=mode.E2*1e-8;

COL='rgbcm';
dox=modePlot.rox;
Ca1=squeeze(Ca(puC,1,:));
%       Ca2=squeeze(Ca(1:PAc:end,2,:));
%      figure(30),
subplot(122)
%hc=figure;
%set(hc,'pos',[256   425   915   525])


for kcur=1:length(puC)
    plot(xfi,Ca1(kcur,:),COL(1),'linewidth',LW(kcur)), hold on,  ResetColor,
end

for kc=2:size(Ca,2)
    Ca3=squeeze(Ca(puC,kc,:));
    for kcur=1:length(puC)
        plot(xfi,Ca3(kcur,:),COL(kc),'linewidth',LW(kcur)), ResetColor,
    end
end
xlim([1e-3 max([2*dox 6])])
xlabel('\rho, \mum')
ylabel('Optical Field Intensity, Arb.Un.')
