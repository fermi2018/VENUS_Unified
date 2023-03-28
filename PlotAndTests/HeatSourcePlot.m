pointHeat=2:mesh.nny-1;

% Heat sources
HeatJoule=reshape(mode.HeatJoule,mesh.nny,mesh.nnx);
HeatRec_RAD=reshape(mode.HeatRec_RAD,mesh.nny,mesh.nnx);
HeatRec_Cap=reshape(mode.HeatRec_Cap,mesh.nny,mesh.nnx);
HeatOptAbs=reshape(mode.HeatOptAbs,mesh.nny,mesh.nnx);
HeatRec_nr=reshape(mode.HeatRec_13,mesh.nny,mesh.nnx);

figure
hold on,box on,grid on
plot(mesh.ygrid(pointHeat)*1e4,HeatJoule(pointHeat,2),'.-')
plot(mesh.ygrid(pointHeat)*1e4,HeatRec_RAD(pointHeat,2),'.-')
plot(mesh.ygrid(pointHeat)*1e4,HeatRec_Cap(pointHeat,2),'.-')
plot(mesh.ygrid(pointHeat)*1e4,HeatOptAbs(pointHeat,2),'.-')
plot(mesh.ygrid(pointHeat)*1e4,HeatRec_nr(pointHeat,2),'.-')
xlim([mesh.ygrid(mesh.ICAV(1))*1e4-0.5 mesh.ygrid(mesh.ICAV(1)+find(diff(mesh.ICAV)>1,1))*1e4+0.5])
set(gca,'yscale','log')
ylabel('Heat Sources, W/cm^3'),xlabel('z, \mum')
legend('Joule','Rad','Ccap','FCA','NR','location','best')

figure
hold on,box on,grid on
plot(mode.ii_dd*1e3,mode.DeltaTmax_Joule,'.-')
plot(mode.ii_dd*1e3,mode.DeltaTmax_RAD,'.-')
plot(mode.ii_dd*1e3,mode.DeltaTmax_Ccap,'.-')
plot(mode.ii_dd*1e3,mode.DeltaTmax_OptAbs,'.-')
plot(mode.ii_dd*1e3,mode.DeltaTmax_srhAu,'.-')
ylabel('\DeltaT, K'),xlabel('Current, mA')
xlim([0 max(mode.ii_dd)*1e3])
legend('Joule','Rad','Ccap','FCA','NR','location','best')


% figure,plot(sum(mode.Pst_dd,1),'o-'),ylabel('Optical power, mW')
% figure,plot(mode.ii_dd*1e3,'o-'),ylabel('Current, mA')
% % keyboard
% % iV=25;
% iV=input('iV?\n');
%
modep=MODEplot{1};

DensityCurrentPlot

HeatJoulep=(squeeze(modep.HeatJoule(iV-1,pointHeat,2)));
HeatRec_RADp=(squeeze(modep.HeatRec_RAD(iV-1,pointHeat,2)));
HeatRec_Capp=(squeeze(modep.HeatRec_Cap(iV-1,pointHeat,2)));
HeatOptAbsp=(squeeze(modep.HeatOptAbs(iV-1,pointHeat,2)));
HeatRec_nrp=(squeeze(modep.HeatRec_13(iV-1,pointHeat,2)));

figure
hold on,box on,grid on
plot(mesh.ygrid(pointHeat)*1e4,HeatJoulep,'.-')
plot(mesh.ygrid(pointHeat)*1e4,HeatRec_RADp,'.-')
plot(mesh.ygrid(pointHeat)*1e4,HeatRec_Capp,'.-')
plot(mesh.ygrid(pointHeat)*1e4,HeatOptAbsp,'.-')
plot(mesh.ygrid(pointHeat)*1e4,HeatRec_nrp,'.-')
xlim([mesh.ygrid(mesh.ICAV(1))*1e4-0.5 mesh.ygrid(mesh.ICAV(1)+find(diff(mesh.ICAV)>1,1))*1e4+0.5])
ylabel('Heat Sources, W/cm^3'),xlabel('Current, mA')
legend('Joule','Rad','Ccap','FCA','NR','location','best')
set(gca,'yscale','log')
title(['Sources @ I=',num2str(mode.ii_dd(iV)*1e3),' mA'])

if mode.quasi1D==1
    fatCorr=mode.fattore_correttivo(iV);
    format short
    T_Joule=max(f_ThermD1ANA_Contributi(mesh,mode.condzTe(iV,:),HeatJoulep*1e12,fatCorr))
    T_Rec_srhAu=max(f_ThermD1ANA_Contributi(mesh,mode.condzTe(iV,:),HeatRec_nrp*1e12,fatCorr))
    T_Rec_Cap=max(f_ThermD1ANA_Contributi(mesh,mode.condzTe(iV,:),HeatRec_Capp*1e12,fatCorr))
    T_Rec_RAD=max(f_ThermD1ANA_Contributi(mesh,mode.condzTe(iV,:),HeatRec_RADp*1e12,fatCorr))
    T_OptAbs=max(f_ThermD1ANA_Contributi(mesh,mode.condzTe(iV,:),HeatOptAbsp*1e12,fatCorr))
    format shorte 

    figure
    hold on,box on,grid on
    plot(mode.ii_dd*1e3,mode.fattore_correttivo,'.-')
    xlim([0.1 mode.ii_dd(end)*1e3])
    xlabel('Current, mA'),ylabel('fattore\_correttivo')
else
%     fattore_correttivo=mode.PDissPred(2:end)./mode.PTherm;
% 
%     figure
%     hold on,box on,grid on
%     plot(mode.ii_dd*1e3,fattore_correttivo,'.-')
%     xlim([0.1 mode.ii_dd(end)*1e3])
%     xlabel('Current, mA'),ylabel('fattore\_correttivo')

end

