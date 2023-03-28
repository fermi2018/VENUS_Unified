for ii=1:length(mode.vv_dd)
   Tcont(ii)=squeeze(MODEplot{1}.Temp(ii,end,25))-293;
   T0(ii)=squeeze(MODEplot{1}.Temp(ii,end,1))-293;
end
    
FatQc=Tcont./T0;

figure
grid on,box on,hold on,
plot(mode.vv_dd,FatQc,'.-')
ylabel('FatQcontact'),xlabel('Voltage, V')

figure
grid on,box on,hold on,
plot(mode.ii_dd*1e3,FatQc,'.-')
ylabel('FatQcontact'),xlabel('Current, mA')

figure
grid on,box on,hold on,
plot(mode.DeltaTmax,FatQc,'.-')
ylabel('FatQcontact'),xlabel('\DeltaT_{max}, K')