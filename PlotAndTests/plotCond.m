figure(80),hold on,grid on,plot(mode.vv_dd,mode.ii_dd*1e3,'.-')
xlabel('Voltage, V'),ylabel('Current, mA')

figure(81),hold on,grid on,plot(mode.ii_dd*1e3,sum(mode.Pst_dd,1),'.-')
xlabel('Current, mA'),ylabel('Optical power, mW')

figure(85),hold on,grid on,plot(mode.ii_dd*1e3,mode.DeltaTmax,'.-')
xlabel('Current, mA'),ylabel('\DeltaT_{max}, K')

figure(1111)
subplot(221)
if isfield(mode,'C1')==1 && mode.C1>1
    semilogy(mode.ii_dd*1e3,mode.cond_intermedio,'.-')
else
    semilogy(mode.ii_dd*1e3,mode.cond_prima,'.-')
end
ylabel('Cond (before)')
hold on,grid on,box on
subplot(222)
if isfield(mode,'C1')==1 && mode.C1>1
    semilogy(mode.vv_dd,mode.cond_intermedio,'.-')
else
    semilogy(mode.vv_dd,mode.cond_prima,'.-')
end
hold on,grid on,box on
subplot(223)
semilogy(mode.ii_dd*1e3,mode.cond_dopo,'.-')
xlabel('Current, mA'),ylabel('Cond (after)')
hold on,grid on,box on
subplot(224)
semilogy(mode.vv_dd,mode.cond_dopo,'.-')
hold on,grid on,box on
xlabel('Voltage, V')
