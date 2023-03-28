% load LUT_GBTBTonT % only needed for the Tvet values in the LUT
% coefficienti del fit della corrente ricavatao in Main_3
% load Interp_Fit_AlGaAs_Temp     % Main_3, T dependence

if strfind(mode.strName,'Julian')
    load Interp_btj_GaAs_Julian_log
else
    load Interp_Fit_AlGaAs_Temp_log
end

if exist('T_vet')
    mode.T_NEGF=T_vet;
else
    mode.T_NEGF=293;
end
% mode.cBTBT=c;
mode.cBTBT=clog;
mode.I0_NEGF=I0;