function [mu_dop_elec,mu_dop_hole]=arora_gaas(T,N);
%
% C:\ISE_TCAD\tcad\8.0\lib\dessis\MaterialDB\GaAs.par
%
Ar_mumin	= [2.1360e+03 ,	21.48];	% [cm^2/Vs]
Ar_alm  	= [-7.4570e-01 ,	-1.1240e+00];	% [1]
Ar_mud  	= [6.3310e+03 ,	3.3120e+02];	% [cm^2/Vs]
Ar_ald  	= [-2.6870e+00 ,	-2.3660e+00];	% [1]
Ar_N0   	= [7.3450e+16 ,	5.1360e+17];	% [cm^(-3)]
Ar_alN  	= [3.535 ,	3.69];	% [1]
Ar_a    	= [0.6273 ,	0.8057];	% [1]
Ar_ala  	= [-1.4410e-01 ,	0.0000e+00];	% [1]
%
T0=300; % K 
%
muminA=Ar_mumin(1).*(T/T0).^Ar_alm(1); 
mudA = Ar_mud(1).*(T/T0).^Ar_ald(1);
N00=Ar_N0(1).*(T/T0).^Ar_alN(1); 
AA = Ar_a(1).*(T/T0).^Ar_ala(1);
mu_dop_elec = muminA + mudA./(1.+(N./N00).^AA);
%
muminA=Ar_mumin(2).*(T/T0).^Ar_alm(2); 
mudA = Ar_mud(2).*(T/T0).^Ar_ald(2);
N00=Ar_N0(2).*(T/T0).^Ar_alN(2); 
AA = Ar_a(2).*(T/T0).^Ar_ala(2);
mu_dop_hole = muminA + mudA./(1.+(N./N00).^AA);