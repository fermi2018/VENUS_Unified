% Loading physical constants
qel=mp('1.6021766208e-019'); % Elementary charge (C)
m0=mp('9.10938291e-31'); % Electron mass (kg)
kB=mp('1.3806488e-23'); % Boltzmann constant (J/K)
h=mp('6.626070040e-34'); % Planck constant (J*s)
hbar=h/mp('2/pi'); % Planck constant (J*s)
Clight=mp('2.99792458e+10'); % Speed of light (cm/s)
mu0=mp('4*pi*1e-9'); % Magnetic permeability (H/cm)
eps0=mp('1')/(mu0*Clight^mp('2')); % Dielectric permittivity (F/cm)

mpqelNorm=qel*mp(mode.CarrierNorm);
mpqelNorm2D=qel*mp(mode.CarrierNorm2D);
