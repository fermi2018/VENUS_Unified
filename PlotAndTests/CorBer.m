kB=1.3806488e-23; % Boltzmann constant (J/K)
Vcori=Vcor(end);
mesh.nny=length(mode.y);
Temp=mode.Tqw(Vcori,1);
in_i = 1:mesh.nny-1;
in_j = 2:mesh.nny;

elec=mode.El(Vcori,:);

elec_i = elec(in_i);
elec_j = elec(in_j);
l_ij = diff(mode.y(1:mesh.nny));

phi_i = mode.Phi0(in_i);
phi_j = mode.Phi0(in_j);

VT = kB*Temp/qel;
Mob=2000;
Dn = VT.*Mob;

delta_ij = (phi_i-phi_j)./VT;

% Jprova = qel.*Dn./l_ij.*(elec_j.*bern(-delta_ij)-elec_i.*bern(delta_ij))
Jprova = Dn.*(elec_j.*bern(-delta_ij)-elec_i.*bern(delta_ij))./l_ij



figure,plot(yqw(2:end),Jprova);%,ygrid,mode.Jn_y(in_i))
% figure,plot(ygrid,mode.Jn_y(in_i));%,ygrid,mode.Jn_y(in_i))
xlim([-100,100])