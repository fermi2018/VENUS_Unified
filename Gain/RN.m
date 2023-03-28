n=linspace(0,20,101)*1e18;

tau=100e-9;
B=2e-10;
C=1e-30;

Rn=n/tau+B*n.^2+C*n.^3;

tau_eqIn=1/tau+2*B*n+3*C*n.^2;

figure, plot(n,Rn),

figure, semilogy(n*1e-18,1e9./tau_eqIn)
xlabel('N  10^{18}/cm^3')
ylabel('tau  ns')