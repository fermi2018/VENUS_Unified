%Model based on Gehrsitz
%Model adjusted for small x, close to band gap
%WL: Wavelenght (m)
%x: aluminum fraction 0 ... 1
%T: Temperature in K
%Return: Refractive index
function [n]=gehrsitz_WL_V02(WL,x)

T = 293;

%Constant values
h=6.6260e-34;
c=299792458;
e=1.6020e-019;

E = h*c/(WL*e);
Elsq=(E/1.239856)^2;



%Parameter anfitten; Alle Energien in E^2/µm
%*** A(x,T) ***
% A0=7.3377+5.534e-4*T-0.356e-6*T^2; %FIT 1
A0=5.9613+7.178e-4*T-0.953e-6*T^2; %FIT 2
A=A0-16.159*x+43.511*x.^2-71.317*x.^3+57.535*x.^4-17.451*x.^5;

%*** E0(x,T) ***
%Faktoren für Bandlücke S.7829
Eg0=1.5192;
Edeb=15.9e-3;
Et0=33.6e-3;
S=1.8;
St0=1.1;
kB=8.61708e-5;
% EGaAs in eV!
EGaAs=Eg0+S*Edeb*(1-coth(Edeb/(2*kB*T)))+St0*Et0*(1-coth(Et0/(2*kB*T)));
% EAlGaAs=(EGaAs+1.36*x+0.22*x^2); 
% EGaAs auf 1/µm skalierne und E0^2!
E0=(EGaAs/1.239856+1.1308*x+0.1436*x.^2).^2;



%*** E1(x,T) ***
% E10=3.791-3.779e-4*T-1.121e-6*T^2; %FIT 1
E10=4.7171-3.237e-4*T-1.358e-6*T^2; %FIT 2
E1=E10+11.006*x-3.08*x.^2;
%*** E2 ***
E2=0.724e-3;
%*** E3 ***
E3=1.331e-3;

%*** C0 ***
C0=1./(50.535-150.7*x-60.209*x.^2+797.16*x.^3-1125*x.^4+503.79*x.^5);
%*** C1 ***
C1=21.5647+113.74*x-122.5*x.^2+108.401*x.^3-47.318*x.^4;
%*** C2 ***
C2=1.55e-3;
%*** C3 ***
C3=2.61e-3;

R=(1-x)*C2/(E2-Elsq)+x*C3/(E3-Elsq);

%prevent of discontinuity close to direct band gap
B=(E0-Elsq);

index_vectorial = B<0.025;
B(index_vectorial) = 0.025;

% if B<0.025
%     B=0.025;
% end

nsq=A+C0./(B)+C1./(E1-Elsq)+R;
n=sqrt(nsq);

