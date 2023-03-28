
function dy = rate_equation(t,y) 

global P

dy = zeros(2,1); tau_s = 3e-9; % carriers lifetime 
tau_p = 1e-12; % photons lifetime 
A = 1e-12; % linear gain costant 
N0 = 1e24; % trasparency carries density 
V = 3.75e-14; % modal volume 
gamma = 1e-5; % gain compression factor 
q = 1.6e-19; % electron charge 
I0 = N0*q*V/tau_s; % trasparency current 
tau_norm = tau_s/tau_p; 
eta = A*tau_p*N0; % efficiency 
I = 2.5*I0+0.1*I0*cos(P*t); % pumping current ( try: from I0 to 3*I0 for example ...and see what happens!) 
dy(1)= I/I0 -y(2)*(y(1) - 1) -y(1); dy(2) = tau_norm*(y(2)*(eta*(y(1) - 1) -1) + gamma*eta*y(1));

