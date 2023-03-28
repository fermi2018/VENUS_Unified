clear 
% close all 
clc 

global P
% https://www.physicsforums.com/threads/laser-rate-equation-using-matlab.755709/

P=1000;
tau_s = 3e-9; 
N0 = 1e24; 
A =1e-12; 
P0 = 1/(A*tau_s); 
TSPAN = [0 30]; 
Y0 =[0 0]; 
[T,Y] = ode23s(@rate_equation,TSPAN,Y0); 
figure(1)
hold on
subplot(2,1,1) 
plot(T*tau_s ,Y(:,1)*N0) 
title('carriers density in high laser level') % carriers density in high laser level 
subplot(2,1,2) 
plot(T*tau_s ,Y(:,2)*P0) 
title('photons density in activer region') % photons density in activer region 

ind=find(T*tau_s>5e-8);
delta = max(Y(ind,2))-min(Y(ind,2))
