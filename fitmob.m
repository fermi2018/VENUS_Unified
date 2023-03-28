% Used in MDPIAS
T=[0 80];
y=[1.2 4];
tv1=linspace(1,80,101);
cot=polyfit(T,y,1);
figure, plot(T,y,'o',tv1,polyval(cot,tv1))
save COT cot

% MDPIAS - right contact 
T=[0 25 80];
y=[1.2 2 5.5];
tv1=linspace(1,80,101);
cot=polyfit(T,y,2);
figure(2222),hold on, plot(T,y,'o',tv1,polyval(cot,tv1))
save COT cot

% Fit VENUS 3? 
T=[0 30 80];
y=[1.2 3 9];
tv1=linspace(1,80,101);
cot=polyfit(T,y,2);
figure(2222),hold on, plot(T,y,'o',tv1,polyval(cot,tv1))
save COT cot

% T=[0 25 80];
% y=[1.2 2 11];
% tv1=linspace(1,80,101);
% cot=polyfit(T,y,2);
% figure, plot(T,y,'o',tv1,polyval(cot,tv1))
% save COT cot

% % Absolute temperature
% T=[0 80 130]+293;
% y=[1.2 4 5];
% tv1=linspace(1,130,101)+293;
% cot=polyfit(T,y,2);
% figure, plot(T,y,'o',tv1,polyval(cot,tv1))
% save COT cot
% 
% T=[0 80 130]+293;
% y=[1.2 4.2 5];
% tv1=linspace(1,130,101)+293;
% cot=polyfit(T,y,2);
% figure, plot(T,y,'o',tv1,polyval(cot,tv1))
% save COT cot