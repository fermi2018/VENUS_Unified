
clear
close all
clc

Tvet=linspace(0,400,1001);
Tstart=300;
DeltaT=50;
tautop=1e-9;
taubottom=1e-12;
% ybottom=1e-9;
% ytop=1e-12;


y=f_RaisedCosine(Tvet,Tstart,DeltaT,tautop,taubottom);

figure
grid on
hold on
box on
plot(Tvet,y)
% set(gca,'yscale','log')
% axis([x(1),x(end),-0.5,1.5])