% Impulso rettangolare. supporto=T; max=1
function [y]= rect(x,T)
n=length(x);
y=zeros(1,n);
i=find(abs(x)<=T/2);
m=length(i);
y(i)=ones(1,m);
return