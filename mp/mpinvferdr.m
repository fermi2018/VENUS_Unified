function [x]=mpinvferdr(y,tol)
%
% y=ferdr(x,1/2), so
% n = NC*ferdr(x,1/2) = NC*y
%
xini = log(y) + y.*(mp('64')+mp('0.05524')*y.*(mp('64')+sqrt(y))).^(mp('-0.25'));
x=xini;
res = Inf;
%
while res>tol
x0 = x; 
x = x0 - (mpferdr(x0,mp('0.5'))-y)./mpferdr(x0,mp('-0.5'));
res = max(abs(mpferdr(x,mp('0.5'))-y)./y);
end
