function [x]=invferdr(y,tol)
%
% y=ferdr(x,1/2), so
% n = NC*ferdr(x,1/2) = NC*y
%
xini = log(y) + y.*(64+0.05524*y.*(64+sqrt(y))).^(-1/4);
x=xini;
res = Inf;
%
while res>tol
x0 = x; 
x = x0 - (ferdr(x0,1/2)-y)./ferdr(x0,-1/2);
res = max(abs(ferdr(x,1/2)-y)./y);
end
