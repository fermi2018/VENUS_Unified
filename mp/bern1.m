function B=bern1(x)
% Evaluates the Bernoulli function
%
%                     x
%            B(x)=----------
%                  exp(x)-1
%
% according to the formula given in Selberherr, p.169.
%
% FB October 03, 2007

B=x./expm1(x);
ind0=find(x==0);
B(ind0)=1;