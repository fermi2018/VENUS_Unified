function B=mpbern(x)
% Evaluates the Bernoulli function
%
%                     x
%            B(x)=----------
%                  exp(x)-1
%
% for the multi-precision code. This does not require any of the Selbelherr
% approximations, that instead result damaging. The only point is to treat
% correctly the eliminable singularity.
%
% Alberto Tibaldi, 13/07/2016

B=x./(exp(x)-mp('1'));
ind0=find(x==mp('0'));
B(ind0)=mp('1');


end