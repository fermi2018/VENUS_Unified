function Bprime=mpdbern(x)
% Evaluates the first derivative of the Bernoulli function
%
%                   exp(x)(1-x)-1
%            B'(x)=---------------
%                    (exp(x)-1)^2
%
% for the multi-precision code. This does not require any of the Selbelherr
% approximations, that instead result damaging. The only point is to treat
% correctly the eliminable singularity.
%
% Alberto Tibaldi, 13/07/2016

Bprime=(exp(x).*(mp('1')-x)-mp('1'))./(exp(x)-mp('1')).^mp('2');
ind0=find(x==mp('0'));
Bprime(ind0)=mp('-0.5');


end