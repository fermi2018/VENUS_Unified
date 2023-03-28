function [n,dn] = mpferdr2D(xnum,xden)
%
% function [n,dn] = ferdr2D(xnum,xden)
%
% This function computes the 2D carrier distribution assuming that the
% energy interval is finite; such choice leads to the introduction of the
% "denominator terms".
%
% The input parameters xnum, xden are expected to be matrices; by this way,
% n, dn are matrices where each row corresponds to a different bound state.
% The variable "n2D" (p2D) is given by the sum of all these contributions.
%
% This function does not contain any information about the effective
% masses, therefore they should be multiplied a posteriori.
%
% n can be either n2D/Nc2D or p2D/Nv2D.
%
% dn can be either dn2D/Nc2D or dp2D/Nv2D; this is the derivative with
% respect to phi/Vt ; so, in the thermodynamic equilibrium solver it is
% already the derivative, whereas it should be divided by Vt to find the
% derivative w.r.t. EFn2D, EFp2D.
%
% Ver. 0.1, 09/03/2016
% Alberto Tibaldi

% "numerator" terms
expnum = exp(xnum);
lognum=log(mp('1')+expnum);
dlognum_phi=expnum./(mp('1')+expnum);

% "denominator" terms
expden = exp(xden);
logden=log(mp('1')+expden);
dlogden_phi=expden./(mp('1')+expden);

% Test debug case: without denominator, this coincides to the unbound
% energy case
% n=lognum;
% dn=dlognum_phi;

% Correct case: difference of logarithms
n=lognum-logden;
dn=dlognum_phi-dlogden_phi;

return