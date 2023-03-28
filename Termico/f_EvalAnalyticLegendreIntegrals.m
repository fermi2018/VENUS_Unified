function [D1,D2]=f_EvalAnalyticLegendreIntegrals(M)
%
% function [D1,D2]=f_EvalAnalyticLegendreIntegrals(M)
% Version 1.0
%
% This function computes the change-of-basis matrices aimed at computing
% analytically the SEM integrals when Legendre polynomials are used as
% basis functions. These expressions are taken from
% Gottlieb-Orszag, Numerical Analysis of Spectral Methods, (A.30) and
% (A.31).
%
% The resulting integrals equal the ones obtained with numerical
% integration and inner product weight equal to 1 (as appropriate for
% Legendre polynomials).
%
% Alberto Tibaldi, 22/01/2016

D1=zeros(M+1,M+1);
for m=0:M-1
    p=m-1;
    while p<M-1
        p=p+2;
        D1(m+1,p+1)=(2*m+1);
    end
end

D2=zeros(M+1,M+1);
for m=0:M-2
    p=m;
    while p<M-1
        p=p+2;
        D2(m+1,p+1)=(p+m+1)*(p-m)*(m+0.5);
    end
end

return