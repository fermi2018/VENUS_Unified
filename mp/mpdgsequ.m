function [R,C]=mpdgsequ(A)
%
A=abs(A);
% smlnum=1e-12;
% bignum=1e20;

smlnum=mp('1e-20');
bignum=mp('1e80');
% smlnum=mp('1e-4900');
% bignum=mp('1e4900');
%
[m,n]=size(A);
r = max(A,[],2);
r = 1./min(max(r,smlnum),bignum); 
R=sparse(1:m,1:m,r);
A=R*A;
c = max(A,[],1);
c = 1./min(max(c,smlnum),bignum);
C=sparse(1:n,1:n,c);