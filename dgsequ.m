function [R,C]=dgsequ(A);
%
A=abs(A);
%
smlnum=1e-12;
bignum=1e20;
%
[m,n]=size(A);
r = max(A,[],2);
r = 1./min(max(r,smlnum),bignum); 
R=sparse(1:m,1:m,r);
A=R*A;
c = max(A,[],1);
c = 1./min(max(c,smlnum),bignum);
C=sparse(1:n,1:n,c);
% C=sparse(1:n,1:n,1,n,n);