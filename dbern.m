function Bprime=dbern(x)
% Evaluates the first derivative of the Bernoulli function
%
%                   exp(x)(1-x)-1
%            B'(x)=---------------
%                    (exp(x)-1)^2
%
% according to the formula given in Selberherr, p.169.
% The parameters x1,...,x5 are evaluated for MATLAB.
%
% FB 21/09/96

% Defines the nodes for the approximation
x1=-36.25; x2=-7.63e-6; x3=-x2; x4=32.92; x5=36.5;

% Evaluates the coefficients
B1=zeros(size(x)); B2=B1; B3=B1; B4=B1; B5=B1;

% B1=1 iff x<=x1
B1=sign(sign(x1-x)+1);

% B2=1 iff x>x1 & x<x2
B2=(1-sign(sign(x1-x)+1)).*(1-sign(sign(x-x2)+1));

% B3=1 iff x>=x2 & x<=x3
B3=(sign(sign(x-x2)+1)).*(sign(sign(x3-x)+1));

% B4=1 iff x>x3 & x<x4
B4=(1-sign(sign(x3-x)+1)).*(1-sign(sign(x-x4)+1));

% B5=1 iff x>=x4 & x<x5
B5=(sign(sign(x-x4)+1)).*(1-sign(sign(x-x5)+1));

% Evaluates the function
Bprime=zeros(size(x))-1.*B1+ ...
 ((exp(x).*(1-x)-1)./(exp(x)-1+1.e-99).^2).*B2- ...
 (1/2).*B3+ ...
 ((exp(-x).*(1-x)-exp(-2*x))./(1-exp(-x)+1.e-99).^2).*B4 ...
 +((1-x).*exp(-x)).*B5;

% kludge: eliminate NaN values
Bprime(x>200)=0; Bprime(x<-200)=-1;