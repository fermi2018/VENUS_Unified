function B=bern(x)
% Evaluates the Bernoulli function
%
%                     x
%            B(x)=----------
%                  exp(x)-1
%
% according to the formula given in Selberherr, p.169.
% The parameters x1,...,x5 are evaluated for MATLAB.
%
% FB October 03, 2007

% Defines the nodes for the approximation
x1=-36.25; x2=-7.63e-6; x3=-x2; x4=32.92; x5=36.5;
%
B=zeros(size(x));
%
B1 = (x<=x1); 
B(B1) = -x(B1);
%
B2= (x>x1)&(x<x2);
B(B2) = +  x(B2)./(exp(x(B2))-1+1.e-99);
%
B3 = (x>=x2)&(x<=x3);
B(B3) = 1-x(B3)/2;
%
B4 = (x>x3)&(x<x4);
B(B4) = x(B4).*exp(-x(B4))./(1-exp(-x(B4))+1.e-99);
%
B5= (x>=x4)&(x<x5);
B(B5) =  x(B5).*exp(-x(B5));