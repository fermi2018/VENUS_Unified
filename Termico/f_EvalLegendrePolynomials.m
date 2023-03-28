
function [Pn,dPndx,d2Pndx2]=f_EvalLegendrePolynomials(n,x)
%
% function [Pn,dPndx,d2Pndx2]=f_EvalLegendrePolynomials(n,x)
% Version 1.0
%
% This function generates the Legendre polynomials of degree n in points x.
% Both n and x are intended to be row vectors, with different dimensions,
% where the "usual" normalization is adopted.
%
% The function is based on recurrence relations for the calculation of both
% functions and first, second derivatives. The values of the function and
% of the derivatives at the [-1,1] interval endpoints are computed
% analytically.
%
% Alberto Tibaldi, 19/01/2016

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finding internal points, x=+1 and x=-1 points (for analytical limits)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
threshold=1e-11; % threshold to identify internal points 
ind=find(abs(x-1)>threshold & abs(x-(-1))>threshold); % internal points
ind_p1=find(abs(x-1)<threshold); % x=+1
ind_m1=find(abs(x-(-1))<threshold); % x=-1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Synthesis of Legendre polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pn=ones(length(n),length(x)); % P_0(x) = 1
Pn(2,:)=x; % P_1(x) = x
if(not(isempty(ind)))
    for indf=2:(length(n)-1)
        Pn(indf+1,ind)=((2*n(indf)+1).*x(ind).*Pn(indf,ind)-n(indf).*Pn(indf-1,ind))./(n(indf)+1);
    end
end
%-- x=+1 limit is already set
Pn_m1=(-1).^n.'; % limit for x=-1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Synthesis of Legendre polynomials first derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,N]=meshgrid(x,n); % define meshgrid to use MATLAB matrix features
dPndx=zeros(size(Pn)); % zero order derivative is 0
if(not(isempty(ind)))
    dPndx(2:end,ind)=(X(2:end,ind).*Pn(2:end,ind)-Pn(1:end-1,ind)).*N(2:end,ind)./(X(2:end,ind).^2-1);
end
dPndx_p1=(n.*(n+1)/2).'; % x=+1 limit
dPndx_m1=((-1).^(n-1).').*dPndx_p1; % x=-1 limit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Synthesis of Legendre polynomials second derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d2Pndx2=zeros(size(Pn));
if(not(isempty(ind)))
    d2Pndx2(2:end,ind)=(2.*X(2:end,ind).*dPndx(2:end,ind)-N(2:end,ind).*(N(2:end,ind)+1).*Pn(2:end,ind))./(1-X(2:end,ind).^2);
end
d2Pndx2_p1=((n.*(n.^3+2.*n.^2-n-2)/8)).'; % x=+1 limit
d2Pndx2_m1=((-1).^(n).*(n.*(n.^3+2.*n.^2-n-2)/8)).'; % x=-1 limit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting limit values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(not(isempty(ind_p1)))
    dPndx(:,ind_p1)=dPndx_p1;
    d2Pndx2(:,ind_p1)=d2Pndx2_p1;
end    
if(not(isempty(ind_m1)))
    Pn(:,ind_m1)=Pn_m1; % function
    dPndx(:,ind_m1)=dPndx_m1;
    d2Pndx2(:,ind_m1)=d2Pndx2_m1;
end

return