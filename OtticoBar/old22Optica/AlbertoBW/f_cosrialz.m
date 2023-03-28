function y=cos_rialz(x,mu,T,alpha)

% definita come sul libro di Lo Presti
if alpha ~=0
y=zeros(size(x));
i=find(abs(x-mu)<T*(1-alpha)/2);
y(i)=1;
i=find(abs(x-mu)<=T*(1+alpha)/2 & abs(x-mu)>=T*(1-alpha)/2);
y(i)=0.5*(1-sin(pi/alpha*(abs(x(i)-mu)/T-.5)));
else y=f_rect(x,mu,T);
end

