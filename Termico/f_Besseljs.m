function j = f_Besseljp(m,x)
%-------------------------------------------
% Funzione per calcolare la derivata prima
% della funzione di Bessel di ordine m tra-
% mite la formula:
% xJm(x)=mJm(x)-xJm+1(x)
% per l'origine si usa lo sviluppo di Taylor 
% arrestato al primo termine
%-------------------------------------------
%flag = 0;
%if x(1) == 0 x(1) = 1;
%   flag = 1;
%end %if x(1) == 0 x(1) = 1;

%[I,J]=find(x==0);
%x(I,J) = ones(size(x(I,J)));
%j = m.*besselj(m,x)./x-besselj(m+1,x);
%if ~isempty(I)  
%   if m ==1
%   	j(I,J) = 1/2;
%   else j(I,J) = 0;
%   end %if m ==1
%end %if I ~= []  

j=(besselj(m-2,x)+besselj(m+2,x)-2*besselj(m,x))/4;
   
return