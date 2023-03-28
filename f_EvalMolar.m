
function [f]=f_EvalMolar(x,y,gvet)

a=gvet(1);
b=gvet(2);
c=gvet(3);

f=a.*x+b.*y+c;

return