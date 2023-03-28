
function [J,dJdx,dJdx2]=f_EvalPiecewiseBesselFunctions(m,k,x)

kx=k*x;
J=besselj(m,kx);
dJdx=diag(k)*f_Besseljp(m,kx);
dJdx2=diag(k.^2)*f_Besseljs(m,kx);

return




