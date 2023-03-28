function y=f_RaisedCosine(Tvec,Tstart,DeltaT,tautop,taubottom)

% function y=f_RaisedCosine(Tvec,Tstart,DeltaT)
%
% This function has been defined, starting from an older function from
% Renato Orta (applied in RCWA codes), to obtain a temperature dependence
% for the capture/escape times.
%
% Input parameters
%
% Tvec: (possibly) vector of temperatures to be considered
% Tstart: starting temperature from which the model changes
% DeltaT: temperature range for the variation of the model
% tautop: value of the capture/escape time for T < Tstart
% tautop: value of the capture/escape time for T > Tstart + DeltaT
%
% Alberto Tibaldi, 05/03/2018

% Writing parameters according to Lo Presti definitions
T=2*Tstart+DeltaT;
alpha=DeltaT/T;
%
% Lo Presti definition
if alpha ~=0
    y=zeros(size(Tvec));
    i=find(abs(Tvec)<T*(1-alpha)/2);
    y(i)=1;
    i=find(abs(Tvec)<=T*(1+alpha)/2 & abs(Tvec)>=T*(1-alpha)/2);
    y(i)=0.5*(1-sin(pi/alpha*(abs(Tvec(i))/T-.5)));
else
    y=f_Rect(Tvec,T);
end

y = (tautop-taubottom)*y+taubottom;

