function [Efp2D]=mpinvferdr_p2D_0(p2D,Ed,Ev0,Nv2D,Vt,LowerBound,UpperBound)
%
% function [Efp2D]=invferdr_p2D(p2D,Ed,Ev0,Nv2D,Vt,LowerBound,UpperBound)
%
% This function is aimed at computing the numerical inverse of ferdr2D_0.
% Such calculation is required, since, if the energy integral is computed
% up to a finite bound or more bound states are used, the inverse cannot be
% written in closed form.
% From the 2D electrons (sum on all the bound states) a unique Fermi level
% EFn2D is computed; this can be used then to compute each contribution to
% the Poisson equation.
% This is done by means of a two-step algorithm: the first step consists of
% finding a guess of the quasi-Fermi level through a low-accuracy bisection
% method. Then, the second step is a high accuracy Newton loop.
%
% p2D (vector): sum of 2D holes of all bound states
%
% Ed (matrix): electrons energy levels; rows: different energy levels;
% columns: different lateral points (to account for different potential
% heights)
%
% Ev0 (matrix): upper energy integration bound (eV)
%
% Nv2D (matrix): 2D Density of States (DoS); it is a matrix since different
% effective masses can be used (1/cm^2)
%
% LowerBound,UpperBound: lower and upper bounds for the bisection algorithm
% for Efp2D (vectors, eV)

%
% Ver. 0.1, 09/03/2016
% Ver. 0.2, 29/03/2016 - introduced automatic guess (with bisection)
% 
% Alberto Tibaldi

% Two-step algorithm: first Dichotomy to achieve a guess, then Newton

tolBisection=mp('1e-2'); % tolerance for solution via bisection method
tolNewton=mp('1e-13'); % tolerance for solution via Newton method
MaxIterNewton=50; % maximum number of iterations for Newton loop

LHS=p2D;
Efp2D_GuessForNewton=(LowerBound+UpperBound)/mp('2');
% Setting to inf to enter in the loop at least once
GuessPreviousStep=inf*ones(size(Efp2D_GuessForNewton),'mp');

while(max(abs(Efp2D_GuessForNewton-GuessPreviousStep))>tolBisection)
    
    GuessPreviousStep=Efp2D_GuessForNewton;
    % Calculation of the upper bound for the current iteration
    xnum=(Ed-ones(size(Ed,1),1,'mp')*UpperBound)./Vt; 
    xden=(Ev0-ones(size(Ed,1),1,'mp')*UpperBound)./Vt;
    [p2D_UB,dum] = mpferdr2D_0(xnum,xden);
    p2D_UB=sum(Nv2D.*p2D_UB,1);
    fU=p2D_UB-LHS;
    %
    % Calculation of the lower bound for the current iteration
    xnum=(Ed-ones(size(Ed,1),1,'mp')*LowerBound)./Vt; 
    xden=(Ev0-ones(size(Ed,1),1,'mp')*LowerBound)./Vt;
    [p2D_LB,dum] = mpferdr2D_0(xnum,xden);
    p2D_LB=sum(Nv2D.*p2D_LB,1);
    fL=p2D_LB-LHS;
    %
    % Calculation of the guess for the current iteration
    xnum=(Ed-ones(size(Ed,1),1,'mp')*Efp2D_GuessForNewton)./Vt; 
    xden=(Ev0-ones(size(Ed,1),1,'mp')*Efp2D_GuessForNewton)./Vt;
    [p2D_C,dum] = mpferdr2D_0(xnum,xden);
    p2D_C=sum(Nv2D.*p2D_C,1);
    fC=p2D_C-LHS;
    %
    % Bisection tests
    TestL=fL.*fC;
    TestU=fU.*fC;    %
    indTestL=find(TestL>mp('0'));
    indTestU=find(TestU>mp('0'));
    %
    % Bisection: reducing the interval for the search
    if(not(isempty(indTestL)))
        LowerBound(indTestL)=Efp2D_GuessForNewton(indTestL);
    end
    if(not(isempty(indTestU)))
        UpperBound(indTestU)=Efp2D_GuessForNewton(indTestU);
    end
    %
    % Upgrading the guess
    Efp2D_GuessForNewton=(LowerBound+UpperBound)/mp('2');
    
end

% Newton loop

% computing the function at the first step
expnum=exp((Ed-ones(size(Ed,1),1,'mp')*Efp2D_GuessForNewton)./Vt);
expden=exp((Ev0-ones(size(Ed,1),1,'mp')*Efp2D_GuessForNewton)./Vt);
lognum=log(mp('1')+expnum);
logden=log(mp('1')+expden);

Upgrade=inf;

NumberOfIterations=0;
% while(max(abs(f./f0))>tolNewton) % check performed on relative function decrease
while(max(abs(Upgrade))>tolNewton && NumberOfIterations<MaxIterNewton) % check performed on the solution
    
    % computing the function at each iteration.
    % "numerator" exponential and terms
    expnum=exp((Ed-ones(size(Ed,1),1,'mp')*Efp2D_GuessForNewton)./Vt);
    lognum=log(mp('1')+expnum);
    ddlognum_Efp2D=expnum./(mp('1')+expnum)./Vt;

    % "denominator" exponential and terms
    expden=exp((Ev0-ones(size(Ed,1),1,'mp')*Efp2D_GuessForNewton)./Vt);
    logden=log(mp('1')+expden);
    ddlogden_Efp2D=expden./(mp('1')+expden)./Vt;
    
    % Test, debug case: no "denominator" terms
    % f=LHS-sum(Nv2D.*lognum,1);
    % df=sum(Nv2D.*ddlognum_Efp2D,1);
    
    % Functional to be minimized: definition (f = LHS - Nv2D*ferdr2D_0)
    % f=LHS-sum(Nv2D.*(lognum-logden),1);
    f=LHS-sum(Nv2D.*(lognum),1);
    % df=sum(Nv2D.*(ddlognum_Efp2D-ddlogden_Efp2D),1);
    df=sum(Nv2D.*(ddlognum_Efp2D),1);
    
    Upgrade=-f./df; % Newton upgrade
    Efp2D_GuessForNewton=Efp2D_GuessForNewton+Upgrade; % upgrading solution
    NumberOfIterations=NumberOfIterations+1;

end

if(NumberOfIterations>=MaxIterNewton)
    disp(['Return from invferdr_p2D since the maximum number of iterations has been reached'])
    disp(['Maximum last upgrade: ',num2str(max(Upgrade))])
end

Efp2D=Efp2D_GuessForNewton;

return