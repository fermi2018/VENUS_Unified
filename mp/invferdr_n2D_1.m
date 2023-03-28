function [Efn2D]=invferdr_n2D_1(n2D,Em,Ec0,Nc2D,Vt,LowerBound,UpperBound)
%
% function [Efn2D]=invferdr_n2D(n2D,Em,Ec0,Nc2D,Vt,LowerBound,UpperBound)
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
% n2D (vector): sum of 2D electrons of all bound states
%
% Em (matrix): electrons energy levels; rows: different energy levels;
% columns: different lateral points (to account for different potential
% heights)
%
% Ec0 (matrix): upper energy integration bound (eV)
%
% Nc2D (matrix): 2D Density of States (DoS); it is a matrix since different
% effective masses can be used (1/cm^2)
%
% LowerBound,UpperBound: lower and upper bounds for the bisection algorithm
% for Efn2D (vectors, eV)
%
% Ver. 0.1, 09/03/2016
% Ver. 0.2, 29/03/2016 - introduced automatic guess (with bisection)
%
% Alberto Tibaldi

% Two-step algorithm: first Dichotomy to achieve a guess, then Newton

tolBisection=1e-2; % tolerance for solution via bisection method
tolNewton=1e-13; % tolerance for solution via Newton method
MaxIterNewton=50; % maximum number of iterations for Newton loop

NBound_C=size(Em,1);

LHS=n2D;
Efn2D_GuessForNewton=(LowerBound+UpperBound)/2;
% Setting to inf to enter in the loop at least once
GuessPreviousStep=inf*ones(size(Efn2D_GuessForNewton));

while(max(abs(Efn2D_GuessForNewton-GuessPreviousStep))>tolBisection)
    
    GuessPreviousStep=Efn2D_GuessForNewton;
    % Calculation of the upper bound for the current iteration
    xnum=(ones(NBound_C,1)*UpperBound-Em)./Vt; 
    xden=(ones(NBound_C,1)*UpperBound-Ec0)./Vt;
    [n2D_UB,dum] = mpferdr2D_0(xnum,xden);
    n2D_UB=sum(Nc2D.*n2D_UB,1);
    fU=n2D_UB-LHS;
    %
    % Calculation of the lower bound for the current iteration
    xnum=(ones(NBound_C,1)*LowerBound-Em)./Vt; 
    xden=(ones(NBound_C,1)*LowerBound-Ec0)./Vt;
    [n2D_LB,dum] = mpferdr2D_0(xnum,xden);
    n2D_LB=sum(Nc2D.*n2D_LB,1);
    fL=n2D_LB-LHS;
    %
    % Calculation of the guess for the current iteration
    xnum=(ones(NBound_C,1)*Efn2D_GuessForNewton-Em)./Vt; 
    xden=(ones(NBound_C,1)*Efn2D_GuessForNewton-Ec0)./Vt;
    [n2D_C,dum] = mpferdr2D_0(xnum,xden);
    n2D_C=sum(Nc2D.*n2D_C,1);
    fC=n2D_C-LHS;
    %
    % Bisection tests
    TestL=fL.*fC;
    TestU=fU.*fC;    %
    indTestL=find(TestL>0);
    indTestU=find(TestU>0);
    %
    % Bisection: reducing the interval for the search
    if(not(isempty(indTestL)))
        LowerBound(indTestL)=Efn2D_GuessForNewton(indTestL);
    end
    if(not(isempty(indTestU)))
        UpperBound(indTestU)=Efn2D_GuessForNewton(indTestU);
    end
    %
    % Upgrading the guess
    Efn2D_GuessForNewton=(LowerBound+UpperBound)/2;
    
end

% Newton loop
% computing the function at the first step
expnum=exp((ones(NBound_C,1)*Efn2D_GuessForNewton-Em)./Vt);
expden=exp((ones(NBound_C,1)*Efn2D_GuessForNewton-Ec0)./Vt);
lognum=log1p(1+expnum);
logden=log1p(1+expden);

% Setting to inf to enter in the loop at least once
Upgrade=inf;

NumberOfIterations=0;
% while(max(abs(f./f0))>tolNewton) % check performed on relative function decrease
while(max(abs(Upgrade))>tolNewton && NumberOfIterations<MaxIterNewton) % check performed on the solution
    
    % computing the function at each iteration.
    % "numerator" exponential and terms
    expnum=exp((ones(NBound_C,1)*Efn2D_GuessForNewton-Em)./Vt);
    lognum=log1p(1+expnum);
    ddlognum_Efn2D=-expnum./(1+expnum)./Vt;

    % "denominator" exponential and terms
    expden=exp((ones(NBound_C,1)*Efn2D_GuessForNewton-Ec0)./Vt);
    logden=log1p(1+expden);
    ddlogden_Efn2D=-expden./(1+expden)./Vt;
    
    % Test, debug case: no "denominator" terms
    % f=LHS-sum(Nc2D.*lognum,1);
    % df=sum(Nc2D.*ddlognum_Efn2D,1);
    
    % Functional to be minimized: definition (f = LHS - Nc2D*ferdr2D_0)
    % f=LHS-sum(Nc2D.*(lognum-logden),1);
    f=LHS-sum(Nc2D.*(lognum),1);
    % df=sum(Nc2D.*(ddlognum_Efn2D-ddlogden_Efn2D),1);
    df=sum(Nc2D.*(ddlognum_Efn2D),1);
    
    Upgrade=-f./df; % Newton upgrade
    Efn2D_GuessForNewton=Efn2D_GuessForNewton+Upgrade; % upgrading solution
    NumberOfIterations=NumberOfIterations+1;
    
end

if(NumberOfIterations>=MaxIterNewton)
    disp(['Return from invferdr_n2D since the maximum number of iterations has been reached'])
    disp(['Maximum last upgrade: ',num2str(max(Upgrade))])
end

Efn2D=Efn2D_GuessForNewton;

return