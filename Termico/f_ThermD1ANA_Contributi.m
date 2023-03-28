function deltaT = f_ThermD1ANA_Contributi(mesh,condzTe,Heat,fattore_correttivo)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of pointers to nodes and equations and related masks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nAir = 0; 
L_air = 0;%200e-7;
nn=mesh.nny-2+nAir;         % pointer to Heat equation
nodeAir = L_air/nAir*ones(1,nAir);
nodeAir = mesh.ygrid(end-1)+cumsum(nodeAir);
node = [mesh.ygrid(2:end-1), nodeAir];  % remove the contacts
%
% Defining indexes for assembling FEM and FD matrices
in=1:nn; % node index
in1=in(1:end-1); in2=in(2:end); % left, right nodes of each element
inr=[in1 in1 in2 in2]; inc=[in1 in2 in1 in2]; in12=[in1 in2];
maskr=true(1,nn); 
mask_in12=maskr(in12); in12=in12(mask_in12);
%
% Defining weighting areas to compute integrals when assembling matrices
mesh.Le = diff(node);
mesh.Lp = zeros(1,nn);
mesh.Lp(1:nn-1) = mesh.Le/2;
mesh.Lp(2:nn) = mesh.Lp(2:nn) + mesh.Le/2;
Lp1=mesh.Lp(in1)/2; Lp1(1)=2*Lp1(1);
Lp2=mesh.Lp(in2)/2; Lp2(end)=2*Lp2(end);
%
TotalHeat = ones(1,nn);
TotalHeat(1:nn) = Heat;
%
iterCount=0;
itermaxTherm=1;
TolTherm=1e-3;
er=inf; 
deltaT_old=0;

while (iterCount<itermaxTherm && er>TolTherm)
    
    % System assembling
    A_1e1e=+condzTe./mesh.Le;
    A_1e2e=-condzTe./mesh.Le;
    A_2e1e=-condzTe./mesh.Le;
    A_2e2e=+condzTe./mesh.Le;
    MM=[A_1e1e, A_1e2e, A_2e1e, A_2e2e];
    
    Amat = sparse(inr,inc,MM,nn,nn);
    
    TotalHeat1 = TotalHeat(in1);
    TotalHeat2 = TotalHeat(in2);
    MM=[TotalHeat1.*Lp1, TotalHeat2.*Lp2];
    
    qvet=sparse(in12,1,MM(mask_in12),nn,1);
    %
    % Enforce the boundary condition
    qvet(1)= 0;
    %
    % Boundary conditions: Heat sink temperature (w.r.t. device)
    Amat(1,:) = 0;
    Amat(1,1) = 1; % line contact on the left
    %
    % Solve the system
    deltaT = Amat\(fattore_correttivo*qvet);
            
    er=abs(1-max(deltaT_old)/max(deltaT));
    if isnan(er)
        er=0;
    end
    
    deltaT_old=deltaT;
    iterCount=iterCount+1;
    
end

% DeltaT = deltaT(1:mesh.nn);

