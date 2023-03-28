function [f,fp,fs,Nfun,MI,MIs,MId]=SEM_zeta(Lv,Nv,zi,kv)
%-- Questo programma risolve il problema agli autovalori della linea di
%   trasmissione, lunga pigreco. Se si mettono entrambe condizioni di
%   Neumann o entrambe di Dirichlet, le radici degli autovalori sono
%   multiple di 1: {1, 2, 3,...}

%-- DISTINGUERE IL CASO A 1 INTERVALLO DA TUTTI GLI ALTRI (LE DUE SONO
%CONDIZIONI AL CONTORNO !!!!)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters of the program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-- General parameters
%Options.Orthogonality=1; % if 1 functions are orthogonal DA IMPLEMENTARE
%Options.BCTHreshold=1e-12; % boundary conditions SVD threshold
%Options.CCTHreshold=1e-12; % boundary conditions SVD threshold
%-- Patch 1


if nargin<4
kv=ones(size(Lv));
end

%'kv', keyboard
icontder=0;
  if icontder==0
   condiz='Continuity';
  else 
   condiz='Derivability';
  end 

for ns=1:length(Lv)
 PatchInfo(ns).Ntot=Nv(ns);
 PatchInfo(ns).L=Lv(ns);
 PatchInfo(ns).BCInfo={condiz,condiz};
 PatchInfo(ns).kFlux=kv(ns); % for flux continuity (kFlux=1 means Derivability)
end

%PatchInfo(2).BCInfo={'Derivability',condiz};

%PatchInfo(1).BCInfo={'Dirichlet','Derivability'};
PatchInfo(1).BCInfo={'Dirichlet',condiz};
PatchInfo(end).BCInfo={condiz,'Neumann'};

if length(Lv)==1
 PatchInfo(1).BCInfo={'Dirichlet','Neumann'};
end

z=zi;
if size(zi,2)<size(zi,1)
 z=zi';
end 

zi=z(1);
zu=Lv(1);
Lvs=cumsum(Lv);
fi=find(z>=zi & z<=zu);
PatchInfo(1).x=z(fi);
zi=zu;
for ns=2:length(Lvs)
 zu=Lvs(ns);
 fi=find(z>zi & z<=zu);
 PatchInfo(ns).x=z(fi);
 zi=zu;
end

%'kv', keyboard
%'centro', keyboard

[PatchInfo]=f_SynthesizeSEMBasisFunctions_N(PatchInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test problem: second order transmission line equation (eigvl. problem)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nfun=PatchInfo(1).Nfun; 

f=[]; 
fp=[]; 
fs=[]; 

[D1,D2]=f_EvalAnalyticLegendreIntegrals(2*PatchInfo(1).Nfun);

for indPatch=1:length(PatchInfo)
    vn=[0:1:PatchInfo(indPatch).Ntot-1];
    Mi=diag(2./(2.*vn+1))./PatchInfo(indPatch).dudx; % Gram matrix
    Ki=(Mi)*D2(1:PatchInfo(indPatch).Ntot,1:PatchInfo(indPatch).Ntot)*PatchInfo(indPatch).dudx^2; % Integrale Der2 * fun
    BoundaryTerm=PatchInfo(indPatch).f_R*PatchInfo(indPatch).dfdx_R.'-PatchInfo(indPatch).f_L*PatchInfo(indPatch).dfdx_L.';
    KKi=(BoundaryTerm-Ki); % Integrale Der1*Der1
    
    PatchInfo(indPatch).M=PatchInfo(indPatch).coeff_mn'*Mi*PatchInfo(indPatch).coeff_mn;
    PatchInfo(indPatch).K=PatchInfo(indPatch).coeff_mn'*Ki*PatchInfo(indPatch).coeff_mn;
    PatchInfo(indPatch).KK=PatchInfo(indPatch).coeff_mn'*KKi*PatchInfo(indPatch).coeff_mn;
    
    M(:,:,indPatch)=PatchInfo(indPatch).M;
    K(:,:,indPatch)=PatchInfo(indPatch).K;
    KK(:,:,indPatch)=PatchInfo(indPatch).KK;
    f=[f, PatchInfo(indPatch).coeff_mn.'*PatchInfo(indPatch).f];
    fp=[fp, PatchInfo(indPatch).coeff_mn.'*PatchInfo(indPatch).dfdx];
    fs=[fs, PatchInfo(indPatch).coeff_mn.'*PatchInfo(indPatch).d2fdx2];
    
end

MI=M;
MIs=K;
MId=KK;

%'cont sem', keyboard