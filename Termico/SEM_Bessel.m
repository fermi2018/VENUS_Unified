function [f,fp,fs,x,wx,Nfun,MI,MId,MI1,fiX]=SEM_Bessel(Rv,Nv,off,m)
%tic

Rmax=sum(Rv);
Rint=cumsum(Rv);
if nargin<4
 m=0;
end 
define_chiprimo_value %-- loading Bessel functions derivatives zeros

OCThr=1e-10; %-- threshold ortonormalità

Nz=5*Nv;
%Nz=[51 61 23];


Rint=[0 Rint];

for iint=1:length(Nv)
%-- Patch 2
 PatchInfo(iint).Ntot=Nv(iint);
 PatchInfo(iint).L=Rv(iint);
 if m==0
  vk=[0, chiprimo_value(m+1,:)];
 else
   vk=chiprimo_value(m+1,:);
 end
 PatchInfo(iint).k=vk([1:PatchInfo(iint).Ntot]+off(iint)).'./Rmax;
 [nodes,weights]=quadad('legen',1,Nz(iint),0); % Gauss-Legendre nodes
 [u,I]=sort(nodes); % sorting nodes ...
 wu=weights(I); % ... and weights ...
 DxDu=PatchInfo(iint).L/2;
 PatchInfo(iint).wx=DxDu*wu;
 PatchInfo(iint).x=PatchInfo(iint).L/2*u+PatchInfo(iint).L/2+Rint(iint);

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Synthesis of the basis functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
xRight=0;
for indPatch=1:length(PatchInfo)

    Ntot=PatchInfo(indPatch).Ntot;
    k=PatchInfo(indPatch).k;
    Lpatch=PatchInfo(indPatch).L;
    x=PatchInfo(indPatch).x;
    x_L=xRight;
    x_R=xRight+Lpatch;
    
    [f,dfdx,dfdx2]=f_EvalPiecewiseBesselFunctions(m,k,x);
    [f_L,dfdx_L,dfdx2_L]=f_EvalPiecewiseBesselFunctions(m,k,x_L);
    [f_R,dfdx_R,dfdx2_R]=f_EvalPiecewiseBesselFunctions(m,k,x_R);
    
    PatchInfo(indPatch).k=k;
    
    PatchInfo(indPatch).f=f;
    PatchInfo(indPatch).dfdx=dfdx;
    PatchInfo(indPatch).dfdx2=dfdx2;
 
    PatchInfo(indPatch).f_L=f_L;
    PatchInfo(indPatch).dfdx_L=dfdx_L;
    PatchInfo(indPatch).dfdx2_L=dfdx2_L;
    
    PatchInfo(indPatch).f_R=f_R;
    PatchInfo(indPatch).dfdx_R=dfdx_R;
    PatchInfo(indPatch).dfdx2_R=dfdx2_R;
    
    xRight=xRight+Lpatch;

    PatchInfo(indPatch).Nfun=PatchInfo(indPatch).Ntot;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orthonormality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for indPatch=1:length(PatchInfo)

    f=PatchInfo(indPatch).f;
    x=PatchInfo(indPatch).x;
    wx=PatchInfo(indPatch).wx;
    Weight=ones(PatchInfo(indPatch).Ntot,1)*wx;
    X=ones(PatchInfo(indPatch).Ntot,1)*x;
    Mpatch=conj(f.*Weight.*X)*f.';
    
    [U,S,V]=svd(Mpatch);
    S=diag(S);
    I=find(S>=OCThr);
    Nfun=PatchInfo(indPatch).Nfun;
    U=U(:,I)./sqrt((ones(Nfun,1)*S(I).'));
    PatchInfo(indPatch).coeff_mn=U;
    PatchInfo(indPatch).Nfun=size(PatchInfo(indPatch).coeff_mn,2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Continuity conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(length(PatchInfo)>1)
Ntot=0;
    for indPatch=1:length(PatchInfo)
        PatchInfo(indPatch).indCC=Ntot+(1:PatchInfo(indPatch).Nfun);
        Ntot=Ntot+PatchInfo(indPatch).Nfun;
    end
    fCC=[];
    
    for indCC=1:(length(PatchInfo)-1) %-- Loop on the interfaces
        fCCRow=zeros(1,Ntot);
        fCC_L=PatchInfo(indCC).f_R;
        fCC_R=PatchInfo(indCC+1).f_L;

        fCCRow(PatchInfo(indCC).indCC)=(PatchInfo(indCC).coeff_mn.'*fCC_L).';
        fCCRow(PatchInfo(indCC+1).indCC)=-(PatchInfo(indCC+1).coeff_mn.'*fCC_R).';
        fCC=[fCC;fCCRow];
        
    end
    
    [U,S,V]=svd(fCC);
    S=diag(S);
    indCC=find(abs(S./S(1))<=1e-12);
    if isempty(indCC)==1
        indCC=length(S)+1;
    else
        indCC=indCC(1);
    end
    coeff_mn_CC=V(:,indCC:end);
    
    CountFun=0;
    for indPatch=1:length(PatchInfo)
        PatchInfo(indPatch).coeff_mn=PatchInfo(indPatch).coeff_mn*coeff_mn_CC(PatchInfo(indPatch).indCC,:);
        PatchInfo(indPatch).Nfun=size(PatchInfo(indPatch).coeff_mn,2); % upgrade Ntot for C.C.
    end
else % 1 patch design
    PatchInfo(end).coeff_mn=eye(PatchInfo(end).Ntot);
    PatchInfo(end).Nfun=size(PatchInfo(end).coeff_mn,2);
end

%-- calcoliamoci gli integrali
%kPlot=k(indk);

coeffBuild=zeros(PatchInfo(1).Nfun,1);
x=[];
wx=[];
f=[];
fp=[];
fs=[];

for indPatch=1:length(PatchInfo)
    
    fl=PatchInfo(indPatch).f;
    dfdx=PatchInfo(indPatch).dfdx;
    dfdx2=PatchInfo(indPatch).dfdx2;
    xl=PatchInfo(indPatch).x;
    wxl=PatchInfo(indPatch).wx;
    Weight=ones(PatchInfo(indPatch).Ntot,1)*wxl;
    X=ones(PatchInfo(indPatch).Ntot,1)*xl;
%    Function=besselj(0,kPlot*x);

    Mpatch=conj(fl.*Weight.*X)*fl.';
    Mpatch1=conj(fl.*Weight./X)*fl.';
    KKpatch=conj(dfdx.*Weight.*X)*dfdx.'; %-- è utile? occhio al peso
%    KKpatch=conj(dfdx.*Weight.*X)*dfdx.'; %-- è utile? occhio al peso
    
  
    M(:,:,indPatch)=PatchInfo(indPatch).coeff_mn'*Mpatch*PatchInfo(indPatch).coeff_mn;
    M1(:,:,indPatch)=PatchInfo(indPatch).coeff_mn'*Mpatch1*PatchInfo(indPatch).coeff_mn;
    KK(:,:,indPatch)=PatchInfo(indPatch).coeff_mn'*KKpatch*PatchInfo(indPatch).coeff_mn;    KK(:,:,indPatch)=PatchInfo(indPatch).coeff_mn'*KKpatch*PatchInfo(indPatch).coeff_mn;

    fiX{indPatch}=length(x)+[1:length(xl)];  
    x=[x,PatchInfo(indPatch).x];
    wx=[wx,PatchInfo(indPatch).wx.*xl];
    f=[f,PatchInfo(indPatch).coeff_mn.'*PatchInfo(indPatch).f];
    fp=[fp,PatchInfo(indPatch).coeff_mn.'*PatchInfo(indPatch).dfdx];
    fs=[fs,PatchInfo(indPatch).coeff_mn.'*PatchInfo(indPatch).dfdx2];
  
%    coeffBuild=coeffBuild+PatchInfo(indPatch).coeff_mn'*conj(f.*Weight.*X)*Function.';
    
end

MI=M;
MI1=M1;
MId=KK;
Nfun=size(f,1);
%toc
%'qui cont new', keyboard