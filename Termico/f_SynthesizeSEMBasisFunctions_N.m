
function [PatchInfo]=f_SynthesizeSEMBasisFunctions(PatchInfo)
xRight=0;
for indPatch=1:length(PatchInfo)
    Lpatch=PatchInfo(indPatch).L;
    x=PatchInfo(indPatch).x;
    PatchInfo(indPatch).dudx=2/Lpatch;
    u=PatchInfo(indPatch).dudx*(x-xRight)-1;
    vn=[0:1:PatchInfo(indPatch).Ntot];
    PatchInfo(indPatch).Ntot=PatchInfo(indPatch).Ntot+1; %-- Ntot+1 poly.
    
    [f,dfdu,d2fdu2]=f_EvalLegendrePolynomials(vn,u); % internal points
    [f_L,dfdu_L,d2fdu2_L]=f_EvalLegendrePolynomials(vn,-1); % left B.C.
    [f_R,dfdu_R,d2fdu2_R]=f_EvalLegendrePolynomials(vn,1); % right B.C.
    
%    'sei qua'
%    keyboard
    PatchInfo(indPatch).f=f;
    PatchInfo(indPatch).dfdx=dfdu*PatchInfo(indPatch).dudx;
    PatchInfo(indPatch).d2fdx2=d2fdu2*PatchInfo(indPatch).dudx^2;
 
    PatchInfo(indPatch).f_L=f_L;
    PatchInfo(indPatch).dfdx_L=dfdu_L*PatchInfo(indPatch).dudx;
    PatchInfo(indPatch).d2fdx2_L=d2fdu2_L*PatchInfo(indPatch).dudx^2;
    
    PatchInfo(indPatch).f_R=f_R;
    PatchInfo(indPatch).dfdx_R=dfdu_R*PatchInfo(indPatch).dudx;
    PatchInfo(indPatch).d2fdx2_R=d2fdu2_R*PatchInfo(indPatch).dudx^2;
    
    xRight=xRight+Lpatch;
    PatchInfo(indPatch).Nfun=PatchInfo(indPatch).Ntot;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary and continuity conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(length(PatchInfo)>1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Left boundary condition
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BCInfo=PatchInfo(1).BCInfo{1};
    if(strcmp(BCInfo,'Dirichlet'))
        fBC=[PatchInfo(1).f_L.']; %-- Dirichlet
    elseif(strcmp(BCInfo,'Neumann'))
        fBC=[PatchInfo(1).dfdx_L.']; %-- Neumann
    else
        error('Wrong left boundary condition')
    end
    [U,S,V]=svd(fBC);
    S=diag(S);
    indBC=find(abs(S./S(1))<=1e-12);
    if isempty(indBC)==1
        indBC=length(S)+1;
    else
        indBC=indBC(1);
    end
    PatchInfo(1).coeff_mn=V(:,indBC:end);
    PatchInfo(1).Nfun=size(PatchInfo(1).coeff_mn,2); % upgrade Ntot for B.C.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Right boundary condition
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BCInfo=PatchInfo(end).BCInfo{2};
    if(strcmp(BCInfo,'Dirichlet'))
        fBC=[PatchInfo(end).f_R.']; %-- Dirichlet
    elseif(strcmp(BCInfo,'Neumann'))
        fBC=[PatchInfo(end).dfdx_R.']; %-- Neumann
    else
        error('Wrong right boundary condition')
    end
    [U,S,V]=svd(fBC);
    S=diag(S);
    indBC=find(abs(S./S(1))<=1e-12);
    if isempty(indBC)==1
        indBC=length(S)+1;
    else
        indBC=indBC(1);
    end
    PatchInfo(end).coeff_mn=V(:,indBC:end);
    PatchInfo(end).Nfun=size(PatchInfo(end).coeff_mn,2); % upgrade Ntot for B.C.
    
    Ntot=0;
    for indPatch=1:length(PatchInfo)
        PatchInfo(indPatch).indCC=Ntot+(1:PatchInfo(indPatch).Nfun);
        Ntot=Ntot+PatchInfo(indPatch).Nfun;
    end
    fCC=[];
    
    for indCC=1:(length(PatchInfo)-1) %-- Loop on the interfaces
        CCInfo_L=PatchInfo(indCC).BCInfo{2};
        CCInfo_R=PatchInfo(indCC+1).BCInfo{1};
        %-- Check if the continuity condition is properly indicated
        if(not(strcmp(CCInfo_L,CCInfo_R)))
            error(['Wrong continuity definitions between patches ',num2str(indCC),' and ',num2str(indCC+1)])
        end
        fCCRow=zeros(1,Ntot);
        if(strcmp(CCInfo_L,'Continuity'))
            fCC_L=PatchInfo(indCC).f_R;
            fCC_R=PatchInfo(indCC+1).f_L;
            if(indCC==1)
                fCC_L=PatchInfo(indCC).coeff_mn.'*fCC_L;
            end
            if indCC==(length(PatchInfo)-1)
                fCC_R=PatchInfo(end).coeff_mn.'*fCC_R;
            end
            fCCRow(PatchInfo(indCC).indCC)=fCC_L.';
            fCCRow(PatchInfo(indCC+1).indCC)=-fCC_R.';
            fCC=[fCC;fCCRow];
            
        elseif(strcmp(CCInfo_L,'Derivability'))
            fCCRowFlux=zeros(1,Ntot);
            
            fCC_L=PatchInfo(indCC).f_R;
            fCC_R=PatchInfo(indCC+1).f_L;
            fCCFlux_L=PatchInfo(indCC).dfdx_R.*PatchInfo(indCC).kFlux;
            fCCFlux_R=PatchInfo(indCC+1).dfdx_L.*PatchInfo(indCC+1).kFlux;
            
            if(indCC==1)
                fCC_L=PatchInfo(indCC).coeff_mn.'*fCC_L;
                fCCFlux_L=PatchInfo(indCC).coeff_mn.'*fCCFlux_L;
            end
            if indCC==(length(PatchInfo)-1)
                fCC_R=PatchInfo(end).coeff_mn.'*fCC_R;
                fCCFlux_R=PatchInfo(end).coeff_mn.'*fCCFlux_R;
            end
            fCCRow(PatchInfo(indCC).indCC)=fCC_L.';
            fCCRow(PatchInfo(indCC+1).indCC)=-fCC_R.';
            fCCRowFlux(PatchInfo(indCC).indCC)=fCCFlux_L.';
            fCCRowFlux(PatchInfo(indCC+1).indCC)=-fCCFlux_R.';
            
            fCC=[fCC;fCCRow;fCCRowFlux];
            
        else
            error(['Wrong continuity definitions between patches ',num2str(indCC),' and ',num2str(indCC+1)])
        end
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
    
    for indPatch=1:length(PatchInfo)
        if indPatch==1 || indPatch==length(PatchInfo)
            PatchInfo(indPatch).coeff_mn=PatchInfo(indPatch).coeff_mn*coeff_mn_CC(PatchInfo(indPatch).indCC,:);
        else
            PatchInfo(indPatch).coeff_mn=coeff_mn_CC(PatchInfo(indPatch).indCC,:);
        end
        PatchInfo(indPatch).Nfun=size(PatchInfo(indPatch).coeff_mn,2); % upgrade Ntot for C.C.
    end
else % 1 patch design
    BCInfo=PatchInfo(end).BCInfo;
    if(strcmp(BCInfo{1},'Dirichlet'))
        fBC=[PatchInfo(end).f_L.']; %-- Dirichlet
    elseif(strcmp(BCInfo{1},'Neumann'))
        fBC=[PatchInfo(end).dfdx_L.']; %-- Neumann
    else
        error('Wrong left boundary condition')
    end
    if(strcmp(BCInfo{2},'Dirichlet'))
        fBC=[fBC;PatchInfo(end).f_R.']; %-- Dirichlet
    elseif(strcmp(BCInfo{2},'Neumann'))
        fBC=[fBC;PatchInfo(end).dfdx_R.']; %-- Neumann
    else
        error('Wrong left boundary condition')
    end
    [U,S,V]=svd(fBC);
    S=diag(S);
    indBC=find(abs(S./S(1))<=1e-12);
    if isempty(indBC)==1
        indBC=length(S)+1;
    else
        indBC=indBC(1);
    end
    PatchInfo(end).coeff_mn=V(:,indBC:end);
    PatchInfo(end).Nfun=size(PatchInfo(end).coeff_mn,2); % upgrade Ntot for B.C.
end

return