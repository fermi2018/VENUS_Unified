function [rvet_inc]=evalMatrices(mesh,mode,uvet)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading physical constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
qel=1.6021766208e-019; % Elementary charge, C
m0=9.10938291e-31; % Electron mass, kg
kB=1.3806488e-23; % Boltzmann constant, J/K
h=6.626070040e-34; % Planck constant, J*s
Clight=2.99792458e+10; % Speed of light, cm/s
mu0=4*pi*1e-9; % Magnetic permeability, H/cm
eps0=1/(mu0*Clight^2); % Dielectric permittivity, F/cm
%
% Temperature (defined on nodes; could be an input, in general)
T=mode.Temp0*ones(1,mesh.nn); % (K)
Vt=kB.*T./qel; % electrical equivalent in temperature, V
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of pointers to nodes and equations and related masks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nn=mesh.nn;         % pointer to electron continuity equations
pp=2*mesh.nn;       % pointer to hole continuity equations
%
% Masks for excluding contact nodes from equations
maskr=true(1,nn); maskr([1,end])=0; % Poisson equation
masks=true(1,nn); masks([1,end])=0; % Continuity equations
maskc=false(1,nn); maskc(mesh.iline)=1; % Current equation
%
% Defining indexes for assembling FEM and FD matrices
in=1:mesh.nn; % node index
in1=in(1:end-1); in2=in(2:end); % left, right nodes of each element
inr=[in1 in1 in2 in2]; % assembling rows, FEM mass/stiffness matrices
inc=[in1 in2 in1 in2]; % assembling columns, FEM mass/stiffness matrices
iic = inr; jjc = inc;  % indexes for current equations
in12=[in1 in2];
%
% Applying masks to indexes. Since mask_inr is used to select the elements
% of the assembled matrices MM, this is used for both inr and inc.
mask_inr=maskr(inr);
mask_in12=maskr(in12); in12=in12(mask_in12);
mask_iic=maskc(iic);
%
inr=inr(mask_inr); inc=inc(mask_inr);
iic=iic(mask_iic); jjc=jjc(mask_iic);
%
% Defining weighting areas to compute integrals when assembling matrices
Lp1=mesh.Lp(in1)/2; Lp1(1)=2*Lp1(1);
Lp2=mesh.Lp(in2)/2; Lp2(end)=2*Lp2(end);
% debug: check Lp1,Lp2: integral of the constant function on the domain
%        this is expected to be equal to the length of the domain
% sum(diag(sparse([in1 in2],[in1 in2],[Lp1 Lp2],nn,nn)),1)
% return
% end debug
%
% Scaling factor
C0=mode.C0;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Off-equilibrium: drift-diffusion analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zero-field mobility (defined on elements)
[mobn0,mobp0]=fmp_EvalDopingMobility(mode,mesh,in1,in2);
%
neq=nn+nn+nn; % Poisson + continuity electrons + continuity holes

if(mode.Quantum==1) % introduce quantum corrections and photon rate eqs.
    
    NQW=length(mesh.iMQW);
    n2Di=[]; p2Di=[];
    n2Ds=[]; p2Ds=[];
    
    for indQW=1:NQW
        n2Di=[n2Di,mode.n2Di{indQW}];
        p2Di=[p2Di,mode.p2Di{indQW}];
    end
    
    vv = neq;       % pointer to first quantum electron equation
    ww=neq+NQW;     % pointer to first quantum hole equation
    neq=neq+2*NQW;  % Poisson + bulk + quantum transport equations
    
    if mode.indV==1
    else
        for indQW=1:NQW
            n2Ds=[n2Ds,mode.n2Dc{indQW}];
            p2Ds=[p2Ds,mode.p2Dc{indQW}];
        end
    end
% 
    if(mode.Optical==1)
        ss=neq;
        NModes=length(mode.Lm);
        neq=neq+NModes;
    end
    
elseif(mode.Quantum==0)
    neq=nn+nn+nn; % Poisson + continuity electrons + continuity holes
end
%
% Adding voltage/current equations
qq = neq;
rr = neq+1;
neq = neq+2;
iic = (rr+1)*ones(size(iic));
%
iv=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEBUG: to let the code start directly from the last voltage, uncomment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conversion of nodal quantities to elements
epsi=(mesh.epsi(in1)+mesh.epsi(in2))/2;
%
uvet=uvet;
v0_dd=mode.Vbias(iv);
if(mode.Quantum==1 && mode.FlatBandQW==0)
    [mesh]=fmp_SynchronizeMacro(geom,mesh,mode);
    % Conversion of nodal quantities to elements
    epsi=(mesh.epsi(in1)+mesh.epsi(in2))/2;
end
%
if(mode.Quantum==1)
    for indQW=1:NQW
        mode.n2Dc{indQW}=abs(uvet(vv+indQW)).*C0; % abs is a kludge :-\
        mode.p2Dc{indQW}=abs(uvet(ww+indQW)).*C0; % abs is a kludge :-\
    end
    if(mode.Optical==1)
        %                 Pst=abs(uvet(ss+1:end)); % abs is a kludge :-\
        NUOVA_Gain
    end
end
%
% loading solution @ previous step
phi=reshape(uvet(1:nn),1,nn);
elec=reshape(uvet(nn+1:pp),1,nn);
hole=reshape(uvet(pp+1:pp+nn),1,nn);
%
efield=-diff(phi)./mesh.Le;
%
if(mode.Quantum==1)
    for indQW=1:NQW
        n2Dc{indQW}=abs(uvet(vv+indQW)); % abs is a kludge :-\
        p2Dc{indQW}=abs(uvet(ww+indQW)); % abs is a kludge :-\
    end
    if(mode.Optical==1)
        Pst=abs(uvet(ss+1)); % abs is a kludge :-\
    end
end
%
Ec = -phi + mesh.Eg/2 - Vt/2.*log(mesh.Nv./mesh.Nc) + mesh.phi_r;
Ev = Ec - mesh.Eg;
if(strcmp(mode.Statistics,'Boltzmann'))
    EFn = Ec + Vt.*log(elec*C0./mesh.Nc);
    elecB=elec;
    delec=elec;
    gamma_n=ones(size(elec));
    dgamma_n=0*ones(size(elec));
    dloggamma_n = 0*ones(size(elec));
    %
    EFp = Ev - Vt.*log(hole*C0./mesh.Nv);
    holeB=hole;
    dhole=hole;
    gamma_p=ones(size(hole));
    dgamma_p=0*ones(size(hole));
    dloggamma_p = 0*ones(size(hole));
    %
elseif(strcmp(mode.Statistics,'Fermi'))
    tmp = abs(elec*C0./mesh.Nc); % kludge: forcing charge to be positive
    EFn = Ec + Vt.*mpinvferdr(tmp,mode.TolConv_invferdr);
    elec=mesh.Nc.*mpferdr((EFn-Ec)./Vt,1/2)/C0; % Recompute with Fermi
    delec=mesh.Nc.*mpferdr((EFn-Ec)./Vt,-1/2)/C0;
    elecB=mesh.Nc.*exp((EFn-Ec)./Vt)/C0; % Boltzmann
    gamma_n=elec./elecB;
    dgamma_n=(1-elec./delec)./elecB;
    ilim=(tmp==0);
    dgamma_n(ilim)=-sqrt(2)./4./mesh.Nc(ilim)/C0;
    dloggamma_n = dgamma_n./gamma_n;
    %
    tmp = abs(hole*C0./mesh.Nv); % kludge: forcing charge to be positive
    EFp = Ev - Vt.*mpinvferdr(tmp,mode.TolConv_invferdr);
    hole=mesh.Nv.*mpferdr((Ev-EFp)./Vt,1/2)/C0; % Recompute with Fermi
    dhole=mesh.Nv.*mpferdr((Ev-EFp)./Vt,-1/2)/C0;
    holeB=mesh.Nv.*exp((Ev-EFp)./Vt)/C0; % Boltzmann
    gamma_p=hole./holeB;
    dgamma_p=(1-hole./dhole)./holeB;
    ilim=(tmp==0);
    dgamma_p(ilim)=-sqrt(2)./4./mesh.Nv(ilim)/C0;
    dloggamma_p = dgamma_p./gamma_p;
else
    error('Please, choose one between Boltzmann and Fermi! ç_ç')
end
dloggamma_n1=dloggamma_n(in1); dloggamma_n2=dloggamma_n(in2);
dloggamma_p1=dloggamma_p(in1); dloggamma_p2=dloggamma_p(in2);
%
phi1=phi(in1); phi2=phi(in2);
elec1=elec(in1); elec2=elec(in2);
hole1=hole(in1); hole2=hole(in2);
%
if(strcmp(mode.Ionization,'Full'))
    dop_dp=mesh.dop_d;
    dop_am=mesh.dop_a;
    ddop_dp=0*ones(size(dop_dp));
    ddop_am=0*ones(size(dop_am));
    %
elseif(strcmp(mode.Ionization,'Incomplete'))
    gD=2; n1 = mesh.Nc.*exp(-mesh.DeltaEd./Vt);
    dop_dp  =   mesh.dop_d./( 1 + gD.*elec*C0./(n1.*gamma_n)); % 1/cm^3
    ddop_dp = - mesh.dop_d./((1 + gD.*elec*C0./(n1.*gamma_n)).^2) ...
        .*gD.*((n1.*gamma_n - elec.*n1.*dgamma_n)./((n1.*gamma_n).^2));
    dop_dp(mesh.dop_d==0)=0;
    ddop_dp(mesh.dop_d==0)=0;
    %
    gA = 4; p1 = mesh.Nv.*exp(-mesh.DeltaEa./Vt);
    dop_am  =    mesh.dop_a./( 1 + gA.*hole*C0./(p1.*gamma_p)); % ionized donor density, 1/cm^3
    ddop_am =  - mesh.dop_a./((1 + gA.*hole*C0./(p1.*gamma_p)).^2) ...
        .*gA.*((p1.*gamma_p - hole.*p1.*dgamma_p)./((p1.*gamma_p).^2));
    dop_am(mesh.dop_a==0)=0;
    ddop_am(mesh.dop_a==0)=0;
    %
else
    error('My friend, ionization can be just Full or Incomplete ç_ç')
end
mode.dop_dp = dop_dp;
mode.dop_am = dop_am;
%
% matrix initialization
Jmat0=sparse(neq,neq);
Kmat0=sparse(neq,neq);
tvet=sparse(neq,1);
%
% Computing field-dependent mobility
[Dn,Dp,dDndphi1,dDndphi2,dDpdphi1,dDpdphi2] = f_EvalFieldMobility(mode,mesh,efield,mobn0,mobp0,Vt,in1,in2);
Dn=Dn;Dp=Dp;dDndphi1=dDndphi1;dDndphi2=dDndphi2;dDpdphi1=dDpdphi1;dDpdphi2=dDpdphi2;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assembling Poisson equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second derivative term
M11=epsi./mesh.Le; M12=-epsi./mesh.Le; M21=M12; M22=epsi./mesh.Le;
MM=[M11 M12 M21 M22]/C0;
Kmat0=Kmat0+sparse(inr,inc,MM(mask_inr),neq,neq);
%
% Doping charge
dop=dop_dp-dop_am;
dop1=dop_dp(in1)-dop_am(in1)/C0;
dop2=dop_dp(in2)-dop_am(in2)/C0;
MM=-qel.*[dop1.*Lp1, dop2.*Lp2];
tvet=tvet+sparse(in12,1,MM(mask_in12),neq,1);
%
% Polarization charge
if(mode.PolCharges~=0)
    if(strcmp(mode.GaFacePosition,'left')==1)
        PolFraction=+mode.PolCharges;
    elseif(strcmp(mode.GaFacePosition,'right')==1)
        PolFraction=-mode.PolCharges;
    end
    rhopol=-PolFraction.*mesh.rhopol.*qel/C0;
    MM=[rhopol(in1).*Lp1 rhopol(in2).*Lp2];
    tvet=tvet+sparse(in12,1,MM(mask_in12),neq,1);
else
    PolFraction=0;
end
%
% Electron charge
MM=qel.*[Lp1, Lp2];
Kmat0=Kmat0+sparse(in12,nn+in12,MM(mask_in12),neq,neq);
%
% Hole charge
MM=-qel.*[Lp1, Lp2];
Kmat0=Kmat0+sparse(in12,pp+in12,MM(mask_in12),neq,neq);
%
% Boundary conditions
Jmat2 = sparse(neq,neq);
Jmat2 = Jmat2 + sparse(qq+1,qq+1,1,neq,neq);
%
% Circuit equation: impedance
Zmat = sparse(mode.Zmat);
Jmat2(qq+1,rr+1) = Zmat;
%
% Enforcing neutrality condition at contacts
ic=find(not(maskr));
Kmat0=Kmat0+sparse(ic,ic,1,neq,neq);
tvet(ic)=-mode.phi_neutr(ic);
%
% Correcting line contact with voltage equation
il=mesh.iline;
Kmat0 = Kmat0 + sparse(il,qq+1,-1,neq,neq);
tvet(qq+1)=-v0_dd;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assembling electron continuity equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Treatment of heterointerfaces
offset_n=mesh.affinity+Vt.*log(mesh.Nc);
offset_n1=offset_n(in1); offset_n2=offset_n(in2);
%
% Fermi statistics correction
fermiCorr_n1=Vt(in1).*log(gamma_n(in1));
fermiCorr_n2=Vt(in2).*log(gamma_n(in2));
%
% Evaluating Bernoulli and its derivatives w.r.t. potential
% Argument of Bernoulli function
delta12=(phi1+offset_n1+fermiCorr_n1)./Vt(in1)-(phi2+offset_n2+fermiCorr_n2)./Vt(in2);
%
% Bernoulli functions
B12=mpbern(delta12); dB12=mpdbern(delta12)./Vt(in1);
B21=mpbern(-delta12); dB21=mpdbern(-delta12)./Vt(in2);
M11=qel.*Dn./mesh.Le.*B12;
M12=-qel.*Dn./mesh.Le.*B21;
M21=-qel.*Dn./mesh.Le.*B12;
M22=qel.*Dn./mesh.Le.*B21;
MM=[M11 M12 M21 M22];
Kmat0=Kmat0+sparse(nn+inr,nn+inc,MM(mask_inr),neq,neq);
Kmat0=Kmat0+sparse(   iic,nn+jjc,-C0*MM(mask_iic),neq,neq);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assembling hole continuity equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Treatment of heterointerfaces
offset_p=mesh.affinity-Vt.*log(mesh.Nv)+mesh.Eg;
offset_p1=offset_p(in1); offset_p2=offset_p(in2);
%
% Fermi statistics correction
fermiCorr_p1=Vt(in1).*log(gamma_p(in1));
fermiCorr_p2=Vt(in2).*log(gamma_p(in2));
%
% Evaluating Bernoulli and its derivatives w.r.t. potential
% Argument of Bernoulli function
delta12=(phi1+offset_p1-fermiCorr_p1)./Vt(in1)-(phi2+offset_p2-fermiCorr_p2)./Vt(in2);
%
% Bernoulli functions
B12=mpbern(delta12); dB12=mpdbern(delta12)./Vt(in2);
B21=mpbern(-delta12); dB21=mpdbern(-delta12)./Vt(in1);
M11=qel.*Dp./mesh.Le.*B21;
M12=-qel.*Dp./mesh.Le.*B12;
M21=-qel.*Dp./mesh.Le.*B21;
M22=qel.*Dp./mesh.Le.*B12;
MM=[M11 M12 M21 M22];
Kmat0=Kmat0+sparse(pp+inr,pp+inc,MM(mask_inr),neq,neq);
Kmat0=Kmat0+sparse(   iic,pp+jjc,C0*MM(mask_iic),neq,neq);
%
% assembling current equations
Kmat0=Kmat0+sparse(rr+1,rr+1,-1,neq,neq);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SRH, radiative, Auger generation-recombination terms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R=0*ones(1,mesh.nn); dRn=0*ones(1,mesh.nn); dRp=0*ones(1,mesh.nn);
iq=2:mesh.nn-1;
% TEMP: for Boltzmann distribution
gamman=ones(1,mesh.nn); gammap=gamman; dGn=0*gamman; dGp=dGn;
n  =  elec(iq); p  =  hole(iq);
np = (C0.*n.*C0.*p-gamman(iq).*gammap(iq).*mesh.nint(iq).^2);
%         mode.RSRH = zeros(1,mesh.nn);
%         mode.Rrad = zeros(1,mesh.nn);
%         mode.RAug = zeros(1,mesh.nn);
for indGR=1:length(mode.GR)
    if strcmp(mode.GR{indGR},'SRH')==1
        % SRH recombination model
        n1 = gamman(iq).*mesh.nint(iq).*exp(mesh.Etrap(iq)./Vt(iq));
        p1 = gammap(iq).*mesh.nint(iq).*exp(-mesh.Etrap(iq)./Vt(iq));
        denSRH = mesh.taup(iq).*(C0*n+n1)+mesh.taun(iq).*(C0*p+p1);
        R(iq) = R(iq) + np./denSRH;
        %
    elseif strcmp(mode.GR{indGR},'rad')==1
        % Radiative recombination model
        R(iq) = R(iq) + mesh.brad(iq).*np;
        %
    elseif strcmp(mode.GR{indGR},'Auger')==1
        % Auger recombination model
        R(iq) = R(iq) + (mesh.Cnnp(iq).*C0.*n+mesh.Cppn(iq).*C0.*p).*np;
        %
    elseif strcmp(mode.GR{indGR},'gopt')==1
        % Optical generation model
        alpha = 7000; % 1/cm
        Eph = 0.5*qel; % J
        Popt = 0.1; % W/cm^2
        Gopt = zeros(1,mesh.nn);
        Gopt(iq) = alpha*Popt/Eph;
        R(iq) = R(iq) - Gopt(iq);
        mode.Gopt = Gopt;
        %
    end
end

if(not(isfield(mode,'CoeffBTBT')))
    mode.CoeffBTBT = 0;
end

if(mode.CoeffBTBT > 1e-50)
%     GBTBT = zeros(1,mesh.nn);
%     dGBTBT_dn = zeros(1,mesh.nn);
%     dGBTBT_dp = zeros(1,mesh.nn);
%     %             Vinterp = mode.EFp(mesh.inTJ(end)) - mode.EFn(mesh.inTJ(1));
%     Vinterp = EFp(mesh.inTJ(end)) - EFn(mesh.inTJ(1));
%     % Vinterp = Ev(mesh.inTJ(end)) - Ec(mesh.inTJ(1));
%     dVinterp_dn = - Vt(mesh.inTJ)/mesh.Nc(mesh.inTJ)./mpferdr((EFn(mesh.inTJ)-Ec(mesh.inTJ))./Vt(mesh.inTJ),-1/2);
%     dVinterp_dp = - Vt(mesh.inTJ)/mesh.Nv(mesh.inTJ)./mpferdr(-(EFp(mesh.inTJ)-Ev(mesh.inTJ))./Vt(mesh.inTJ),-1/2);
%     %             dVinterp_dn = - Vt(mesh.inTJ)/mesh.Nc(mesh.inTJ)./mpferdr((EFn(mesh.inTJ)-Ec(mesh.inTJ))./Vt(mesh.inTJ),-1/2);
%     %             dVinterp_dp = - Vt(mesh.inTJ)/mesh.Nv(mesh.inTJ)./mpferdr(-(EFp(mesh.inTJ)-Ev(mesh.inTJ))./Vt(mesh.inTJ),-1/2);
%     [GBTBT_interp,dGBTBT_dn_interp,dGBTBT_dp_interp] = f_InterpGBTBT(mode.strName,Vinterp,dVinterp_dn,dVinterp_dp);
%     GBTBT(mesh.inTJ) = GBTBT_interp(2:end);
%     dGBTBT_dn(mesh.inTJ) = dGBTBT_dn_interp(2:end);
%     dGBTBT_dp(mesh.inTJ) = dGBTBT_dp_interp(2:end);
%     %             dGBTBT_dn(mesh.inTJ(1)) = dGBTBT_dn_interp(2);
%     %             dGBTBT_dp(mesh.inTJ(end)) = dGBTBT_dp_interp(end);
%     
%     
%     mode.GBTBT = GBTBT*mode.CoeffBTBT;
    %
%     R = R + GBTBT*mode.CoeffBTBT;
    
    R = R + mode.GBTBT;
    %             dRn = dRn + dGBTBT_dn*mode.CoeffBTBT;
    %             dRp = dRp + dGBTBT_dp*mode.CoeffBTBT;
else
    mode.GBTBT = 0*ones(1,mesh.nn);
end
%
% Assembling in electron continuity equation
% Recombination term
MM = qel.*[R(in1).*Lp1 R(in2).*Lp2]/C0;
tvet = tvet + sparse(nn+in12,1,MM(mask_in12),neq,1);
%
% Assembling in hole continuity equation
% Recombination term
MM = qel.*[R(in1).*Lp1 R(in2).*Lp2]/C0;
tvet = tvet + sparse(pp+in12,1,MM(mask_in12),neq,1);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assembling quantum corrections and related jacobians
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(mode.Quantum==1)
    % Basic parameters
    NQW=length(mesh.iMQW);
    mode.N2D=0; mode.P2D=0;
    
    if(mode.Optical==1)
        Gmod=zeros(NModes,1);
        vph=Clight./mode.nindexQW;
    end
    for indQW=1:NQW
        % Population in each QW
        n2D=n2Dc{indQW};
        p2D=p2Dc{indQW};
        iQW = mesh.iMQW{indQW};  % central mesh index of the QW
        IQW = mesh.in{mesh.iMQWlayer{indQW}}; % all the QW mesh indeces
        WQW = mesh.WMQW{indQW}; % QWs width
        TQW = T(iQW);   % QW temperature
        VtQW = Vt(iQW);
        %
        % Quantities coming from macroQW computed in
        % fmp_SynchronizeMacro.m
        tauscatn = mesh.tauscatnMQW{indQW}; % QW electron capture lifetime
        tauscatp = mesh.tauscatpMQW{indQW}; % QW hole capture lifetime
        meffnQW = mesh.meffnMQW{indQW}; % QW electron effective mass
        meffpQW = mesh.meffpMQW{indQW}; % QW hole effective mass
        %
        % Energy differences between barrier and QW bottom level
        DeltaEc = mesh.DeltaEcQW{indQW}(1); % Conduction band energy
        DeltaEv = mesh.DeltaEvQW{indQW}(1); % Valence band energy
        
        if mode.pointEig==0
            lambda_C=mesh.lambda_C{indQW};  % QW CB eigenvalues (from Schrodinger)
            lambda_V=mesh.lambda_V{indQW};  % QW VB eigenvalues (from Schrodinger)
        else
            lambda_C=0;
            lambda_V=0;
        end
        
        Psi_C2=mesh.Psi_C2{indQW};  % QW CB eigenvector (from Schrodinger)
        Psi_V2=mesh.Psi_V2{indQW};  % QW VB eigenvector (from Schrodinger)
        Nc2D = mesh.Nc2D{indQW};    % QW available conduction band states
        Nv2D = mesh.Nv2D{indQW};    % QW available valence band states
        %
        Ec0 = Ec(iQW);  % bulk (barrier) CB energy level at QW mesh point
        Ev0 = Ev(iQW);  % bulk (barrier) VB energy level at QW mesh point
        Em = Ec0 - abs(DeltaEc) + lambda_C; % eigenvalue referred to bottom of QW
        Ed = Ev0 + abs(DeltaEv) - lambda_V; % eigenvalue referred to top of QB
        N2 = 4*pi*m0/h^2*qel.*sum(meffnQW.*(Ec0-Em),1)/10000;
        P2 = 4*pi*m0/h^2*qel.*sum(meffpQW.*(Ed-Ev0),1)/10000;
        %
        if(n2D>=N2 || isnan(n2D))
            n2D=0.99*N2;
        end
        if(p2D>=P2 || isnan(p2D))
            p2D=0.99*P2;
        end
        % Compute QW state filling
        denN = N2;
        denP = P2;
        rho_n{indQW}=n2D*C0./denN;
        rho_p{indQW}=p2D*C0./denP;
        drho_n{indQW}=C0./denN;
        drho_p{indQW}=C0./denP;
        %
        % Compute bound electron quasi Fermi level with inverse
        % Fermi integral from n2D, bulk and bound energy levels,
        % and bulk quasi fermi levels
        [EFn2D]=mpinvferdr_n2D(n2D*C0,Em,Ec0,Nc2D,VtQW,EFp(iQW)-20,EFn(iQW)+20);
        xnum=(ones(length(lambda_C),1)*EFn2D-Em)./VtQW; xden=(ones(length(lambda_C),1)*EFn2D-Ec0)./VtQW;
        [n2Dm,dn2D_2] = mpferdr2D(xnum,xden); % compute m-th level population with Fermi integral
        n2Dm=Nc2D.*n2Dm/C0;
        dn2Dm_EFn2D=Nc2D.*dn2D_2./C0./VtQW;
        dn2D_EFn2D=sum(dn2Dm_EFn2D,1);
        %
        N2D = reshape(Psi_C2.'*n2Dm,1,nn);  % Bound population computed by means of the eigenfucntion
        dN2D_EFn2D = reshape(Psi_C2.'*dn2Dm_EFn2D,1,nn);
        dn2D_EFn2Dvec = reshape((ones(size(Psi_C2))).'*dn2Dm_EFn2D,1,nn);
        dN2D_n2D = dN2D_EFn2D./dn2D_EFn2Dvec;   % derivative with respect to n2D
        %
        % Compute bound holes quasi Fermi level with inverse
        % Fermi integral from p2D, bulk and bound energy levels,
        % and bulk quasi fermi levels
        [EFp2D]=mpinvferdr_p2D(p2D*C0,Ed,Ev0,Nv2D,VtQW,EFp(iQW)-20,EFn(iQW)+20);
        xnum=(Ed-ones(length(lambda_V),1)*EFp2D)./VtQW; xden=(Ev0-ones(length(lambda_V),1)*EFp2D)./VtQW;
        [p2Dm,dp2D_2] = mpferdr2D(xnum,xden);
        p2Dm=Nv2D.*p2Dm/C0;
        dp2Dm_EFp2D=-Nv2D.*dp2D_2./C0./VtQW;
        dp2D_EFp2D=sum(dp2Dm_EFp2D,1);
        %
        P2D = reshape(Psi_V2.'*p2Dm,1,nn);
        dP2D_EFp2D = reshape(Psi_V2.'*dp2Dm_EFp2D,1,nn);
        dp2D_EFp2Dvec = reshape((ones(size(Psi_V2))).'*dp2Dm_EFp2D,1,nn);
        dP2D_p2D = dP2D_EFp2D./dp2D_EFp2Dvec;
        if strcmp(mode.QuantumAssem,'VENUS')==1
            %
            % Assembling QW electron charge in Poisson equation
            % known term (similar to doping)
            MM = qel.*[N2D(in1).*Lp1, N2D(in2).*Lp2];
            tvet = tvet + sparse(in12,1,MM(mask_in12),neq,1);
            %
            % Assembling QW hole charge in Poisson equation
            % known term (similar to doping)
            MM = - qel.*[P2D(in1).*Lp1, P2D(in2).*Lp2];
            tvet = tvet + sparse(in12,1,MM(mask_in12),neq,1);
            %
        else if strcmp(mode.QuantumAssem,'NUSOD')==1
                % NUSOD way (potential is inserted in Kmat0)
                %
                dN2D_n2D = Psi_C2(1,:);
                MM = qel*[dN2D_n2D(in1).*Lp1, dN2D_n2D(in2).*Lp2];
                Kmat0 = Kmat0 + sparse(in12,vv+indQW,MM(mask_in12),neq,neq);
                %
                dP2D_p2D = Psi_V2(1,:);
                MM = -qel*[dP2D_p2D(in1).*Lp1, dP2D_p2D(in2).*Lp2];
                Kmat0 = Kmat0 + sparse(in12,ww+indQW,MM(mask_in12),neq,neq);
                %
                N2D = dN2D_n2D(1,:)*n2D;
                P2D = dP2D_p2D(1,:)*p2D;
            end
        end
        %
        mode.N2D = mode.N2D + N2D*C0;
        mode.P2D = mode.P2D + P2D*C0;
        %
        if strcmp(mode.QuantumAssem,'VENUS')==1
            % 2-D electron capture term ------------------------
            Ccapn = - (1-exp((EFn2D-EFn(iQW))/VtQW)).*(1-n2D.*C0/N2)./tauscatn.*elec(iQW).*C0;
        else
            % NUSOD way
            taucapn = tauscatn;
            tauescn = taucapn.*mode.rn_ec{indQW}(1);
            %
            Ccapn = - (1-n2D*C0./N2)*elec(iQW)*C0./taucapn + n2D*C0./(WQW*tauescn);
            %
        end
        %
        % Assembling 2-D electron capture term ------------------------
        MM    = qel.*(Ccapn.*WQW)/C0;
        tvet  = tvet + sparse(vv+indQW,1,MM,neq,1);
        %
        % Assembling 3-D electron capture term ------------------------
        MM    = - qel.*mesh.Lp(IQW).*Ccapn/C0;
        tvet  = tvet + sparse(nn+IQW,1,MM,neq,1);
        %
        if strcmp(mode.QuantumAssem,'VENUS')==1
            % 2-D hole capture term ---------------------------------------
            Ccapp = - (1-exp((EFp(iQW)-EFp2D)/VtQW)).*(1-p2D*C0/P2)./tauscatp.*hole(iQW)*C0;
        else
            % NUSOD way
            taucapp = tauscatp;
            tauescp = taucapp.*mode.rp_ec{indQW}(1);
            
            Ccapp = - (1-p2D*C0./P2)*hole(iQW)*C0./taucapp + p2D*C0./(WQW*tauescp);
        end        
        %
        % Assembling 2-D hole capture term ------------------------
        MM    = qel.*(Ccapp.*WQW)/C0;
        tvet  = tvet + sparse(ww+indQW,1,MM,neq,1);
        %
        % Assembling 3-D hole capture term ------------------------
        MM    = - qel.*mesh.Lp(IQW).*Ccapp/C0;
        tvet  = tvet + sparse(pp+IQW,1,MM,neq,1);
        %
        % Recombination of bound carriers
        R2D = 0; dR2D_n2D=0; dR2D_p2D=0;
        n2Di=mode.n2Di{indQW}; p2Di=mode.p2Di{indQW};
        %
        for indGR=1:length(mode.GRQuantum)
            %
            np2D=C0*n2D.*C0*p2D-n2Di.*p2Di;
            dnp2D_n=C0.*C0.*p2D; dnp2D_p=C0.*C0.*n2D;
            if strcmp(mode.GRQuantum{indGR},'SRH')==1
                % SRH recombination model
                % Shockley-Read-Hall recombination term
                denSRH2D = mesh.taup(iQW).*(C0.*n2D+n2Di)+mesh.taun(iQW).*(C0.*p2D+p2Di);
                RSRH2D = np2D./denSRH2D;
                dRSRH_n2D = (C0.*C0.*p2D.*denSRH2D - mesh.taup(iQW).*C0.*np2D)./(denSRH2D.^2);
                dRSRH_p2D = (C0.*C0.*n2D.*denSRH2D - mesh.taun(iQW).*C0.*np2D)./(denSRH2D.^2);
                %
                R2D = R2D + RSRH2D;
                dR2D_n2D = dR2D_n2D + dRSRH_n2D;
                dR2D_p2D = dR2D_p2D + dRSRH_p2D;
                %
            elseif strcmp(mode.GRQuantum{indGR},'rad')==1
                
                % Radiative recombination term
                Brad2D = mesh.brad(iQW)/WQW; % Brad in cm^2/s (/WQW!)
                Rrad2D = Brad2D.*np2D; % cm^2/s
                dRrad2D_n2D = Brad2D.*dnp2D_n;
                dRrad2D_p2D = Brad2D.*dnp2D_p;
                %
                R2D = R2D + Rrad2D;
                dR2D_n2D = dR2D_n2D + dRrad2D_n2D;
                dR2D_p2D = dR2D_p2D + dRrad2D_p2D;
                %
            elseif strcmp(mode.GRQuantum{indGR},'Auger')==1
                
                % Auger recombination model
                Cnnp2D=mesh.Cnnp(iQW)./(WQW.^2);
                Cppn2D=mesh.Cppn(iQW)./(WQW.^2);
                RAuger2D = (Cnnp2D.*C0.*n2D+Cppn2D.*C0.*p2D).*np2D;
                dRAuger_n2D = 2*Cnnp2D.*np2D.*C0 + Cppn2D.*p2D.*dnp2D_n.*C0;
                dRAuger_p2D = 2*Cppn2D.*np2D.*C0 + Cnnp2D.*n2D.*dnp2D_p.*C0;
                R2D = R2D + RAuger2D;
                dR2D_n2D = dR2D_n2D + dRAuger_n2D;
                dR2D_p2D = dR2D_p2D + dRAuger_p2D;
            end
        end
        %
        mode.RSRH2D{iv,indQW} = RSRH2D;
        mode.Rrad2D{iv,indQW} = Rrad2D;
        mode.RAuger2D{iv,indQW} = RAuger2D;
        %
        % 2D electron equation recombination term ---------------------------------
        MM=qel.*R2D/C0;
        tvet = tvet + sparse(vv+indQW,1,MM,neq,1);
        %
        % 2D hole equation recombination term -------------------------------------
        MM=qel.*R2D/C0;
        tvet = tvet + sparse(ww+indQW,1,MM,neq,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Stimulated emission: coupling with optical simulator
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(mode.Optical==1)
            nQW=n2D;
            pQW=p2D;
            %
            Gamma_z=mode.Gamma_z(indQW);
            %
            %-- Stimulated emission: rate equation
            for indMode=1:NModes
                
                [g,dgE,dgH,rsp,drspE,drspH] = f_InterpGain_lin(double(nQW)*C0,double(pQW)*C0,indQW,indMode);
                %
                g=vph*g*mode.nlG;
                dgE=vph*dgE*mode.nlG*C0;
                dgH=vph*dgH*mode.nlG*C0;
                %
                rsp=rsp;
                drspE=drspE*C0;
                drspH=drspH*C0;
                %
                gm=g;
                dgm_n=dgE;
                dgm_p=dgH;
                Rsp=rsp;
                Lm=mode.Lm;
                
                %%% Saving control parameters which are needed to
                %%% add extra bias points, helping the simulation
                %%% when optical threshold is reached
                mode.gm = gm;
                mode.Lm = Lm;
                mode.vph = vph;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                %-- Spontaneous emission power
                Si=sum(Gamma_z*mode.fPES(indMode).*Rsp);
                tvet(ss+indMode)=Si;
                %
                %-- Multiplicative term times stimulated modal power
                MM=sum(Gamma_z.*gm-Lm/NQW);
                Kmat0=Kmat0+sparse(ss+indMode,ss+indMode,MM,neq,neq);
                %
                Gmod(indMode)=Gmod(indMode)+sum(Gamma_z.*gm/(vph)); % for display only
                %
                %-- Stimulated emission: recombination term
                Rst=g.*Pst(indMode).*mode.fPdif(indMode)./mode.Area*WQW;
                MM=qel.*Rst/C0;
                tvet=tvet+sparse(vv+indQW,1,MM,neq,1);
                tvet=tvet+sparse(ww+indQW,1,MM,neq,1);
                %
            end
        end
    end
    
end%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assembling system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jmat=Kmat0+Jmat0+Jmat2;
rvet=(Kmat0+Jmat2)*uvet+tvet;
rvet_inc = rvet;
%

end
