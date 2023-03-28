function [Kmat0,Jmat0,Jmat2,uvet,rvet,mode,tvet]=assem_mp(geom,mesh,mode,uvet,v0_dd)
%
% define constants
%
mpConstants
smp_LoadConstants
Vt=mp(mode.Vt); % electrical equivalent in temperature, V
% =============================================================================================100
% Loading solution database
% =============================================================================================100
nmodes=mode.nmodes; % number of optical modes (from VELM)
if(isfield(mesh,'nnxQW'))
    nnQW=mesh.nnxQW{1};          % number of lateral points in the QW (x-directed)
else
    nnQW=0;
end
NQW=mesh.NMQW;                          % number of quantum wells
nt=mesh.nt;                             % number of triangles
nm=geom.nm;                             % number of active contacts
nc=geom.nc;                             % number of contacts
nl=mode.ntrap;                          % number of trap levels
%
if mode.EqPermutazione==0
    %
    nn=mesh.nn;                   % -> elec eqs.
    pp=nn+mode.nflg*nn;           % -> hole eqs.
    tt=pp+mode.pflg*nn;           % -> trap eqs.
    qq=tt+mode.tflg*nl*nn;        % -> circuit eqs.
    rr=qq+nm;                     % -> current eqs.
    vv=rr+nm;                     % -> n2D eqs.
    ww=vv+length(mesh.xgrid)*NQW;               % -> p2D eqs.
    ss=ww+length(mesh.xgrid)*NQW;               % -> optical rate eqs.
else
    nn=mesh.nn;                             % -> elec eqs.
    pp=nn+mode.nflg*nn;                     % -> hole eqs.
    tt=pp+mode.pflg*nn;                     % -> trap eqs.
    vv=tt+mode.tflg*nl*nn;                  % -> n2D eqs
    ww=vv+length(mesh.xgrid)*NQW;           % -> p2D eqs
    ss=ww+length(mesh.xgrid)*NQW;           % -> optical rate eqs.
    qq=ss+mode.nmodes;
    rr=qq+nm;
end
%
% total number of equations/unknowns
neq=nn+mode.nflg*nn+mode.pflg*nn+mode.tflg*nl*nn+2*nm+mode.oflg*not(mode.firstrun)*2*NQW*length(mesh.xgrid)+mode.oflg*mode.nmodes;
%
% =============================================================================================100
% Loading and computing geometry parameters
% =============================================================================================100
% geometrical data
triangle=mesh.triangle; node=mp(mesh.node); contact=mesh.contact; iq=mesh.iq;
% indexes to the nodes 1, 2, 3 of each triangle
in1=triangle(1,:); in2=triangle(2,:); in3=triangle(3,:);
%
% x and y positions of node triangles
x1=node(1,in1); x2=node(1,in2); x3=node(1,in3);
y1=node(2,in1); y2=node(2,in2); y3=node(2,in3);
% distances from node triangles (FEM computations)
r3x=x2-x1; r1x=x3-x2; r2x=x1-x3;
r3y=-y2+y1; r1y=-y3+y2; r2y=-y1+y3;
% FEM basis functions
f23xx=r2x.*r3x;
f23yy=r2y.*r3y;
f23xy=r2x.*r3y;
f32xy=r3x.*r2y;
% triangle side lengths
l1=sqrt(r1x.^two+r1y.^two); l2=sqrt(r2x.^two+r2y.^two); l3=sqrt(r3x.^two+r3y.^two);
xm1=half*(x2+x3); xm2=half*(x3+x1); xm3=half*(x1+x2);
ym1=half*(y2+y3); ym2=half*(y3+y1); ym3=half*(y1+y2);
% split points
xcc=half*(x2.*f32xy+x1.*f32xy-y2.*f23yy-x3.*f23xy-x1.*f23xy+y3.*f23yy)./(f32xy-f23xy);
ycc=half*(x2.*f23xx-y2.*f23xy-y1.*f23xy-x3.*f23xx+y3.*f32xy+y1.*f32xy)./(f32xy-f23xy);
% triangle mid-points
s1=sqrt((xm1-xcc).^two+(ym1-ycc).^two)*zero;
s2=sqrt((xm2-xcc).^two+(ym2-ycc).^two);
s3=sqrt((xm3-xcc).^two+(ym3-ycc).^two);
Se=abs(r2x.*r1y-r2y.*r1x)./two; % triangle area
% split areas (for each triangle, the area associated to node 1, 2, 3)
Se1=half*(xm3.*ycc-xcc.*ym3+x1.*ym3-y1.*xm3+xcc.*ym2-xm2.*ycc-x1.*ym2+y1.*xm2);
Se2=half*(xm1.*ycc-xcc.*ym1+x2.*ym1-y2.*xm1+xcc.*ym3-xm3.*ycc-x2.*ym3+y2.*xm3);
Se3=half*(xm2.*ycc-xcc.*ym2+x3.*ym2-y3.*xm2+xcc.*ym1-xm1.*ycc-x3.*ym1+y3.*xm1);
%
% marking triangles outside semiconductor regions
semiconductor=find(geom.semiconductor);
it=not(ismember(triangle(4,:),semiconductor));
% by setting to zero Se1, Se2, Se3, we ensure that continuity equations
% and trap equations are assembled only in the semiconductor regions!
Se1(it)=zero; Se2(it)=zero; Se3(it)=zero;
if(mode.AugerExpulsion3D==1)
    SeCAV1=zero*ones(1,mesh.nt);
    SeCAV2=SeCAV1;
    SeCAV3=SeCAV1;
    SeCAV1(mesh.ITrCAV)=Se1(mesh.ITrCAV);
    SeCAV2(mesh.ITrCAV)=Se2(mesh.ITrCAV);
    SeCAV3(mesh.ITrCAV)=Se3(mesh.ITrCAV);
end
%
% =============================================================================================100
% Cylindrical coordinates correction
% =============================================================================================100
G1 = one; G2 = one; G3 = one;
if((isfield(mode,'symmetry'))&&(strcmp(mode.symmetry,'Cylindrical-Y')))
    xc = (x1+x2+x3)/three;
    Se1=twopi*xc.*Se1;
    Se2=twopi*xc.*Se2;
    Se3=twopi*xc.*Se3;
    if(mode.AugerExpulsion3D==1)
        SeCAV1=twopi*xc.*SeCAV1;
        SeCAV2=twopi*xc.*SeCAV2;
        SeCAV3=twopi*xc.*SeCAV3;
    end
    G1 = twopi*xc.*mpGji(x2,x3,mp(mode.tolGji));
    G2 = twopi*xc.*mpGji(x3,x1,mp(mode.tolGji));
    G3 = twopi*xc.*mpGji(x1,x2,mp(mode.tolGji));  
end
%
if((isfield(mode,'symmetry'))&&(strcmp(mode.symmetry,'Cylindrical-X')))
    yc = (y1+y2+y3)/three;
    Se1=twopi*yc.*Se1;
    Se2=twopi*yc.*Se2;
    Se3=twopi*yc.*Se3;
    if(mode.AugerExpulsion3D==1)
        SeCAV1=twopi*yc.*SeCAV1;
        SeCAV2=twopi*yc.*SeCAV2;
        SeCAV3=twopi*yc.*SeCAV3;
    end
    G1 = twopi*yc.*Gji(y2,y3,mode.tolGji);
    G2 = twopi*yc.*Gji(y3,y1,mode.tolGji);
    G3 = twopi*yc.*Gji(y1,y2,mode.tolGji); 
end
% =============================================================================================100
% here we define the logical masks for dirichlet boundary conditions;
% these masks are used when assembling the model equations:
%
%    1) maskr -> Poisson equation
%    2) masks -> elec/hole continuity equations, trap equations
%    3) maskt -> current equations
%
maskr=true(1,nn); masks=true(1,nn); maskt=false(1,nn);
for ic=1:nc
    ii=(contact==ic); % loop on contacts
    switch geom.contact_type(ic)
        case 1 % ohmic contact
            maskr(ii)=0; % Poisson equation
            masks(ii)=0; % elec/hole equations
        case 2 % Schottky contact
            maskr(ii)=0; % Poisson equation
        case 3 % Schottky contact, vsurf=Inf
            maskr(ii)=0; % Poisson equation
            masks(ii)=0; % elec/hole equations
        otherwise, error('ls_solve -> contact type unknown!'), 
    end
    if(ic<=nm)
        maskt(ii)=1; 
    end
end
    %
    masks(not(iq))=0; % continuity equations are eliminated outside semiconductors;
    maskt(not(iq))=0; % current equations are eliminated outside semiconductors;
    % =============================================================================================100
    % here we define the pointers to the matrix entries
    ijr=[in1 in2 in3]; ijs=ijr;
    iir=[in1 in1 in1 in2 in2 in2 in3 in3 in3]; iis=iir; iiu=iir;
    jjr=[in1 in2 in3 in1 in2 in3 in1 in2 in3]; jjs=jjr; jjt=jjr;
    %
    mask_iir=maskr(iir);
    mask_iis=masks(iir);
    mask_iit=maskt(iir);
    %
    mask_ijr=maskr(ijr);
    mask_ijs=masks(ijr);
    %
    % the pointers IIT and JJT are necessary for the current equations;
    in=zeros(1,nn);
    for ic=1:nm
        ii=find(contact==1); % contacts
        in(ii)=rr+(ic-1)+1;
    end
    iit=in(iir);
    %
    iir=iir(mask_iir); jjr=jjr(mask_iir);
    iis=iis(mask_iis); jjs=jjs(mask_iis);
    iit=iit(mask_iit); jjt=jjt(mask_iit);
    %
    ijr=ijr(mask_ijr);
    ijs=ijs(mask_ijs);
    %
    % =============================================================================================100
    % Loading material parameters and previous solution
    % =============================================================================================100
    % loading material parameters
    if isfield(mode,'Fat_Dop')
        FatDop=mp(mode.Fat_Dop);
    else
        FatDop=one;
    end
    dop_a=FatDop*mp(mesh.dop_a); dop_d=FatDop*mp(mesh.dop_d); epsxx=mp(mesh.epsxx);
    mobn0=mp(mesh.mobn0_t); mobp0=mp(mesh.mobp0_t); Nc=mp(mesh.Nc); Nv=mp(mesh.Nv);
    nint=mp(mesh.nint);
    Eg=mp(mesh.Eg);
    affinity=mp(mesh.affinity);
    
    if isfield(mode,'ireno')==1
        ireno=mode.ireno;
        if ireno==1
            Reno=2.7*10^(-8)*mode.dop_am.^(1/3);
            Eg=mesh.Eg-Reno;
            affinity=mesh.affinity-Reno*.4;
            %
            'EG',keyboard
        end
    end
    % loading neutrality solution
    phi_neutr=mp(mode.phi_neutr);
    %
    % loading previous solution
    phi=mp(uvet(1:nn)); % electric potential
    if mode.IdriveON==1
        i0_dd=mp(v0_dd);    % current driving
        Ymat=mp(mode.Zmat);
        
        v_ls=mp(uvet(qq+1));  % contact voltages
        i_ls=-mp(v_ls)*(Ymat.')+i0_dd;  % contact currents
        uvet(rr+1)=i_ls;
    else
        i_ls=mp(uvet(rr+1)); % contact currents
        v_ls=-i_ls*mp(mode.Zmat.')+mp(v0_dd); % contact voltages
        uvet(qq+1)=v_ls;
    end
    if(mode.nflg), elec=abs(mp(uvet(nn+(1:nn)))); end % electron concentration
    if(mode.pflg), hole=abs(mp(uvet(pp+(1:nn)))); end % hole concentration
    if(mode.oflg)
        Pst=abs(mp(uvet(ss+(1:nmodes))));
        for indQW=1:NQW
            n2Dc{indQW} = abs(mp(uvet(vv+(indQW-1)*length(mesh.xgrid)+(1:nnQW))));
            p2Dc{indQW} = abs(mp(uvet(ww+(indQW-1)*length(mesh.xgrid)+(1:nnQW))));
        end
    end
    % =============================================================================================100
    % these vectors are necessary for the treatment of heterointerfaces
    offset_n=zero*ones(1,nn); offset_p=zero*ones(1,nn);
    offset_n(iq)=affinity(iq)+Vt(iq).*log(Nc(iq));
    offset_p(iq)=affinity(iq)-Vt(iq).*log(Nv(iq))+Eg(iq);
    offset_n1=offset_n(in1); offset_n2=offset_n(in2); offset_n3=offset_n(in3);
    offset_p1=offset_p(in1); offset_p2=offset_p(in2); offset_p3=offset_p(in3);
    Vt1 = Vt(in1); Vt2 = Vt(in2); Vt3 = Vt(in3);
    % =============================================================================================100
    % band diagram, eV
    mode.ecb = NaN*zero*ones(1,nn);
    mode.evb = NaN*zero*ones(1,nn);
    mode.EFn = zero*ones(1,nn);
    mode.EFp = zero*ones(1,nn);
    mode.Efi = NaN*zero*ones(1,nn);
    %
    mode.ecb(iq) = - phi(iq) + Eg(iq)/two - Vt(iq)./two.*log(Nv(iq)./Nc(iq)) + mp(mesh.phi_r(iq));
    mode.evb(iq) = mode.ecb(iq) - Eg(iq);
    mode.Efi(iq) = (mode.ecb(iq) + mode.evb(iq))/two + Vt(iq)./two.*log(Nv(iq)./Nc(iq));
    if(mode.nflg), mode.EFn(iq) = mode.ecb(iq) + Vt(iq).*log(elec(iq)./Nc(iq)); end
    if(mode.pflg), mode.EFp(iq) = mode.evb(iq) - Vt(iq).*log(hole(iq)./Nv(iq)); end
    %
    if((isfield(mode,'stats'))&&(strcmp(mode.stats,'Fermi'))) % Fermi ---------
        %
        if(mode.nflg)
            elec = abs(elec); % kludge: forcing electron charge to be positive
            tmp = elec./Nc;
            mode.EFn = mode.ecb + Vt.*mpinvferdr(tmp,mp(mode.tolinvferdr)); 
        end
        nF  = Nc.*mpferdr((mode.EFn - mode.ecb)./Vt,half); % 1/cm^3;
        dnF = Nc.*mpferdr((mode.EFn - mode.ecb)./Vt,-half);
        nB =  Nc.*exp  ((mode.EFn - mode.ecb)./Vt); % 1/cm^3;
        gamman = nF./nB; gamman(not(iq)) = one;
        dgamman = (one - nF./dnF)./nB;
        if(not(mode.firstrun))
            ii = tmp<mp('1e-9'); dgamman(ii) = - sqrt(two)/mp('4')./Nc(ii); % non-degenerate limit
        end
        dlGn = dgamman./gamman; dlGn(not(iq)) = zero;
        dlGn1 = dlGn(in1); dlGn2 = dlGn(in2); dlGn3 = dlGn(in3);
        %
        if(mode.pflg)
            hole = abs(hole); % kludge: forcing hole charge to be positive
            tmp = hole./Nv;
            mode.EFp = mode.evb - Vt.*mpinvferdr(tmp,mp(mode.tolinvferdr)); 
        end
        pF  = Nv.*mpferdr((- mode.EFp + mode.evb)./Vt,half); % 1/cm^3;
        dpF = Nv.*mpferdr((- mode.EFp + mode.evb)./Vt,-half);
        pB =  Nv.*exp  ((- mode.EFp + mode.evb)./Vt); % 1/cm^3;
        gammap = pF./pB; gammap(not(iq)) = one;
        dgammap = (one - pF./dpF)./pB;
        if(not(mode.firstrun))
            ii = tmp<mp('1e-9'); dgammap(ii) = - sqrt(two)/mp('4')./Nv(ii); % non-degenerate limit
        end
        dlGp = dgammap./gammap; dlGp(not(iq)) = zero;
        dlGp1 = dlGp(in1); dlGp2 = dlGp(in2); dlGp3 = dlGp(in3);
        %
    else % Boltzmann ----------------------------------------------------------
        %
        if(mode.nflg)
            elec = abs(elec); % kludge: forcing electron charge to be positive
            tmp = elec./Nc;
            mode.EFn = mode.ecb + Vt.*log(tmp); 
        end
        nB =  Nc.*exp((mode.EFn - mode.ecb)./Vt); % 1/cm^3;
        nF  = nB; % 1/cm^3;
        dnF = nF;
        gamman = ones(size(nB),'mp');
        dgamman = zero*ones(size(nB));
        dlGn = dgamman;
        dlGn1 = dlGn(in1); dlGn2 = dlGn(in2); dlGn3 = dlGn(in3);
        %
        if(mode.pflg)
            hole = abs(hole); % kludge: forcing hole charge to be positive
            tmp = hole./Nv;
            mode.EFp = mode.evb - Vt.*log(tmp); 
        end
        pB =  Nv.*exp((- mode.EFp + mode.evb)./Vt); % 1/cm^3;
        pF  = pB; % 1/cm^3;
        dpF = pF;
        gammap = ones(size(pB),'mp');
        dgammap = zero*ones(size(pB));
        dlGp = one./gammap.*dgammap; dlGp(not(iq)) = zero;
        dlGp1 = dlGp(in1); dlGp2 = dlGp(in2); dlGp3 = dlGp(in3); 
    end % ------------
    %
    % =============================================================================================100
    phi=reshape(phi,1,nn);
    phi1=phi(in1); phi2=phi(in2); phi3=phi(in3);
    if(mode.nflg)
        elec=reshape(elec,1,nn);
        elec1=elec(in1); elec2=elec(in2); elec3=elec(in3);
    end
    if(mode.pflg)
        hole=reshape(hole,1,nn);
        hole1=hole(in1); hole2=hole(in2); hole3=hole(in3);
    end
    % =============================================================================================100
    % Mobility models
    [Dn, Dp, dDndphi1, dDndphi2, dDndphi3, dDpdphi1, dDpdphi2, dDpdphi3] = ...
        mobility(mode,mobn0,mobp0,Se,r1x,r2x,r3x,r1y,r2y,r3y);
    % the mobility function should be called for each material
    % here we set to zero the mobility outside semiconductors
    Dn(it)=zero; Dp(it)=zero;
    % compute mobility from Einstein relation
    mode.mobn=mp(Dn)./mp(mode.Vt_tr); mode.mobp=mp(Dp)./mp(mode.Vt_tr); % save field-dependent mobilities
    % =============================================================================================100
    % GR models
    % =============================================================================================100
    R=zero*ones(1,nn); dRn=zero*ones(1,nn); dRp=zero*ones(1,nn);
    RSRH=zero*ones(1,nn);
    Rrad=zero*ones(1,nn);
    RAuger=zero*ones(1,nn);
    GBTBT=zero*ones(1,nn);
    
    RAuger_n3D = zero*ones(1,nn); RAuger_p3D = zero*ones(1,nn);
    dRAuger_n3D_n3D = zero*ones(1,nn); dRAuger_n3D_p3D = zero*ones(1,nn);
    dRAuger_p3D_n3D = zero*ones(1,nn); dRAuger_p3D_p3D = zero*ones(1,nn);
    %
    % add current equations
    Jmat0=mp(sparse(neq,neq),'mp');
    tvet=mp(sparse(neq,1),'mp');
    
    if(not(mode.firstrun)) % At equilibrium, R=0
        for indGR=1:length(mode.GR)
            if(not(mode.nflg && mode.pflg))
                error('ls_solve -> select both elec and hole eqs. with GR models')
            end
            n  =  elec(iq); p  =  hole(iq);
            np = (n.*p-gamman(iq).*gammap(iq).*nint(iq).^two);
            dnp_n = (p-dgamman(iq).*gammap(iq).*nint(iq).^two);
            dnp_p = (n-dgammap(iq).*gamman(iq).*nint(iq).^two);
            if strcmp(mode.GR{indGR},'SRH')==1
                % SRH recombination model
                n1 = gamman(iq).*nint(iq).*exp(mp(mesh.Etrap(iq))./Vt(iq));
                p1 = gammap(iq).*nint(iq).*exp(-mp(mesh.Etrap(iq))./Vt(iq));
                denSRH = mp(mesh.taup(iq)).*(n+n1)+mp(mesh.taun(iq)).*(p+p1);
                RSRH(iq) = np./denSRH;
                R(iq) = R(iq) + RSRH(iq);
                dRn(iq) = dRn(iq) + ((p-dgamman(iq).*gammap(iq).*nint(iq).^two).*denSRH - ...
                    mp(mesh.taup(iq)).*(one+dgamman(iq).*n1).*np)./(denSRH.^two);
                dRp(iq) = dRp(iq) + ((n-dgammap(iq).*gamman(iq).*nint(iq).^two).*denSRH - ...
                    mp(mesh.taun(iq)).*(one+dgammap(iq).*p1).*np)./(denSRH.^two);
            elseif strcmp(mode.GR{indGR},'rad')==1
                % radiative recombination model
                Rrad(iq) = mp(mesh.brad(iq)).*np;
                R(iq) = R(iq) + Rrad(iq);
                dRn(iq)=dRn(iq)+mp(mesh.brad(iq)).*dnp_n;
                dRp(iq)=dRp(iq)+mp(mesh.brad(iq)).*dnp_p;
            elseif strcmp(mode.GR{indGR},'Auger')==1
                % Auger recombination model
                %
                RAuger_n3D(iq) = mp(mesh.Cnnp(iq)).*n.*np;
                RAuger_p3D(iq) = mp(mesh.Cppn(iq)).*p.*np;
                
                dRAuger_n3D_n3D(iq) = two*mp(mesh.Cnnp(iq)).*np;
                dRAuger_n3D_p3D(iq) = mp(mesh.Cnnp(iq)).*dnp_p;
                
                dRAuger_p3D_n3D(iq) = mp(mesh.Cppn(iq)).*dnp_n;
                dRAuger_p3D_p3D(iq) = two*mp(mesh.Cppn(iq)).*np;
                %
                RAuger(iq) = RAuger_n3D(iq) + RAuger_p3D(iq);
                %
                R(iq) = R(iq) + RAuger(iq);
                dRn(iq) = dRn(iq) + dRAuger_n3D_n3D(iq) + dRAuger_p3D_n3D(iq);
                dRp(iq) = dRp(iq) + dRAuger_n3D_p3D(iq) + dRAuger_p3D_p3D(iq);
            elseif (strcmp(mode.GR{indGR},'BTBT')==1) 
                %
                % Load the NEGF LUT
                if (abs(v0_dd)==zero)
                    LUT_GBTBT
                end
                %
                % BTBT generation rate, from NEGF
                %
                for indexBTJ = 1:mode.nBTJ
                    %
                    GBTBT=zeros(1,nn);dRnTJ=zeros(1,mesh.nn); dRpTJ=dRnTJ; dRphin=dRnTJ; dRphip=dRpTJ;
                    %
                    iTJfake=indexBTJ;
                    %
                    IBTJ=sort(unique([mesh.IBTJ{indexBTJ,:}]));
                    ITJp_vet=mesh.iLBTJ{iTJfake,end};                     % right BTJ node, (no passivation)
                    ITJn_vet=mesh.iRBTJ{iTJfake,1};                     % left BTJ node
                    iTJp=ITJp_vet(1);
                    iTJn=ITJn_vet(1);
                    %
                    iRBTJP=ITJn_vet;                                % RIGHT BTJ nodes, all the columns
                    iLBTJP=ITJp_vet;                                % LEFT BTJ nodes, all the columns
                    %
                    IIBTJP=repmat(1:mesh.nnx,mesh.nny,1);   IIBTJP=IIBTJP(:).';
                    %
                    IRBTJP=iRBTJP(IIBTJP);  % needed to assemble the GBTBT derivatives only in the right and left nodes of the BTJ (as the QW: INMQWP)
                    ILBTJP=iLBTJP(IIBTJP);
                    %
                    SeTJ1=zero*ones(1,mesh.nt); SeTJ2=SeTJ1; SeTJ3=SeTJ1;
                    SeTJ1(unique([mesh.ITrBTJ{indexBTJ,:}]))=Se1(unique([mesh.ITrBTJ{indexBTJ,:}]));
                    SeTJ2(unique([mesh.ITrBTJ{indexBTJ,:}]))=Se2(unique([mesh.ITrBTJ{indexBTJ,:}]));
                    SeTJ3(unique([mesh.ITrBTJ{indexBTJ,:}]))=Se3(unique([mesh.ITrBTJ{indexBTJ,:}]));
                    %
                    % Rule of thumb:
                    % - BTJ is pn: n --> mesh.inTJ(end); p --> mesh.inTJ(1)
                    % - BTJ is np (2020Tibaldi_PRAPP): n --> mesh.inTJ(1); p --> mesh.inTJ(end)
                    %
                    Vint=mode.EFp(ITJp_vet)-mode.EFn(ITJn_vet);
                    Vint(isnan(Vint))=zero;
                    mode.VTJ(indexBTJ,mode.ind_v0)=Vint(1);
                    % Initialization of Vinterp derivatives for Jacobian matrix
                    % dn -> right (end); dp -> left (1); 1 -> right; 2 -> left
                    % Vt,ecb, evb, EFn, EFp should ALL be taken from "real"
                    % nodes, and not only from the central section!
                    % iTJend -> IBTJright; iTJ1 -> IBTJleft
                    dVint_dn = - Vt(ITJn_vet)./Nc(ITJn_vet)./mpferdr((mode.EFn(ITJn_vet)-mode.ecb(ITJn_vet))./Vt(ITJn_vet),-half);
                    dVint_dp = - Vt(ITJp_vet)./Nv(ITJp_vet)./mpferdr(-(mode.EFp(ITJp_vet)-mode.evb(ITJp_vet))./Vt(ITJp_vet),-half);
                    dVint_dn(isnan(dVint_dn))=zero; dVint_dp(isnan(dVint_dp))=zero;
                    % La derivata rispetto al potenziale si calcola usando EFn = Ec
                    % + altro,
                    % dphi1 --> dRphi1, TJ(end); dphin --> dRphi2, TJ(1)
                    dVint_dphip = -one; % derivata rispetto a EFp (-1 per colpa di Ev = -q*phi + puffi)
                    dVint_dphin = +one; % derivata rispetto a EFn (-1 per colpa di Ev = -q*phi + puffi, altro -1 perché ha segno -)
                    
                    % TIBALDI JACOBIAN GBTBT
                    % lunghezza della regione dove spalmo il tasso
                    LTJ = mp(mesh.ygrid(iTJn))-mp(mesh.ygrid(iTJp));
                    
                    pDeg=size(mode.cBTBT,2)-1;
                    
                    VV=zero*ones(length(ITJp_vet),pDeg+1);
                    dVV=zero*ones(length(ITJp_vet),pDeg);
                    %
                    [~,iRagTJ]=min(abs(mp(mesh.xgrid)*mp('1e4')-mp(mode.rAperture)));
                    %
%                     for iRad=1:length(ITJp_vet)
                    for iRad=1:iRagTJ
                        
                        [~,iT]=min(abs(mp(mode.T_NEGF)-mean(mp(mesh.T(ITJp_vet(iRad):ITJn_vet(iRad))))));
                        cT=mp(mode.cBTBT(iT,:)); % fitting coefficient at a given temperature
                        %             cT=c;       % fitting of NEGF w/out temperature
                        %
                        % coefficienti della derivata del fit della corrente (3 values,
                        % polyfit 3)
                        dc =cT(1:pDeg).*[pDeg:-1:1];
                        
                        VV(iRad,:) = Vint(iRad).^[pDeg:-1:0];
                        dVV(iRad,:) = Vint(iRad).^[pDeg-1:-1:0];
                        
                    end
                                        
                    I_interp=mp('10').^(VV*cT')-mp(mode.I0_NEGF(iT));
                    II=1/LTJ*(-I_interp/qel);
                    II(iRagTJ+1:end)=zero;
                    GBTBT(IBTJ)=II(IIBTJP(IBTJ))/CarrierNorm;
                    
                    ITJ=-VV*cT'*mode.AreaOx;
                    RTJ=abs(Vint)'./abs(ITJ);
                    mode.sigmaTJ=1./double(RTJ.*repmat(LTJ,mesh.nnx,1));

                    JTJ=-VV*cT';
                    mode.HeatTJ=abs(double(JTJ.*Vint'./repmat(LTJ,mesh.nnx,1)));
                    
                    dI_interp=log(mp('10'))*mp('10').^(VV*cT').*(dVV*dc');
                    dII=one/LTJ*(-dI_interp/qel)'/CarrierNorm;
                    
                    dRnTJ(IBTJ) = dII(IIBTJP(IBTJ)).*dVint_dn(IIBTJP(IBTJ));
                    dRpTJ(IBTJ) = dII(IIBTJP(IBTJ)).*dVint_dp(IIBTJP(IBTJ));
                    dRphin(IBTJ) = dII(IIBTJP(IBTJ)).*dVint_dphin;
                    dRphip(IBTJ) = dII(IIBTJP(IBTJ)).*dVint_dphip;
                    %
                    % GBTBT derivatives assembled
                    %
                    % Assembling in electron continuity equation
                    % Recombination term
                    MM = qel.*[SeTJ1.*GBTBT(in1) SeTJ2.*GBTBT(in2) SeTJ3.*GBTBT(in3)];
                    tvet = tvet + sparse(ijs+nn,1,MM(mask_ijs),neq,1,'mp');
                    %
                    % Derivative of recombination term w.r.t. electrons
                    MM = qel.*[SeTJ1.*dRnTJ(in1) SeTJ2.*dRnTJ(in2) SeTJ3.*dRnTJ(in3)];
                    Jmat0 = Jmat0 + sparse(ijs+nn,nn+IRBTJP(ijs),MM(mask_ijs),neq,neq,'mp');
                    % Derivative of recombination term w.r.t. holes
                    MM = qel.*[SeTJ1.*dRpTJ(in1) SeTJ2.*dRpTJ(in2) SeTJ3.*dRpTJ(in3)];
                    Jmat0 = Jmat0 + sparse(ijs+nn,pp+ILBTJP(ijs),MM(mask_ijs),neq,neq,'mp');
                    % Derivative of recombination term w.r.t. phi (left node)
                    MM = qel.*[SeTJ1.*dRphin(in1) SeTJ2.*dRphin(in2) SeTJ3.*dRphin(in3)];
                    Jmat0 = Jmat0 + sparse(ijs+nn,IRBTJP(ijs),MM(mask_ijs),neq,neq,'mp');
                    % Derivative of recombination term w.r.t. phi (left node)
                    MM = qel.*[SeTJ1.*dRphip(in1) SeTJ2.*dRphip(in2) SeTJ3.*dRphip(in3)];
                    Jmat0 = Jmat0 + sparse(ijs+nn,ILBTJP(ijs),MM(mask_ijs),neq,neq,'mp');
                    %
                    % Assembling in hole continuity equation
                    % Recombination term
                    MM = qel.*[SeTJ1.*GBTBT(in1) SeTJ2.*GBTBT(in2) SeTJ3.*GBTBT(in3)];
                    tvet = tvet + sparse(ijs+pp,1,MM(mask_ijs),neq,1,'mp');
                    %
                    % Derivative of recombination term w.r.t. holes
                    MM = qel.*[SeTJ1.*dRpTJ(in1) SeTJ2.*dRpTJ(in2) SeTJ3.*dRpTJ(in3)];
                    Jmat0 = Jmat0 + sparse(ijs+pp,pp+ILBTJP(ijs),MM(mask_ijs),neq,neq,'mp');
                    % Derivative of recombination term w.r.t. electrons
                    MM = qel.*[SeTJ1.*dRnTJ(in1) SeTJ2.*dRnTJ(in2) SeTJ3.*dRnTJ(in3)];
                    Jmat0 = Jmat0 + sparse(ijs+pp,nn+IRBTJP(ijs),MM(mask_ijs),neq,neq,'mp');
                    % Derivative of recombination term w.r.t. phi (left node)
                    MM = qel.*[SeTJ1.*dRphin(in1) SeTJ2.*dRphin(in2) SeTJ3.*dRphin(in3)];
                    Jmat0 = Jmat0 + sparse(ijs+pp,IRBTJP(ijs),MM(mask_ijs),neq,neq,'mp');
                    % Derivative of recombination term w.r.t. phi (left node)
                    MM = qel.*[SeTJ1.*dRphip(in1) SeTJ2.*dRphip(in2) SeTJ3.*dRphip(in3)];
                    Jmat0 = Jmat0 + sparse(ijs+pp,ILBTJP(ijs),MM(mask_ijs),neq,neq,'mp');
                end
            end
        end
        mode.RSRH = double(RSRH);
        mode.Rrad = double(Rrad);
        mode.RAuger = double(RAuger);
        mode.GBTBT = double(GBTBT);
    end
    % =============================================================================================100
    % assem Poisson equation
    % =============================================================================================100
    % #############################################################################################100
    if((isfield(mode,'ionization'))&&(strcmp(mode.ionization,'Incomplete')))
        gD = two; n1 = Nc.*exp(-mp(mesh.DeltaEd)./Vt);
        dop_dp  =    dop_d./( one + gD.*nF./(n1.*gamman)); % ionized donor density, 1/cm^3
        ddop_dp =  - dop_d./((one + gD.*nF./(n1.*gamman)).^two) ...
            .*gD.*((n1.*gamman - nF.*n1.*dgamman)./((n1.*gamman).^two));
        dop_dp(dop_d==0)=zero;
        ddop_dp(dop_d==0)=zero;
        %
        gA = mp('4'); p1 = Nv.*exp(-mp(mesh.DeltaEa)./Vt);
        dop_am  =    dop_a./( one + gA.*pF./(p1.*gammap)); % ionized donor density, 1/cm^3
        ddop_am =  - dop_a./((one + gA.*pF./(p1.*gammap)).^two) ...
            .*gA.*((p1.*gammap - pF.*p1.*dgammap)./((p1.*gammap).^two));
        dop_am(dop_a==0)=zero;
        ddop_am(dop_a==0)=zero;
        %
    else 
        dop_dp = dop_d;
        dop_am = dop_a;
    end
    %
    % save ionized dopant profiles
    mode.dop_dp = double(dop_dp);
    mode.dop_am = double(dop_am);
    % =============================================================================================100
    it=ismember(triangle(4,:),find(geom.reg2contact)); % triangles inside conductors
    M11 = ( s3./l3.*G3 + s2./l2.*G2).*epsxx; M11(it)=zero;
    M12 = (-s3./l3.*G3             ).*epsxx; M12(it)=zero;
    M13 = (-s2./l2.*G2             ).*epsxx; M13(it)=zero;
    M22 = ( s3./l3.*G3 + s1./l1.*G1).*epsxx; M22(it)=zero;
    M23 = (-s1./l1.*G1             ).*epsxx; M23(it)=zero;
    M33 = ( s1./l1.*G1 + s2./l2.*G2).*epsxx; M33(it)=zero;
    MM=[M11 M12 M13 M12 M22 M23 M13 M23 M33];
    %
    Kmat0=sparse(iir,jjr,MM(mask_iir),neq,neq,'mp');
    % =============================================================================================100
    %
    % doping charge in Poisson rhs
    MM=mpqelNorm.*[Se1.*dop_dp(in1) Se2.*dop_dp(in2) Se3.*dop_dp(in3)];
    tvet=tvet+sparse(ijr,1,-MM(mask_ijr),neq,1,'mp');
    %
    MM=mpqelNorm.*[Se1.*dop_am(in1) Se2.*dop_am(in2) Se3.*dop_am(in3)];
    tvet=tvet+sparse(ijr,1,MM(mask_ijr),neq,1,'mp');
    %
    % incomplete ionization
    if((isfield(mode,'ionization'))&&(strcmp(mode.ionization,'Incomplete')))
        if(mode.nflg)
            MM=mpqelNorm.*[ddop_dp(in1).*Se1 ddop_dp(in2).*Se2 ddop_dp(in3).*Se3];
            Jmat0=Jmat0+sparse(ijr,ijr+nn,-MM(mask_ijr),neq,neq,'mp'); 
        end
        if(mode.pflg)
            MM=mpqelNorm.*[ddop_am(in1).*Se1 ddop_am(in2).*Se2 ddop_am(in3).*Se3];
            Jmat0=Jmat0+sparse(ijr,ijr+pp,MM(mask_ijr),neq,neq,'mp'); 
        end
    end
    % =============================================================================================100
    % elec charge in Poisson rhs
    if(mode.nflg), MM=mpqelNorm.*[Se1 Se2 Se3];
        Kmat0 = Kmat0 + sparse(ijr,ijr+nn,MM(mask_ijr),neq,neq,'mp'); 
    end
    % =============================================================================================100
    % hole charge in Poisson rhs
    if(mode.pflg), MM=-mpqelNorm.*[Se1 Se2 Se3];
        Kmat0 = Kmat0 + sparse(ijr,ijr+pp,MM(mask_ijr),neq,neq,'mp'); 
    end
    % #############################################################################################100
    % =============================================================================================100
    % assem electron continuity equation
    % =============================================================================================100
    if(mode.nflg) % ##############################################################################100
        %
        % Bernoulli functions
        delta12 = (phi1 + offset_n1 + Vt1.*log(gamman(in1)))./Vt1 - (phi2 + offset_n2 + Vt2.*log(gamman(in2)))./Vt2;
        delta23 = (phi2 + offset_n2 + Vt2.*log(gamman(in2)))./Vt2 - (phi3 + offset_n3 + Vt3.*log(gamman(in3)))./Vt3;
        delta31 = (phi3 + offset_n3 + Vt3.*log(gamman(in3)))./Vt3 - (phi1 + offset_n1 + Vt1.*log(gamman(in1)))./Vt1;
        %
        B12=G3.*mpbern(delta12);  dB12=G3.*mpdbern(delta12)./Vt1;
        B21=G3.*mpbern(-delta12); dB21=G3.*mpdbern(-delta12)./Vt2;
        B23=G1.*mpbern(delta23);  dB23=G1.*mpdbern(delta23)./Vt2;
        B32=G1.*mpbern(-delta23); dB32=G1.*mpdbern(-delta23)./Vt3;
        B31=G2.*mpbern(delta31);  dB31=G2.*mpdbern(delta31)./Vt3;
        B13=G2.*mpbern(-delta31); dB13=G2.*mpdbern(-delta31)./Vt1;
        % =============================================================================================100
        dDndphi1=mp(dDndphi1).*Vt1; dDndphi2=mp(dDndphi2).*Vt2; dDndphi3=mp(dDndphi3).*Vt3;
        % =============================================================================================100
        M11=Dn.*(s3./l3.*B12+s2./l2.*B13);
        M12=-Dn.*s3./l3.*B21;
        M13=-Dn.*s2./l2.*B31;
        M21=-Dn.*s3./l3.*B12;
        M22=Dn.*(s1./l1.*B23+s3./l3.*B21);
        M23=-Dn.*s1./l1.*B32;
        M31=-Dn.*s2./l2.*B13;
        M32=-Dn.*s1./l1.*B23;
        M33=Dn.*(s2./l2.*B31+s1./l1.*B32);
        MM=qel.*[M11 M12 M13 M21 M22 M23 M31 M32 M33];
        %
        Kmat0=Kmat0+sparse(iis+nn,jjs+nn,MM(mask_iis),neq,neq,'mp'); % Assembly electron eqs.
        Kmat0=Kmat0+sparse(iit   ,jjt+nn,-MM(mask_iit),neq,neq,'mp'); % Assembly current eq.
        % derivatives w.r.t. potential ================================================================100
        M11=Dn.*((s3./l3.*dB12+s2./l2.*dB13).*elec1+s3./l3.*dB21.*elec2+s2./l2.*dB31.*elec3);
        M12=Dn.*(-s3./l3.*dB12.*elec1-s3./l3.*dB21.*elec2);
        M13=Dn.*(-s2./l2.*dB13.*elec1-s2./l2.*dB31.*elec3);
        M21=Dn.*(-s3./l3.*dB12.*elec1-s3./l3.*dB21.*elec2);
        M22=Dn.*(s3./l3.*dB12.*elec1+(s1./l1.*dB23+s3./l3.*dB21).*elec2+s1./l1.*dB32.*elec3);
        M23=Dn.*(-s1./l1.*dB23.*elec2-s1./l1.*dB32.*elec3);
        M31=Dn.*(-s2./l2.*dB13.*elec1-s2./l2.*dB31.*elec3);
        M32=Dn.*(-s1./l1.*dB23.*elec2-s1./l1.*dB32.*elec3);
        M33=Dn.*(s2./l2.*dB13.*elec1+s1./l1.*dB23.*elec2+(s2./l2.*dB31+s1./l1.*dB32).*elec3);
        MM=qel.*[M11 M12 M13 M21 M22 M23 M31 M32 M33];
        %
        Jmat0=Jmat0+sparse(iis+nn,jjs,MM(mask_iis),neq,neq,'mp'); % Assembly electron eqs.
        Jmat0=Jmat0+sparse(iit   ,jjt,-MM(mask_iit),neq,neq,'mp'); % Assembly current eq.
        % Fermi statistics ============================================================================100
        if((isfield(mode,'stats'))&&(strcmp(mode.stats,'Fermi')))
            M11=Dn.*((s3./l3.*dB12+s2./l2.*dB13).*elec1+s3./l3.*dB21.*elec2+s2./l2.*dB31.*elec3);
            M12=Dn.*(-s3./l3.*dB12.*elec1-s3./l3.*dB21.*elec2);
            M13=Dn.*(-s2./l2.*dB13.*elec1-s2./l2.*dB31.*elec3);
            M21=Dn.*(-s3./l3.*dB12.*elec1-s3./l3.*dB21.*elec2);
            M22=Dn.*(s3./l3.*dB12.*elec1+(s1./l1.*dB23+s3./l3.*dB21).*elec2+s1./l1.*dB32.*elec3);
            M23=Dn.*(-s1./l1.*dB23.*elec2-s1./l1.*dB32.*elec3);
            M31=Dn.*(-s2./l2.*dB13.*elec1-s2./l2.*dB31.*elec3);
            M32=Dn.*(-s1./l1.*dB23.*elec2-s1./l1.*dB32.*elec3);
            M33=Dn.*(s2./l2.*dB13.*elec1+s1./l1.*dB23.*elec2+(s2./l2.*dB31+s1./l1.*dB32).*elec3);
            MM=qel.*[Vt1.*M11.*dlGn1 Vt2.*M12.*dlGn2 Vt3.*M13.*dlGn3 ...
                Vt1.*M21.*dlGn1 Vt2.*M22.*dlGn2 Vt3.*M23.*dlGn3 ...
                Vt1.*M31.*dlGn1 Vt2.*M32.*dlGn2 Vt3.*M33.*dlGn3];
            %
            Jmat0=Jmat0+sparse(iis+nn,jjs+nn,MM(mask_iis),neq,neq,'mp'); % Assembly electron eqs.
            Jmat0=Jmat0+sparse(iit   ,jjt+nn,-MM(mask_iit),neq,neq,'mp'); % Assembly current eq.
        end
        
        % =============================================================================================100
        M11=((s3./l3.*B12+s2./l2.*B13).*elec1-s3./l3.*B21.*elec2-s2./l2.*B31.*elec3).*dDndphi1;
        M12=((s3./l3.*B12+s2./l2.*B13).*elec1-s3./l3.*B21.*elec2-s2./l2.*B31.*elec3).*dDndphi2;
        M13=((s3./l3.*B12+s2./l2.*B13).*elec1-s3./l3.*B21.*elec2-s2./l2.*B31.*elec3).*dDndphi3;
        M21=(-s3./l3.*B12.*elec1+(s1./l1.*B23+s3./l3.*B21).*elec2-s1./l1.*B32.*elec3).*dDndphi1;
        M22=(-s3./l3.*B12.*elec1+(s1./l1.*B23+s3./l3.*B21).*elec2-s1./l1.*B32.*elec3).*dDndphi2;
        M23=(-s3./l3.*B12.*elec1+(s1./l1.*B23+s3./l3.*B21).*elec2-s1./l1.*B32.*elec3).*dDndphi3;
        M31=(-s2./l2.*B13.*elec1-s1./l1.*B23.*elec2+(s2./l2.*B31+s1./l1.*B32).*elec3).*dDndphi1;
        M32=(-s2./l2.*B13.*elec1-s1./l1.*B23.*elec2+(s2./l2.*B31+s1./l1.*B32).*elec3).*dDndphi2;
        M33=(-s2./l2.*B13.*elec1-s1./l1.*B23.*elec2+(s2./l2.*B31+s1./l1.*B32).*elec3).*dDndphi3;
        MM=qel.*[M11 M12 M13 M21 M22 M23 M31 M32 M33];
        %
        Jmat0=Jmat0+sparse(iis+nn,jjs,MM(mask_iis),neq,neq,'mp'); % Assembly electron eqs.
        Jmat0=Jmat0+sparse(iit   ,jjt,-MM(mask_iit),neq,neq,'mp'); % Assembly current eq.
        % =============================================================================================100
        % direct recombination
        MM=qel.*[Se1.*R(in1) Se2.*R(in2) Se3.*R(in3)];
        tvet=tvet+sparse(ijs+nn,1,MM(mask_ijs),neq,1,'mp');
        % direct recombination, derivatives
        MM=qel.*[Se1.*dRn(in1) Se2.*dRn(in2) Se3.*dRn(in3)];
        Jmat0=Jmat0+sparse(ijs+nn,ijs+nn,MM(mask_ijs),neq,neq,'mp');
        if(mode.pflg)
            MM=qel.*[Se1.*dRp(in1) Se2.*dRp(in2) Se3.*dRp(in3)];
            Jmat0=Jmat0+sparse(ijs+nn,ijs+pp,MM(mask_ijs),neq,neq,'mp');
        end
        if(mode.AugerExpulsion3D==1)
            RAugerExpulsion=zero*ones(1,nn);
            dRAugerExpulsion_n3D=zero*ones(1,nn);
            dRAugerExpulsion_p3D=zero*ones(1,nn);
            RAugerExpulsion(mesh.ICAV)=RAuger_n3D(mesh.ICAV);
            dRAugerExpulsion_n3D(mesh.ICAV)=dRAuger_n3D_n3D(mesh.ICAV);
            dRAugerExpulsion_p3D(mesh.ICAV)=dRAuger_n3D_p3D(mesh.ICAV);
            %
            MM=qel.*[SeCAV1.*RAugerExpulsion(in1) SeCAV2.*RAugerExpulsion(in2) SeCAV3.*RAugerExpulsion(in3)];
            tvet=tvet+sparse(ijs+nn,1,MM(mask_ijs),neq,1,'mp');
            MM=qel.*[SeCAV1.*dRAugerExpulsion_n3D(in1) SeCAV2.*dRAugerExpulsion_n3D(in2) SeCAV3.*dRAugerExpulsion_n3D(in3)];
            Jmat0=Jmat0+sparse(ijs+nn,ijs+nn,MM(mask_ijs),neq,neq,'mp');
            if(mode.pflg)
                MM=qel.*[SeCAV1.*dRAugerExpulsion_p3D(in1) SeCAV2.*dRAugerExpulsion_p3D(in2) SeCAV3.*dRAugerExpulsion_p3D(in3)];
                Jmat0=Jmat0+sparse(ijs+nn,ijs+pp,MM(mask_ijs),neq,neq,'mp');
            end
        end
        % =============================================================================================100
        % Schottky contacts
        % =============================================================================================100
        for ic=1:geom.nc
            ii=find((mesh.contact==ic) & (mesh.iq)); % semiconductor nodes on contact ic
            nq=length(ii);
            if((geom.contact_type(ic)==2) && (nq>0))
                jj=1:nq;
                ll(jj)=mp(mesh.contact_seg(ic,ii));
                vsurfn(jj)=mp(mesh.vsurfn(ii));
                ii=find((contact==ic) & iq);
                n_eq(1:nq)=nint(ii).*exp(phi(ii)-v_ls((ic-one)+one)./Vt(ii));
                tvet(ii+nn)=tvet(ii+nn)-(vsurfn.*n_eq.*ll).';
                Kmat0=Kmat0+sparse(ii+nn,ii+nn,vsurfn.*ll,neq,neq,'mp');
            end
        end
        % =============================================================================================100
    end % #########################################################################################100
    % =============================================================================================100
    % assem hole continuity equation
    % =============================================================================================100
    if(mode.pflg) % ##############################################################################100
        %
        % Bernoulli functions
        delta12 = (phi1 + offset_p1 - Vt1.*log(gammap(in1)))./Vt1 - (phi2 + offset_p2 - Vt2.*log(gammap(in2)))./Vt2;
        delta23 = (phi2 + offset_p2 - Vt2.*log(gammap(in2)))./Vt2 - (phi3 + offset_p3 - Vt3.*log(gammap(in3)))./Vt3;
        delta31 = (phi3 + offset_p3 - Vt3.*log(gammap(in3)))./Vt3 - (phi1 + offset_p1 - Vt1.*log(gammap(in1)))./Vt1;
        %
        B12=G3.*mpbern(delta12);  dB12=G3.*mpdbern(delta12)./Vt2;
        B21=G3.*mpbern(-delta12); dB21=G3.*mpdbern(-delta12)./Vt1;
        B23=G1.*mpbern(delta23);  dB23=G1.*mpdbern(delta23)./Vt3;
        B32=G1.*mpbern(-delta23); dB32=G1.*mpdbern(-delta23)./Vt2;
        B31=G2.*mpbern(delta31);  dB31=G2.*mpdbern(delta31)./Vt1;
        B13=G2.*mpbern(-delta31); dB13=G2.*mpdbern(-delta31)./Vt3;
        % =============================================================================================100
        dDpdphi1=mp(dDpdphi1).*Vt1; dDpdphi2=mp(dDpdphi2).*Vt2; dDpdphi3=mp(dDpdphi3).*Vt3;
        % =============================================================================================100
        M11=Dp.*(s3./l3.*B21+s2./l2.*B31);
        M12=-Dp.*s3./l3.*B12;
        M13=-Dp.*s2./l2.*B13;
        M21=-Dp.*s3./l3.*B21;
        M22=Dp.*(s1./l1.*B32+s3./l3.*B12);
        M23=-Dp.*s1./l1.*B23;
        M31=-Dp.*s2./l2.*B31;
        M32=-Dp.*s1./l1.*B32;
        M33=Dp.*(s2./l2.*B13+s1./l1.*B23);
        MM=qel.*[M11 M12 M13 M21 M22 M23 M31 M32 M33];
        %
        Kmat0=Kmat0+sparse(iis+pp,jjs+pp,MM(mask_iis),neq,neq,'mp'); % Assembly hole eqs.
        Kmat0=Kmat0+sparse(iit   ,jjt+pp,MM(mask_iit),neq,neq,'mp'); % Assembly current eq.
        % derivatives w.r.t. potential ================================================================100
        M11=Dp.*((-s3./l3.*dB21-s2./l2.*dB31).*hole1-s3./l3.*dB12.*hole2-s2./l2.*hole3.*dB13);
        M12=Dp.*(s3./l3.*dB21.*hole1+s3./l3.*dB12.*hole2);
        M13=Dp.*(s2./l2.*dB31.*hole1+s2./l2.*hole3.*dB13);
        M21=Dp.*(s3./l3.*dB21.*hole1+s3./l3.*dB12.*hole2);
        M22=Dp.*(-s3./l3.*dB21.*hole1+(-s1./l1.*dB32-s3./l3.*dB12).*hole2-s1./l1.*hole3.*dB23);
        M23=Dp.*(s1./l1.*dB32.*hole2+s1./l1.*hole3.*dB23);
        M31=Dp.*(s2./l2.*dB31.*hole1+s2./l2.*hole3.*dB13);
        M32=Dp.*(s1./l1.*dB32.*hole2+s1./l1.*hole3.*dB23);
        M33=Dp.*(-s2./l2.*dB31.*hole1-s1./l1.*dB32.*hole2+(-s2./l2.*dB13-s1./l1.*dB23).*hole3);
        MM=qel.*[M11 M12 M13 M21 M22 M23 M31 M32 M33];
        %
        Jmat0=Jmat0+sparse(iis+pp,jjs,MM(mask_iis),neq,neq,'mp');
        Jmat0=Jmat0+sparse(iit   ,jjt,MM(mask_iit),neq,neq,'mp');
        % Fermi statistics ============================================================================100
        if((isfield(mode,'stats'))&&(strcmp(mode.stats,'Fermi')))
            M11=Dp.*((-s3./l3.*dB21-s2./l2.*dB31).*hole1-s3./l3.*dB12.*hole2-s2./l2.*hole3.*dB13);
            M12=Dp.*(s3./l3.*dB21.*hole1+s3./l3.*dB12.*hole2);
            M13=Dp.*(s2./l2.*dB31.*hole1+s2./l2.*hole3.*dB13);
            M21=Dp.*(s3./l3.*dB21.*hole1+s3./l3.*dB12.*hole2);
            M22=Dp.*(-s3./l3.*dB21.*hole1+(-s1./l1.*dB32-s3./l3.*dB12).*hole2-s1./l1.*hole3.*dB23);
            M23=Dp.*(s1./l1.*dB32.*hole2+s1./l1.*hole3.*dB23);
            M31=Dp.*(s2./l2.*dB31.*hole1+s2./l2.*hole3.*dB13);
            M32=Dp.*(s1./l1.*dB32.*hole2+s1./l1.*hole3.*dB23);
            M33=Dp.*(-s2./l2.*dB31.*hole1-s1./l1.*dB32.*hole2+(-s2./l2.*dB13-s1./l1.*dB23).*hole3);
            MM=-qel.*[Vt1.*M11.*dlGp1 Vt2.*M12.*dlGp2 Vt3.*M13.*dlGp3 ...
                Vt1.*M21.*dlGp1 Vt2.*M22.*dlGp2 Vt3.*M23.*dlGp3 ...
                Vt1.*M31.*dlGp1 Vt2.*M32.*dlGp2 Vt3.*M33.*dlGp3];
            %
            Jmat0=Jmat0+sparse(iis+pp,jjs+pp,MM(mask_iis),neq,neq,'mp'); % Assembly hole eqs.
            Jmat0=Jmat0+sparse(iit   ,jjt+pp,MM(mask_iit),neq,neq,'mp'); % Assembly current eq.
        end
        % =============================================================================================100
        M11=((s3./l3.*B21+s2./l2.*B31).*hole1-s3./l3.*B12.*hole2-s2./l2.*hole3.*B13).*dDpdphi1;
        M12=((s3./l3.*B21+s2./l2.*B31).*hole1-s3./l3.*B12.*hole2-s2./l2.*hole3.*B13).*dDpdphi2;
        M13=((s3./l3.*B21+s2./l2.*B31).*hole1-s3./l3.*B12.*hole2-s2./l2.*hole3.*B13).*dDpdphi3;
        M21=(-s3./l3.*B21.*hole1+(s1./l1.*B32+s3./l3.*B12).*hole2-s1./l1.*hole3.*B23).*dDpdphi1;
        M22=(-s3./l3.*B21.*hole1+(s1./l1.*B32+s3./l3.*B12).*hole2-s1./l1.*hole3.*B23).*dDpdphi2;
        M23=(-s3./l3.*B21.*hole1+(s1./l1.*B32+s3./l3.*B12).*hole2-s1./l1.*hole3.*B23).*dDpdphi3;
        M31=(-s2./l2.*B31.*hole1-s1./l1.*B32.*hole2+(s2./l2.*B13+s1./l1.*B23).*hole3).*dDpdphi1;
        M32=(-s2./l2.*B31.*hole1-s1./l1.*B32.*hole2+(s2./l2.*B13+s1./l1.*B23).*hole3).*dDpdphi2;
        M33=(-s2./l2.*B31.*hole1-s1./l1.*B32.*hole2+(s2./l2.*B13+s1./l1.*B23).*hole3).*dDpdphi3;
        MM=qel.*[M11 M12 M13 M21 M22 M23 M31 M32 M33];
        %
        Jmat0=Jmat0+sparse(iis+pp,jjs,MM(mask_iis),neq,neq,'mp'); % Assembly hole eqs.
        Jmat0=Jmat0+sparse(iit   ,jjt,MM(mask_iit),neq,neq,'mp'); % Assembly current eq.
        % =============================================================================================100
        % direct recombination, derivatives
        MM=qel.*[Se1.*R(in1) Se2.*R(in2) Se3.*R(in3)];
        tvet=tvet+sparse(ijs+pp,1,MM(mask_ijs),neq,1,'mp');
        % direct recombination, derivatives
        MM=qel.*[Se1.*dRp(in1) Se2.*dRp(in2) Se3.*dRp(in3)];
        Jmat0=Jmat0+sparse(ijs+pp,ijs+pp,MM(mask_ijs),neq,neq,'mp');
        if(mode.nflg)
            MM=qel.*[Se1.*dRn(in1) Se2.*dRn(in2) Se3.*dRn(in3)];
            Jmat0=Jmat0+sparse(ijs+pp,ijs+nn,MM(mask_ijs),neq,neq,'mp');
        end
        if(mode.AugerExpulsion3D==1)
            RAugerExpulsion=zero*ones(1,nn);
            dRAugerExpulsion_p3D=zero*ones(1,nn);
            dRAugerExpulsion_n3D=zero*ones(1,nn);
            RAugerExpulsion(mesh.ICAV)=RAuger_p3D(mesh.ICAV);
            dRAugerExpulsion_p3D(mesh.ICAV)=dRAuger_p3D_p3D(mesh.ICAV);
            dRAugerExpulsion_n3D(mesh.ICAV)=dRAuger_p3D_n3D(mesh.ICAV);
            %
            MM=qel.*[SeCAV1.*RAugerExpulsion(in1) SeCAV2.*RAugerExpulsion(in2) SeCAV3.*RAugerExpulsion(in3)];
            tvet=tvet+sparse(ijs+pp,1,MM(mask_ijs),neq,1,'mp');
            MM=qel.*[SeCAV1.*dRAugerExpulsion_p3D(in1) SeCAV2.*dRAugerExpulsion_p3D(in2) SeCAV3.*dRAugerExpulsion_p3D(in3)];
            Jmat0=Jmat0+sparse(ijs+pp,ijs+pp,MM(mask_ijs),neq,neq,'mp');
            if(mode.nflg)
                MM=qel.*[SeCAV1.*dRAugerExpulsion_n3D(in1) SeCAV2.*dRAugerExpulsion_n3D(in2) SeCAV3.*dRAugerExpulsion_n3D(in3)];
                Jmat0=Jmat0+sparse(ijs+pp,ijs+nn,MM(mask_ijs),neq,neq,'mp');
            end
        end
        % =============================================================================================100
        % Schottky contacts
        for ic=1:geom.nc
            ii=find((mesh.contact==ic)&(mesh.iq)); % semiconductor nodes on contact ic
            nq=length(ii);
            if((geom.contact_type(ic)==2)&&(nq>0))
                ll(1:nq) = mp(mesh.contact_seg(ic,ii));
                vsurfp(1:nq) = mp(mesh.vsurfp(ii));
                ii=find((contact==ic)&iq);
                jj=ii;
                p_eq(1:nq) = nint(ii).*exp(-phi(ii)+v_ls((ic-one)+one)./Vt(ii));
                tvet(ii+pp)=tvet(ii+pp)-(vsurfp.*p_eq.*ll).';
                Kmat0=Kmat0+sparse(ii+pp,ii+pp,vsurfp.*ll,neq,neq,'mp');
            end
        end
        % =============================================================================================100
    end % #########################################################################################100
    % =============================================================================================100
    % =============================================================================================100
    % Add circuit equations: V-Z*I-V0=0 (V-driven); I-Y*V-I0=0 (I-driven)
    % =============================================================================================100
    Zmat=mp(sparse(mode.Zmat));
    Jmat2=mp(sparse(neq,neq));
    
    if mode.IdriveON==1   % current driving
        Jmat2=Jmat2+sparse(qq+1,rr+1,one,neq,neq,'mp'); % I
        %
        for ic=1:nm
            ii=(ic-1)+1;
            for jc=1:nm
                jj=(jc-1)+1;
                Jmat2(qq+ii,qq+jj)=Ymat(ii,jj);
            end
        end % -Y*V
        tvet(qq+(1:nm))=-i0_dd; % -V0
    else    % voltage driving
        Jmat2=Jmat2+sparse(qq+1,qq+1,one,neq,neq,'mp'); % V
        %
        for ic=1:nm
            ii=(ic-1)+1;
            for jc=1:nm
                jj=(jc-1)+1;
                Jmat2(qq+ii,rr+jj)=Zmat(ii,jj);
            end
        end % -Z*I
        tvet(qq+(1:nm))=-v0_dd; % -V0
    end
    %
    % =============================================================================================100
    % Add current equations
    % =============================================================================================100
    % Using -qel instead of -1 because of equations normalization by qel
    Kmat0=Kmat0+sparse(rr+(1:nm),rr+(1:nm),-one,neq,neq,'mp');
    %
    % =============================================================================================100
    % Assembling quantum corrections and related jacobians
    % =============================================================================================100
    if(mode.oflg)
        iTappo=mode.iTappo;
        N2Dtot=zero; P2Dtot=zero; Psp=zero;
        mode.gmod=zero*ones(nmodes,1);
        mode.matgain=0;
        %        mode.DeltaN=0;
        mode.RSRHQW=zeros(1,mesh.nn);
        mode.RradQW=zeros(1,mesh.nn);
        mode.RAugerQW=zeros(1,mesh.nn);
        mode.RLeakageQW=zeros(1,mesh.nn);
        mode.RradQW_SBE=zeros(1,mesh.nn);
        mode.Ccapn3D=zeros(1,mesh.nn);
        mode.Ccapp3D=zeros(1,mesh.nn);
        IntCcapp=[zero zero];
        IntCcapn=[zero zero];
        GaSum=zero;
        IntSRH=zero;
        IntAug=zero;
        IntStim=zero;
        IntRad=zero;
        IntLea=zero;
        IntRecQW3Dn=zero; % integral over all the QW (3D) recombinations
        IntRecQW3Dp=zero; % integral over all the QW (3D) recombinations
        ImatQW = sparse(neq,neq);
        tvetQW = sparse(neq,1);
        MMgain=zero;
        for indQW=1:NQW
            %
            % Basic parameters: geometry
            WQW = mp(mesh.vWMQW{indQW}); %-- cm
            IQW = unique(mesh.IMQW{indQW}); % questo unique metterlo in rectmesh
            inQW = mesh.inMQW{indQW}; % index to QW nodes
            inQWP = mesh.inMQWP{indQW}; % including passivation
            VV = vv+(indQW-1)*length(mesh.xgrid);
            WW = ww+(indQW-1)*length(mesh.xgrid);
            %
            iiQW1=1:nnQW-1; iiQW2=2:nnQW;
            maskQW=true(1,nnQW);
            % maskQW([1,end])=zero;
            % maskQW([1,2,end-1,end])=zero;
            %
            SeQW1=zero*ones(1,mesh.nt); SeQW2=SeQW1; SeQW3=SeQW1;
            SeQW1(mesh.ITrMQW{indQW})=Se1(mesh.ITrMQW{indQW});
            SeQW2(mesh.ITrMQW{indQW})=Se2(mesh.ITrMQW{indQW});
            SeQW3(mesh.ITrMQW{indQW})=Se3(mesh.ITrMQW{indQW});
            ijrQW=[iiQW1 iiQW2]; inQW1=inQW(iiQW1); inQW2=inQW(iiQW2);
            %
            IIQWP = repmat(1:mesh.nnx,mesh.nny,1); IIQWP=IIQWP(:).'; % for Poisson
            
            INQWP=inQWP(IIQWP); % for Ccap recombinations in 3D continuity equations
            %
            xQW = node(1,inQW); yQW = node(2,inQW);
            xcQW = (xQW(2:end) + xQW(1:end-1))/two; % edge centers, cm
            ycQW = (yQW(2:end) + yQW(1:end-1))/two; % edge centers, cm
            LeQW = xQW(2:end) - xQW(1:end-1); % edge length, cm
            Lp = zero*ones(1,nnQW); % box length, cm
            Lp(1:(nnQW-1)) = LeQW/two; Lp(2:nnQW) = Lp(2:nnQW) + LeQW/two;
            %
            Lp1=Lp(iiQW1)/two; Lp1(1)=two*Lp1(1);
            Lp2=Lp(iiQW2)/two; Lp2(end)=two*Lp2(end);
            Gcyl=1;
            if((isfield(mode,'symmetry'))&&(strcmp(mode.symmetry,'Cylindrical-Y')))
                Lp1=twopi.*Lp1.*xcQW;
                Lp2=twopi.*Lp2.*xcQW;
                Gcyl=twopi.*xcQW;
                yQW=yQW(1);
            end
            %
            if((isfield(mode,'symmetry'))&&(strcmp(mode.symmetry,'Cylindrical-X')))
                Lp1=twopi.*Lp1.*ycQW;
                Lp2=twopi.*Lp2.*ycQW;
                Gcyl=twopi.*ycQW;
                xQW=xQW(1);
            end
            % % debug: check Lp1,Lp2: integral of the constant function on the domain
            % %        this is expected to be equal to the length (z) of the QW
            % sum(diag(sparse([inQW1 inQW2],[inQW1 inQW2],[Lp1 Lp2],nn,nn)),1)
            %
            Vt2D = Vt(inQW);
            %
            tauscatn = mp(mesh.tauscatnMQW{indQW})/mode.FatAdiab; % capture time, s
            tauscatp = mp(mesh.tauscatpMQW{indQW})/mode.FatAdiab; % capture time, s
            meffn2D = mp(mesh.meffnMQW{indQW});
            meffp2D = mp(mesh.meffpMQW{indQW});
            DeltaEc = mp(mesh.DeltaEcQW{indQW});
            DeltaEv = mp(mesh.DeltaEvQW{indQW});
            lambda_C=mp(mesh.lambda_C{indQW});
            lambda_V=mp(mesh.lambda_V{indQW});
            Psi_C2=mp(mesh.Psi_C2{indQW});
            Psi_V2=mp(mesh.Psi_V2{indQW});
            NBound_C=length(lambda_C); % number of conduction band bound states
            NBound_V=length(lambda_V); % number of valence band bound states
            Ec0 = ones(NBound_C,1,'mp')*mode.ecb(inQW);
            Ev0 = ones(NBound_V,1,'mp')*mode.evb(inQW);
            Em = Ec0 - DeltaEc + lambda_C*ones(1,nnQW,'mp'); % eigenvalue referred to bottom of QW
            Ed = Ev0 + DeltaEv - lambda_V*ones(1,nnQW,'mp'); % eigenvalue referred to top of QB
            %        'TAROCCO N2 @@@@@@@@'
            %        N2 =.85* 4*pi/h^2*m0*qel*(sum((meffn2D*ones(1,mesh.nnxQW)).*(Ec0-Em),1)/10000);
            if isfield(mode,'Fat_N2')
             Fat_N2=mp(mode.Fat_N2);            
             Fat_P2=mp(mode.Fat_P2);            
            else
             Fat_N2=one;
             Fat_P2=one;
            end
            
            N2 = Fat_N2 * mp('4*pi')/h^two*m0*qel*(sum((meffn2D*ones(1,mesh.nnxQW{indQW},'mp')).*(Ec0-Em),1)/mp('10000'))/CarrierNorm2D;
            P2 = Fat_P2 * mp('4*pi')/h^two*m0*qel*(sum((meffp2D*ones(1,mesh.nnxQW{indQW},'mp')).*(Ed-Ev0),1)/mp('10000'))/CarrierNorm2D;
            
            Nc2D = mp(mesh.Nc2D{indQW});%*ones(1,nnQW);
            Nv2D = mp(mesh.Nv2D{indQW});%*ones(1,nnQW);
            MVt2D_C=ones(NBound_C,1,'mp')*Vt2D;
            MVt2D_V=ones(NBound_V,1,'mp')*Vt2D;
            %
            %           Dn2Dd = mode.FAT_idiffusioneQW.*mesh.mobnQW{indQW}.*Vt2D(iiQW1); % 2D diffusivity
            %           Dp2Dd = mode.FAT_idiffusioneQW.*mesh.mobpQW{indQW}.*Vt2D(iiQW1);  % 2D diffusivity
            Dn2Dd = mp(mode.FAT_idiffusioneQW_E).*mp(mesh.mobnQW{indQW}).*Vt2D(iiQW1); % 2D diffusivity
            Dp2Dd = mp(mode.FAT_idiffusioneQW_H).*mp(mesh.mobpQW{indQW}).*Vt2D(iiQW1);  % 2D diffusivity
            
            Damb=one./(one./Dn2Dd+one./Dp2Dd);
            %
            if(mode.iambipolarQW==1)
                Dn2D = Damb;
                Dp2D = Damb;
            else
                Dn2D = Dn2Dd;
                Dp2D = Dp2Dd;
            end
            %
            EFn3D = ones(NBound_C,1,'mp')*mode.EFn(inQW);
            EFp3D = ones(NBound_V,1,'mp')*mode.EFp(inQW);
            %
            iivQW=[iiQW1 iiQW1 iiQW2 iiQW2]; jjvQW=[iiQW1 iiQW2 iiQW1 iiQW2];
            invQW=[inQW1 inQW1 inQW2 inQW2]; jnvQW=[inQW1 inQW2 inQW1 inQW2];
            iivQWcur=[iiQW1 iiQW1]; jjvQWcur=[iiQW1 iiQW2];
            mask_iivQW=maskQW(iivQW);
            mask_iivQWcur=maskQW(iivQWcur);
            %
            n2D=n2Dc{indQW}; p2D=p2Dc{indQW};
            n2Di=mp(mode.n2Di{indQW}); p2Di=mp(mode.p2Di{indQW});
            %
            % checking if n2D>N2 or p2D>P2
            % indNaN=find(isnan(n2D));
            % if(not(isempty(indNaN))),n2D(indNaN)=mean(n2D); end

            if(iTappo==1 | iTappo==3 | iTappo==6)
                indN2=find(n2D>=N2);if(not(isempty(indN2))),n2D(indN2)=mp('0.99')*N2(indN2);end
                % indNaN=find(isnan(p2D));
                % if(not(isempty(indNaN))),p2D(indNaN)=mean(p2D); end
                indP2=find(p2D>=P2);if(not(isempty(indP2))),p2D(indP2)=mp('0.99')*P2(indP2);end
            end
            %
            % Evaluation of 2D quasi-Fermi levels
            % 2D electrons
            
            if(iTappo==1 | iTappo==3 | iTappo==6 | iTappo==7)
                EFn2D=mpinvferdr_n2D(n2D,Em,Ec0,Nc2D,MVt2D_C,EFp3D(1,:)-mp('20'),EFn3D(1,:)+mp('20'));
            elseif(iTappo==0 | iTappo==2 | iTappo==4  | iTappo==8)
                EFn2D=mpinvferdr_n2D_0(n2D,Em,Ec0,Nc2D,MVt2D_C,EFp3D(1,:)-mp('20'),EFn3D(1,:)+mp('20'));
            end
            
            xnum=(ones(NBound_C,1,'mp')*EFn2D-Em)./MVt2D_C; xden=(ones(NBound_C,1,'mp')*EFn2D-Ec0)./MVt2D_C;
            if(iTappo==1 | iTappo==3 | iTappo==6 | iTappo==7)
                [n2Dm,dn2D_2] = mpferdr2D(xnum,xden);
            elseif(iTappo==0 | iTappo==2 | iTappo==4  | iTappo==8)
                [n2Dm,dn2D_2] = mpferdr2D_0(xnum,xden);
            end
            n2Dm=Nc2D.*n2Dm;
            dn2Dm_EFn2D=Nc2D.*dn2D_2./MVt2D_C;
            dn2D_EFn2D=sum(dn2Dm_EFn2D,1);
            %
            % 2D holes
            if(iTappo==1 | iTappo==3  | iTappo==6 | iTappo==7)
                EFp2D=mpinvferdr_p2D(p2D,Ed,Ev0,Nv2D,MVt2D_V,EFp3D(1,:)-mp('20'),EFn3D(1,:)+mp('20'));
            elseif(iTappo==0 | iTappo==2 | iTappo==4  | iTappo==8)
                EFp2D=mpinvferdr_p2D_0(p2D,Ed,Ev0,Nv2D,MVt2D_V,EFp3D(1,:)-mp('20'),EFn3D(1,:)+mp('20'));
            end
            xnum=(Ed-ones(NBound_V,1,'mp')*EFp2D)./MVt2D_V; xden=(Ev0-ones(NBound_V,1,'mp')*EFp2D)./MVt2D_V;
            if(iTappo==1 | iTappo==3  | iTappo==6 | iTappo==7)
                [p2Dm,dp2D_2] = mpferdr2D(xnum,xden);
            elseif(iTappo==0 | iTappo==2 | iTappo==4  | iTappo==8)
                [p2Dm,dp2D_2] = mpferdr2D_0(xnum,xden);
            end
            p2Dm=Nv2D.*p2Dm;
            dp2Dm_EFp2D=-Nv2D.*dp2D_2./MVt2D_V;
            dp2D_EFp2D=sum(dp2Dm_EFp2D,1);
            %
            % Assembing 2D electrons charge in Poisson equation
            n2DmPois = zero*ones(size(Nc2D,1),length(mesh.xgrid));
            dn2Dm_EFn2DPois = n2DmPois;
            n2DmPois(:,1:size(n2Dm,2))=n2Dm;
            dn2Dm_EFn2DPois(:,1:size(dn2Dm_EFn2D,2))=dn2Dm_EFn2D;
            N2D = reshape(Psi_C2.'*n2DmPois,1,nn);
            N2Dtot = N2Dtot + N2D;
            dN2D_EFn2D = reshape(Psi_C2.'*dn2Dm_EFn2DPois,1,nn);
            dn2D_EFn2Dvec = reshape(ones(NBound_C,mesh.nny,'mp').'*dn2Dm_EFn2DPois,1,nn);
            dN2D_n2D = dN2D_EFn2D./dn2D_EFn2Dvec;
            dN2D_n2D(isnan(dN2D_n2D)) = zero;
            % qelNorm (not qelNorm2D) because N2D is multiplied by Psi
            % (1/cm)
            MM = mpqelNorm2D.*[Se1.*N2D(in1) Se2.*N2D(in2) Se3.*N2D(in3)];
            tvet = tvet + sparse(ijr,1,MM(mask_ijr),neq,1,'mp');
            MM = mpqelNorm2D.*[Se1.*dN2D_n2D(in1) Se2.*dN2D_n2D(in2) Se3.*dN2D_n2D(in3)];
            Jmat0 = Jmat0 + sparse(ijr,VV+IIQWP(ijr),MM(mask_ijr),neq,neq,'mp');
            %
            % Assembing 2D holes charge in Poisson equation
            p2DmPois = zero*ones(size(Nv2D,1),length(mesh.xgrid),'mp');
            dp2Dm_EFp2DPois = p2DmPois;
            p2DmPois(:,1:size(p2Dm,2))=p2Dm;
            dp2Dm_EFp2DPois(:,1:size(dp2Dm_EFp2D,2))=dp2Dm_EFp2D;
            P2D = reshape(Psi_V2.'*p2DmPois,1,nn);
            P2Dtot = P2Dtot + P2D;
            dP2D_EFp2D = reshape(Psi_V2.'*dp2Dm_EFp2DPois,1,nn);
            dp2D_EFp2Dvec = reshape(ones(NBound_V,mesh.nny,'mp').'*dp2Dm_EFp2DPois,1,nn);
            dP2D_p2D = dP2D_EFp2D./dp2D_EFp2Dvec;
            dP2D_p2D(isnan(dP2D_p2D)) = zero;
            MM = - mpqelNorm2D.*[Se1.*P2D(in1) Se2.*P2D(in2) Se3.*P2D(in3)];
            tvet = tvet + sparse(ijr,1,MM(mask_ijr),neq,1,'mp');
            MM = - mpqelNorm2D.*[Se1.*dP2D_p2D(in1) Se2.*dP2D_p2D(in2) Se3.*dP2D_p2D(in3)];
            Jmat0 = Jmat0 + sparse(ijr,WW+IIQWP(ijr),MM(mask_ijr),neq,neq,'mp');
            %
            % 2-D electron continuity equation ----------------------------------------
            if mode.idiffusioneQW==2 | mode.idiffusioneQW==3
                if(mode.idiffusioneQW==3)
                    EwC=Ec0(1,:)-DeltaEc;
                    nB2D=mean(Nc2D,1).*exp((EFn2D-EwC)./Vt2D);
                    % nB2D=sum(Nc2D,1).*exp((EFn2D-EwC)./Vt2D);
                    % nB2D=Nc2D(1).*exp((EFn2D-EwC)./Vt2D);
                    gamman2D=n2D./nB2D;
                    lnGn2D=log(gamman2D);
                    dlnGn2D_n2D=((one-n2D./(Vt2D.*dn2D_EFn2D))./nB2D)./gamman2D;
                    delta12 = (phi(inQW1)+Vt2D(iiQW1).*lnGn2D(iiQW1))./Vt2D(iiQW1) - (phi(inQW2)+Vt2D(iiQW2).*lnGn2D(iiQW2))./Vt2D(iiQW2);
                elseif(mode.idiffusioneQW==2)
                    delta12 = phi(inQW1)./Vt2D(iiQW1) - phi(inQW2)./Vt2D(iiQW2);
                end
                n2D1=n2D(iiQW1); n2D2=n2D(iiQW2);
                B12=mpbern(delta12); dB12=mpdbern(delta12)./Vt2D(iiQW1);
                B21=mpbern(-delta12); dB21=mpdbern(-delta12)./Vt2D(iiQW2);
            elseif mode.idiffusioneQW==1
                B12 = one; B21 = one; % diffusion-only equations
            end
            %
            if mode.idiffusioneQW~=0
                M11 =   Gcyl.*Dn2D./LeQW.*B12;
                M12 = - Gcyl.*Dn2D./LeQW.*B21;
                M21 = - Gcyl.*Dn2D./LeQW.*B12;
                M22 =   Gcyl.*Dn2D./LeQW.*B21;
                MM  =   qel.*[M11 M12 M21 M22];
                Kmat0 = Kmat0 + sparse(VV+iivQW(mask_iivQW),VV+jjvQW(mask_iivQW),MM(mask_iivQW),neq,neq,'mp');
                MM  =   qel.*[M11 M12];
                ImatQW = ImatQW + sparse(VV+iivQWcur(mask_iivQWcur),VV+jjvQWcur(mask_iivQWcur),-MM(mask_iivQWcur),neq,neq,'mp');
            end
            %
            if mode.idiffusioneQW==2 | mode.idiffusioneQW==3
                % Derivative w.r.t. potential
                M11 =   Gcyl.*Dn2D./LeQW.*(dB12.*n2D1+dB21.*n2D2);
                M12 = - Gcyl.*Dn2D./LeQW.*(dB21.*n2D2+dB12.*n2D1);
                M21 = - Gcyl.*Dn2D./LeQW.*(dB21.*n2D2+dB12.*n2D1);
                M22 =   Gcyl.*Dn2D./LeQW.*(dB21.*n2D2+dB12.*n2D1);
                MM  =   qel.*[M11 M12 M21 M22];
                Jmat0=Jmat0+sparse(VV+iivQW(mask_iivQW),jnvQW(mask_iivQW),MM(mask_iivQW),neq,neq,'mp');
                %
                % Derivative Fermi factor
                if(mode.idiffusioneQW==3)
                    MM  =   qel.*[Vt2D(iiQW1).*M11.*dlnGn2D_n2D(iiQW1) ...
                        Vt2D(iiQW2).*M12.*dlnGn2D_n2D(iiQW2) ...
                        Vt2D(iiQW1).*M21.*dlnGn2D_n2D(iiQW1) ...
                        Vt2D(iiQW2).*M22.*dlnGn2D_n2D(iiQW2)];
                    Jmat0=Jmat0+sparse(VV+iivQW(mask_iivQW),VV+jjvQW(mask_iivQW),MM(mask_iivQW),neq,neq,'mp');
                end
            end
            
            %
            % 2-D hole continuity equation --------------------------------------------
            if mode.idiffusioneQW==2 | mode.idiffusioneQW==3
                if(mode.idiffusioneQW==3)
                    EwV=Ev0(1,:)+DeltaEv;
                    pB2D=mean(Nv2D,1).*exp((EwV-EFp2D)./Vt2D); % arbitrary :-(
                    % pB2D=sum(Nv2D,1).*exp((EwV-EFp2D)./Vt2D); % arbitrary :-(
                    % pB2D=Nv2D(1).*exp((EwV-EFp2D)./Vt2D); % arbitrary :-(
                    gammap2D=p2D./pB2D;
                    lnGp2D=log(gammap2D);
                    dlnGp2D_p2D=((one+p2D./(Vt2D.*dp2D_EFp2D))./pB2D)./gammap2D;
                    delta12 = (phi(inQW1)-Vt2D(iiQW1).*lnGp2D(iiQW1))./Vt2D(iiQW1) - (phi(inQW2)-Vt2D(iiQW2).*lnGp2D(iiQW2))./Vt2D(iiQW2);
                    %
                elseif(mode.idiffusioneQW==2)
                    
                    delta12 = phi(inQW1)./Vt2D(iiQW1) - phi(inQW2)./Vt2D(iiQW2);
                end
                p2D1=p2D(iiQW1); p2D2=p2D(iiQW2);
                B12=mpbern(delta12); dB12=mpdbern(delta12)./Vt2D(iiQW2);
                B21=mpbern(-delta12); dB21=mpdbern(-delta12)./Vt2D(iiQW1);
            elseif mode.idiffusioneQW==1
                B12 = one; B21 = one;  % diffusion-only equations
            end
            %
            if mode.idiffusioneQW~=0
                M11 =   Gcyl.*Dp2D./LeQW.*B21;
                M12 = - Gcyl.*Dp2D./LeQW.*B12;
                M21 = - Gcyl.*Dp2D./LeQW.*B21;
                M22 =   Gcyl.*Dp2D./LeQW.*B12;
                MM  =   qel.*[M11 M12 M21 M22];
                Kmat0 = Kmat0 + sparse(WW+iivQW(mask_iivQW),WW+jjvQW(mask_iivQW),MM(mask_iivQW),neq,neq,'mp');
                MM  =   qel.*[M11 M12];
                ImatQW = ImatQW + sparse(WW+iivQWcur(mask_iivQWcur),WW+jjvQWcur(mask_iivQWcur),MM(mask_iivQWcur),neq,neq,'mp');
            end
            %
            % Derivative w.r.t. potential
            if mode.idiffusioneQW==2 | mode.idiffusioneQW==3
                M11 =   Gcyl.*Dp2D./LeQW.*(-dB21.*p2D1-dB12.*p2D2);
                M12 = - Gcyl.*Dp2D./LeQW.*(-dB12.*p2D2-dB21.*p2D1);
                M21 = - Gcyl.*Dp2D./LeQW.*(-dB21.*p2D1-dB12.*p2D2);
                M22 =   Gcyl.*Dp2D./LeQW.*(-dB12.*p2D2-dB21.*p2D1);
                MM  =   qel.*[M11 M12 M21 M22];
                Jmat0 = Jmat0+sparse(WW+iivQW(mask_iivQW),jnvQW(mask_iivQW),MM(mask_iivQW),neq,neq,'mp');
                %
                % Derivative Fermi factor
                if(mode.idiffusioneQW==3)
                    MM  = - qel.*[Vt2D(iiQW1).*M11.*dlnGp2D_p2D(iiQW1) ...
                        Vt2D(iiQW2).*M12.*dlnGp2D_p2D(iiQW2)
                        Vt2D(iiQW1).*M21.*dlnGp2D_p2D(iiQW1) ...
                        Vt2D(iiQW2).*M22.*dlnGp2D_p2D(iiQW2)];
                    Jmat0 = Jmat0+sparse(WW+iivQW(mask_iivQW),WW+jjvQW(mask_iivQW),MM(mask_iivQW),neq,neq,'mp');
                end
            end
            %
            % 2D capture terms
            % 2D electrons capture terms
            %        tauscatnT=tauscatn.*exp(-(TQW-mode.T0)/mode.T_tauscat);
            %        tauscatpT=tauscatp.*exp(-(TQW-mode.T0)/mode.T_tauscat);
            tauscatnT=tauscatn.*ones(1,nnQW,'mp');
            tauscatpT=tauscatp.*ones(1,nnQW,'mp');
            
            if(mode.taucarrier~=1) % no Vallone
                if isfield(mode,'T_tauscat')
                    TQW=mesh.T(inQW);
                    T_tauscat=mode.T_tauscat;
                    Tcap_EXP=mode.Tcap_EXP;
                    %                'qui ED', keyboard
                    if T_tauscat~=0
                        %                    TQWm=mean(TQW);
                        %                    tauscatnT=tauscatn.*exp((TQW-mode.T0)/mode.T_tauscat);
                        %                    tauscatpT=tauscatp.*exp((TQW-mode.T0)/mode.T_tauscat);
                        
%                        if Tcap_EXP==Inf
%                            
%                            Fat_cap=exp((TQW-mode.T0)/mode.T_tauscat);
%                        elseif Tcap_EXP==-Inf
%                            
%                            DeT=mode.T_tauscat;
%                            tauRat=mode.tauRat;
%                            
%                        else
                            %                     Fat_cap=(1+(TQW-mode.T0)/mode.T_tauscat).^Tcap_EXP;
                            Fat_cap=exp((TQW-mp(mode.T300))/mp(mode.T_tauscat)).^Tcap_EXP;
%                        end
                        
                        if mode.leakage>0
                            Tstart=mode.Tstart;
                            taubottom=mode.tL_high;
                            DeT=mode.leakage;
                            tautop=mode.tL_low;
                            Fat_capleak=f_RaisedCosine(TQW,Tstart,DeT,taubottom,tautop);
                        end
                        Fcap=one;
                        %                    tauscatnT=tauscatn.*Fcap;
                        %                    tauscatpT=tauscatp.*Fcap;
                        tauscatnT=tauscatn.*Fat_cap;
                        tauscatpT=tauscatp.*Fat_cap;
                        %
                        % Se decidiamo di usare il modello di Vallone
                        
                        %                    tauscatnT=tauscatn.*ones(size(Fat_cap));
                        %                    tauscatpT=tauscatp.*ones(size(Fat_cap));
                    end
                    tauRat=mode.tauRat;
                    tauescn=tauscatn./Fat_cap*tauRat.*ones(1,nnQW,'mp');
                    tauescp=tauscatp./Fat_cap*tauRat.*ones(1,nnQW,'mp');
                    mode.FatCap=Fat_cap;
                    
                end
                
            else % Vallone
                if(isfield(mode,'Fat_cap_e'))
                    tauscatnT=tauscatn.*mode.Fat_cap_e;
                    tauscatpT=tauscatp.*mode.Fat_cap_h;
                end
                tauRat=mode.tauRat;
                tauescn=tauscatn*tauRat.*ones(1,nnQW,'mp');
                tauescp=tauscatp*tauRat.*ones(1,nnQW,'mp');
                mode.FatCap=mode.Fat_cap_e;
                
            end            
            
            tauleakn=tauescn;
            tauleakp=tauescp;
            %                    tauescn=tauscatn.*tauRat.*ones(size(Fat_cap));
            %                    tauescp=tauscatp.*tauRat.*ones(size(Fat_cap));
            %                    tauescp=tauscatn.*Fat_cap
            
            dN3D_Efn3D =    dnF(inQW)./Vt2D;
            
            if(iTappo==1) % esattamente il termine degli svizzeri
                FCapn=exp((EFn2D-EFn3D(1,:))./Vt2D);
                Ccapn = - (one-FCapn).*(one-n2D./N2)./tauscatnT.*nF(inQW);
                dCcapn_N3D =  - (one-FCapn).*(one-n2D./N2)./tauscatnT + ...
                    - FCapn.*(one-n2D./N2).*nF(inQW)./tauscatnT./dN3D_Efn3D./Vt2D; %th der. inversa
                dCcapn_n2D =  + FCapn.*(one-n2D./N2).*nF(inQW)./tauscatnT./dn2D_EFn2D./Vt2D + ...
                    + (one-FCapn)./N2./tauscatnT.*nF(inQW);
            elseif(iTappo==0) % termine degli svizzeri senza il tappo
                FCapn=exp((EFn2D-EFn3D(1,:))./Vt2D);
                Ccapn = - (one-FCapn)./tauscatnT.*nF(inQW);
                dCcapn_N3D =  - (one-FCapn)./tauscatnT + ...
                    - FCapn.*nF(inQW)./tauscatnT./dN3D_Efn3D./Vt2D; %th der. inversa
                dCcapn_n2D =  + FCapn.*nF(inQW)./tauscatnT./dn2D_EFn2D./Vt2D;
            elseif(iTappo==2) % termine Debernardi
                Ccapn = - ((nF(inQW)-mode.n3Di(inQW))./tauscatnT - (n2D-n2Di)./(WQW*tauescn));
                dCcapn_N3D =  - one./tauscatnT;
                dCcapn_n2D =  + one./(WQW*tauescn);
            elseif(iTappo==4) % termine degli svizzeri con leakage
                FCapn=exp((EFn2D-EFn3D(1,:))./Vt2D);
                Ccapn = - (one-FCapn).*(one-n2D./N2)./tauscatnT.*nF(inQW) + (n2D-n2Di)./(WQW*tauleakn);
                dCcapn_N3D =  - (one-FCapn).*(one-n2D./N2)./tauscatnT + ...
                    - FCapn.*(one-n2D./N2).*nF(inQW)./tauscatnT./dN3D_Efn3D./Vt2D; %th der. inversa
                dCcapn_n2D =  + FCapn.*(one-n2D./N2).*nF(inQW)./tauscatnT./dn2D_EFn2D./Vt2D + ...
                    + (one-FCapn)./N2./tauscatnT.*nF(inQW) + one./(WQW*tauleakn);
            elseif(iTappo==3) % svizzeri senza tappo e con leakage
                FCapn=exp((EFn2D-EFn3D(1,:))./Vt2D);
                Ccapn = - (one-FCapn)./tauscatnT.*nF(inQW) + (n2D-n2Di)./(WQW*tauleakn);
                dCcapn_N3D =  - (one-FCapn)./tauscatnT + ...
                    - FCapn.*nF(inQW)./tauscatnT./dN3D_Efn3D./Vt2D; %th der. inversa
                dCcapn_n2D =  + FCapn.*nF(inQW)./tauscatnT./dn2D_EFn2D./Vt2D + one./(WQW*tauleakn);
            elseif(iTappo==6) % cattura-fuga alla Cappelluti/Gioannini maniera
                taucapn = tauscatnT;
                tauescn = taucapn.*mode.rn{indQW};
                %
                Ccapn = - (one-n2D./N2).*nF(inQW)./taucapn + n2D./(WQW*tauescn);
                dCcapn_N3D =  - (one-n2D./N2)./taucapn;
                dCcapn_n2D =  + one./N2./taucapn.*nF(inQW) + one./(WQW*tauescn);
            elseif(iTappo==7 | iTappo==8) % tappo esponenziale
                FCapn=exp((EFn2D-EFn3D(1,:))./Vt2D);
                eD=exp(-n2D./N2);
                Ccapn = - (one-FCapn).*eD./tauscatnT.*nF(inQW);
                dCcapn_N3D =  - (one-FCapn).*eD./tauscatnT + ...
                    - FCapn.*eD.*nF(inQW)./tauscatnT./dN3D_Efn3D./Vt2D; %th der. inversa
                dCcapn_n2D =  + FCapn.*eD.*nF(inQW)./tauscatnT./dn2D_EFn2D./Vt2D + ...
                    + (one-FCapn).*eD./N2./tauscatnT.*nF(inQW);
            end
            
            CcapnQW=Ccapn*CarrierNorm/CarrierNorm2D;
            dCcapn_N3DQW=dCcapn_N3D*CarrierNorm/CarrierNorm2D;
            dCcapn_n2DQW=dCcapn_n2D*CarrierNorm/CarrierNorm2D;
            
            CcapnM1=CcapnQW;
            CcapnM2=CcapnQW;
            
    %        Xlimi=mesh.xgrid2(3)*1.5;
      %      fii=find(mesh.xgrid<=Xlimi);
       %     fie=find(mesh.xgrid>Xlimi);
            fii=1:length(mesh.xgrid);
            fie=1:length(mesh.xgrid);
            CcapnM1(fie)=zero;
            CcapnM2(fii)=zero;
            MM1 = qel.*WQW.*[Lp1.*CcapnM1(iiQW1) Lp2.*CcapnM1(iiQW2)];
            MM2 = qel.*WQW.*[Lp1.*CcapnM2(iiQW1) Lp2.*CcapnM2(iiQW2)];
            
            %             MM = qel.*WQW.*[Lp1.*Ccapn(iiQW1) Lp2.*Ccapn(iiQW2)];
            IntCcapn=IntCcapn+[sum(MM1) sum(MM2)];  % C/s
            
            mode.IntCcapn=IntCcapn;
            
            %             xdx=[zero Lp1];
            %  figure, plot(1e4*mesh.xgrid,[JpQW_x+JnQW_x; JpQW_y+JnQW_y])
            %             Jx=JpQW_x+JnQW_x;
            %             Jy=JpQW_y+JnQW_y;
            %
            %             JnQW_x=mode.Jn_x(inQW);
            %             JnQW_y=mode.Jn_y(inQW);
            %             JpQW_x=mode.Jp_x(inQW);
            %             JpQW_y=mode.Jp_y(inQW);
            
            %
            % assembling in 2D electrons equations
            MM = qel.*WQW.*[Lp1.*CcapnQW(iiQW1) Lp2.*CcapnQW(iiQW2)];
            tvet = tvet + sparse(VV+ijrQW,1,MM,neq,1,'mp');
            MM = qel.*WQW.*[Lp1.*dCcapn_N3DQW(iiQW1) Lp2.*dCcapn_N3DQW(iiQW2)];
            Jmat0 = Jmat0 + sparse(VV+ijrQW,nn+inQW(ijrQW),MM,neq,neq,'mp');
            MM = qel.*WQW.*[Lp1.*dCcapn_n2DQW(iiQW1) Lp2.*dCcapn_n2DQW(iiQW2)];
            Jmat0 = Jmat0 + sparse(VV+ijrQW,VV+ijrQW,MM,neq,neq,'mp');
            %
            % assembling in 3D electrons equations
            Ccapn3D = zero*ones(1,mesh.nn);
            dCcapn3D_N3D = zero*ones(1,mesh.nn);
            dCcapn3D_n2D = zero*ones(1,mesh.nn);
            %
            Ccapn3D(IQW)=Ccapn(IIQWP(IQW));
            dCcapn3D_N3D(IQW)=dCcapn_N3D(IIQWP(IQW));
            dCcapn3D_n2D(IQW)=dCcapn_n2D(IIQWP(IQW));
            %
            MM = - qel.*[SeQW1.*Ccapn3D(in1) SeQW2.*Ccapn3D(in2) SeQW3.*Ccapn3D(in3)];
            tvet = tvet + sparse(nn+ijr,1,MM(mask_ijr),neq,1,'mp');
            MM = - qel.*[SeQW1.*dCcapn3D_N3D(in1) SeQW2.*dCcapn3D_N3D(in2) SeQW3.*dCcapn3D_N3D(in3)];
            Jmat0 = Jmat0 + sparse(nn+ijr,nn+INQWP(ijr),MM(mask_ijr),neq,neq,'mp');
            MM = - qel.*[SeQW1.*dCcapn3D_n2D(in1) SeQW2.*dCcapn3D_n2D(in2) SeQW3.*dCcapn3D_n2D(in3)];
            Jmat0 = Jmat0 + sparse(nn+ijr,VV+IIQWP(ijr),MM(mask_ijr),neq,neq,'mp');
            
            
            %             MM = [SeQW1(in1) SeQW2(in2) SeQW3(in3)];
            %             MM = [SeQW1 SeQW2 SeQW3];
            %             figure(34891),hold on
            %             plot(node(1,:),node(2,:),'b*')
            %             plot(node(1,IQW),node(2,IQW),'ro')
            RQW=zero*ones(1,nn);
            RQW(IQW)=R(IQW);
            %            MM = qel.*[SeQW1.*(Ccapn3D(in1)+RQW(in1)) SeQW2.*(Ccapn3D(in2)+RQW(in2)) SeQW3.*(Ccapn3D(in3)+RQW(in3))];
            %             MM = qel.*[SeQW1.*(Ccapn3D(in1)) SeQW2.*(Ccapn3D(in2)) SeQW3.*(Ccapn3D(in3))];
            %             Int_Anal=pi*(mesh.xgrid(end)-mesh.xgrid(1))^2*WQW*indQW
            %            IntRecQW3Dn=IntRecQW3Dn+sum(sparse(ijr,1,MM(mask_ijr),nn,1));
            
            
            VaIn=RQW;
            MM = qel.*[SeQW1.*VaIn(in1) SeQW2.*VaIn(in2) SeQW3.*VaIn(in3)];
            IntRecQW3Dn=IntRecQW3Dn+sum(sparse(ijr,1,MM(mask_ijr),nn,1,'mp'),1);
            
            %             ErrRelInt=abs(Int_Anal-IntegralQW)./abs(Int_Anal)
            
            mode.IntRecN=double(IntRecQW3Dn);
            
            dP3D_Efp3D =  - dpF(inQW)./Vt2D;
            
            if(iTappo==1) % esattamente il termine degli svizzeri
                FCapp=exp(-(EFp2D-EFp3D(1,:))./Vt2D);
                Ccapp = - (one-FCapp).*(one-p2D./P2)./tauscatpT.*pF(inQW);
                dCcapp_P3D =  - (one-FCapp).*(one-p2D./P2)./tauscatpT + ...
                    + FCapp.*(one-p2D./P2).*pF(inQW)./tauscatpT./dP3D_Efp3D./Vt2D; %th der. inversa
                dCcapp_p2D =  - FCapp.*(one-p2D./P2).*pF(inQW)./tauscatpT./dp2D_EFp2D./Vt2D + ...
                    + (one-FCapp)./P2./tauscatpT.*pF(inQW);
            elseif(iTappo==0) % termine degli svizzeri senza il tappo
                FCapp=exp(-(EFp2D-EFp3D(1,:))./Vt2D);
                Ccapp = - (one-FCapp)./tauscatpT.*pF(inQW);
                dCcapp_P3D =  - (one-FCapp)./tauscatpT + ...
                    + FCapp.*pF(inQW)./tauscatpT./dP3D_Efp3D./Vt2D; %th der. inversa
                dCcapp_p2D =  - FCapp.*pF(inQW)./tauscatpT./dp2D_EFp2D./Vt2D;
            elseif(iTappo==2) % termine Debernardi
                Ccapp = - ((pF(inQW)-mode.p3Di(inQW))./tauscatpT - (p2D-p2Di)./(WQW*tauescp));
                dCcapp_P3D =  - one./tauscatpT;
                dCcapp_p2D =  + one./(WQW*tauescp);
            elseif(iTappo==4) % termine degli svizzeri con leakage
                FCapp=exp(-(EFp2D-EFp3D(1,:))./Vt2D);
                Ccapp = - (one-FCapp).*(one-p2D./P2)./tauscatpT.*pF(inQW) + (p2D-p2Di)./(WQW*tauleakp);
                dCcapp_P3D =  - (one-FCapp).*(one-p2D./P2)./tauscatpT + ...
                    + FCapp.*(one-p2D./P2).*pF(inQW)./tauscatpT./dP3D_Efp3D./Vt2D; %th der. inversa
                dCcapp_p2D =  - FCapp.*(one-p2D./P2).*pF(inQW)./tauscatpT./dp2D_EFp2D./Vt2D + ...
                    + (one-FCapp)./P2./tauscatpT.*pF(inQW) + one./(WQW*tauleakp);
            elseif(iTappo==3) % svizzeri senza tappo e con leakage
                FCapp=exp(-(EFp2D-EFp3D(1,:))./Vt2D);
                Ccapp = - (one-FCapp)./tauscatpT.*pF(inQW) + (p2D-p2Di)./(WQW*tauleakp);
                dCcapp_P3D =  - (one-FCapp)./tauscatpT + ...
                    + FCapp.*pF(inQW)./tauscatpT./dP3D_Efp3D./Vt2D; %th der. inversa
                dCcapp_p2D =  - FCapp.*pF(inQW)./tauscatpT./dp2D_EFp2D./Vt2D + one./(WQW*tauleakp);
            elseif(iTappo==6) % cattura-fuga alla Cappelluti/Gioannini maniera
                taucapp = tauscatpT;
                tauescp = taucapp.*mode.rp{indQW};
                %
                Ccapp = - (one-p2D./P2).*pF(inQW)./taucapp + p2D./(WQW.*tauescp);
                dCcapp_P3D =  - (one-p2D./P2)./taucapp;
                dCcapp_p2D =  + one./P2./taucapp.*pF(inQW) + one./(WQW*tauescp);
            elseif(iTappo==7 | iTappo==8) % mio
                pD=exp(-p2D./P2);
                FCapp=exp(-(EFp2D-EFp3D(1,:))./Vt2D);
                Ccapp = - (one-FCapp).*pD./tauscatpT.*pF(inQW);
                dCcapp_P3D =  - (one-FCapp).*pD./tauscatpT + ...
                    + FCapp.*pD.*pF(inQW)./tauscatpT./dP3D_Efp3D./Vt2D; %th der. inversa
                dCcapp_p2D =  - FCapp.*pD.*pF(inQW)./tauscatpT./dp2D_EFp2D./Vt2D + ...
                    + (one-FCapp).*pD./P2./tauscatpT.*pF(inQW);                
                
            end
            
            %load stopCap
            %if istopCap==1
            % 'Stop Ccap', keyboard
            %end
            CcappQW=Ccapp*CarrierNorm/CarrierNorm2D;
            dCcapp_P3DQW=dCcapp_P3D*CarrierNorm/CarrierNorm2D;
            dCcapp_p2DQW=dCcapp_p2D*CarrierNorm/CarrierNorm2D;
			
            MM = qel.*WQW.*[Lp1.*CcappQW(iiQW1) Lp2.*CcappQW(iiQW2)];
            %             IntCcapp=IntCcapp+sum(MM);  % C/s
            
            CcappM1=CcappQW;
            CcappM2=CcappQW;
            CcappM1(fie)=zero;
            CcappM2(fii)=zero;
            MM1 = qel.*WQW.*[Lp1.*CcappM1(iiQW1) Lp2.*CcappM1(iiQW2)];
            MM2 = qel.*WQW.*[Lp1.*CcappM2(iiQW1) Lp2.*CcappM2(iiQW2)];
            
            IntCcapp=IntCcapp+[sum(MM1) sum(MM2)];  % C/s
            
            
            
            mode.IntCcapp=double(IntCcapp);
            
            
            %
            % assembling in 2D holes equations
            MM = qel.*WQW.*[Lp1.*CcappQW(iiQW1) Lp2.*CcappQW(iiQW2)];
            tvet = tvet + sparse(WW+ijrQW,1,MM,neq,1,'mp');
            MM = qel.*WQW.*[Lp1.*dCcapp_P3DQW(iiQW1) Lp2.*dCcapp_P3DQW(iiQW2)];
            Jmat0 = Jmat0 + sparse(WW+ijrQW,pp+inQW(ijrQW),MM,neq,neq,'mp');
            MM = qel.*WQW.*[Lp1.*dCcapp_p2DQW(iiQW1) Lp2.*dCcapp_p2DQW(iiQW2)];
            Jmat0 = Jmat0 + sparse(WW+ijrQW,WW+ijrQW,MM,neq,neq,'mp');
            %
            % assembling in 3D holes equations
            Ccapp3D = zero*ones(1,mesh.nn);
            dCcapp3D_N3D = zero*ones(1,mesh.nn);
            dCcapp3D_n2D = zero*ones(1,mesh.nn);
            %
            Ccapp3D(IQW)=Ccapp(IIQWP(IQW));
            dCcapp3D_N3D(IQW)=dCcapp_P3D(IIQWP(IQW));
            dCcapp3D_n2D(IQW)=dCcapp_p2D(IIQWP(IQW));
            %
            MM = - qel.*[SeQW1.*Ccapp3D(in1) SeQW2.*Ccapp3D(in2) SeQW3.*Ccapp3D(in3)];
            tvet = tvet + sparse(pp+ijr,1,MM(mask_ijr),neq,1,'mp');
            MM = - qel.*[SeQW1.*dCcapp3D_N3D(in1) SeQW2.*dCcapp3D_N3D(in2) SeQW3.*dCcapp3D_N3D(in3)];
            Jmat0 = Jmat0 + sparse(pp+ijr,pp+INQWP(ijr),MM(mask_ijr),neq,neq,'mp');
            MM = - qel.*[SeQW1.*dCcapp3D_n2D(in1) SeQW2.*dCcapp3D_n2D(in2) SeQW3.*dCcapp3D_n2D(in3)];
            Jmat0 = Jmat0 + sparse(pp+ijr,WW+IIQWP(ijr),MM(mask_ijr),neq,neq,'mp');
            
            %            VaIn=Ccapp3D+RQW;
            %            MM = qel.*[SeQW1.*(Ccapp3D(in1)+RQW(in1)) SeQW2.*(Ccapp3D(in2)+RQW(in2)) SeQW3.*(Ccapp3D(in3)+RQW(in3))];
            %            MM = qel.*[SeQW1.*VaIn(in1) SeQW2.*VaIn(in2) SeQW3.*VaIn(in3)];
            
            VaIn=RQW;
            MM = qel.*[SeQW1.*VaIn(in1) SeQW2.*VaIn(in2) SeQW3.*VaIn(in3)];
            
            IntRecQW3Dp=IntRecQW3Dp+sum(sparse(ijr,1,MM(mask_ijr),nn,1,'mp'),1);
            
            mode.IntRecP=double(IntRecQW3Dp);
            
            %
            % Recombination of bound carriers
            R2D = zero; dR2D_n2D=zero; dR2D_p2D=zero;
            np2D=n2D.*p2D-n2Di.*p2Di;
            %
            % Shockley-Read-Hall recombination term
            taunSRH=mp(mesh.taunQW(inQW));
            taupSRH=mp(mesh.taupQW(inQW));
            denSRH2D = taupSRH.*(n2D+n2Di)+taunSRH.*(p2D+p2Di);
            RSRH2D = np2D./denSRH2D;
            dRSRH_n2D = (p2D.*denSRH2D - taupSRH.*np2D)./(denSRH2D.^two);
            dRSRH_p2D = (n2D.*denSRH2D - taunSRH.*np2D)./(denSRH2D.^two);
            R2D = R2D + RSRH2D;
            dR2D_n2D = dR2D_n2D + dRSRH_n2D;
            dR2D_p2D = dR2D_p2D + dRSRH_p2D;
            
            VI=RSRH2D;
            MM = qel.*[Lp1.*VI(iiQW1) Lp2.*VI(iiQW2)];
            IntSRH=IntSRH+sum(MM);
            
            %
            % Radiative recombination term
            % Semiconductor Bloch Equations (SBE) model
            %            nQW=n2D./WQW;
            %            pQW=p2D./WQW;
            %            np=(nQW+pQW)/2-(n2Di+p2Di)/2/WQW;
            %            IntRecQW3Dn
            
            np_E=n2D-n2Di;
            np_H=p2D-p2Di;
            
            np_E=n2D;
            np_H=p2D;
            
            [Rsp,dRspE,dRspH]=f_InterpRsp_lin(double(np_E),double(np_H),indQW);
            
            Rsp=mp(Rsp)*mp(mode.fat_gain);
            dRspE=mp(dRspE)*mp(mode.fat_gain);
            dRspH=mp(dRspH)*mp(mode.fat_gain);
%             if CarrierNorm>1
%                 Rrad2D_SBE = Rsp; % from 1/(s*cm^3) to 1/(s*cm^2)
%                 Rrad2D = Rsp; % from 1/(s*cm^3) to 1/(s*cm^2)
%                 dRrad2D_n2D = dRspE; % 1/cm and cm simplify
%                 dRrad2D_p2D = dRspH; % 1/cm and cm simplify
%             else
                Rrad2D_SBE = Rsp*WQW; % from 1/(s*cm^3) to 1/(s*cm^2)
                Rrad2D = Rsp*WQW; % from 1/(s*cm^3) to 1/(s*cm^2)
                dRrad2D_n2D = dRspE*WQW; % 1/cm and cm simplify
                dRrad2D_p2D = dRspH*WQW; % 1/cm and cm simplify
%             end
            
            VI=Rrad2D;
            MM = qel.*[Lp1.*VI(iiQW1) Lp2.*VI(iiQW2)];
            IntRad=IntRad+sum(MM);
            
            
            MM = [Lp1.*Rrad2D(iiQW1) Lp2.*Rrad2D(iiQW2)];
            frsp=mp(mode.frsp);
            Psp=Psp+frsp*mp('1000')*sum(sparse(MM),1)*h*(Clight*mp('1e-2')/mean(mp('1e-9')*mp(mode.vlambda))); % milliwatt
            
            % "Brad" model
            %         Brad2D = mesh.brad(inQW)./WQW; % Brad in cm^2/s (/WQW!)
            %         Rrad2D = Brad2D.*np2D; % cm^(-2)/s
            %         dRrad2D_n2D = Brad2D.*p2D;
            %         dRrad2D_p2D = Brad2D.*n2D;
            %
            R2D = R2D + Rrad2D;
            dR2D_n2D = dR2D_n2D + dRrad2D_n2D;
            dR2D_p2D = dR2D_p2D + dRrad2D_p2D;
            %
            % Auger recombination model
%             if CarrierNorm>1
%                 Cnnp2D=mode.FatAuger23D*mesh.Cnnp(inQW);
%                 Cppn2D=mode.FatAuger23D*mesh.Cppn(inQW);
%             else
                Cnnp2D=mp(mode.FatAuger23D)*mp(mesh.Cnnp(inQW))./(WQW.^two);
                Cppn2D=mp(mode.FatAuger23D)*mp(mesh.Cppn(inQW))./(WQW.^two);
%             end
            RAuger_n2D = Cnnp2D.*n2D.*np2D;
            RAuger_p2D = Cppn2D.*p2D.*np2D;
            dRAuger_n2D_n2D = two*Cnnp2D.*np2D; % derivative of RAuger_n2D w.r.t. n2D
            dRAuger_n2D_p2D = Cnnp2D.*n2D.*n2D;  % derivative of RAuger_n2D w.r.t. p2D
            dRAuger_p2D_n2D = Cppn2D.*p2D.*p2D; % derivative of RAuger_p2D w.r.t. n2D
            dRAuger_p2D_p2D = two*Cppn2D.*np2D;  % derivative of RAuger_p2D w.r.t. p2D
            RAuger2D = RAuger_n2D + RAuger_p2D;
            dRAuger_n2D = dRAuger_n2D_n2D + dRAuger_p2D_n2D;
            dRAuger_p2D = dRAuger_n2D_p2D + dRAuger_p2D_p2D;
            
            %             dRAuger_n2D = 2*Cnnp2D.*np2D + Cppn2D.*p2D.*p2D;
            %             dRAuger_p2D = 2*Cppn2D.*np2D + Cnnp2D.*n2D.*n2D;
            
            VI=RAuger2D;
            MM = qel.*[Lp1.*VI(iiQW1) Lp2.*VI(iiQW2)];
            IntAug=IntAug+sum(MM);
            
            R2D = R2D + RAuger2D;
            dR2D_n2D = dR2D_n2D + dRAuger_n2D;
            dR2D_p2D = dR2D_p2D + dRAuger_p2D;
            
            % Auger expulsion model, from 2013Deppner_SPIE
            if(mode.AugerExpulsion==1)
                % Assembling n2D expulsion term (n2D equation only)
                MM = qel.*[Lp1.*RAuger_n2D(iiQW1) Lp2.*RAuger_n2D(iiQW2)];
                tvet = tvet + sparse(VV+ijrQW,1,MM,neq,1,'mp');
                % Assembling n2D expulsion jacobian w.r.t. n2D
                MM = qel.*[Lp1.*dRAuger_n2D_n2D(iiQW1) Lp2.*dRAuger_n2D_n2D(iiQW2)];
                Jmat0 = Jmat0 + sparse(VV+ijrQW,VV+ijrQW,MM,neq,neq,'mp');
                % Assembling n2D expulsion jacobian w.r.t. p2D
                MM = qel.*[Lp1.*dRAuger_n2D_p2D(iiQW1) Lp2.*dRAuger_n2D_p2D(iiQW2)];
                Jmat0 = Jmat0 + sparse(VV+ijrQW,WW+ijrQW,MM,neq,neq,'mp');
                %
                % Assembling p2D expulsion term (p2D equation only)
                MM = qel.*[Lp1.*RAuger_p2D(iiQW1) Lp2.*RAuger_p2D(iiQW2)];
                tvet = tvet + sparse(WW+ijrQW,1,MM,neq,1,'mp');
                % Assembling p2D expulsion jacobian w.r.t. p2D
                MM = qel.*[Lp1.*dRAuger_p2D_p2D(iiQW1) Lp2.*dRAuger_p2D_p2D(iiQW2)];
                Jmat0 = Jmat0 + sparse(WW+ijrQW,WW+ijrQW,MM,neq,neq,'mp');
                % Assembling p2D expulsion jacobian w.r.t. n2D
                MM = qel.*[Lp1.*dRAuger_p2D_n2D(iiQW1) Lp2.*dRAuger_p2D_n2D(iiQW2)];
                Jmat0 = Jmat0 + sparse(WW+ijrQW,VV+ijrQW,MM,neq,neq,'mp');
                %
                % Auger regeneration: Gaussian profile
                Deltaz_E=mp('190e-7');
                Deltaz_H=mp('170e-7');
                
                B=mp(mode.Auger_broad);
                f_reg=mp(mode.Fat_regeneration);
                
                GaussEnvelopE=exp(-two*((mp(mesh.ygrid)-(yQW+Deltaz_E))./(B*mp('1e-7'))).^two);
                GaussEnvelopH=exp(-two*((mp(mesh.ygrid)-(yQW-Deltaz_H))./(B*mp('1e-7'))).^two);
                GaussEnvelop=exp(-two*((mp(mesh.ygrid)-yQW)./(B*mp('1e-7'))).^two);
                GaSum=GaSum+GaussEnvelop;
                %                pausak
                % Norm2=sqrt(pi/2)*B*1e-7;
                Norm2=trapz(mp(mesh.ygrid),GaussEnvelop);
                Norm2E=trapz(mp(mesh.ygrid),GaussEnvelopE);
                Norm2H=trapz(mp(mesh.ygrid),GaussEnvelopH);
                GaussEnvelop=GaussEnvelop./Norm2;
                GaussEnvelopE=GaussEnvelopE/Norm2E;
                GaussEnvelopH=GaussEnvelopH/Norm2H;
                %
                % Assembling n2D regeneration term (n3D equation only)
                GAuger_n2D      = - f_reg*RAuger_n2D;
                GAuger_n3D      = reshape(GaussEnvelopE.'*GAuger_n2D,1,nn);
                MM              = qel.*[Se1.*GAuger_n3D(in1) Se2.*GAuger_n3D(in2) Se3.*GAuger_n3D(in3)];
                tvet            = tvet + sparse(nn+ijs,1,MM(mask_ijs),neq,1,'mp');
                % Assembling n3D regeneration jacobian w.r.t. n2D
                dGAuger_n2D_n2D = - f_reg*dRAuger_n2D_n2D;
                dGAuger_n3D_n2D = reshape(GaussEnvelopE.'*dGAuger_n2D_n2D,1,nn);
                MM              = qel.*[Se1.*dGAuger_n3D_n2D(in1) Se2.*dGAuger_n3D_n2D(in2) Se3.*dGAuger_n3D_n2D(in3)];
                Jmat0           = Jmat0 + sparse(nn+ijs,VV+IIQWP(ijs),MM(mask_ijs),neq,neq,'mp');
                % Assembling n3D regeneration jacobian w.r.t. p2D
                dGAuger_n2D_p2D = - f_reg*dRAuger_n2D_p2D;
                dGAuger_n3D_p2D = reshape(GaussEnvelopE.'*dGAuger_n2D_p2D,1,nn);
                MM              = qel.*[Se1.*dGAuger_n3D_p2D(in1) Se2.*dGAuger_n3D_p2D(in2) Se3.*dGAuger_n3D_p2D(in3)];
                Jmat0           = Jmat0 + sparse(nn+ijs,WW+IIQWP(ijs),MM(mask_ijs),neq,neq,'mp');
                %
                % Assembling p2D regeneration term (p3D equation only)
                GAuger_p2D      = - f_reg*RAuger_p2D;
                GAuger_p3D      = reshape(GaussEnvelopH.'*GAuger_p2D,1,nn);
                MM              = qel.*[Se1.*GAuger_p3D(in1) Se2.*GAuger_p3D(in2) Se3.*GAuger_p3D(in3)];
                tvet            = tvet + sparse(pp+ijs,1,MM(mask_ijs),neq,1,'mp');
                % Assembling p3D regeneration jacobian w.r.t. p2D
                dGAuger_p2D_p2D = - f_reg*dRAuger_p2D_p2D;
                dGAuger_p3D_p2D = reshape(GaussEnvelopH.'*dGAuger_p2D_p2D,1,nn);
                MM              = qel.*[Se1.*dGAuger_p3D_p2D(in1) Se2.*dGAuger_p3D_p2D(in2) Se3.*dGAuger_p3D_p2D(in3)];
                Jmat0           = Jmat0 + sparse(pp+ijs,WW+IIQWP(ijs),MM(mask_ijs),neq,neq,'mp');
                % Assembling p3D regeneration jacobian w.r.t. n2D
                dGAuger_p2D_n2D = - f_reg*dRAuger_p2D_n2D;
                dGAuger_p3D_n2D = reshape(GaussEnvelopH.'*dGAuger_p2D_n2D,1,nn);
                MM              = qel.*[Se1.*dGAuger_p3D_n2D(in1) Se2.*dGAuger_p3D_n2D(in2) Se3.*dGAuger_p3D_n2D(in3)];
                Jmat0           = Jmat0 + sparse(pp+ijs,VV+IIQWP(ijs),MM(mask_ijs),neq,neq,'mp'); 
                
            end
            
            %
            % Additional SRH-like leakage terms
            if(mode.leakage>0)
                % tauleakn
                % tauleakp
                %                tauleakn=1e-8;
                %                tauleakp=1e-8;
                tauleakn=Fat_capleak;
                tauleakp=Fat_capleak;
                %                 RLeakage2D=(n2D-n2Di)./tauleakn + (p2D-p2Di)./tauleakp;
                %                 dRLeakage_n2D = 1./tauleakn;
                %                 dRLeakage_p2D = 1./tauleakp;
                %            np2D=n2D.*p2D-n2Di.*p2Di;
                denLeakage2D = tauleakp.*(n2D+n2Di)+tauleakn.*(p2D+p2Di);
                RLeakage2D = np2D./denLeakage2D;
                VI=RLeakage2D;
                MM = qel.*[Lp1.*VI(iiQW1) Lp2.*VI(iiQW2)];
                IntLea=IntLea+sum(MM);
                
                dRLeakage_n2D = (p2D.*denLeakage2D - tauleakp.*np2D)./(denLeakage2D.^two);
                dRLeakage_p2D = (n2D.*denLeakage2D - tauleakn.*np2D)./(denLeakage2D.^two);
                
                R2D = R2D + RLeakage2D;
                dR2D_n2D = dR2D_n2D + dRLeakage_n2D;
                dR2D_p2D = dR2D_p2D + dRLeakage_p2D;
            else
                RLeakage2D=zero*ones(1,mesh.nnx);
            end
            %
            % 2D electron equation recombination term ---------------------------------
            MM = qel.*[Lp1.*R2D(iiQW1) Lp2.*R2D(iiQW2)];
            tvet = tvet + sparse(VV+ijrQW,1,MM,neq,1,'mp');
            MM = qel.*[Lp1.*dR2D_p2D(iiQW1) Lp2.*dR2D_p2D(iiQW2)];
            Jmat0 = Jmat0 + sparse(VV+ijrQW,WW+ijrQW,MM,neq,neq,'mp');
            MM = qel.*[Lp1.*dR2D_n2D(iiQW1) Lp2.*dR2D_n2D(iiQW2)];
            Jmat0 = Jmat0 + sparse(VV+ijrQW,VV+ijrQW,MM,neq,neq,'mp');
            %
            % 2D hole equation recombination term -------------------------------------
            MM = qel.*[Lp1.*R2D(iiQW1) Lp2.*R2D(iiQW2)];
            tvet = tvet + sparse(WW+ijrQW,1,MM,neq,1,'mp');
            MM = qel.*[Lp1.*dR2D_p2D(iiQW1) Lp2.*dR2D_p2D(iiQW2)];
            Jmat0 = Jmat0 + sparse(WW+ijrQW,WW+ijrQW,MM,neq,neq,'mp');
            MM = qel.*[Lp1.*dR2D_n2D(iiQW1) Lp2.*dR2D_n2D(iiQW2)];
            Jmat0 = Jmat0 + sparse(WW+ijrQW,VV+ijrQW,MM,neq,neq,'mp');
            %
            mode.RSRHQW(IQW)=mode.RSRHQW(IQW)+double(RSRH2D(IIQWP(IQW))./WQW);
            mode.RradQW(IQW)=mode.RradQW(IQW)+double(Rrad2D(IIQWP(IQW))./WQW);
            mode.RradQW_SBE(IQW)=mode.RradQW_SBE(IQW)+double(Rrad2D_SBE(IIQWP(IQW))./WQW);
            mode.RAugerQW(IQW)=mode.RAugerQW(IQW)+double(RAuger2D(IIQWP(IQW))./WQW);
            mode.RLeakageQW(IQW)=mode.RLeakageQW(IQW)+double(RLeakage2D(IIQWP(IQW))./WQW);
            %
            %x=mesh.xgrid;
            %		Cn=sum(reshape(mode.Ccapn3D,mesh.nny,mesh.nnx));
            %		Cp=sum(reshape(mode.Ccapp3D,mesh.nny,mesh.nnx));
            % Icn=qel*twopi*trapz(x,x.*Cn)*mp('1000')*WQW;
            % Icp=qel*twopi*trapz(x,x.*Cp)*mp('1000')*WQW;
            %
            mode.Em{indQW}=double(Em); mode.n2D{indQW}=double(n2D); mode.EFn2D{indQW,mode.ind_v0}=double(EFn2D);
            mode.Ed{indQW}=double(Ed); mode.p2D{indQW}=double(p2D); mode.EFp2D{indQW,mode.ind_v0}=double(EFp2D);
            mode.N2{indQW}=double(N2); mode.P2{indQW}=double(P2);
            mode.Ccapn{indQW,mode.ind_v0}=double(Ccapn); mode.Ccapp{indQW,mode.ind_v0}=double(Ccapp);
            %
            mode.Ccapn3D=double(mode.Ccapn3D)-double(Ccapn3D); mode.Ccapp3D=double(mode.Ccapp3D)-double(Ccapp3D);
            %
            % =========================================================================
            % Assembling photon rate equations: stimulated emission
            % =========================================================================
            vph=Clight./mp(mode.nindexQW); % phase velocity in GaAs quantum well
            TQW=mp(mesh.T(inQW));
            diffMQW=find(abs(diff(cell2mat(mesh.MQWcell)))>2);
            %
            for indMode=1:nmodes
                
                Gamma_z=mp(mode.Gamma_z(indMode,indQW));

                if isempty(diffMQW)
                    % 1 MQW region
                    Campo_attivo=mode.Gamma_z(indMode,:)./mode.Gamma_z(indMode,ceil(end/2));
                else
                    % more than one MQW region
                    Campo_attivo=mode.Gamma_z(indMode,:)./mode.Gamma_z(indMode,ceil(diffMQW(1)/2));
                end
                Fattore_attivo=Campo_attivo(indQW);
                
                [g,dgE,dgH,rsp,drspE,drspH] = f_InterpGain_lin(double(n2D),double(p2D),indQW,indMode);
                                
                %
                if isfield(mode,'fat_gainG')
                 Fat_soloGain=mp(mode.fat_gainG);
                else
                 Fat_soloGain=one;
                end
                
                if isfield(mode,'nlG')
                 nlG=mp(mode.nlG);
                else
                 nlG=one;
                end                
                
                fat_gain=mp(mode.fat_gain);
                g=vph*fat_gain*Fat_soloGain*mp(g).*nlG;

                mode.g{indQW,indMode}=double(g); % saving just for the last mode..
                
                dgE=vph*fat_gain*Fat_soloGain*mp(dgE).*nlG;
                dgH=vph*fat_gain*Fat_soloGain*mp(dgH).*nlG;
                
                rsp=mp(rsp)*fat_gain;
                drspE=mp(drspE)*fat_gain;
                drspH=mp(drspH)*fat_gain;
                %
                %-- Stimulated emission: rate equation
                mode.matgain=double(mode.matgain)+g/vph; % saving material gain for VELM
                %            mode.DeltaN=mode.DeltaN+DeltaN;
                %
                E2 = mp(mode.E2(indMode,1:nnQW)); % electric field intensity (normalized)
                %
                FLos=(mode.FLos);
                Lm = FLos*mp(mode.Lm(indMode)); % losses, 1/s                
                %
                gE = g.*E2; % gain-field product
                MM = [Lp1.*gE(iiQW1) Lp2.*gE(iiQW2)];
                
                
                
                gm = sum(diag(sparse(ijrQW,ijrQW,MM,nnQW,nnQW,'mp')),1);
                dgm_nE = dgE.*E2;
                dgm_nH = dgH.*E2;
                %
                rspE = rsp.*E2;
                MM = [Lp1.*rspE(iiQW1) Lp2.*rspE(iiQW2)];
                rspMod = sum(diag(sparse(ijrQW,ijrQW,MM,nnQW,nnQW,'mp')),1);
                dRsp_nE = drspE.*E2;
                dRsp_nH = drspH.*E2;
                %
                fPES=mp(mode.fPES);
                % Assembling spontaneous emission power
                Sq=Gamma_z.*fPES(indMode).*rspMod;
                mode.Sq(mode.ind_v0)=double(Sq);
                tvet(ss+indMode)=tvet(ss+indMode)+Sq;
                %
                % Assembling derivatives of photon equation w.r.t. carriers
                PJacobE = Gamma_z.*dgm_nE*Pst(indMode) + dRsp_nE.*Gamma_z.*fPES(indMode);
                MM = [Lp1.*PJacobE(iiQW1) Lp2.*PJacobE(iiQW2)];
                Jmat0 = Jmat0+sparse(ss+indMode,VV+ijrQW,MM,neq,neq,'mp');
                PJacobH = Gamma_z.*dgm_nH*Pst(indMode) + dRsp_nH.*Gamma_z.*fPES(indMode);
                MM = [Lp1.*PJacobH(iiQW1) Lp2.*PJacobH(iiQW2)];
                Jmat0 = Jmat0+sparse(ss+indMode,WW+ijrQW,MM,neq,neq,'mp');
                %
                % Assembling derivative w.r.t. Pst of photon rate equation (linear)
                MM = (Gamma_z.*gm-Lm/NQW);    % la somma su tutti i QW deve dare Lm
                MMgain=MM+MMgain;
                %
                Kmat0=Kmat0+sparse(ss+indMode,ss+indMode,MM,neq,neq,'mp');
                %
                % for display only
                mode.gmod(indMode)=double(mode.gmod(indMode)+Gamma_z.*gm./vph);
                mode.lmod(indMode)=double(Lm./vph);
                %
                % Equazione Portatori 2D
                %-- Stimulated emission: recombination term
                Rst = gE.*Pst(indMode).*Fattore_attivo*mp(mode.fPdif(indMode))/CarrierNorm2D*WQW;
                % 2D electron equation recombination term ---------------------------------
                MM = qel.*[Lp1.*Rst(iiQW1) Lp2.*Rst(iiQW2)];
                
                IntStim=IntStim+sum(MM);
                
                tvet = tvet + sparse(VV+ijrQW,1,MM,neq,1,'mp');
                tvet = tvet + sparse(WW+ijrQW,1,MM,neq,1,'mp');
                %
                % Derivative of electrons equation Rst with respect to electrons
                dRst_n2D = Fattore_attivo*dgm_nE.*Pst(indMode).*mode.fPdif(indMode)/CarrierNorm2D*WQW;
                MM = qel.*[Lp1.*dRst_n2D(iiQW1) Lp2.*dRst_n2D(iiQW2)];
                Jmat0 = Jmat0+sparse(VV+ijrQW,VV+ijrQW,MM,neq,neq,'mp');
                %
                % Derivative of holes equation Rst with respect to holes
                dRst_p2D = Fattore_attivo*dgm_nH.*Pst(indMode).*mode.fPdif(indMode)/CarrierNorm2D*WQW;
                MM = qel.*[Lp1.*dRst_p2D(iiQW1) Lp2.*dRst_p2D(iiQW2)];
                Jmat0 = Jmat0+sparse(WW+ijrQW,WW+ijrQW,MM,neq,neq,'mp');
                %
                % Derivative of electrons equation Rst with respect to holes
                MM = qel.*[Lp1.*dRst_p2D(iiQW1) Lp2.*dRst_p2D(iiQW2)];
                Jmat0 = Jmat0+sparse(VV+ijrQW,WW+ijrQW,MM,neq,neq,'mp');
                %
                % Derivative of holes equation with Rst respect to electrons
                MM = qel.*[Lp1.*dRst_n2D(iiQW1) Lp2.*dRst_n2D(iiQW2)];
                Jmat0 = Jmat0+sparse(WW+ijrQW,VV+ijrQW,MM,neq,neq,'mp');
                %
                % Derivatives of electrons and holes Rst equation w.r.t. Pst
                %Rst = gE.*Pst(indMode).*mode.fPdif(indMode)*WQW;
                dRst_Pi = Fattore_attivo*gE.*mp(mode.fPdif(indMode))*WQW/CarrierNorm2D;
                MM = qel.*[Lp1.*dRst_Pi(iiQW1) Lp2.*dRst_Pi(iiQW2)];
                Jmat0 = Jmat0 + sparse(VV+ijrQW,ss+indMode,MM,neq,neq,'mp');
                Jmat0 = Jmat0 + sparse(WW+ijrQW,ss+indMode,MM,neq,neq,'mp');
                
            end   %modes
            
            tvetQW(VV+[1:mesh.nnxQW{indQW}]) = tvet(VV+[1:mesh.nnxQW{indQW}]);
            tvetQW(WW+[1:mesh.nnxQW{indQW}]) = tvet(WW+[1:mesh.nnxQW{indQW}]);
            
        end  % MQW
        
        mode.MMgain(mode.ind_v0)=MMgain;

        curQW = ImatQW*uvet'+tvetQW;
        nnQW=mesh.nnxQW{1};
        iiQW=1:nnQW;
        for indQW = 1:NQW
            WQW = double(mesh.vWMQW{indQW}); %-- cm
            VV = vv+(indQW-1)*length(mesh.xgrid);
            WW = ww+(indQW-1)*length(mesh.xgrid);
            mode.JnQW{indQW} = double(curQW(VV+iiQW)/WQW);
            mode.JpQW{indQW} = double(curQW(WW+iiQW)/WQW);
            
        end
%         figure,plot(mesh.xgrid,mode.JnQW{1},mesh.xgrid,mode.Jn_x(mesh.inMQW{indQW}),'--'),title('Electron currents'),legend('2D','3D')
%         figure,plot(mesh.xgrid,mode.JpQW{1},mesh.xgrid,mode.Jp_x(mesh.inMQW{indQW}),'--'),title('Hole currents'),legend('2D','3D')

        Scheck=mode.gmod'./mode.lmod; % stimulated emission check
        %
        % Debug prints
        for indMode=1:nmodes
            fprintf('Pst = %.5e\n',Pst(indMode).')
            fprintf('Gamma_z*Gm/Lm = %.8e\n',full(Scheck(indMode)))
            disp([' '])
        end
%         load stopCap
%         if istopCap==1
%             'Stop Ccap', keyboard
%         end
        mode.N2D=double(N2Dtot); mode.P2D=double(P2Dtot);
        %
        % saving variables
        mode.Pst=double(Pst); mode.Scheck=double(Scheck);
        %mode.nQW{mode.ind_v0}=nQW; mode.pQW{mode.ind_v0}=pQW;
        mode.matgain=mode.matgain/(NQW*nmodes); mode.Psp=double(Psp);
    
        % assembling additional equations for passivation
        if(mesh.nnx-mesh.nnxQW{1}>zero)
            iiP = mesh.nnxQW{1}+1:mesh.nnx;
            for indQW=1:mesh.NMQW
                VV = vv+(indQW-1)*length(mesh.xgrid);
                WW = ww+(indQW-1)*length(mesh.xgrid);
                Kmat0 = Kmat0 + sparse(VV+iiP,VV+iiP,1,neq,neq,'mp');
                Kmat0 = Kmat0 + sparse(WW+iiP,WW+iiP,1,neq,neq,'mp');
            end
        end
        
    end
    %mode.DeltaN=mode.DeltaN/(NQW*nmodes);
    %
    % =============================================================================================100
    % Apply boundary conditions, Poisson equation: phi-V-phi_r=0
    % =============================================================================================100
    ii=find(not(maskr));
    Kmat0=Kmat0+sparse(ii,ii,one,neq,neq,'mp'); % phi
    %
    for ic=1:nm
        ii=find(contact==ic);
        Kmat0=Kmat0+sparse(ii,qq+(ic-1)+1,-one,neq,neq,'mp');
    end % - V
    %
    for ic=1:nc; ii=find(contact==ic); % loop on contacts
        switch geom.contact_type(ic)
            case 1 % ohmic contact
                tvet(ii)=-phi_neutr(ii); % -phi_r
            case 2 % Schottky contact
                tvet(ii)=(mp(geom.workfun(ic))+mp(mesh.phi_rr)); % -phi_r
            case 3 % Schottky contact, vsurf=Inf
                tvet(ii)=(mp(geom.workfun(ic))+mp(mesh.phi_rr)); % -phi_r
            otherwise, error('contact type unknown!')
        end
    end
        % =============================================================================================100
        % Apply boundary conditions, continuity equations
        % =============================================================================================100
        ii=find(not(masks));
        if(mode.nflg), Jmat0=Jmat0+sparse(ii+nn,ii+nn,one,neq,neq,'mp'); end
        if(mode.pflg), Jmat0=Jmat0+sparse(ii+pp,ii+pp,one,neq,neq,'mp'); end
        %
        % =============================================================================================100
        % Compute residual
        % =============================================================================================100
%         Jmat0=double(Jmat0); Kmat0=double(Kmat0); uvet=double(uvet); 
%         tvet=double(tvet); Jmat2=double(Jmat2);
        
        rvet = (Kmat0 + Jmat2) * uvet.' + tvet;
        
        if mode.oflg
            Intot=sum([IntAug IntSRH IntRad IntStim IntLea])
            IntCa=sum(IntCcapp)/WQW
            mode.IntRec=double([IntSRH IntRad IntAug IntStim IntLea]);
        end
        
        if(mode.AugerExpulsion==1)
            mode.GaSum=GaSum;
            %                  figure(489330), plot(mesh.ygrid*1e4,GaSum),
            %                  xlim([354 356])
            %                  keyboard
        end
        