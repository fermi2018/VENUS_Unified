%========================================================================76
% test pin heterojunction
%========================================================================76
clear
clear global
close all
colordef white
dbstop if error
addpath('termico')
addpath('rumenta')
addpath('generageom')
addpath('rumenta/new17Optica')
dbstop if error

irest=2; % restart flag; 0 restarts from 0, 2 restarts from ended simulation and need additional Vadd points


rad_setting='_Feb21';  
Last_Workspace='Mer';  % Deltalam=-3  Rel=3
VELMinput='VelmR';

isperF='old';     % file da caricare per Standard
%isperF='new';

CONTRIBUTI_TERMICI=1;
IPLOT=0;


ivelMstop=0;
vv_double=0;

flagCallVELM=1;
mode1=0;
iVELM=0;  % se 1, calcolo sempre VELM

irel=1; % if 1, relief; if 0, standard VCSEL
nqw=3;

eval(['settings_vari',rad_setting])
 'RIGENERO STRUTTURA'
 structureName
 Sub_Str_from_VELM
%end

geom.QWorientation='X'; % orient. quantum well (parallel to x (rho): VCSEL)
% geom.QWorientation='Y'; % orient. quantum well (parallel to y (z): microROD)






mode.iderGLUT=1;   % = 0 fa vecchio fit
mode.iderGLUTnew=2;   % = 0 fa vecchio fit
eval(['settings_vari',rad_setting])
%carica LUT e la mette nei common
Glut4Dinc(mode)


% mode.symmetry='Cylindrical-X'; % rotation around x axis
%========================================================================76
% solve at equilibrium
% mode.Zmat=8.5; % linear network embedding nonlinear device
mode.Zmat=0; % linear network embedding nonlinear device

% mesh generation
%'me', keyboard
mesh=rectmesh(geom); % create tensorial mesh
%'me', keyboard
mode.T0=300;
T=mode.T0*ones(1,mesh.nn);


mesh.DeltaT=T-mode.T0;

mode.DeltaTmax=0;

fil_str=[structureName,'.str'];


% Temporaneo: da sostituire con ci? che sar? preso dal codice Density
% Matrix

%dd=[[0:0.4:1.2], [1.3:0.1:1.6], [1.65:0.05:1.75], [1.76:0.02:2.4]]; % col bulk, per res
dd=[[0:0.4:1.2], [1.3:0.1:1.6], [1.65:0.02:2.8]]; % col bulk, per res
dd=[[0:0.4:1.2], [1.3:0.1:1.6], [1.65:0.02:1.96]]; % col bulk, per res
%dd=[0:0.4:1.2]; % col bulk, per res

%'v0dd', keyboard
%dd=[[0:0.1:1.3], [1.35:0.05:1.8]];
%dd=sort([dd dd]);
mode.v0_dd=dd;

eval(['settings_vari',rad_setting])

mode.vlambda=860*ones(NUMERO_MODI,1); % inizializzazione
mode.matgain=0;

Dn=mode.matgain;
mode.DeltaN=0;

mesh=loadmacro(geom,mesh,mode,T);

inoplot=1;
if inoplot==0
node='off'; triangle='off'; color='on'; vpath='off'; arrows='off';
sd=1:geom.nd; scale=[]; cmap=[]; grido='on'; cbar='off';
figure
plot_tri_mesh(geom,mesh,[],sd,color,triangle,node,grido,cbar,cmap,scale,vpath,arrows)
drawnow

figure(1001)
set(gcf,'Position',[144 514 1313 452])
subplot(1,2,1)
axis on
grid on
hold on
box on
plot(1e7*mesh.node(2,:),mesh.xmol,'b','LineWidth',2)
xlabel('z (nm)')
ylabel('Molar fraction')
subplot(1,2,2)
axis on
grid on
hold on
box on
plot(1e7*mesh.node(2,:),mesh.dop_d,'b','LineWidth',2)
plot(1e7*mesh.node(2,:),mesh.dop_a,'r','LineWidth',2)
xlabel('z (nm)')
ylabel('Doping profile cm^{-3}')
legend('N_D','N_A')
set(gca,'yscale','log')
drawnow
%
disp([' ']);
disp(['Verifying molar fraction and doping properties']);
disp([' ']);
% pausak
%
% return
end

mode=p_solve1(geom,mesh,mode); % first run: thermodynamic equilibrium
mode.n3Di=mode.elec;
mode.p3Di=mode.hole;

% nmodes=length(mode.Lm); % number of optical modes (VELM)
% mode.nmodes=4;
mode.nmodes=NUMERO_MODI;
mode.lambda=mode.vlambda;

% To initialize matgain, we use carriers in the central well
%indCenter=floor((mesh.NMQW+1)/2);
%nQW=median(mode.n2D{indCenter});
%pQW=median(mode.p2D{indCenter});
%TQW=300*ones(size(nQW));
%[g_eq] = f_InterpGain4D(nQW,pQW,mode.Deltalam+mode.vlambda(1)*ones(size(nQW)),TQW);
g_eq=0;
mode.matgain=g_eq;
Tprec=0;
PTherm=0;

% initializations
mode.iconv=0;
semiconductor=find(geom.semiconductor);
it=not(ismember(mesh.triangle(4,:),semiconductor));

% compute electric field (in time domain), V/um
efieldx_t=-mode.phi*mesh.gradx*1e-4; efieldx_t(it)=0;
efieldy_t=-mode.phi*mesh.grady*1e-4; efieldy_t(it)=0;
mode.efield_x=reshape(pdeprtni(mesh.node,mesh.triangle(1:4,:),efieldx_t),1,mesh.nn);
mode.efield_y=reshape(pdeprtni(mesh.node,mesh.triangle(1:4,:),efieldy_t),1,mesh.nn);
mode.Deltan=0; % Drude Delta refractive index (real part)
mode.alpha=0;  % Drude Delta refractive index (imag part)



% mode_equil=mode;
% mesh_equil=mesh;
% save('ThEq.mat','mode_equil','mesh_equil');

% return


% =============================================================================================100
% solution vector uvet: {phi,elec,hole,elec_t,v_dd,i_dd}
% vector indexes:
%
%    in    -> pointer to Poisson equation
%    in+nn -> pointer to electron continuity equation
%    in+pp -> pointer to hole continuity equation
%    in+tt -> pointer to trap equations
%    in+qq -> pointer to circuit equations
%    in+rr -> pointer to current equations
%    iv+vv -> pointer to 2d electron continuity equations
%    iv+ww -> pointer to 2d hole continuity equations
%    in+ss -> pointer to photon rate equations
%
nm=geom.nm;               % number of active contacts
if(isfield(mesh,'nnxQW'))
    nnQW=mesh.nnxQW;         % number of lateral points in the QW
else
    nnQW=0;
end
NQW=mesh.NMQW;                % number of quantum wells
%
nl=mode.ntrap;                % number of trap levels
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
neq=nn+mode.nflg*nn+mode.pflg*nn+mode.tflg*nl*nn+2*nm+mode.oflg*not(mode.firstrun)*2*NQW*nnQW+mode.oflg*mode.nmodes;

% We first compute the solution vector with all the information we have at this point;
% notice that here the solution is in frequency domain
uvet=zeros(1,neq);
uvet(1:nn)=mode.phi;
if(mode.nflg), uvet(nn+(1:nn))=mode.elec; end
if(mode.pflg), uvet(pp+(1:nn))=mode.hole; end
if(isfield(mode,'v_dd')), uvet(qq+1)=mode.v_dd; end
if(isfield(mode,'i_dd')), uvet(rr+1)=mode.i_dd; end
if(isfield(mode,'n2D'))
    for indQW=1:NQW
        uvet(vv+(indQW-1)*nnQW+(1:nnQW))=mode.n2D{indQW};
    end
end
if(isfield(mode,'p2D'))
    for indQW=1:NQW
        uvet(ww+(indQW-1)*nnQW+(1:nnQW))=mode.p2D{indQW};
    end
end

if(mode.oflg & isfield(mode,'Pst')),
    uvet(ss+(1:mode.nmodes))=mode.Pst;
end
%

s_LoadConstants

mode_old=mode;
indv=1;

DeltaTempVELMLast=0;
npMaxLast=0;
Pst=0;
indVELM=0;

if irest>0
        save LSW Last_Workspace irest rad_setting
%       !copy Last_Workspace.mat Last_Workspace_Old.mat        
%'qui', keyboard
     if irest==2
        eval(['!copy ', Last_Workspace,'.mat  ',  Last_Workspace,'_OLD.mat']) 
     else   
        eval(['!copy ', Last_Workspace,'_stop.mat  ',  Last_Workspace,'_OLDstop.mat']) 
     end   
        clear
        load LSW
        %        clear global
        close all
     if irest==2
        eval(['load ',Last_Workspace])
     else   
        eval(['load ',Last_Workspace,'_stop'])
     end        

        load LSW        
%        load Last_Workspace
        indv=indv+1;
    eval(['settings_vari',rad_setting])
    if irest==2
            mode.v0_dd=[mode.v0_dd,Vadd];

    end

end
%            'verifica', keyboard
pippoflag=0;
if(pippoflag==1)
    disp('hai messo pippoflag, ne sei proprio convinto?')
    keyboard
end
%'ferna', keyboard

if irel==1
    load('Meas_ReliefVCSEL.mat')
    % load('Meas_Relief.mat')
else
 if strcmp(isperF,'old')
    load('Meas_StandardVCSEL_ebr3.mat')  %bookchapter
  else  
    load('Meas_StandardVCSEL_fbu3.mat')   % altro
  end 
end

Rmeas=diff(Vmeas)./diff(Imeas/1000);
Rmeas=[Rmeas(1:end-2)];
Imeas_res=Imeas(1:end-3);
%figure,plot(Imeas,Rmeas,'.')

%'qui contr', keyboard

mesh.fCondTer=fCondTer;
mesh.fCondTerZ=fCondTerZ;
mesh.fatt_dndT=fatt_dndT;
ContaBisez=0;

NVbias=length(mode.v0_dd);
while(indv<=NVbias) % outer bias-loop %$$$$$$$$$$$$$$$$$$$$$$$$
    
    elapsed_time=cputime;
    %    mode_old=mode;
    %    mesh_old=mesh;
    uvet_old=uvet;
    v0_dd = mode.v0_dd(indv)
    if(mode.Tflg & v0_dd>mode.minthermalvoltage & flagCallVELM) % OCCHIO!!! flagCallVELM
        Tprec=mesh.Tprec;
        % Recombination-only heating
%        'termico', keyboard
PPst=sum(modePlot.Pst_dd,1)+modePlot.Psp_dd;        
PElec=modePlot.ii_dd*1000.*modePlot.vv_dd;        
PDiss=PElec(end)-PPst(end);
mesh.PDiss=PDiss;
%'sto', keyboard
        mode.iplotTerm=1;
        %        mode.iTfig=1;
%                'DeltaT', keyboard
        [DeltaT,Tprec,PTherm,T_Contributi]=f_ThermicFun(Tprec,mesh,mode,StrTT,IPLOT);
        DeltaTJoule=T_Contributi{1}; 
        DeltaTRec_srhAu=T_Contributi{2}; 
        DeltaTRec_Ccap=T_Contributi{3};         
        DeltaTRec_RAD=T_Contributi{4};         
        DeltaTOptAbs=T_Contributi{5};   
        DeltaT=reshape(DeltaT,1,mesh.nn);
        %        'DeltaT', keyboard
    elseif(flagCallVELM==1)
        DeltaT=zeros(1,mesh.nn);
        DeltaTJoule=DeltaT; DeltaTOptAbs=DeltaT;
        DeltaTRec_srhAu=DeltaT; DeltaTRec_Ccap=DeltaT; 
        DeltaTRec_RAD=DeltaT; 
        Tprec=0; % initalization of Tprec=0 for thermic self-consistency
    end
    mesh.DeltaT=DeltaT;
    mesh.DeltaTvelm=DeltaT;
%    'qui delta', keyboard
    DeltaTmax=max(DeltaT);
    DeltaTmax_Joule=max(max(DeltaTJoule)); 
    DeltaTmax_srhAu=max(max(DeltaTRec_srhAu));
    DeltaTmax_Ccap=max(max(DeltaTRec_Ccap));
    DeltaTmax_RAD=max(max(DeltaTRec_RAD));
    DeltaTmax_OptAbs=max(max(DeltaTOptAbs));
  
    
    % giusto per non buttare via quel plot, faccio cos?.
    if(mode.oflg)
        n2Dtmp=0; p2Dtmp=0; n3Dtmp=0; p3Dtmp=0;
        for indQW=1:NQW
            n2Dtmp=n2Dtmp+mode.n2D{indQW};
            p2Dtmp=p2Dtmp+mode.p2D{indQW};
            n3Dtmp=n3Dtmp+mode.elec(mesh.inMQW{indQW});
            p3Dtmp=p3Dtmp+mode.hole(mesh.inMQW{indQW});
        end
        
        npMax=max(n2Dtmp+p2Dtmp)/2;
        nMax=max(n2Dtmp);
        pMax=max(p2Dtmp);
        n3Max=max(n3Dtmp);
        p3Max=max(p3Dtmp);
        %     nMax=max(mode.n2D);
        %     pMax=max(mode.p2D);
        %     n3Max=max(mode.elec(mesh.inQW));
        %     p3Max=max(mode.hole(mesh.inQW));
        mode.npMaxVet(indv)=npMax;
        mode.nMaxVet(indv)=nMax;
        mode.pMaxVet(indv)=pMax;
        mode.n3MaxVet(indv)=n3Max;
        mode.p3MaxVet(indv)=p3Max;
        
        % calcolo guadagno e derivate
        %    'qui derivate', keyboard
        NUOVA_Gain
    end
    
    load stopDD
    if istop==1
      if indv>1
        savePlot
        eval(['save ',Last_Workspace,'_Plot',rad_setting,' modePlot'])            
        eval(['save ',Last_Workspace,'stop'])
      end
        'pausa simulazione', keyboard
    end
    
    
    if(mode.oflg)
        %     flagPort=(flagCallVELM.*abs(DeltaTmax-DeltaTempVELMLast)>=mode.mintempVELM);
        DELTAT=abs(DeltaTmax-DeltaTempVELMLast)
        flagTemp=(DELTAT>=mode.mintempVELM);
        flagPort=(abs(npMax-npMaxLast)>=mode.minPorVELM);
        flagPot=(Pst>=.1);
        %        flagPot=(Pst*1000>=.1);
        %        if isfield(mode,'Scheck')
        %          flagTemp=(flagTemp &        max(full(mode.Scheck))<.8)
        %        end
        %flagTemp=0
        %        flagPot=0
        %       if isfield(mode,'Gmod')
        %        iGmo=length(find(mode.Gmod>0));
        %        if iGmo>2 & flagPot
        %         iVELM=1;
        %        end
        %       end
        %        if( (indv==1 | flagTemp | flagPort | flagPot | iVELM) & flagCallVELM==1)
        HARD=1
        flagPot=0;
        FLAG_velm0=( indv==1 | (flagCallVELM & flagTemp) | iVELM | flagPot) & HARD;
        flagvv=0;
        if indv>1
            if length(VELMInfo)>0
                vlast=mode.v0_dd(VELMInfo(indVELM).indVoltage);
                flagvv=vlast>=v0_dd;
            end
        end
        
        FLAG_velm=((FLAG_velm0 & ~flagvv & flagCallVELM) & HARD)
        
        %keyboard
        if FLAG_velm
            ' Calcolo VELM !!!!!! '
            ' Calcolo VELM !!!!!! '
            ' Calcolo VELM !!!!!! '
            ' Calcolo VELM !!!!!! '
            ' Calcolo VELM !!!!!! '
            if ivelMstop==1
                keyboard
            end
            
            %         iVELM=0;
            
            DeltaTempVELMLast=DeltaTmax;
            
            
            %            disp(['prima di VELM']),            keyboard
%            if(pippoflag==1)
%                load('Workspace_prima_VELM.mat')
%                disp('Ho caricato il workspace')
%                keyboard
%                pippoflag=0;
%            else
%                save('Workspace_prima_VELM.mat')
%                disp('Salvato workspace prima di VELM')
%            end
            
            %if(v0_dd>=2.31 & v0_dd<=2.39)
            % fileName=['Ws_PrimaVELM_V_',num2str(v0_dd),'.mat'];
            % save(fileName)
            % end
            
            colordef black
            %            'prima di Velm', keyboard
            [velm] = f_CallVELM(mesh,mode,mode1,ParVet,VelmOptions,fil_str);
            %            clear global
            %            mode.verbVELM=0;
            colordef white
            
            mode.Lm=velm.Lm(1:mode.nmodes); % prendo solo il primo modo (per ora)
            mode.vlambda=velm.vlambda;
            mode.NQW=velm.NQW; % prendo solo il primo modo (per ora)
            mode.Gamma_z=velm.Gamma_z;
            mode.nindexQW=velm.nindexQW;
            mode.fPES=velm.fPES*1000;
            vph=Clight./mode.nindexQW;
            mode.fPdif=TAROCCO*velm.fPdif/vph/1000;
            % mode.fPdif=velm.fPdif/2.3; % ATTENZIONE CHE HO MESSO UN /2 PER FAR TORNARE LE COSEEEEEEEEEEEEEEEEEEEEE
            mode.E2=velm.E2; % original
            %            'vagawt['
            %            keyboard
            
            % % setting non-zero field to avoid problems
            % indZero=find(mode.E2(:,1)<1e-5);
            % mode.E2(indZero,1)=10;
            % mode.E2=mode.E2+10;
            
            mode.Inviluppo_SW=velm.Inviluppo_SW;
            
            indVELM=indVELM+1;
            
            VELMInfo(indVELM).Lm=velm.Lm(1:mode.nmodes); % prendo solo il primo modo (per ora)
            VELMInfo(indVELM).vlambda=velm.vlambda;
            VELMInfo(indVELM).NQW=velm.NQW; % prendo solo il primo modo (per ora)
            VELMInfo(indVELM).Gamma_z=velm.Gamma_z;
            VELMInfo(indVELM).nindexQW=velm.nindexQW;
            VELMInfo(indVELM).nindexRef=velm.nindexRef;
            VELMInfo(indVELM).fPES=velm.fPES;
            VELMInfo(indVELM).fPdif=velm.fPdif;
            VELMInfo(indVELM).E2=velm.E2;
            VELMInfo(indVELM).Inviluppo_SW=velm.Inviluppo_SW;
            VELMInfo(indVELM).SW=velm.SW;
            
            VELMInfo(indVELM).v0_dd=v0_dd;
            VELMInfo(indVELM).DeltaTmax=DeltaTmax;
            VELMInfo(indVELM).indVoltage=indv;
            
            
            
            mode1.LamVelm(indVELM)=velm.vlambda(1);
            mode1.TmVelm(indVELM)=DeltaTmax;
            mode1.NmVelm(indVELM)=nMax;
            mode1.PmVelm(indVELM)=pMax;
            mode1.IndVelm(indVELM)=indv;
            mode1.LmVelm(indVELM,:)=velm.Lm(1);
            mode1.E2Velm(indVELM,:)=velm.E2(1,:)';
            npMaxLast=npMax;
            % VELMinput
            salva_per_VELM
            
            flagVELMWasCalled=1;
            
            colordef white
            
        end
        
    end
    
    T=mode.T0+DeltaT;
    
    Vt=kB.*T./qel;
    
    mesh=loadmacro(geom,mesh,mode,T);
    
    mode.Vt=Vt;
    Vt_tr=pdeintrp(mesh.node,mesh.triangle(1:4,:),Vt.'); % T on triangles
    mode.Vt_tr=Vt_tr;
    
    mode=f_ComputeNeutrality(geom,mesh,mode);
    
    if(mode.Tflg==1)
        %        mode.DeltaTmax(indv)=DeltaTmax;
        mesh.Tprec=Tprec;
    end
    
    mode_old=mode;
    mesh_old=mesh;
    
    iterMag1=0;
    t0_Newton=1;
    
    fprintf('\n Drift-diffusion analysis\n')
    mode.iconv(indv) = 0; % flag, 1 convergence, 2 convergence failure
    iter = 0; % number of outer iterations
    iter_t = 0; % number of inner iterations
    mode.res(indv) = Inf; % residual error
    % t = 1; % damping parameter % QUAAAAAA
    %
    
    
    
    
    
    %    if(v0_dd>=2.05 & v0_dd<2.1)
    if(v0_dd>=20.05)
        
        
        % effetti=[1 1 1 1 1]
        effetti=[0 0 0 0 0]
        
        if mode.Pst_dd(end)>2
            effetti=[0 0 0 0 0];
        end
        
        
        if effetti(1)==1
            mode.fPES=VELMInfo(13).fPES*1e3;
        end
        if effetti(2)==1
            mode.Gamma_z=VELMInfo(13).Gamma_z;
        end
        if effetti(3)==1
            mode.fPdif=VELMInfo(13).fPdif/vph/1e3;
        end
        if effetti(4)==1
            mode.E2=VELMInfo(13).E2;
        end
        if effetti(5)==1
            mode.vlambda=VELMInfo(13).vlambda;
        end
        'fPES Gamma_z  fPdif E2 lam'
        
        
        %     'prima di while iter',keyboard
    end
    
    
    
    clear Res
    
    
    while(iter<mode.maxiter && iter_t<mode.maxiter_t), % inner while-loop
        % assem jacobian matrix in time domain
        mode.ind_v0=indv;

        %      'prima di assem', keyboard
        [Kmat0,Jmat0,Jmat2,uvet,rvet,mode]=assem(geom,mesh,mode,uvet,v0_dd);
        if(mode.oflg)
            if full(mode.Scheck)>1
                iterMag1=iterMag1+1;
                if iterMag1>1
                    %         t0_Newton=.5
                    % break
                    
                end
            end
        end
        % solve matrix equation
        Jmat = Kmat0 + Jmat0 + Jmat2;
        [R,C] = dgsequ(Jmat);
        
        
        
        % la merda suprema peccato terrificante quando uno lo fa un
        % programmatore muore
        %   R(ss+1,:)=R(ss+1,:)/1e15;
        %   R(vv+[1:mesh.nnx],:)=R(vv+[1:mesh.nnx],:)/1e15;
        %   R(ww+[1:mesh.nnx],:)=R(ww+[1:mesh.nnx],:)/1e15;
        
        
        
        % che paura
        R(ss+1:NUMERO_MODI,:)=R(ss+1:NUMERO_MODI,:)*1e30;
        % R(vv+[1:mesh.nnx],:)=R(vv+[1:mesh.nnx],:)/1e15;
        % R(ww+[1:mesh.nnx],:)=R(ww+[1:mesh.nnx],:)/1e15;
        
        
        res = norm(R*rvet);
        Res(iter+1)=res;
        
        %        'residuo', keyboard
        %    if full(mode.Scheck)>1
        %        disp('dopo calcolo residuo')
        %        keyboard
        %        end
        
        %==============================================================================================100
        if(isnan(res)),
            disp('================================================')
            disp('=         Residual error is NaN! :-(           =')
            disp('================================================')
            break
        end
        % if(res<mode.res(indv)) % uncomment for Bank-Rose
        if(iter==0 && mode.report), fprintf('Iteration     Residual\n')
            fprintf(sprintf('%4i%20.4e\n',iter,res)), end
        if(iter>=1 && mode.report), fprintf('Iteration     Residual            t\n')
            fprintf(sprintf('%4i%20.4e%15.2e\n',iter,res,t)), end
        mode.res(indv)=res; % save residual
        
        if(res<mode.tolconv), mode.iconv(indv)=1; mode.res(indv) = res; break, end
        iter = iter + 1;
        if(res<mode.tolconv*1e3)
            if res>Res(iter)
                mode.iconv(indv)=1; mode.res(indv) = res; break,
            end
        end
        
        % uvet2 = uvet; % save solution, uncomment for Bank-Rose
        iter_t=0; t=t0_Newton; % restore damping parameter % QUAAAAAA
        %=========================================
        % solve matrix equation
        Jmat = R*Jmat*C;
        [Lmat,Umat,pm,qm] = lu(Jmat);
        deltau = - C*(qm*(Umat\(Lmat\(pm*(R*rvet)))));
        uvet = uvet + t*deltau.'; % update
        %============================================
        % else % uncomment for Bank-Rose
        %     t=t/2; % uncomment for Bank-Rose
        %     iter_t = iter_t + 1; % uncomment for Bank-Rose
        %     uvet = uvet2 + t*deltau.'; % update, uncomment for Bank-Rose
        % end % uncomment for Bank-Rose
    end % outer while-loop
    if(not(mode.iconv(indv)) & not(isnan(res))),
        disp('================================================')
        disp('=        Convergence failure!!!! :-(           =')
        disp('================================================')
        Scheck=mode.gmod./mode.lmod';
        indWrongScheck=find(Scheck>=1);
        if length(indWrongScheck)>1
            disp('================================================')
            disp('Too much code breaks due to Scheck > 1 condition')
            disp(['Restarting code with Pst=',num2str(mode.ScheckMultiplicationFactor),'*Pst'])
            disp('================================================')
            
            uvet(ss+[indWrongScheck])=mode.ScheckMultiplicationFactor*uvet(ss+[indWrongScheck]);
            Pst=uvet(ss+[indWrongScheck]);
        end
    end
    if ((not(mode.iconv(indv)) | isnan(res)) & ContaBisez<ContaBisezMax)
        ContaBisez=ContaBisez+1
        %      keyboard
        uvet=uvet_old;
        mode=mode_old;
        mesh=mesh_old;
        flagConv=0;
        if(mode.oflg)
            if(flagVELMWasCalled==1),
                VELMInfo(indVELM).indVoltage=indv-1;
                flagVELMWasCalled=0;
            end
            flagCallVELM=0;
        end
        indv=indv-1;
        Vb1=mode.v0_dd(1:indv);
        VbNew=(mode.v0_dd(indv)+mode.v0_dd(indv+1))/2;
        Vb2=mode.v0_dd(indv+1:end);
        % aggiunge tensione
        mode.v0_dd=[Vb1,VbNew,Vb2];
        NVbias=NVbias+1;
    else
        flagConv=1;
        flagCallVELM=1;
        ContaBisez=0;
    end
    if (ContaBisez>=ContaBisezMax)
        break
    end
    %
    % =====================================================================
    % save DD solution
    % =====================================================================
    if(flagConv==1) % if residual is not NaN
        mode.phi=uvet(1:nn);
        iq = mesh.iq;
        mode.elec=uvet(nn+(1:nn));
        if((isfield(mode,'stats'))&&(strcmp(mode.stats,'Fermi')))
            tmp = mode.elec(iq)./mesh.Nc(iq);
            nferinv = invferdr(tmp,mode.tolinvferdr);
            mode.EFn(iq) = mode.ecb(iq) + Vt(iq).*nferinv; 
        end
        %
        mode.hole=uvet(pp+(1:nn));
        if((isfield(mode,'stats'))&&(strcmp(mode.stats,'Fermi')))
            tmp = mode.hole(iq)./mesh.Nv(iq);
            nferinv = invferdr(tmp,mode.tolinvferdr);
            mode.EFp(iq) = mode.evb(iq) - Vt(iq).*nferinv; 
        end
        %
        fprintf('Elapsed time: %g, s\n',(cputime-elapsed_time))
        mode.elapsed_time=(cputime-elapsed_time)/60;
        
        v_dd=uvet(qq+1); % saving voltage
        i_dd=uvet(rr+1); % saving current
        
        if(mode.oflg==0)
            mode.ii_dd(indv)=i_dd;
            mode.vv_dd(indv)=v_dd;
            mode.vv0_dd(indv)=v0_dd;
            % compute electric field (in time domain), V/um
            efieldx_t=-mode.phi*mesh.gradx*1e-4; efieldx_t(it)=0;
            efieldy_t=-mode.phi*mesh.grady*1e-4; efieldy_t(it)=0;
            mode.efield_x=reshape(pdeprtni(mesh.node,mesh.triangle(1:4,:),efieldx_t),1,mesh.nn);
            mode.efield_y=reshape(pdeprtni(mesh.node,mesh.triangle(1:4,:),efieldy_t),1,mesh.nn);
            mode=f_EvalHeatingTerms(geom,mesh,mode);
        else
            Pst=uvet(ss+[1:mode.nmodes]);
            if(max(mode.Scheck)>=1)
                Scheck=mode.Scheck;
                mode=mode_old;
                mesh=mesh_old;
                uvet=uvet_old;
                if(flagVELMWasCalled==1),
                    VELMInfo(indVELM).indVoltage=indv-1;
                    flagVELMWasCalled=0;
                end
                indv=indv-1;
                Vb1=mode.v0_dd(1:indv);
                VbNew=(mode.v0_dd(indv)+mode.v0_dd(indv+1))/2;
                VbNew=v0_dd;
                Vb2=mode.v0_dd(indv+1:end);
                % aggiunge tensione
                mode.v0_dd=[Vb1,VbNew,Vb2];
                NVbias=NVbias+1;
                flagCallVELM=0;
                
                %%%%%%%%
                
                % indv=indv-1;
                % Vb1=mode.v0_dd(1:indv);
                % VbNew=(mode.v0_dd(indv)+mode.v0_dd(indv+1))/2;
                % Vb2=mode.v0_dd(indv+1:end);
                % mode.v0_dd=[Vb1,VbNew*[1 1],Vb2];
                % NVbias=NVbias+1;
                % flagCallVELM=0;
                % countScheck=countScheck+1;
                %%%%%%%
                
                countScheck=countScheck+1;
                if(countScheck>mode.maxScheckRepeat)
                    indWrongScheck=find(Scheck>=1);
                    disp('================================================')
                    disp('Too much code breaks due to Scheck > 1 condition')
                    disp(['Restarting code with Pst=',num2str(mode.ScheckMultiplicationFactor),'*Pst'])
                    disp('================================================')
                    %                   keyboard
                    uvet(ss+[indWrongScheck])=mode.ScheckMultiplicationFactor*uvet(ss+[indWrongScheck]);
                    Pst=uvet(ss+[indWrongScheck]);
                end
            else
                countScheck=0;
                flagCallVELM=1;
                mode.vv0_dd(indv)=v0_dd;
                mode.vv_dd(indv)=v_dd;
                mode.ii_dd(indv)=i_dd;
                mode.Pst_dd(:,indv)=Pst(:);
                mode.Psp_dd(:,indv)=mode.Psp;
                mode.lambda(:,indv)=mode.vlambda;
                mode.DeltaTmax(indv)=DeltaTmax;
                mode.DeltaTmax_Joule(indv)=DeltaTmax_Joule;
                mode.DeltaTmax_srhAu(indv)=DeltaTmax_srhAu;
                mode.DeltaTmax_Ccap(indv)=DeltaTmax_Ccap;
                mode.DeltaTmax_RAD(indv)=DeltaTmax_RAD;
                mode.DeltaTmax_OptAbs(indv)=DeltaTmax_OptAbs;
                
                mode.Gmod(:,indv)=full(mode.gmod);
                mode.Lmod(:,indv)=mode.lmod;
                mode.nQW{indv}=mode.n2D;
                mode.pQW{indv}=mode.p2D;
                mode.PTherm(:,indv)=PTherm;
% da togliere
               mode.IntCcapN(indv)=mode.IntCcapn;
               mode.IntCcapP(indv)=mode.IntCcapp;  
% calcolo leakage
               leakage               
                
                % compute electric field (in time domain), V/um
                efieldx_t=-mode.phi*mesh.gradx*1e-4; efieldx_t(it)=0;
                efieldy_t=-mode.phi*mesh.grady*1e-4; efieldy_t(it)=0;
                mode.efield_x=reshape(pdeprtni(mesh.node,mesh.triangle(1:4,:),efieldx_t),1,mesh.nn);
                mode.efield_y=reshape(pdeprtni(mesh.node,mesh.triangle(1:4,:),efieldy_t),1,mesh.nn);
                mode=f_EvalHeatingTerms(geom,mesh,mode);
            end
        end
    end
    
    vind=[];
    for indplot=1:length(VELMInfo)
        vind=[vind,VELMInfo(indplot).indVoltage];
    end
    
    PPst=sum(mode.Pst_dd,1)+mode.Psp_dd;

    if(IPLOT==1)
       DD_PLOT
    end
    %     modevet(indv)=mode;
    %     meshvet(indv)=mesh;
    
    %     if(indv>1)
    %         fileName=['Workspace_',num2str(mode.v0_dd(indv-1)),'V.mat'];
    %         delete(fileName);
    %     end
    
    %     fileName=['Workspace_',num2str(mode.v0_dd(indv)),'V.mat'];
    %     if(mode.oflg==1)
    %         save(fileName,'modevet','meshvet','VELMInfo','geom','-v7.3');
    %     else
    %         save(fileName,'modevet','meshvet','geom','-v7.3');
    %     end
    % salvo per velm
    
    if(flagConv)
%        eval(['save ',Last_Workspace])
        savePlot
        eval(['save ',Last_Workspace,'_Plot',rad_setting,' modePlot'])        
% save('Last_plot.mat','modePlot')        
    end
    
    indv=indv+1;
    
end
if irest==0
      indv=indv-1;
      eval(['save ',Last_Workspace])
end      
%      save Last_Workspace

PElec=mode.ii_dd*1000.*mode.vv_dd;
PDiss=PElec-PPst;
PDissmeas=Vmeas.*Imeas-Lmeas;
figure,
hold on
grid on
plot(Imeas,PDissmeas,'ro')
plot(mode.ii_dd.*1000,PDiss,'b','LineWidth',2)
plot(mode.ii_dd*1000,mode.PTherm,'k','LineWidth',2);
xlabel('Current, mA')
ylabel('Dissipated power, mW')
legend('Measurements','Electric - Optic','Thermic')

return

close all
figure(43892)
hold on
A=J_Y(3,:);
A=A-A(end);
Norm=A(1);
A=A./A(1);
B=J_Y(18,:)
B=B-B(end);
Norm=B(1);
B=B./B(1);
plot(mesh.xgrid,A,mesh.xgrid,B)






load Meas_ReliefVCSEL

Resistence=diff(mode.vv0_dd)./diff(mode.ii_dd);
Resistence=[Resistence(1),Resistence];
figure
hold on
grid on
plot(mode.ii_dd*1000,Resistence,Imeasinterp,Resmeasinterp)
xlabel('Current, mA')
ylabel('Resistance, \Omega')
xlim([0.5,mode.ii_dd(end)*1000])
ylim([80,240])
legend('Simulation','Measurements','Location','Best')

ResFake=8.5;

figure
hold on
grid on
plot(mode.vv_dd+ResFake*mode.ii_dd,mode.ii_dd*1000,Vmeas,Imeas)
xlabel('Voltage, V')
xlabel('Current, mA')
legend('Simulation','Measurements','Location','Best')

for indv=1:length(mode.nQW)
    nQWmax(indv)=max(mode.nQW{indv});
    pQWmax(indv)=max(mode.pQW{indv});
end
figure
hold on
grid on
box on
plot(mode.ii_dd*1000,nQWmax/mesh.WQW,mode.ii_dd*1000,pQWmax/mesh.WQW)
legend('electrons 2D','holes 2D')
xlabel('Current, mA')

%np=(nQWmax+pQWmax)/2/mesh.WQW;
%vph=3e10/3.6;

%[g] = f_InterpGainEH(np,848,mode.T0+mode.DeltaTmax);
%g=g*mode.fat_gain;
%figure(10)
%hold on
%grid on
%box on
%plot(mode.ii_dd*1000,g/vph)
%xlabel('Current, mA')

JN_X=reshape(mode.Jn_x,mesh.nny,mesh.nnx);
JN_Y=reshape(mode.Jn_y,mesh.nny,mesh.nnx);
JP_X=reshape(mode.Jp_x,mesh.nny,mesh.nnx);
JP_Y=reshape(mode.Jp_y,mesh.nny,mesh.nnx);


J_X=JN_X+JP_X;
J_Y=JN_Y+JP_Y;

mesh.X=reshape(mesh.node(1,:),mesh.nny,mesh.nnx);
mesh.Y=reshape(mesh.node(2,:),mesh.nny,mesh.nnx);
 step_x=1;
 step_y=1;
 sX=1:step_x:size(mesh.X,2);
 sY=30:step_y:size(mesh.X,1);
 X0=mesh.X(1,:)*1e4;
 Y0=mesh.Y(:,1)*1e4;
 dx=diff(X0);
 sogx=.3;
 su=0;
 sX=[];
 for k=1:length(dx)
  su=su+dx(k);
  if su>sogx
   su=0;
   sX=[sX k];
  end
 end 
 figure, plot(X0,'.')
 hold on
 plot(sX,X0(sX),'ro'), pausak
 
  dx=diff(Y0);
  sogx=.3;
  su=0;
  sY=[];
  for k=1:length(dx)
   su=su+dx(k);
   if su>sogx
    su=0;
    sY=[sY k];
   end
  end 
  figure, plot(Y0,'.')
  hold on
  plot(sY,Y0(sY),'ro'), pausak
 
 
 SX=mesh.X(sY,sX)*1e4;
 SY=mesh.Y(sY,sX)*1e4;

 figure,quiver(SX,SY,J_X(sY,sX),J_Y(sY,sX),.1), ylim([350 360]), xlim([0 10])



figure,contourf(mesh.X*1e4,mesh.Y*1e4,(J_X)),title('J_\rho'),shading interp,xlabel('\rho'),ylabel('z')%,set(gca,'zscale','log')
figure,contourf(mesh.X*1e4,mesh.Y*1e4,(J_Y)),title('J_z'),shading interp,xlabel('\rho'),ylabel('z')%,set(gca,'zscale','log')

figure,[hC,hC]=contourf(mesh.X,mesh.Y,J_X,-1e5:1e3:1e5),title('J_\rho')
set(hC,'linestyle','none')
caxis([-5e3,5e3])

figure,[hC,hC]=contourf(mesh.X,mesh.Y,J_Y,-1e5:1e3:1e5),title('J_\rho')
set(hC,'linestyle','none')
caxis([-5e3,5e3])

x=mesh.xgrid;
for indy=1:mesh.nny
    
    jn_y=JN_Y(indy,:)+JP_Y(indy,:);
    curr(indy)=trapz(x,2.*pi.*x.*jn_y);
    
    
end
figure,plot(mode.vv0_dd,1000*mode.ii_dd),xlabel('V'),title('current vs V')

figure,plot(1e7*mesh.ygrid,curr*1000),xlabel('z'),title('current vs z')

figure,plot(1e7*mesh.node(2,:),mode.HeatJoule_y),xlabel('z'),title('Joule heat z')

figure,plot(1e7*mesh.node(2,:),mode.HeatJoule_x),xlabel('z'),title('Joule heat \rho')

figure,plot(1e7*mesh.node(2,:),mode.HeatRec_nr),xlabel('z'),title('NR recomb heat')

el=mode.elec; el=el/max(el).*5;
ho=mode.hole; ho=ho/max(ho).*5;
% el=mesh.mobn_n.*mode.elec; el=el/max(el).*5;
% ho=mesh.mobp_n.*mode.hole; ho=ho/max(ho).*5;

figure,plot(1e7*mesh.node(2,:),1./mode.sigma,1e7*mesh.node(2,:),el,'r--',1e7*mesh.node(2,:),ho,'g--'),xlabel('z'),title('\rho = 1/\sigma')

ind=1:mesh.nny;
ind=1:mesh.nn;
figure,hold on,
plot(1e7*mesh.node(2,ind),mode.elec(ind),'b',1e7*mesh.node(2,ind),mode.hole(ind),'r','Linewidth',2),
% plot(1e7*mesh.node(2,ind),mode.N2D(ind),'k',1e7*mesh.node(2,ind),mode.P2D(ind),'m','Linewidth',2),
plot(1e7*mesh.node(2,ind),mode.elec(ind)+mode.N2D(ind),'k',1e7*mesh.node(2,ind),mode.hole(ind)+mode.P2D(ind),'m','Linewidth',2),
set(gca,'yscale','log')

figure,plot((mesh.node(2,end)-mesh.node(2,:))*1e7,mode.sigma)
SIGMA=reshape(mode.sigma,mesh.nny,mesh.nnx);


figure,[hC,hC]=contourf(mesh.X,mesh.Y,log10(1./SIGMA),-5:0.1:1),title('\sigma')
set(hC,'linestyle','none')
caxis([-3,0.5])


JOULE=reshape(mode.HeatJoule,mesh.nny,mesh.nnx);
figure,[hC,hC]=contourf(mesh.X,mesh.Y,log10(JOULE),-16:0.1:-3),title('Joule')
set(hC,'linestyle','none')
caxis([-9,-3])


NONRAD=reshape(mode.HeatRec_nr,mesh.nny,mesh.nnx);
figure,[hC,hC]=contourf(mesh.X,mesh.Y,log10(abs(NONRAD)),-30:0.1:-3),title('Non-rad sum')
set(hC,'linestyle','none')
caxis([-9,-3])

NONRAD_QW=reshape(mode.HeatRec_nr_QW,mesh.nny,mesh.nnx);
figure,[hC,hC]=contourf(mesh.X,mesh.Y,log10(abs(NONRAD_QW)),-30:0.1:-3),title('Non-rad QW')
set(hC,'linestyle','none')
caxis([-9,-3])

NONRAD_bulk=reshape(mode.HeatRec_nr_bulk,mesh.nny,mesh.nnx);
figure,[hC,hC]=contourf(mesh.X,mesh.Y,log10(abs(NONRAD_bulk)),-30:0.1:-3),title('Non-rad bulk')
set(hC,'linestyle','none')
caxis([-9,-3])



% save('Workspace_2_5_V.mat')

% vlam=[]; vind=[];
% for indplot=1:length(VELMInfo)
%     vind=[vind,VELMInfo(indplot).indVoltage];
%     vlam=[vlam,VELMInfo(indplot).vlambda];
% end
%
%
% voltind=1:10:136;
% for indplot=1:length(voltind)
%
%     figure(289)
%     hold on
%     grid on
%     box on
%     ind=voltind(indplot);
%     np=(mode.nQW{ind}+mode.pQW{ind})./2./mesh.WQW;
%     plot(1e7*mesh.xgrid,np)
%
% end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Saving output data
if((isfield(mode,'symmetry'))&&(strcmp(mode.symmetry,'Cylindrical-X')||strcmp(mode.symmetry,'Cylindrical-Y')))
    save(['output_' structureName '_cylindrical.mat'],'mesh','mode','geom')
else
    save(['output_' structureName '.mat'],'mesh','mode','geom')
end


return

figure(2)
set(gcf,'Position',[425 164 814 646])
axis on
grid on
hold on
box on
plot(mesh.node(2,:)*1e7,mode.ecb,'b','LineWidth',2)
plot(mesh.node(2,:)*1e7,mode.evb,'b','LineWidth',2)
plot(mesh.node(2,:)*1e7,mode.EFn,'k-.','LineWidth',2)
plot(mesh.node(2,:)*1e7,mode.EFp,'k-.','LineWidth',2)
xlabel('z (nm)')
ylabel('Band diagram (eV)')
set(gca,'FontSize',14)

figure(3)
set(gcf,'Position',[369 320 1031 420])
subplot(1,2,1)
axis on
grid on
hold on
box on
plot(mesh.node(2,:)*1e7,mode.elec,'LineWidth',2)
if(mode.oflg)
    plot(mesh.node(2,:)*1e7,mode.elec+mode.N2D,'k-','LineWidth',2)
end
xlabel('z (nm)')
ylabel('n(z) (1/cm^3)')
set(gca,'FontSize',14)
set(gca,'YScale','log')
subplot(1,2,2)
axis on
grid on
hold on
box on
plot(mesh.node(2,:)*1e7,mode.hole,'LineWidth',2)
if(mode.oflg)
    plot(mesh.node(2,:)*1e7,mode.hole+mode.P2D,'k-','LineWidth',2)
end
xlabel('z (nm)')
ylabel('p(z) (1/cm^3)')
set(gca,'FontSize',14)
set(gca,'YScale','log')

