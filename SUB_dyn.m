
CurDyn=[];
flagCallVELM=1;
ContaBisezMax=3;

if IPLOT==1
    portfig=figure(kpar);
    set(portfig,'pos',[1546         611         364         373])
end

if IPAR==41
    if kpar==1
        save iplot IPLOT
    else
        load iplot
    end
end

%========================================================================76
% solve at equilibrium

%'rect', keyboard
% mesh generation

if irest==0
    
    mesh=rectmesh(geom); % create tensorial mesh
    if mode.quasi1D==1
        %         mesh.puntatore=mesh.nny+(2:mesh.nny-1);
        mesh.puntatore=(2:mesh.nny-1);
    end
    
    T=mode.T0*ones(1,mesh.nn);
    
    mesh.DeltaT=T-mode.T0;
    % mesh.DeltaT=0;
    
    mesh.dz_Term=ParMore.dz_Term;
    %'dopo mesh', keyboard
    
    mode.DeltaTmax=0;
    
    fil_str=[structureName,'.str'];
    mode.matgain=0;
    mode.vlambda=850*ones(NUMERO_MODI,1); % inizializzazione
    
    % Temporaneo: da sostituire con ciò che sarà preso dal codice Density
    % Matrix
    Dn=mode.matgain;
    mode.DeltaN=0;
    
    mode.v0_dd=dd;
    
    s_LoadConstants
    %'load', keyboard
    TKT=T;
    if isfield(mode,'KT')
        if mode.KT==0
            TKT=mode.T0*ones(size(T));
        end
    end
    
    Vt=kB.*TKT./qel;
    mode.Vt=Vt;
    
    %'load', keyboard
    
    
    
    % Material dependent quantities are computed and extracted
    
    mesh=loadmacro(geom,mesh,mode,T);
    
    %'macro', keyboard
    
    mesh.Tgain=mesh.T;
    
    
    
    if abs(IPLOT)>=1
        node='off'; triangle='off'; color='on'; vpath='off'; arrows='off';
        sd=1:geom.nd; scale=[]; cmap=[]; grido='on'; cbar='off';
        %sd=1:geom.nd; scale=[]; cmap=[]; grido='off'; cbar='off';
        figure
        plot_tri_meshPLOT(geom,mesh,[],sd,color,triangle,node,grido,cbar,cmap,scale,vpath,arrows)
        title([' Grigliato: Punti X Y   ', num2str(length(mesh.xgrid)),'  x  ',num2str(length(mesh.ygrid))])
        axis normal
        drawnow
        pausak
        node='off'; triangle='off'; color='on'; vpath='off'; arrows='off';
        %sd=1:geom.nd; scale=[]; cmap=[]; grido='on'; cbar='off';
        sd=1:geom.nd; scale=[]; cmap=[]; grido='off'; cbar='off';
        figure
        plot_tri_meshPLOT(geom,mesh,[],sd,color,triangle,node,grido,cbar,cmap,scale,vpath,arrows)
        title([' Grigliato: Punti XxY   ', num2str(length(mesh.brad))])
        
        drawnow
        pausak
        
        
        figure(1001)
        set(gcf,'Position',[144  179 1107  787])
        subplot(2,2,1)
        axis on
        grid on
        hold on
        box on
        plot(1e7*mesh.node(2,1:mesh.nny),mesh.xmol(1:mesh.nny),'b','LineWidth',2)
        xlabel('z (nm)')
        ylabel('Molar fraction')
        subplot(2,2,2)
        axis on
        grid on
        hold on
        box on
        plot(1e7*mesh.node(2,1:mesh.nny),mesh.dop_d(1:mesh.nny),'b','LineWidth',2)
        plot(1e7*mesh.node(2,1:mesh.nny),mesh.dop_a(1:mesh.nny),'r','LineWidth',2)
        xlabel('z (nm)')
        ylabel('Doping profile cm^{-3}')
        legend('N_D','N_A','location','Best')
        set(gca,'yscale','log')
        drawnow
        %
        disp([' ']);
        disp(['Verifying molar fraction and doping properties']);
        disp([' ']);
        subplot(2,2,3)
        axis on
        grid on
        hold on
        box on
        % Eg=mesh.Eg;
        % plot(1e7*mesh.ygrid,(Eg(1:mesh.nny)),'b','LineWidth',2)
        %Ec=-mesh.affinity;
        %Ec = mesh.Eg/2 - mesh.T*kB/qel./2.*log(mesh.Nv./mesh.Nc) + mesh.phi_r;
        Ec = mesh.Eg/2  + mesh.phi_r;
        Ev=Ec-mesh.Eg;
        plot(1e7*mesh.ygrid,(Ec(1:mesh.nny)),'b','LineWidth',2)
        plot(1e7*mesh.ygrid,(Ev(1:mesh.nny)),'r','LineWidth',2)
        xlabel('z (nm)')
        ylabel('Cold Band diagram, eV')
        legend('E_c','E_v','location','Best')
        %set(gca,'yscale','log')
        drawnow
        
        pausak
        if IPLOT==2
            IPLOT=0;
        end
    end %IPLOT
    
    %
    % return
    
    %'solve1', keyboard
    
    mode=p_solve1(geom,mesh,mode); % first run: thermodynamic equilibrium
    mode.n3Di=mode.elec;
    mode.p3Di=mode.hole;
    
    % Carriers for VELM
    mode.elecABS=mode.elec*1e-18*mode.CarrierNorm;
    mode.holeABS=mode.hole*1e-18*mode.CarrierNorm;
    
    % nmodes=length(mode.Lm); % number of optical modes (VELM)
    % mode.nmodes=4;
    %'NUM sui', keyboard
    mode.nmodes=NUMERO_MODI;
    mode.lambda=850*ones(mode.nmodes,1);
    
    % To initialize matgain, we use carriers in the central well
    %indCenter=floor((mesh.NMQW+1)/2);
    %nQW=median(mode.n2D{indCenter});
    %pQW=median(mode.p2D{indCenter});
    %TQW=300*ones(size(nQW));
    %[g_eq] = f_InterpGainEH(nQW,pQW,mode.Deltalam+mode.vlambda(1)*ones(size(nQW)),TQW);
    %mode.matgain=g_eq;
    g_eq=0;
    mode.matgain=g_eq;
    
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
        nnQW=mesh.nnxQW{1};         % number of lateral points in the QW
    else
        nnQW=0;
    end
    NQW=mesh.NMQW;                % number of quantum wells
    %
    nl=mode.ntrap;                % number of trap levels
    %
    nn=mesh.nn;                   % -> elec eqs.
    pp=nn+mode.nflg*nn;           % -> hole eqs.
    tt=pp+mode.pflg*nn;           % -> trap eqs.
    qq=tt+mode.tflg*nl*nn;        % -> circuit eqs.
    rr=qq+nm;                     % -> current eqs.
    vv=rr+nm;                     % -> n2D eqs.
    ww=vv+length(mesh.xgrid)*NQW;               % -> p2D eqs.
    ss=ww+length(mesh.xgrid)*NQW;               % -> optical rate eqs.
    %
    % total number of equations/unknowns
    neq=nn+mode.nflg*nn+mode.pflg*nn+mode.tflg*nl*nn+2*nm+mode.oflg*not(mode.firstrun)*2*NQW*length(mesh.xgrid)+mode.oflg*mode.nmodes;
    
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
            uvet(vv+(indQW-1)*length(mesh.xgrid)+(1:nnQW))=mode.n2D{indQW};
        end
    end
    if(isfield(mode,'p2D'))
        for indQW=1:NQW
            uvet(ww+(indQW-1)*length(mesh.xgrid)+(1:nnQW))=mode.p2D{indQW};
        end
    end
    
    if(mode.oflg & isfield(mode,'Pst')),
        uvet(ss+(1:mode.nmodes))=mode.Pst;
    end
    %
    
    s_LoadConstants
    
    PTherm=0;
    
    mode_old=mode;
    indv=1;
    
    DeltaTempVELMLast=0;
    npMaxLast=0;
    Pst=0;
    indVELM=0;
else
    Glut4Dinc(mode)
end %irest

%'irest', keyboard
if irest>0
    if kpar==1
        save LSW Last_Workspace irest rad_setting IREST nomeSav iDyn kpar kpvet IPvet NP PMAT rad_setting rad_settingV
    else
        save LSW Last_Workspace irest rad_setting IREST nomeSav iDyn kpar kpvet IPvet NP PMAT rad_setting  rad_settingV MODEplot
    end
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
    %'ferma', keyboard
    
    %     if irest==2
    %         eval(['load ',Last_Workspace])
    %     else
    %         eval(['load ',Last_Workspace,'_stop'])
    %     end
    
    load LSW
    fis= strfind(Last_Workspace,'\');
    Last_Workspac=Last_Workspace(fis(end)+1:end);
    %'ver', keyboard
    %DirName=structure(1:fis(end));
    eval(['load ',nomeSW,'Contributi_',nomeSav])
    %        load Contributi
    ASSEGNO_Effetti
    eval(['settings_vari',rad_setting])
    if IPAR>0
        ppar=PMAT{IPAR}(kpar);
        stringa_parametro=[Svar{IPAR},'=',num2str(ppar)]
        eval([stringa_parametro,';'])
    end
    ASSEGNO_mode
    %        'Effetti cuora', keyboard
    le_VELMinfo=length(VELMInfo);
    %'ver', keyboard
    %        load Last_Workspace
    %clear MODEplot
    indv=indv+1;
    
    if irest==2
        mode.v0_dd=[mode.v0_dd,Vadd];
        
    end
    
end


pippoflag=0;
if(pippoflag==1)
    disp('hai messo pippoflag, ne sei proprio convinto?')
    keyboard
end


%figure,plot(Imeas,Rmeas,'.')


%settings_vari_sub

if IOLDsw==0
    mesh.fCondTer=fCondTer;
    mesh.fCondTerZ=fCondTer*fCondTerZ;
else
    mesh.fCondTer=fCondTer;
    mesh.fCondTerZ=fCondTerZ;
end

mesh.dndT=fatt_dndT*dndT;
mesh.dndT1D=dndT1D;

ContaBisez=0;
DTM0=mode.T300;
DeltaTmax=0;
I_mA=0;
if ~exist('Imassimo')
    Imassimo=1000;
end
NVbias=length(mode.v0_dd);

CondPotBreak=0;
clear PtotV

mode.deltaT=0;
mode.IdriveON=0;

while(indv<=NVbias && DeltaTmax<DTM0 && I_mA<Imassimo && CondPotBreak==0 ) % outer bias-loop %$$$$$$$$$$$$$$$$$$$$$$$$
    
    elapsed_time=cputime;
    %    mode_old=mode;
    %    mesh_old=mesh;
    uvet_old=uvet;
    v0_dd = mode.v0_dd(indv)
    
    if(mode.Tflg && v0_dd>mode.minthermalvoltage && flagCallVELM) % OCCHIO!!! flagCallVELM
        Tprec=mesh.Tprec;
        % Recombination-only heating
        mode.iplotTerm=0;
        PPst=sum(modePlot.Pst_dd,1)+modePlot.Psp_dd;
        PElec=modePlot.ii_dd*1000.*modePlot.vv_dd;
        PDiss=PElec(end)-PPst(end);
        mode.PDiss(indv)=PDiss;
mode.PDissPred(indv)=PDiss;
        mode.iplotTerm=1;
        %        mode.iTfig=1;
        %           'DeltaT', keyboard
        if isfield(mode,'IsoThermal')
            IsoThermal=mode.IsoThermal;
        else
            IsoThermal=0;
        end
        if IsoThermal==0
            if isfield(mode,'quasi1D') && mode.quasi1D==1
                [DeltaTold,PTherm,T_Contributi]=f_ThermD1ANA(mesh,mode,StrTT);
                Tprec=0;
                
                StimaTempWU
                DeltaT=DeltaTold*Fat_STIMA_Temp;
                
                mode.FatMob=polyval(cot,max(max(DeltaT)));
            else
                [DeltaTold,Tprec,PTherm,T_Contributi,K,TotalHeat]=f_ThermicFun(Tprec,mesh,mode,StrTT,IPLOT);
                mode.FatMob=1;
                
                StimaTempWU
                DeltaT=DeltaTold*Fat_STIMA_Temp;
            end
        else
            if mode.quasi1D==1
                mesh.DeltaT=0;
            else
                DeltaT=zeros(1,mesh.nn);
            end
            PTherm=0;
            for tcj=1:5
                T_Contributi{tcj}=0;
            end
        end
        
        %[DeltaT,Tprec,PTherm,T_Contributi]=f_ThermicFun(Tprec,mesh,mode,StrTT,IPLOT);
        DeltaTJoule=T_Contributi{1};
        DeltaTRec_srhAu=T_Contributi{2};
        DeltaTRec_Ccap=T_Contributi{3};
        DeltaTRec_RAD=T_Contributi{4};
        DeltaTOptAbs=T_Contributi{5};
        DeltaT=C_Temp*reshape(DeltaT,1,mesh.nn);
    elseif(flagCallVELM==1)
        DeltaT=zeros(1,mesh.nn);
        DeltaTJoule=DeltaT; DeltaTOptAbs=DeltaT;
        DeltaTRec_srhAu=DeltaT; DeltaTRec_Ccap=DeltaT;
        DeltaTRec_RAD=DeltaT;
        Tprec=0; % initalization of Tprec=0 for thermic self-consistency
        
    end
    if isfield(mode,'epsNLg')
        nlG=mode.epsNLg;
        mode.nlG=1;
        if isfield(mode,'E2')
            %       'qui per nlG', keyboard
            UU = (mode.Pst.*mode.fPdif)*mode.E2;
            uudu=1./(1+mode.epsNLg*UU);
            mode.nlG=uudu(1:mesh.nnxQW{1});
        end
    else
        mode.nlG=1;
    end
    
    
    mesh.DeltaT=DeltaT;
    
    
    if effetti(6)==0
        mesh.DeltaTvelm=DeltaT;
    else
        'TEMPERATURA Fissa in VELM'
    end
    
    DeltaTmax=max(DeltaT);
    if DeltaTmax-mode.DeltaTmax(end)>DTM0/4
        'incremento di temperatura troppo grande, esco '
        break
    end
    DeltaTmax_Joule=max(max(DeltaTJoule));
    DeltaTmax_srhAu=max(max(DeltaTRec_srhAu));
    DeltaTmax_Ccap=max(max(DeltaTRec_Ccap));
    DeltaTmax_RAD=max(max(DeltaTRec_RAD));
    DeltaTmax_OptAbs=max(max(DeltaTOptAbs));
    
    % giusto per non buttare via quel plot, faccio così.
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
        
        if iold==2
            % calcolo guadagno e derivate
            %    'qui derivate', keyboard
            NUOVA_Gain
            % carico dati tau di Vallone
            global CO coMh Tn Th
            if mode.taucarrier==1
                load saTaue
            end
        end
    end
    
    load stopDD
    if istop==1
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
            
            
            if effetti(7)==1
                VelmOptions.ianti_gui=0;  % tolgo anti_guiding
            end
            %if v0_dd==2.34
            %            save VE1
            %             if v0_dd>=2
            %                        'prima di Velm', keyboard
            %             end
            colordef black
            %                        'prima di Velm', keyboard
            iLoadVELM=0; % 1, check if VELM results are already computed; 0, always enters in VELM
            %
            if (iLoadVELM==1 && exist(['dati/VELM/velm_',strSave,'_quasi1D=',num2str(mode.quasi1D),'.mat'])==2 && v0_dd==0)
                % File exists.
                load(['velm_',strSave,'_quasi1D=',num2str(mode.quasi1D),'.mat'])
            else
                % File does not exist.
                if isfield(mode,'quasi1D') && mode.quasi1D==1
                    [velm] = f_CallVELM_1D(mesh,mode,mode1,ParVet,VelmOptions,fil_str,mode.quasi1D);
                else
                    [velm] = f_CallVELM(mesh,mode,mode1,ParVet,VelmOptions,fil_str);
                end
                % stores VELM output to re-use it for successive simulation
                if v0_dd==0
                    save(['dati/VELM/velm_',strSave,'_quasi1D=',num2str(mode.quasi1D)],'velm')
                end
            end
            
            %            clear global
            %            mode.verbVELM=0;
            colordef white
            if isfield(mode,'quasi1D') && mode.quasi1D==1
                velm.E2=ones(mode.nmodes,mesh.nnx)/(mode.AreaOx);
                velm.SW=0;
            end
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
            %'salve lambda', keyboard
            mode1.TmVelm(indVELM)=DeltaTmax;
            
            mode1.NmVelm(indVELM)=nMax;
            mode1.PmVelm(indVELM)=pMax;
            mode1.IndVelm(indVELM)=indv;
            mode1.LmVelm(indVELM,:)=velm.Lm(1);
            mode1.E2Velm(indVELM,:)=velm.E2(1,:)';
            npMaxLast=npMax;
            %	                        'salva per Velm', keyboard
            salva_per_VELM
            
            flagVELMWasCalled=1;
            
            colordef white
        elseif (isfield(mode,'quasi1D')==1 && mode.quasi1D==1) && indVELM>=3
            fitStart=2; fitDegree=2;
            [mode]=velmFitting(indVELM,fitStart,fitDegree,VELMInfo,mode,v0_dd);
            %             'Interp',keyboard
        end
    end
    
    T=mode.T0+DeltaT;
    mesh.Tgain=T;
    %    'qui ver T', keyboard
    T=mode.T0+DeltaT*mode.C_Temp_DD;
    
    TKT=T;
    if isfield(mode,'KT')
        if mode.KT==0
            TKT=mode.T0*ones(size(T));
        end
    end
    
    Vt=kB.*TKT./qel;
    %    mode.Vt=Vt;
    
    mesh=loadmacro(geom,mesh,mode,T);
    
    mode.Vt=Vt;
    Vt_tr=pdeintrp(mesh.node,mesh.triangle(1:4,:),Vt.'); % T on triangles
    mode.Vt_tr=Vt_tr;
    
    mode_old=mode;
    mesh_old=mesh;
    
    mode=f_ComputeNeutrality(geom,mesh,mode);
    
    if(mode.Tflg==1)
        %        mode.DeltaTmax(indv)=DeltaTmax;
        mesh.Tprec=Tprec;
    end
    
    % calcolo carrier-dependent tau Vallone
    if(mode.taucarrier==1)
        indQW=1; % primo QW
        n2D=mode.n2D{indQW}/mesh.WMQW{indQW};
        p2D=mode.p2D{indQW}/mesh.WMQW{indQW};
        n3D=mode.elec(mesh.inMQW{indQW});
        p3D=mode.hole(mesh.inMQW{indQW});
        if(max(n2D) > 1e8)
            mode.Fat_cap_e=fit_tauEnl(n2D,n3D);
            mode.Fat_cap_h=fit_tauH(p2D,p3D);
            %            'VALLONE', keyboard
        end
    end
    
    
    iterMag1=0;
    t0_Newton=1;
    
    fprintf('\n Drift-diffusion analysis\n')
    mode.iconv(indv) = 0; % flag, 1 convergence, 2 convergence failure
    iter = 0; % number of outer iterations
    iter_t = 0; % number of inner iterations
    mode.res(indv) = Inf; % residual error
    % t = 1; % damping parameter % QUAAAAAA
    %
    
    %             'keyboard le_VEl', keyboard
    
    %    if(v0_dd>=2.05 & v0_dd<2.1)
    if effetti(1)==1
        mode.fPES=VELMInfo(le_VELMinfo).fPES*1e3;
    end
    if effetti(2)==1
        mode.Gamma_z=VELMInfo(le_VELMinfo).Gamma_z;
    end
    if effetti(3)==1
        %            mode.fPdif=TAROCCO*velm.fPdif/vph/1000;
        mode.fPdif=TAROCCO*VELMInfo(le_VELMinfo).fPdif/vph/1e3;
    end
    if effetti(4)==1
        mode.E2=VELMInfo(le_VELMinfo).E2;
    end
    if effetti(5)==1
        mode.vlambda=VELMInfo(le_VELMinfo).vlambda;
    end
    'fPES Gamma_z  fPdif E2 lam'
    
    
    %                 'prima di while iter',keyboard
    
    clear Res
    
    
    while(iter<mode.maxiter && iter_t<mode.maxiter_t), % inner while-loop
        
        indCountMultiplicationFactor=1; %%% temporaneo
        
        
        % assem jacobian matrix in time domain
        mode.ind_v0=indv;
        
        %        figure, plot(mode.vv0_dd,mode.Pst_dd)
        %      'prima di assem', keyboard
        
        % prima di assem
        if iDyn==0
            %             mode_num=mode;
            %             uvet_num=uvet;
            
            %             if (mode.NumJac==1 && round(v0_dd,3)==round(mode.Vnum,3))
            %                 save WorkSpace_Jacobian_TJ_2columns.mat
            %                 keyboard
            %             end
            %             [Kmat0,Jmat0,Jmat2,uvet,rvet,mode]=assem(geom,mesh,mode,uvet,v0_dd);
            [Kmat0,Jmat0,Jmat2,uvet,rvet,mode]=assem_GBT(geom,mesh,mode,uvet,v0_dd);
        else
            [Kmat0,Jmat0,Jmat1,Jmat2,uvet,rvet,mode]=assem_FM(geom,mesh,mode,uvet,v0_dd);
        end
        
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
        %  R(ss+[1:NUMERO_MODI],:)=R(ss+[1:NUMERO_MODI],:)*1000;
        %   R(vv+[1:mesh.nnx],:)=R(vv+[1:mesh.nnx],:)/1e15;
        %   R(ww+[1:mesh.nnx],:)=R(ww+[1:mesh.nnx],:)/1e15;
        
        %'% aumento la sensibilità delle equazioni di potenza ottica        '
        
        %        R(ss+[1:NUMERO_MODI],:)=R(ss+[1:NUMERO_MODI],:)*1e5;
        
        
        res = norm(R*rvet);
        Res(iter+1)=res;
        %    if full(mode.Scheck)>1
        %        disp('dopo calcolo residuo')
        %        keyboard
        %        end
        
        %==============================================================================================100
        if(isnan(res))
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
        % % % %      if(res<mode.tolconv*1e4)
        % % % %       if res>Res(iter)
        % % % %       mode.iconv(indv)=1; mode.res(indv) = res; break,
        % % % %      end
        % % % %      end
        
        % uvet2 = uvet; % save solution, uncomment for Bank-Rose
        iter_t=0; t=t0_Newton; % restore damping parameter % QUAAAAAA
        %=========================================
        % solve matrix equation
        Jmat = R*Jmat*C;
        [Lmat,Umat,pm,qm] = lu(Jmat);
        deltau = - C*(qm*(Umat\(Lmat\(pm*(R*rvet)))));
        iter = iter + 1;
        uvet = uvet + t*deltau.'; % update
        %============================================
        % else % uncomment for Bank-Rose
        %     t=t/2; % uncomment for Bank-Rose
        %     iter_t = iter_t + 1; % uncomment for Bank-Rose
        %     uvet = uvet2 + t*deltau.'; % update, uncomment for Bank-Rose
        % end % uncomment for Bank-Rose
        tmpVar
    end % outer while-loop
    if(not(mode.iconv(indv)) & not(isnan(res))),
        disp('================================================')
        disp('=        Convergence failure!!!! :-(           =')
        disp('================================================')
        
        'parte cancellata'
        Scheck=mode.gmod./mode.lmod';
        indWrongScheck=find(Scheck>=1);
        if length(indWrongScheck)>0
            disp('================================================')
            disp('Too much code breaks due to Scheck > 1 condition')
            disp(['Restarting code with Pst=',num2str(mode.ScheckMultiplicationFactor),'*Pst'])
            disp('================================================')
            
            uvet(ss+[indWrongScheck])=mode.ScheckMultiplicationFactor*uvet(ss+[indWrongScheck]);
            Pst=uvet(ss+[indWrongScheck]);
        end
        'fine parte cancellata'
        %keyboard
    end
    if ((not(mode.iconv(indv)) | isnan(res)) & ContaBisez<ContaBisezMax*length(mode.ScheckMultiplicationFactor))
        ContaBisez=ContaBisez+1
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
        if(indCountMultiplicationFactor==1)
            Vb1=mode.v0_dd(1:indv);
            VbNew=(mode.v0_dd(indv)+mode.v0_dd(indv+1))/2;
            Vb2=mode.v0_dd(indv+1:end);
            % aggiunge tensione
            mode.v0_dd=[Vb1,VbNew,Vb2];
            NVbias=NVbias+1;
        elseif(indCountMultiplicationFactor>=length(mode.ScheckMultiplicationFactor))
            indCountMultiplicationFactor=1;
            Vb1=mode.v0_dd(1:indv);
            VbNew=(mode.v0_dd(indv)+mode.v0_dd(indv+1))/2;
            Vb2=mode.v0_dd(indv+1:end);
            % aggiunge tensione
            mode.v0_dd=[Vb1,VbNew,Vb2];
            NVbias=NVbias+1;
        else
            disp('che gli facciamo fare qui ???')
            
            
        end
    else
        flagConv=1;
        flagCallVELM=1;
        ContaBisez=0;
    end
    if (ContaBisez>=ContaBisezMax)
        'ContaBisez',break
    end
    %
    % =====================================================================
    % save DD solution
    % =====================================================================
    if(flagConv==1) % if residual is not NaN
        %         mode.tolconv=mode.tolconv*2;
        mode.phi=uvet(1:nn);
        iq = mesh.iq;
        mode.elec=abs(uvet(nn+(1:nn)));
        if((isfield(mode,'stats'))&&(strcmp(mode.stats,'Fermi')))
            tmp = mode.elec(iq)./mesh.Nc(iq);
            nferinv = invferdr(tmp,mode.tolinvferdr);
            mode.EFn(iq) = mode.ecb(iq) + Vt(iq).*nferinv; end
        %
        mode.hole=abs(uvet(pp+(1:nn)));
        if((isfield(mode,'stats'))&&(strcmp(mode.stats,'Fermi'))),
            tmp = mode.hole(iq)./mesh.Nv(iq);
            nferinv = invferdr(tmp,mode.tolinvferdr);
            mode.EFp(iq) = mode.evb(iq) - Vt(iq).*nferinv; end
        %
        fprintf('Elapsed time: %g, min\n',(cputime-elapsed_time)/60)
        mode.elapsed_time=(cputime-elapsed_time)/60;
        
        v_dd=uvet(qq+1); % saving voltage
        i_dd=uvet(rr+1)*mode.CarrierNorm; % saving current
        
        if mode.oflg==1 && mode.quasi1D==1
            PPst=mode.Pst+mode.Psp;   % mW/cm2
            if indv>1
                PElec=i_dd*v_dd*1e3; %mW/cm2
                %             PElec=mode.I_dd*mode.Vbias*1e3; %mW/cm2
                PDiss=PElec-PPst;
            else
                PElec=0;
                PPst=0;
                PDiss=0;
            end
            
            mode.PDiss(indv)=PDiss;
            mode.PDissPred=PDiss;
            
            
            if (mode.Pst>0.01)
                xstart=-2;
                if indv<=length(mode.v0_dd)
                    xVolt=mode.v0_dd(indv+[xstart:0]);
                    yPdiss=mode.PDiss(indv+[xstart:0]);
                    
                    coeffPdiss=polyfit(xVolt,log10(yPdiss),1);
                    if indv<length(mode.v0_dd)
                        Pdissfit=10.^polyval(coeffPdiss,mode.v0_dd(indv+1));
                    else
                        Vbias=mode.v0_dd(end)+diff(mode.v0_dd([end-1,end]));
                        Pdissfit=10.^polyval(coeffPdiss,Vbias);
                    end
                    
                    mode.PDissPred(indv)=Pdissfit;
                elseif (indv>length(mode.v0_dd) && isfield(mode,'Ibias'))
                    if  indv<length(mode.i0_dd)
                        xCurr=mode.i0_dd(indv+[xstart:0]);
                        yPdiss=mode.PDiss(indv+[xstart:0]);
                        
                        coeffPdiss=polyfit(xCurr,log10(yPdiss),1);
                        
                        Pdissfit=10.^polyval(coeffPdiss,mode.i0_dd(indv+1));
                        
                        mode.PDissPred(indv)=Pdissfit;
                    end
                end
            end
        end
        
        if(mode.oflg==0)
            mode.ii_dd(indv)=i_dd;
            mode.vv_dd(indv)=v_dd;
            mode.vv0_dd(indv)=v0_dd;
            % compute electric field (in time domain), V/um
            efieldx_t=-mode.phi*mesh.gradx*1e-4; efieldx_t(it)=0;
            efieldy_t=-mode.phi*mesh.grady*1e-4; efieldy_t(it)=0;
            mode.efield_x=reshape(pdeprtni(mesh.node,mesh.triangle(1:4,:),efieldx_t),1,mesh.nn);
            mode.efield_y=reshape(pdeprtni(mesh.node,mesh.triangle(1:4,:),efieldy_t),1,mesh.nn);
            if(mode.Tflg==1)
                mode=f_EvalHeatingTerms(geom,mesh,mode);
            end
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
                flagConv=0;
                % Vb1=mode.v0_dd(1:indv);
                % VbNew=(mode.v0_dd(indv)+mode.v0_dd(indv+1))/2;
                % VbNew=v0_dd;
                % Vb2=mode.v0_dd(indv+1:end);
                % aggiunge tensione
                % mode.v0_dd=[Vb1,VbNew,Vb2];
                % NVbias=NVbias+1;
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
                    disp(['Restarting code with Pst=',num2str(mode.ScheckMultiplicationFactor(indCountMultiplicationFactor)),'*Pst'])
                    disp('================================================')
                    uvet(ss+[indWrongScheck])=mode.ScheckMultiplicationFactor(indCountMultiplicationFactor)*uvet(ss+[indWrongScheck]);
                    Pst=uvet(ss+[indWrongScheck]);
                    indCountMultiplicationFactor=indCountMultiplicationFactor+1;
                    %                     keyboard
                end
            else
                countScheck=0;
                flagCallVELM=1;
                mode.vv0_dd(indv)=v0_dd;
                mode.vv_dd(indv)=v_dd;
                mode.ii_dd(indv)=i_dd;
                I_mA=i_dd*1000;
                Ptotal=sum(Pst(:));
                PtotV(indv)=Ptotal;
                if indv>10
                    diPot=diff(PtotV);
                    if diPot(end)<0
                        if Ptotal<PotMin & max(PtotV)>PotMin
                            CondPotBreak=1;
                        end
                    end
                end
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
                mode.IntCcapN(indv,:)=mode.IntCcapn;
                mode.IntCcapP(indv,:)=mode.IntCcapp;
                mode.IntRecN(indv,:)=mode.IntRecN;
                mode.IntRecP(indv,:)=mode.IntRecP;
                %               mode.IntCcapN(indv,:)=mode.IntCcapn;
                %               mode.IntCcapP(indv,:)=mode.IntCcapp;
                mode.IntREc(indv,:)=mode.IntRec;
                
                
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
    
    if IPLOT==1 && mode.oflg==1
        figure(portfig)
        len=1:length(mode.ii_dd);
        plot(1000*mode.ii_dd,1e-12*mode.nMaxVet(len)/mesh.NMQW,'o-',...
            1000*mode.ii_dd,1e-12*mode.pMaxVet(len)/mesh.NMQW,'+-')
        hold on
        plot(1000*mode.ii_dd,1e-18*mode.n3MaxVet(len),'--',...
            1000*mode.ii_dd,1e-18*mode.p3MaxVet(len),'--')
        xlabel('Current (mA)')
        ylabel('N-P density (PER well) 1e12/cm^2)')
    end
    %            figure(portfig)
    %            plot(1000*mode.ii_dd,1e-12*mode.nMaxVet/3,'o-',...
    %            1000*mode.ii_dd,1e-12*mode.pMaxVet/3,'+-')
    %            xlabel('Current (mA)')
    %            ylabel('N-P density 1e12/cm^2)')
    %       figure
    %               plot(1000*mode.ii_dd,1e-12*mode.nMaxVet(1:end-1),'o-',...
    %           1000*mode.ii_dd,1e-12*mode.pMaxVet(1:end-1),'+-')
    
    if(mode.oflg)
        vind=[];
        for indplot=1:length(VELMInfo)
            vind=[vind,VELMInfo(indplot).indVoltage];
        end
        if IPLOT==1
            figure(1234+kfig),clf
            set(gcf,'Position',[ 58   547   600   376])
            subplot(1,2,1)
            hold on
            grid on
            box on
            plot(Vmeas,Imeas,'r','LineWidth',2)
            plot(mode.vv0_dd,mode.ii_dd*1000,'b.','LineWidth',2)
            plot(mode.vv0_dd(vind),mode.ii_dd(vind)*1000,'bo','LineWidth',2)
            %        axis([0 mode.vv0_dd(vind(end))+.2 0 mode.ii_dd(vind(end))+1])
            axis([0 mode.vv0_dd(end)+.2 0 1000*mode.ii_dd(end)+1])
            xlabel('Voltage, V')
            ylabel('Current, mA')
            title(['Vcurrent = ',num2str(v0_dd),' V'])
            %        legend('Measurements','Simulation','Location','Best')
            
            subplot(1,2,2)
            hold on
            grid on
            box on
            plot(Imeas,Lmeas,'r','LineWidth',2)
            PPst=sum(mode.Pst_dd,1)+mode.Psp_dd;
            
            plot(mode.ii_dd*1000,mode.Pst_dd,'.','LineWidth',2)
            %        axis([0 mode.ii_dd(vind(end))+1 0 mode.Pst_dd(vind(end))+.5])
            axis([0 1000*mode.ii_dd(end)+1 0 max(sum(mode.Pst_dd,1))*1.1+.1])
            %        legend('Measurements','Simulation','Location','Best')
            plot(mode.ii_dd(vind)*1000,PPst(vind),'bo','LineWidth',2)
            
            
            xlabel('Current, mA')
            ylabel('Optical power, mW')
            title(['Voltage ',num2str(indv),' of ',num2str(NVbias)])
            drawnow
            
            figure(1235+kfig),clf
            set(gcf,'Position',[  634    73   579   365])
            subplot(1,2,1)
            hold on
            grid on
            box on
            plot(mode.vv0_dd,mode.DeltaTmax,'k','LineWidth',2)
            plot(mode.vv0_dd,mode.DeltaTmax_Rec,'r--','LineWidth',2)
            plot(mode.vv0_dd,mode.DeltaTmax_Joule,'c--','LineWidth',2)
            plot(mode.vv0_dd,mode.DeltaTmax_OptAbs,'g--','LineWidth',2)
            plot(mode.vv0_dd(vind),mode.DeltaTmax(vind),'ko','LineWidth',2)
            xlabel('Voltage, V')
            ylabel('Temperature rise, K')
            %        legend('Total','Non-rad rec.','Joule','Opt. abs.','Location','Best')
            subplot(1,2,2)
            hold on
            grid on
            box on
            plot(mode.vv0_dd,mode.Gmod,'k','LineWidth',2)
            plot(mode.vv0_dd,mode.Lmod,'r--','LineWidth',2)
            if(mode.nmodes>1)
                plot(mode.vv0_dd(vind),mode.Gmod(:,vind)','ko','LineWidth',2)
                plot(mode.vv0_dd(vind),mode.Lmod(:,vind)','ro','LineWidth',2)
            else
                plot(mode.vv0_dd(vind),mode.Gmod(vind),'ko','LineWidth',2)
                plot(mode.vv0_dd(vind),mode.Lmod(vind),'ro','LineWidth',2)
            end
            xlabel('Voltage, V')
            ylabel('Modal gain vs losses, cm^{-1}')
            
            gipos=find(mode.Gmod(1,:)>0);
            if length(gipos)>1
                xli=mode.vv0_dd(gipos([1 end]) ) +[-.1,.1];
                axis([xli 0 max(max(mode.Lmod))] )
            end
            %        legend('Gain','Losses','Location','Best')
            %keyboard
            drawnow
            
            
            if(indv>1)
                figure(1236+kfig),clf
                set(gcf,'Position',[631   547   585   365])
                subplot(1,2,1)
                hold on
                grid on
                box on
                Resistence=diff(mode.vv0_dd)./diff(mode.ii_dd);
                Resistence=[Resistence(1),Resistence];
                % plot(mode.ii_dd*1000,Resistence,Imeasinterp,Resmeasinterp)
                plot(Imeasinterp,Resmeasinterp,'r-','LineWidth',2)
                plot(Imeas_res,Rmeas,'.')
                lres=1:length(Resistence);
                plot(mode.ii_dd(lres)*1000,Resistence(lres),'go','LineWidth',2)
                ylim([40,120])
                %            xlim([0.5,1000*mode.ii_dd(end)+.1])
                xlabel('Current, mA')
                ylabel('Differential resistence, \Omega')
                %            legend('Measurements','Simulation','Location','Best')
                subplot(1,2,2)
                hold on
                grid on
                box on
                plot(Cur,LAM,'ro-','LineWidth',2)
                plot(mode.ii_dd*1000,mode.lambda,'b.','LineWidth',2)
                %            legend('Measurements','Simulation','Location','Best')
                if(mode.nmodes>1)
                    plot(mode.ii_dd(vind)*1000,mode.lambda(:,vind),'bo','LineWidth',2)
                else
                    plot(mode.ii_dd(vind)*1000,mode.lambda(vind),'bo','LineWidth',2)
                end
                %            xlim([0.,1000*mode.ii_dd(end)+.1])
                xlabel('Current, mA')
                ylabel('Modal wavelength, nm')
                drawnow
            end
        end  %IPLOT
    else
        if IPLOT==1
            figure(1234+kfig),clf
            set(gcf,'Position',[268 533 1096 420])
            %             subplot(1,2,1)
            hold on
            grid on
            box on
            %             plot(Vmeas,Imeas,'r','LineWidth',2)
            plot(mode.vv0_dd,mode.ii_dd*1000,'b.','LineWidth',2)
            xlabel('Voltage, V')
            ylabel('Current, mA')
            title(['Vcurrent = ',num2str(v0_dd),' V'])
            
            keyboard
            %        legend('Measurements','Simulation','Location','Best')
            %             subplot(1,2,2)
            %             if(indv>1)
            %                 hold on
            %                 grid on
            %                 box on
            %                 Resistence=diff(mode.vv0_dd)./diff(mode.ii_dd);
            %                 Resistence=[Resistence(1),Resistence];
            %                 plot(Imeasinterp,Resmeasinterp,'r-','LineWidth',2)
            %                 plot(mode.ii_dd*1000,Resistence,'b.','LineWidth',2)
            %                 ylim([75,120])
            %                 xlabel('Current, mA')
            %                 ylabel('Differential resistence, \Omega')
            %                 %            legend('Measurements','Simulation','Location','Best')
            %                 drawnow
            %             end
        end %IPLOT
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
    
    
    
    %    save('Last_Workspace.mat')
    %        'verifica save', keyboard
    if(flagConv)
        %        eval(['save ',Last_Workspace])
        savePlotD
        %        'verifica save', keyboard
        eval(['save ',Last_Workspace,'_Plot',rad_setting,' modePlot'])
        MODEplot{kpar}=modePlot;
        eval(['save ',nomeSave,num2str(IPAR),'.mat MODEplot'])
        % save('Last_plot.mat','modePlot')
        if exist('V0_save')
            if v0_dd==V0_save
                save Riparto
            end
        end
    end
    %    if v0_dd==1.70
    %     'FERMA per CONFRONTO',
    %     keyboard
    %    end
    %     if input(' Plottare quantità e confrontare con D1ANA? ENTER: NO, 1: Sì ') == 1
    %         Plot_intermediate
    %         Plot_Cmp2D1ANA
    %
    %     end
    
    indv=indv+1;
    
    
    %  'prima save dentro save', keyboard
    %eval(['save ',nomeSW,'Contributi_',nomeSav,' ',' MODEplot -append '])
    
    %  save Contributi MODEplot -append
    
end


%% CURRENT DRIVING SECTION
if mode.Idrive==1
    mode.IdriveON=1;
    Iadd=[mode.ii_dd(end)+mode.Istep*1e-3:mode.Istep*1e-3:Imassimo*1e-3];
    mode.i0_dd=[mode.ii_dd Iadd];
    NVbiasVdrive=NVbias;
    NVbias=length(mode.i0_dd);
    
    
    while(indv<=NVbias && DeltaTmax<DTM0 && I_mA<Imassimo && CondPotBreak==0 ) % outer bias-loop %$$$$$$$$$$$$$$$$$$$$$$$$
        
        elapsed_time=cputime;
        %    mode_old=mode;
        %    mesh_old=mesh;
        uvet_old=uvet;
        %     v0_dd = mode.v0_dd(indv)
        i0_dd = mode.i0_dd(indv)
        
        if(mode.Tflg && v_dd>mode.minthermalvoltage && flagCallVELM) % OCCHIO!!! flagCallVELM
            Tprec=mesh.Tprec;
            % Recombination-only heating
            mode.iplotTerm=0;
            PPst=sum(modePlot.Pst_dd,1)+modePlot.Psp_dd;
            PElec=modePlot.ii_dd*1000.*modePlot.vv_dd;
            PDiss=PElec(end)-PPst(end);
            mode.PDiss(indv)=PDiss;
            mode.iplotTerm=1;
            %        mode.iTfig=1;
            %           'DeltaT', keyboard
            if isfield(mode,'IsoThermal')
                IsoThermal=mode.IsoThermal;
            else
                IsoThermal=0;
            end
            if IsoThermal==0
                if isfield(mode,'quasi1D') && mode.quasi1D==1
                    [DeltaTold,PTherm,T_Contributi]=f_ThermD1ANA(mesh,mode,StrTT);
                    Tprec=0;
                    
                    StimaTempWU
                    DeltaT=DeltaTold*Fat_STIMA_Temp;
                    
                    mode.FatMob=polyval(cot,max(max(DeltaT)));
                else
                    [DeltaTold,Tprec,PTherm,T_Contributi,K,TotalHeat]=f_ThermicFun(Tprec,mesh,mode,StrTT,IPLOT);
                    mode.FatMob=1;
                    
                    StimaTempWU
                    DeltaT=DeltaTold*Fat_STIMA_Temp;
                end
            else
                if mode.quasi1D==1
                    mesh.DeltaT=0;
                else
                    DeltaT=zeros(1,mesh.nn);
                end
                PTherm=0;
                for tcj=1:5
                    T_Contributi{tcj}=0;
                end
            end
            
            %[DeltaT,Tprec,PTherm,T_Contributi]=f_ThermicFun(Tprec,mesh,mode,StrTT,IPLOT);
            DeltaTJoule=T_Contributi{1};
            DeltaTRec_srhAu=T_Contributi{2};
            DeltaTRec_Ccap=T_Contributi{3};
            DeltaTRec_RAD=T_Contributi{4};
            DeltaTOptAbs=T_Contributi{5};
            DeltaT=C_Temp*reshape(DeltaT,1,mesh.nn);
        elseif(flagCallVELM==1)
            DeltaT=zeros(1,mesh.nn);
            DeltaTJoule=DeltaT; DeltaTOptAbs=DeltaT;
            DeltaTRec_srhAu=DeltaT; DeltaTRec_Ccap=DeltaT;
            DeltaTRec_RAD=DeltaT;
            Tprec=0; % initalization of Tprec=0 for thermic self-consistency
            
        end
        if isfield(mode,'epsNLg')
            nlG=mode.epsNLg;
            mode.nlG=1;
            if isfield(mode,'E2')
                %       'qui per nlG', keyboard
                UU = (mode.Pst.*mode.fPdif)*mode.E2;
                uudu=1./(1+mode.epsNLg*UU);
                mode.nlG=uudu(1:mesh.nnxQW{1});
            end
        else
            mode.nlG=1;
        end
        
        
        mesh.DeltaT=DeltaT;
        
        
        if effetti(6)==0
            mesh.DeltaTvelm=DeltaT;
        else
            'TEMPERATURA Fissa in VELM'
        end
        
        DeltaTmax=max(DeltaT);
        if DeltaTmax-mode.DeltaTmax(end)>DTM0/4
            'incremento di temperatura troppo grande, esco '
            break
        end
        DeltaTmax_Joule=max(max(DeltaTJoule));
        DeltaTmax_srhAu=max(max(DeltaTRec_srhAu));
        DeltaTmax_Ccap=max(max(DeltaTRec_Ccap));
        DeltaTmax_RAD=max(max(DeltaTRec_RAD));
        DeltaTmax_OptAbs=max(max(DeltaTOptAbs));
        
        % giusto per non buttare via quel plot, faccio così.
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
            
            if iold==2
                % calcolo guadagno e derivate
                %    'qui derivate', keyboard
                NUOVA_Gain
                % carico dati tau di Vallone
                global CO coMh Tn Th
                if mode.taucarrier==1
                    load saTaue
                end
            end
        end
        
        load stopDD
        if istop==1
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
            
            flagPot=0;
            FLAG_velm0=( indv==1 | (flagCallVELM & flagTemp) | iVELM | flagPot) & HARD;
            flagvv=0;
            if indv>1
                if length(VELMInfo)>0
                    vlast=mode.i0_dd(VELMInfo(indVELM).indVoltage);
                    flagvv=vlast>=v_dd;
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
                
                
                if effetti(7)==1
                    VelmOptions.ianti_gui=0;  % tolgo anti_guiding
                end
                colordef black
                %                        'prima di Velm', keyboard
                iLoadVELM=0; % 1, check if VELM results are already computed; 0, always enters in VELM
                %
                if (iLoadVELM==1 && exist(['dati/VELM/velm_',strSave,'_quasi1D=',num2str(mode.quasi1D),'.mat'])==2 && v0_dd==0)
                    % File exists.
                    load(['velm_',strSave,'_quasi1D=',num2str(mode.quasi1D),'.mat'])
                else
                    % File does not exist.
                    if isfield(mode,'quasi1D') && mode.quasi1D==1
                        [velm] = f_CallVELM_1D(mesh,mode,mode1,ParVet,VelmOptions,fil_str,mode.quasi1D);
                    else
                        [velm] = f_CallVELM(mesh,mode,mode1,ParVet,VelmOptions,fil_str);
                    end
                    % stores VELM output to re-use it for successive simulation
                    if v0_dd==0
                        save(['dati/VELM/velm_',strSave,'_quasi1D=',num2str(mode.quasi1D)],'velm')
                    end
                end
                %            clear global
                %            mode.verbVELM=0;
                colordef white
                if isfield(mode,'quasi1D') && mode.quasi1D==1
                    velm.E2=ones(mode.nmodes,mesh.nnx)/(mode.AreaOx);
                    velm.SW=0;
                end
                mode.Lm=velm.Lm(1:mode.nmodes); % prendo solo il primo modo (per ora)
                mode.vlambda=velm.vlambda;
                mode.NQW=velm.NQW; % prendo solo il primo modo (per ora)
                mode.Gamma_z=velm.Gamma_z;
                mode.nindexQW=velm.nindexQW;
                mode.fPES=velm.fPES*1000;
                vph=Clight./mode.nindexQW;
                mode.fPdif=TAROCCO*velm.fPdif/vph/1000;
                mode.E2=velm.E2; % original
                %            'vagawt['
                %            keyboard
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
                %'salve lambda', keyboard
                mode1.TmVelm(indVELM)=DeltaTmax;
                
                mode1.NmVelm(indVELM)=nMax;
                mode1.PmVelm(indVELM)=pMax;
                mode1.IndVelm(indVELM)=indv;
                mode1.LmVelm(indVELM,:)=velm.Lm(1);
                mode1.E2Velm(indVELM,:)=velm.E2(1,:)';
                npMaxLast=npMax;
                %	                        'salva per Velm', keyboard
                salva_per_VELM
                
                flagVELMWasCalled=1;
                
                colordef white
            elseif (isfield(mode,'quasi1D')==1 && mode.quasi1D==1) && indVELM>=4
                fitStart=3; fitDegree=2;
                [mode]=velmFitting(indVELM,fitStart,fitDegree,VELMInfo,mode,i0_dd);
                %             'Interp',keyboard
            end
        end
        
        T=mode.T0+DeltaT;
        mesh.Tgain=T;
        T=mode.T0+DeltaT*mode.C_Temp_DD;
        
        TKT=T;
        if isfield(mode,'KT')
            if mode.KT==0
                TKT=mode.T0*ones(size(T));
            end
        end
        
        Vt=kB.*TKT./qel;
        
        mesh=loadmacro(geom,mesh,mode,T);
        
        mode.Vt=Vt;
        Vt_tr=pdeintrp(mesh.node,mesh.triangle(1:4,:),Vt.'); % T on triangles
        mode.Vt_tr=Vt_tr;
        
        mode_old=mode;
        mesh_old=mesh;
        
        mode=f_ComputeNeutrality(geom,mesh,mode);
        
        if(mode.Tflg==1)
            %        mode.DeltaTmax(indv)=DeltaTmax;
            mesh.Tprec=Tprec;
        end
        
        % calcolo carrier-dependent tau Vallone
        if(mode.taucarrier==1)
            indQW=1; % primo QW
            n2D=mode.n2D{indQW}/mesh.WMQW{indQW};
            p2D=mode.p2D{indQW}/mesh.WMQW{indQW};
            n3D=mode.elec(mesh.inMQW{indQW});
            p3D=mode.hole(mesh.inMQW{indQW});
            if(max(n2D) > 1e8)
                mode.Fat_cap_e=fit_tauEnl(n2D,n3D);
                mode.Fat_cap_h=fit_tauH(p2D,p3D);
                %            'VALLONE', keyboard
            end
        end
        
        
        iterMag1=0;
        t0_Newton=1;
        
        fprintf('\n Drift-diffusion analysis\n')
        mode.iconv(indv) = 0; % flag, 1 convergence, 2 convergence failure
        iter = 0; % number of outer iterations
        iter_t = 0; % number of inner iterations
        mode.res(indv) = Inf; % residual error
        % t = 1; % damping parameter % QUAAAAAA
        %
        
        %             'keyboard le_VEl', keyboard
        
        if effetti(1)==1
            mode.fPES=VELMInfo(le_VELMinfo).fPES*1e3;
        end
        if effetti(2)==1
            mode.Gamma_z=VELMInfo(le_VELMinfo).Gamma_z;
        end
        if effetti(3)==1
            %            mode.fPdif=TAROCCO*velm.fPdif/vph/1000;
            mode.fPdif=TAROCCO*VELMInfo(le_VELMinfo).fPdif/vph/1e3;
        end
        if effetti(4)==1
            mode.E2=VELMInfo(le_VELMinfo).E2;
        end
        if effetti(5)==1
            mode.vlambda=VELMInfo(le_VELMinfo).vlambda;
        end
        'fPES Gamma_z  fPdif E2 lam'
        
        clear Res
        
        
        while(iter<mode.maxiter && iter_t<mode.maxiter_t), % inner while-loop
            
            indCountMultiplicationFactor=1; %%% temporaneo
            
            
            % assem jacobian matrix in time domain
            mode.ind_v0=indv;
            
            % prima di assem
            if iDyn==0
                [Kmat0,Jmat0,Jmat2,uvet,rvet,mode]=assem_GBT(geom,mesh,mode,uvet,i0_dd);
            else
                [Kmat0,Jmat0,Jmat1,Jmat2,uvet,rvet,mode]=assem_FM(geom,mesh,mode,uvet,i0_dd);
            end
            
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
            
            res = norm(R*rvet);
            Res(iter+1)=res;
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
            % % % %      if(res<mode.tolconv*1e4)
            % % % %       if res>Res(iter)
            % % % %       mode.iconv(indv)=1; mode.res(indv) = res; break,
            % % % %      end
            % % % %      end
            
            % uvet2 = uvet; % save solution, uncomment for Bank-Rose
            iter_t=0; t=t0_Newton; % restore damping parameter % QUAAAAAA
            %=========================================
            % solve matrix equation
            Jmat = R*Jmat*C;
            [Lmat,Umat,pm,qm] = lu(Jmat);
            deltau = - C*(qm*(Umat\(Lmat\(pm*(R*rvet)))));
            iter = iter + 1;
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
            
            'parte cancellata'
            Scheck=mode.gmod./mode.lmod';
            indWrongScheck=find(Scheck>=1);
            if length(indWrongScheck)>0
                disp('================================================')
                disp('Too much code breaks due to Scheck > 1 condition')
                disp(['Restarting code with Pst=',num2str(mode.ScheckMultiplicationFactor),'*Pst'])
                disp('================================================')
                
                uvet(ss+[indWrongScheck])=mode.ScheckMultiplicationFactor*uvet(ss+[indWrongScheck]);
                Pst=uvet(ss+[indWrongScheck]);
            end
            'fine parte cancellata'
            %keyboard
        end
        if ((not(mode.iconv(indv)) | isnan(res)) & ContaBisez<ContaBisezMax*length(mode.ScheckMultiplicationFactor))
            ContaBisez=ContaBisez+1
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
            if(indCountMultiplicationFactor==1)
                Ib1=mode.i0_dd(1:indv);
                IbNew=(mode.i0_dd(indv)+mode.i0_dd(indv+1))/2;
                Ib2=mode.i0_dd(indv+1:end);
                % aggiunge tensione
                mode.i0_dd=[Ib1,IbNew,Ib2];
                NVbias=NVbias+1;
            elseif(indCountMultiplicationFactor>=length(mode.ScheckMultiplicationFactor))
                indCountMultiplicationFactor=1;
                Ib1=mode.i0_dd(1:indv);
                IbNew=(mode.i0_dd(indv)+mode.i0_dd(indv+1))/2;
                Ib2=mode.v0_dd(indv+1:end);
                % aggiunge tensione
                mode.i0_dd=[Ib1,IbNew,Ib2];
                NVbias=NVbias+1;
            else
                disp('che gli facciamo fare qui ???')
                
                
            end
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
            %         mode.tolconv=mode.tolconv*2;
            mode.phi=uvet(1:nn);
            iq = mesh.iq;
            mode.elec=abs(uvet(nn+(1:nn)));
            if((isfield(mode,'stats'))&&(strcmp(mode.stats,'Fermi')))
                tmp = mode.elec(iq)./mesh.Nc(iq);
                nferinv = invferdr(tmp,mode.tolinvferdr);
                mode.EFn(iq) = mode.ecb(iq) + Vt(iq).*nferinv; end
            %
            mode.hole=abs(uvet(pp+(1:nn)));
            if((isfield(mode,'stats'))&&(strcmp(mode.stats,'Fermi'))),
                tmp = mode.hole(iq)./mesh.Nv(iq);
                nferinv = invferdr(tmp,mode.tolinvferdr);
                mode.EFp(iq) = mode.evb(iq) - Vt(iq).*nferinv; end
            %
            fprintf('Elapsed time: %g, min\n',(cputime-elapsed_time)/60)
            mode.elapsed_time=(cputime-elapsed_time)/60;
            
            v_dd=uvet(qq+1); % saving voltage
            i_dd=uvet(rr+1); % saving current
            
            if mode.oflg==1 && mode.quasi1D==1
                PPst=mode.Pst+mode.Psp;   % mW/cm2
                if indv>1
                    PElec=i_dd*v_dd*1e3; %mW/cm2
                    %             PElec=mode.I_dd*mode.Vbias*1e3; %mW/cm2
                    PDiss=PElec-PPst;
                else
                    PElec=0;
                    PPst=0;
                    PDiss=0;
                end
                
                mode.PDiss(indv)=PDiss;
                mode.PDissPred=PDiss;
                
                
                if (mode.Pst>0.01)
                    xstart=-2;
                    if indv<=length(mode.v0_dd)
                        xVolt=mode.v0_dd(indv+[xstart:0]);
                        yPdiss=mode.PDiss(indv+[xstart:0]);
                        
                        coeffPdiss=polyfit(xVolt,log10(yPdiss),1);
                        if indv<length(mode.v0_dd)
                            Pdissfit=10.^polyval(coeffPdiss,mode.v0_dd(indv+1));
                        else
                            Vbias=mode.v0_dd(end)+diff(mode.v0_dd([end-1,end]));
                            Pdissfit=10.^polyval(coeffPdiss,Vbias);
                        end
                        
                        mode.PDissPred(indv)=Pdissfit;
                    elseif (indv>length(mode.v0_dd) && isfield(mode,'Ibias'))
                        if  indv<length(mode.i0_dd)
                            xCurr=mode.i0_dd(indv+[xstart:0]);
                            yPdiss=mode.PDiss(indv+[xstart:0]);
                            
                            coeffPdiss=polyfit(xCurr,log10(yPdiss),1);
                            
                            Pdissfit=10.^polyval(coeffPdiss,mode.i0_dd(indv+1));
                            
                            mode.PDissPred(indv)=Pdissfit;
                        end
                    end
                end
            end
            
            if(mode.oflg==0)
                mode.ii_dd(indv)=i_dd;
                mode.vv_dd(indv)=v_dd;
                mode.vv0_dd(indv)=v0_dd;
                % compute electric field (in time domain), V/um
                efieldx_t=-mode.phi*mesh.gradx*1e-4; efieldx_t(it)=0;
                efieldy_t=-mode.phi*mesh.grady*1e-4; efieldy_t(it)=0;
                mode.efield_x=reshape(pdeprtni(mesh.node,mesh.triangle(1:4,:),efieldx_t),1,mesh.nn);
                mode.efield_y=reshape(pdeprtni(mesh.node,mesh.triangle(1:4,:),efieldy_t),1,mesh.nn);
                if(mode.Tflg==1)
                    mode=f_EvalHeatingTerms(geom,mesh,mode);
                end
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
                    flagConv=0;
                    flagCallVELM=0;
                    
                    countScheck=countScheck+1;
                    if(countScheck>mode.maxScheckRepeat)
                        indWrongScheck=find(Scheck>=1);
                        disp('================================================')
                        disp('Too much code breaks due to Scheck > 1 condition')
                        disp(['Restarting code with Pst=',num2str(mode.ScheckMultiplicationFactor(indCountMultiplicationFactor)),'*Pst'])
                        disp('================================================')
                        uvet(ss+[indWrongScheck])=mode.ScheckMultiplicationFactor(indCountMultiplicationFactor)*uvet(ss+[indWrongScheck]);
                        Pst=uvet(ss+[indWrongScheck]);
                        indCountMultiplicationFactor=indCountMultiplicationFactor+1;
                        %                     keyboard
                    end
                else
                    countScheck=0;
                    flagCallVELM=1;
                    mode.vv0_dd(indv)=v0_dd;
                    mode.vv_dd(indv)=v_dd;
                    mode.ii_dd(indv)=i_dd;
                    I_mA=i_dd*1000;
                    Ptotal=sum(Pst(:));
                    PtotV(indv)=Ptotal;
                    if indv>10
                        diPot=diff(PtotV);
                        if diPot(end)<0
                            if Ptotal<PotMin & max(PtotV)>PotMin
                                CondPotBreak=1;
                            end
                        end
                    end
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
                    mode.IntCcapN(indv,:)=mode.IntCcapn;
                    mode.IntCcapP(indv,:)=mode.IntCcapp;
                    mode.IntRecN(indv,:)=mode.IntRecN;
                    mode.IntRecP(indv,:)=mode.IntRecP;
                    %               mode.IntCcapN(indv,:)=mode.IntCcapn;
                    %               mode.IntCcapP(indv,:)=mode.IntCcapp;
                    mode.IntREc(indv,:)=mode.IntRec;
                    
                    
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
        
        if IPLOT==1 && mode.oflg==1
            figure(portfig)
            len=1:length(mode.ii_dd);
            plot(1000*mode.ii_dd,1e-12*mode.nMaxVet(len)/mesh.NMQW,'o-',...
                1000*mode.ii_dd,1e-12*mode.pMaxVet(len)/mesh.NMQW,'+-')
            hold on
            plot(1000*mode.ii_dd,1e-18*mode.n3MaxVet(len),'--',...
                1000*mode.ii_dd,1e-18*mode.p3MaxVet(len),'--')
            xlabel('Current (mA)')
            ylabel('N-P density (PER well) 1e12/cm^2)')
        end
        %            figure(portfig)
        %            plot(1000*mode.ii_dd,1e-12*mode.nMaxVet/3,'o-',...
        %            1000*mode.ii_dd,1e-12*mode.pMaxVet/3,'+-')
        %            xlabel('Current (mA)')
        %            ylabel('N-P density 1e12/cm^2)')
        %       figure
        %               plot(1000*mode.ii_dd,1e-12*mode.nMaxVet(1:end-1),'o-',...
        %           1000*mode.ii_dd,1e-12*mode.pMaxVet(1:end-1),'+-')
        
        if(mode.oflg)
            vind=[];
            for indplot=1:length(VELMInfo)
                vind=[vind,VELMInfo(indplot).indVoltage];
            end
            if IPLOT==1
                figure(1234+kfig),clf
                set(gcf,'Position',[ 58   547   600   376])
                subplot(1,2,1)
                hold on
                grid on
                box on
                plot(Vmeas,Imeas,'r','LineWidth',2)
                plot(mode.vv0_dd,mode.ii_dd*1000,'b.','LineWidth',2)
                plot(mode.vv0_dd(vind),mode.ii_dd(vind)*1000,'bo','LineWidth',2)
                %        axis([0 mode.vv0_dd(vind(end))+.2 0 mode.ii_dd(vind(end))+1])
                axis([0 mode.vv0_dd(end)+.2 0 1000*mode.ii_dd(end)+1])
                xlabel('Voltage, V')
                ylabel('Current, mA')
                title(['Vcurrent = ',num2str(v0_dd),' V'])
                %        legend('Measurements','Simulation','Location','Best')
                
                subplot(1,2,2)
                hold on
                grid on
                box on
                plot(Imeas,Lmeas,'r','LineWidth',2)
                PPst=sum(mode.Pst_dd,1)+mode.Psp_dd;
                
                plot(mode.ii_dd*1000,mode.Pst_dd,'.','LineWidth',2)
                %        axis([0 mode.ii_dd(vind(end))+1 0 mode.Pst_dd(vind(end))+.5])
                axis([0 1000*mode.ii_dd(end)+1 0 max(sum(mode.Pst_dd,1))*1.1+.1])
                %        legend('Measurements','Simulation','Location','Best')
                plot(mode.ii_dd(vind)*1000,PPst(vind),'bo','LineWidth',2)
                
                
                xlabel('Current, mA')
                ylabel('Optical power, mW')
                title(['Voltage ',num2str(indv),' of ',num2str(NVbias)])
                drawnow
                
                figure(1235+kfig),clf
                set(gcf,'Position',[  634    73   579   365])
                subplot(1,2,1)
                hold on
                grid on
                box on
                plot(mode.vv0_dd,mode.DeltaTmax,'k','LineWidth',2)
                plot(mode.vv0_dd,mode.DeltaTmax_Rec,'r--','LineWidth',2)
                plot(mode.vv0_dd,mode.DeltaTmax_Joule,'c--','LineWidth',2)
                plot(mode.vv0_dd,mode.DeltaTmax_OptAbs,'g--','LineWidth',2)
                plot(mode.vv0_dd(vind),mode.DeltaTmax(vind),'ko','LineWidth',2)
                xlabel('Voltage, V')
                ylabel('Temperature rise, K')
                %        legend('Total','Non-rad rec.','Joule','Opt. abs.','Location','Best')
                subplot(1,2,2)
                hold on
                grid on
                box on
                plot(mode.vv0_dd,mode.Gmod,'k','LineWidth',2)
                plot(mode.vv0_dd,mode.Lmod,'r--','LineWidth',2)
                if(mode.nmodes>1)
                    plot(mode.vv0_dd(vind),mode.Gmod(:,vind)','ko','LineWidth',2)
                    plot(mode.vv0_dd(vind),mode.Lmod(:,vind)','ro','LineWidth',2)
                else
                    plot(mode.vv0_dd(vind),mode.Gmod(vind),'ko','LineWidth',2)
                    plot(mode.vv0_dd(vind),mode.Lmod(vind),'ro','LineWidth',2)
                end
                xlabel('Voltage, V')
                ylabel('Modal gain vs losses, cm^{-1}')
                
                gipos=find(mode.Gmod(1,:)>0);
                if length(gipos)>1
                    xli=mode.vv0_dd(gipos([1 end]) ) +[-.1,.1];
                    axis([xli 0 max(max(mode.Lmod))] )
                end
                %        legend('Gain','Losses','Location','Best')
                %keyboard
                drawnow
                
                
                if(indv>1)
                    figure(1236+kfig),clf
                    set(gcf,'Position',[631   547   585   365])
                    subplot(1,2,1)
                    hold on
                    grid on
                    box on
                    Resistence=diff(mode.vv0_dd)./diff(mode.ii_dd);
                    Resistence=[Resistence(1),Resistence];
                    % plot(mode.ii_dd*1000,Resistence,Imeasinterp,Resmeasinterp)
                    plot(Imeasinterp,Resmeasinterp,'r-','LineWidth',2)
                    plot(Imeas_res,Rmeas,'.')
                    lres=1:length(Resistence);
                    plot(mode.ii_dd(lres)*1000,Resistence(lres),'go','LineWidth',2)
                    ylim([40,120])
                    %            xlim([0.5,1000*mode.ii_dd(end)+.1])
                    xlabel('Current, mA')
                    ylabel('Differential resistence, \Omega')
                    %            legend('Measurements','Simulation','Location','Best')
                    subplot(1,2,2)
                    hold on
                    grid on
                    box on
                    plot(Cur,LAM,'ro-','LineWidth',2)
                    plot(mode.ii_dd*1000,mode.lambda,'b.','LineWidth',2)
                    %            legend('Measurements','Simulation','Location','Best')
                    if(mode.nmodes>1)
                        plot(mode.ii_dd(vind)*1000,mode.lambda(:,vind),'bo','LineWidth',2)
                    else
                        plot(mode.ii_dd(vind)*1000,mode.lambda(vind),'bo','LineWidth',2)
                    end
                    %            xlim([0.,1000*mode.ii_dd(end)+.1])
                    xlabel('Current, mA')
                    ylabel('Modal wavelength, nm')
                    drawnow
                end
            end  %IPLOT
        else
            if IPLOT==1
                figure(1234+kfig),clf
                set(gcf,'Position',[268 533 1096 420])
                %             subplot(1,2,1)
                hold on
                grid on
                box on
                %             plot(Vmeas,Imeas,'r','LineWidth',2)
                plot(mode.vv0_dd,mode.ii_dd*1000,'b.','LineWidth',2)
                xlabel('Voltage, V')
                ylabel('Current, mA')
                title(['Vcurrent = ',num2str(v0_dd),' V'])
                
                keyboard
                %        legend('Measurements','Simulation','Location','Best')
                %             subplot(1,2,2)
                %             if(indv>1)
                %                 hold on
                %                 grid on
                %                 box on
                %                 Resistence=diff(mode.vv0_dd)./diff(mode.ii_dd);
                %                 Resistence=[Resistence(1),Resistence];
                %                 plot(Imeasinterp,Resmeasinterp,'r-','LineWidth',2)
                %                 plot(mode.ii_dd*1000,Resistence,'b.','LineWidth',2)
                %                 ylim([75,120])
                %                 xlabel('Current, mA')
                %                 ylabel('Differential resistence, \Omega')
                %                 %            legend('Measurements','Simulation','Location','Best')
                %                 drawnow
                %             end
            end %IPLOT
        end
        
        
        if(flagConv)
            %        eval(['save ',Last_Workspace])
            savePlotD
            %        'verifica save', keyboard
            eval(['save ',Last_Workspace,'_Plot',rad_setting,' modePlot'])
            MODEplot{kpar}=modePlot;
            eval(['save ',nomeSave,num2str(IPAR),'.mat MODEplot'])
            % save('Last_plot.mat','modePlot')
            if exist('V0_save')
                if v0_dd==V0_save
                    save Riparto
                end
            end
        end
        
        indv=indv+1;
        
        
        %  'prima save dentro save', keyboard
        %eval(['save ',nomeSW,'Contributi_',nomeSav,' ',' MODEplot -append '])
        
        %  save Contributi MODEplot -append
        
    end
end

indv=indv-1;

%'qui irest', keyboard
if irest==0
    %     save sadu MODEplot kpvet kpar IPvet PMAT NP rad_setting  rad_settingV
    %     clear MODEplot kpvet kpar IPvet PMAT NP rad_setting  rad_settingV
    %     eval(['save ',Last_Workspace])
    %     load sadu
    % SAVE MESH AND MODE AS IN D1ANA
    eval(['save ',Last_Workspace,'_',strSave,'_',nomeSav,'.mat geom mesh mode'])
else
    save ultimo2
end
Resistence=diff(mode.vv0_dd)./diff(mode.ii_dd);
ind=5:length(Resistence);
