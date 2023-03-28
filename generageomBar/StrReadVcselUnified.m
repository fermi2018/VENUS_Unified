flgStop=0;
if flgStop==1
'UNIFIED', keyboard
end
Buffer=ParMore.Buffer;
Mesa=ParMore.mesa;
Contact=ParMore.contact;

repdd=rep;  % repetitions of each layer
isop=1;


% P side - eq p parameters extraction
if isfield(mireq,'p')==1
    mir_p=mireq.p;
    Np_pdd=mir_p.P_f-mir_p.P_i;    % number of p-doped pairs substituted by the eq layer
else
    isop=-1;
    mir_p.P_i=0;
    mir_p.P_f=0;
    Np_pdd=0;
    ParMore.mireq.p.x_eq=[  ];
    ParMore.mireq.p.Dop_eq=0;
    ParMore.mireq.p.Nmesh=[];
    ParMore.mireq.p.Mateq='niente';
end
%
% N side (bottom DBR) - eq n parameters extraction
if isfield(mireq,'n')==1
    mir_n=mireq.n;
    Np_ndd=mir_n.P_f-mir_n.P_i;
else
    ParMore.mireq.n.x_eq=[];
    ParMore.mireq.n.Dop_eq=0;
    ParMore.mireq.n.Nmesh=[];
    ParMore.mireq.n.Mateq='niente';
    mir_n.P_i=0;
    mir_n.P_f=0;
    Np_ndd=0;
end
%
% N side (top DBR)- eq n parameters extraction
if isfield(mireq,'n2')==1
    mir_n2=mireq.n2;
    Np_n2dd=mir_n2.P_f-mir_n2.P_i;
else
    ParMore.mireq.n2.x_eq=[];
    ParMore.mireq.n2.Dop_eq=0;
    ParMore.mireq.n2.Nmesh=[];
    ParMore.mireq.n2.Mateq='niente';
    mir_n2.P_i=0;
    mir_n2.P_f=0;
    Np_n2dd=0;
end

%
% P side equivalent medium

repe=rep;   % repetitions of each layer

pun=[1:length(rep)]';   % indeces from 1 to number of layers in the str file


if isfield(mireq,'p')

fiP=find(Na(:,1)>0);    % indeces for Na>1
fi=find(rep(fiP,1)>1);  % indeces for Na>1 && #repetitions>1
difi=diff(fi);          % understand the position of "holes" in various repetitions
fdis=find(difi>1);      % understand the position of "holes" in various repetitions
if length(fdis)>0
    pip=0;
    for kk=1:length(fdis)
        fiS{kk}=fiP(fi(pip+1:fdis(kk)));
        pip=fdis(kk);
        np(kk)=rep(fiS{kk}(1));
    end
    kk=kk+1;
    fiS{kk}=fiP(fi(pip+1:end)); % cells containing the #repetitions of each "block"
    np(kk)=rep(fiS{kk}(1));     % # of repetition of each block identified in fiS
    [~,npMa]=max(np);  % idenntifies the max number of repetitions and its position in np
    fi=fiS{npMa};       % indeces of the str layers that repeat the most among the p-doped ones
else
    fi=fiP(fi);
end


if length(fi)>0
    ret=rep(fi,:); % matrix containing #repetitions and #layers involved
    Np_mirp=ret(1,1);  % #repetition of the block
else
    Np_mirp=0;
end

%'qui verp', keyboard
if mir_p.P_i>Np_mirp || mir_p.P_f>Np_mirp
    disp('ERRORE strato equivalente p-side: n. paia equivalenti maggiore di quelle reali'), keyboard
end
if mir_p.P_i<0 || mir_p.P_f<0
    disp('ERRORE strato equivalente p-side: n. paia minore di zero'), keyboard
end

% Identifies the eq layers
if mir_p.P_i<mir_p.P_f
    if length(fi)>0
        repP=rep(fi,:);
        repPs=repP;
        repef=repP;
        repP(:,1)=mir_p.P_i;    % first column containing #repetitions is replaced by the "starting" layer
        repPs(:,1)=Np_pdd;      % equal to ret
        repPr=repP;
        repef(:,1)=rep(fi(1),1)-mir_p.P_f;  % p-doped layers identified are replaced with 0
        NL_pdd=repPs(1,2);
        repe(fi,1)=rep(fi(1),1)-mir_p.P_f;  % p-doped layers identified are replaced with 0
        %'qui int', keyboard
        if fi(1)>1
            fiPr=(1:fi(1)-1)';  % indeces from 1 to first repeated layer
            repP=repe(fiPr,:);  % layers repetitions BEFORE the actual repetitions
            repe=repe(fi(1):end,:); % repe excluding the layers identified in repP
        else
            fiPr=[];
            repP=[];
        end
        ret=rep(fi,:);  % description of repeated layers (=repPs)
    end
    %'qui', keyboard
    
    if mir_p.P_i>0 & mir_p.P_f<Np_mirp
        repe=[repP; repef; -repPs; repPr; repe(length(repPr)+1:end,:)];
        pun=[fiPr; fi; fi; pun(fi(1):end)];
    elseif mir_p.P_i==0 & mir_p.P_f>0 & mir_p.P_f<Np_mirp
        % repep=
        repe=[repP; -repPs;  repe];
        pun=[fiPr; fi; [fi(1):pun(end)]'];
    elseif mir_p.P_i>0 & mir_p.P_f==Np_mirp
        repe=[repP; -repPs; repe];
        %pun=[fi; fi; pun];
        %fi0=find(repe(:,1)~=0);
        %repe=repe(fi0,:);
        %pun=pun(fi0);
        pun=[fiPr; fi; [fi(1):pun(end)]'];
        fi0=find(abs(repe(:,1))==0);
        repe(fi0,:)=repPr;
    elseif mir_p.P_i==0 & mir_p.P_f==Np_mirp & Np_mirp~=0
        %'ver qui', keyboard
        % fi=find(rep(:,1)==Np_mirp);
        repe=rep;
        repe(fi,:)=-repPs;  % replaces the positive #repetitions and #layers with NEGATIVE
        % repe=[repP; -repPs;  repe];
        % pun=[fiPr; fi; pun];
        % fi0=find(repe(:,1)~=0);
        % repe=repe(fi0,:);
        % pun=pun(fi0);
        pun=[1:length(repe)]';  % indeces from 1 to #layers in str file
    end
end  %no equivalente

end % Na equivalente

%%
% N side (bottom DBR) equivalent medium (NEW)
rep=repe;

fiP=find(Nd(pun,1)>0);  % indeces for Nd>1
fid=find(abs(rep(fiP,1))>1);    % indeces for Nd>1 && #repetitions>1
difi=diff(fid);          % understand the position of "holes" in various repetitions
fdis=find(difi>1);      % understand the position of "holes" in various repetitions

if length(fdis)>0
    pip=0;
    for kk=1:length(fdis)
        fiS{kk}=fiP(fid(pip+1:fdis(kk)));
        pip=fdis(kk);
        np(kk)=rep(fiS{kk}(1));
    end
    kk=kk+1;
    fiS{kk}=fiP(fid(pip+1:end)); % cells containing the #repetitions of each "block"
    np(kk)=rep(fiS{kk}(1));     % # of repetition of each block identified in fiS
    [~,npMa]=max(np);  % idenntifies the max number of repetitions and its position in np
    fi=fiS{npMa};       % indeces of the str layers that repeat the most among the p-doped ones
        
    fim=fi(1);  % first index of the n-side repetitions
    reu=repe(fi,:); % takes from repe the #reps ans #layers (similar for ret for p-side)
    fiu(1)=fi(end)+1;   % index of first layer AFTER the p-side repetitions
    fiu(2)=length(rep); % index of the last layer of the str file
    % fiu is similar to fiPr for the p-side
    repef=repe(fiu(1):fiu(2),:); % takes the repe part AFTER the repetitions (similar to repP)
    
    fius=fiu;
    
    repP=rep(fi,:);
    repPs=rep(fi,:);
    % These two are similar to the repP and repPs already used for n-side
    repP(:,1)=mir_n.P_i;    % contains the #layers (n) "untouched" by the eq medium and #layers
    repPs(:,1)=Np_ndd;      % contains the #layers (n) "replaced" by the eq medium and #layers
    pumis=pun(fi);  % indeces of n-side repetitions (equal to fi)
    pumlast=pun(fi(end)+1:length(pun)); % indeces of layer AFTER the repetitions
   
    Np_mirn=reu(1,1);
    
    if mir_n.P_i>Np_mirn ||  mir_n.P_f>Np_mirn
        fprintf('ERRORE strato equivalente n-side: n. paia equivalenti maggiore di quelle reali\n'), keyboard
    end
    if mir_n.P_i<0 || mir_n.P_f<0
        fprintf('ERRORE strato equivalente n-side: n. paia minore di zero\n'), keyboard
    end
    %if mir_n.P_i>0 &  mir_n.P_f<Np_mirn
    iul=0;
    %'qui', keyboard
    if mir_n.P_i<mir_n.P_f
        
        repe(fi,1)=rep(fi(1),1)-mir_n.P_f; % places 0s where n-side repetitions are present
        
        if mir_n.P_i>0
            irep=2;
            pum0=0:repP(1,2)-1; % indeces from 0 to #layers-1 involved by the repetition
            reped=[repP; -repPs];
            pum=0:2*repP(1,2)-1; % indeces from 0 to 2*(#layers)-1 involved by the repetition
            % substitutes in repe the #repetitions in repeated layers and in layers till the end with reped
            % repe also becomes longer
            repe(fim+pum,:)=reped;
            pumad=[pumis; pumis];
            if mir_n.P_f<Np_mirn
                fiu=fiu+2*repP(1,2);
                repe(fiu(1):fiu(2),:)=repef;
                repe(fim+pum0+2*reu(1,2),:)=reu;
                repe(fim+pum0+2*reu(1,2),1)=reu(1,1)-mir_n.P_f;
            else
                pumad=[pumis];
                pum=pum0;
                fiu=fiu+repP(1,2); % fiu is augmented by the number of layers "excluded" by the eq medium
                % Now repe is made longer to include the 33 layers the will
                % become "equivalent" and the remaining 3 excluded by this
                % procedure: the #layer added depends on the #layers
                % involved by the repetitions
                repe(fiu(1):fiu(2),:)=repef; % Adds the features of the layers AFTER the reps, repe is longer
            end
        elseif mir_n.P_i==0 & mir_n.P_f<=Np_mirn
            %else
            irep=1;
            pumad=[pumis];
            
            pum=0:repP(1,2)-1;
            reped=[-repPs];
            repe(fiu(1):fiu(2),:)=repef;
            repe(fim+pum,:)=reped;
            if reu(1,1)-mir_n.P_f>0
                fiu=fiu+repP(1,2);
                repe(fim+pum+reu(1,2),1)=reu(1,1)-mir_n.P_f;
            else
                iul=1;
            end
            %'qui',keyboard
            %elseif mir_n.P_i==0 & mir_n.P_f==Np_mirn
            % iul=1;
            
        end
        if iul==0
            pun(fi(end)+pum+1)=pumad;
            pun=[pun; pumlast]; % pun is enlarged accordingly to "new" repe
        end
        
    end % strato equivalente
    % Na, Nd, d
elseif length(fid)>0
    fi=fiP(fid);
    
    fim=fi(1);  % first index of the n-side repetitions
    reu=repe(fi,:); % takes from repe the #reps ans #layers (similar for ret for p-side)
    fiu(1)=fi(end)+1;   % index of first layer AFTER the p-side repetitions
    fiu(2)=length(rep); % index of the last layer of the str file
    % fiu is similar to fiPr for the p-side
    repef=repe(fiu(1):fiu(2),:); % takes the repe part AFTER the repetitions (similar to repP)
    
    fius=fiu;
    
    repP=rep(fi,:);
    repPs=rep(fi,:);
    % These two are similar to the repP and repPs already used for n-side
    repP(:,1)=mir_n.P_i;    % contains the #layers (n) "untouched" by the eq medium and #layers
    repPs(:,1)=Np_ndd;      % contains the #layers (n) "replaced" by the eq medium and #layers
    pumis=pun(fi);  % indeces of n-side repetitions (equal to fi)
    pumlast=pun(fi(end)+1:length(pun)); % indeces of layer AFTER the repetitions
   
    Np_mirn=reu(1,1);
    
    if mir_n.P_i>Np_mirn ||  mir_n.P_f>Np_mirn
        fprintf('ERRORE strato equivalente n-side: n. paia equivalenti maggiore di quelle reali\n'), keyboard
    end
    if mir_n.P_i<0 || mir_n.P_f<0
        fprintf('ERRORE strato equivalente n-side: n. paia minore di zero\n'), keyboard
    end
    %if mir_n.P_i>0 &  mir_n.P_f<Np_mirn
    iul=0;
    %'qui', keyboard
    if mir_n.P_i<mir_n.P_f
        
        repe(fi,1)=rep(fi(1),1)-mir_n.P_f; % places 0s where n-side repetitions are present
        
        if mir_n.P_i>0
            irep=2;
            pum0=0:repP(1,2)-1; % indeces from 0 to #layers-1 involved by the repetition
            reped=[repP; -repPs];
            pum=0:2*repP(1,2)-1; % indeces from 0 to 2*(#layers)-1 involved by the repetition
            % substitutes in repe the #repetitions in repeated layers and in layers till the end with reped
            % repe also becomes longer
            repe(fim+pum,:)=reped;
            pumad=[pumis; pumis];
            if mir_n.P_f<Np_mirn
                fiu=fiu+2*repP(1,2);
                repe(fiu(1):fiu(2),:)=repef;
                repe(fim+pum0+2*reu(1,2),:)=reu;
                repe(fim+pum0+2*reu(1,2),1)=reu(1,1)-mir_n.P_f;
            else
                pumad=[pumis];
                pum=pum0;
                fiu=fiu+repP(1,2); % fiu is augmented by the number of layers "excluded" by the eq medium
                % Now repe is made longer to include the 33 layers the will
                % become "equivalent" and the remaining 3 excluded by this
                % procedure: the #layer added depends on the #layers
                % involved by the repetitions
                repe(fiu(1):fiu(2),:)=repef; % Adds the features of the layers AFTER the reps, repe is longer
            end
        elseif mir_n.P_i==0 & mir_n.P_f<=Np_mirn
            %else
            irep=1;
            pumad=[pumis];
            
            pum=0:repP(1,2)-1;
            reped=[-repPs];
            repe(fiu(1):fiu(2),:)=repef;
            repe(fim+pum,:)=reped;
            if reu(1,1)-mir_n.P_f>0
                fiu=fiu+repP(1,2);
                repe(fim+pum+reu(1,2),1)=reu(1,1)-mir_n.P_f;
            else
                iul=1;
            end
            %'qui',keyboard
            %elseif mir_n.P_i==0 & mir_n.P_f==Np_mirn
            % iul=1;
            
        end
        if iul==0
            pun(fi(end)+pum+1)=pumad;
            pun=[pun; pumlast]; % pun is enlarged accordingly to "new" repe
        end
        
    end % strato equivalente
    % Na, Nd, d
end  %strato n

%%
% N side (TOP DBR) - equivalent medium

if isfield(mireq,'n2')
%     pun=[1:length(repe)]';   % indeces from 1 to number of layers in the str file
    
    fiP=find(Nd(:,1)>0);    % indeces for Nd>1
    fi=find(rep(fiP,1)>1);  % indeces for Nd>1 && #repetitions>1
    difi=diff(fi);          % understand the position of "holes" in various repetitions
    fdis=find(difi>1);      % understand the position of "holes" in various repetitions
    if length(fdis)>0
        pip=0;
        for kk=1:length(fdis)
            fiS{kk}=fiP(fi(pip+1:fdis(kk)));
            pip=fdis(kk);
            np(kk)=repe(fiS{kk}(1));
        end
        kk=kk+1;
        fiS{kk}=fiP(fi(pip+1:end)); % cells containing the #repetitions of each "block"
        np(kk)=repe(fiS{kk}(1));     % # of repetition of each block identified in fiS
        [~,npMa]=max(np);  % idenntifies the max number of repetitions and its position in np
        fi=fiS{npMa};       % indeces of the str layers that repeat the most among the p-doped ones
    else
        fi=fiP(fi);
    end
    
    
    if length(fi)>0
        ret=repe(fi,:); % matrix containing #repetitions and #layers involved
        Np_mirn2=ret(1,1);  % #repetition of the block
    else
        Np_mirn2=0;
    end
    
    %'qui verp', keyboard
    if mir_n2.P_i>Np_mirn2 || mir_n2.P_f>Np_mirn2
        fprintf('ERRORE strato equivalente n-side: n. paia equivalenti maggiore di quelle reali\n'), keyboard
    end
    if mir_n2.P_i<0 || mir_n2.P_f<0
        fprintf('ERRORE strato equivalente n-side: n. paia minore di zero\n'), keyboard
    end
    
    % Identifies the eq layers
    if mir_n2.P_i<mir_n2.P_f
        if length(fi)>0
            repP=rep(fi,:);
            repPs=rep(fi,:);
            repef=rep(fi,:);
            repP(:,1)=mir_n2.P_i;    % first column containing #repetitions is replaced by the "starting" layer
            repPs(:,1)=Np_n2dd;      % equal to ret; contains the #layers of the repetition substituted by mireq
            repef(:,1)=rep(fi(1),1)-mir_n2.P_f;  % p-doped layers identified are replaced with 0
            NL_n2dd=repPs(1,2);
            repe(fi,1)=rep(fi(1),1)-mir_n2.P_f;  % p-doped layers identified are replaced with 0
            %'qui int', keyboard
            if fi(1)>1
                fiPr=(1:fi(1)-1)';  % indeces from 1 to first repeated layer
                repPr=repe(fiPr,:);  % layers repetitions BEFORE the actual repetitions
                repe=repe(fi(1):end,:); % repe excluding the layers identified in repP
            else
                fiPr=[];
                repPr=[];
            end
            ret=rep(fi,:);  % description of repeated layers (=repPs)
        end
        %'qui', keyboard
        
        if mir_n2.P_i>0 && mir_n2.P_f<Np_mirn2
            repe=[repPr; repef; -repPs; repP; repe(length(repP)+1:end,:)];
            pun=[fiPr; fi; fi; pun(fi(1):end)];
        elseif mir_n2.P_i==0 && mir_n2.P_f>0 && mir_n2.P_f<Np_mirn2
            % repep=
            repe=[repPr; -repPs; repP; repe];
            pun=[fiPr; fi; [fi(1):pun(end)]'];
        elseif mir_n2.P_i>0 && mir_n2.P_f==Np_mirn2
            repe=[repPr; -repPs; repP; repe(length(repP)+1:end,:)];
            pun=[1:length(repPr), fi', fi:pumad(end), pumad', pumlast']';
            fi0=find(repe(:,1)~=0);
            repe=repe(fi0,:);
            pun=pun(fi0);
        elseif mir_n2.P_i==0 && mir_n2.P_f==Np_mirn2 && Np_mirn2~=0
            %'ver qui', keyboard
            % fi=find(rep(:,1)==Np_mirp);
            repe=rep;
            repe(fi,:)=-repPs;  % replaces the positive #repetitions and #layers with NEGATIVE
            % repe=[repP; -repPs;  repe];
            % pun=[fiPr; fi; pun];
            % fi0=find(repe(:,1)~=0);
            % repe=repe(fi0,:);
            % pun=pun(fi0);
            pun=[1:length(repe)]';  % indeces from 1 to #layers in str file
        end
    end  %no equivalente
end
%%
%if Buffer.thdd>0
d(end)=Buffer.thdd;
%else
% d=d(1:end-1);
%end

fprintf('After mireq stacking - Check "pun"\n')%, keyboard

rep=repe; % repe contains now the upgraded repe: equivalent layers are labeled with NEGATIVE number of reps
% Mesh contains the number of mesh points of each str layer;
% In this way the number of mesh points in the reps "excluded" by the eq
% medium (in case mireq.P_i>0) are assigned to a correct number of mesh
% points (therefore mesh_yd is longer than mesh). Similar procedure is
% applied to ragv, x, xmat, Nadd, Nddd
mesh_yd=mesh(pun); 
ragv=rag(pun,:);
x=x(pun,:);
xmat=ParDD.xm(pun,:);
firag=find(ragv.*xmat<0);
raot=ragv(firag);

xdd=x;

klv=1:length(rep);  % indeces of the updated repe
Nadd=Na(pun,:);
Nddd=Nd(pun,:);
% Old values are overwritten by the updated versions
d=d(pun);
Na=Nadd;
Nd=Nddd;
mesh=mesh_yd;
%

layerStacking_NOlitho

% 'qui basta impilamento', keyboard
%fine

klv=1:length(NA1_i);

for kl=klv
    PM=Pmat{kl};
    material(kl,1:length(PM))=PM;
    if kl<10
        Dch(kl,1:9)=['Layer_i_',num2str(kl)];
    elseif kl<100 & kl>9
        Dch(kl,1:10)=['Layer_i_',num2str(kl)];
    elseif kl<1000 & kl>99
        Dch(kl,1:11)=['Layer_i_',num2str(kl)];
    end
end

fi=find(NA1_i(:,1)>0);
for kl=fi'
    if kl<10
        Dch(kl,1:9)=['Layer_p_',num2str(kl)];
    elseif kl<100 & kl>9
        Dch(kl,1:10)=['Layer_p_',num2str(kl)];
    elseif kl<1000 & kl>99
        Dch(kl,1:11)=['Layer_p_',num2str(kl)];
    end
end

fi=find(ND1_i(:,1)>0);
for kl=fi'
    if kl<10
        Dch(kl,1:9)=['Layer_n_',num2str(kl)];
    elseif kl<100 & kl>9
        Dch(kl,1:10)=['Layer_n_',num2str(kl)];
    elseif kl<1000 & kl>99
        Dch(kl,1:11)=['Layer_n_',num2str(kl)];
    end
end

reprep=[];
k=1;
kc=1;
fiCA=[];
while k<=length(rep)
    
    rei=rep(k,1);
    nei=rep(k,2);
    if nei==0
        nei=1;
    end
    if length(find(fiCav-k)==0)==1
        fiCAa=1;
    else
        fiCAa=0;
    end
    
    
    if rei<0
        reprep=[reprep; -1];
    elseif rei==1i
        reprep=[reprep; rei];
    else
        reprep=[reprep; repmat(ones(nei,1),rei,1)];
        fiCAa=fiCAa*repmat(ones(nei,1),rei,1);
    end
    fiCA=[fiCA; fiCAa];
    k=k+abs(nei);
    %  pausak
end
%'reprep', keyboard
% looks for imaginary reps, corresponding to QWs
lab=Dch;
fiq=find(reprep==1i);
finq=find(reprep>=1);
ic=0;

if length(fiq)>0
    for kl=fiq'
        ic=ic+1;
        %  'ci', keyboard
        lab(kl,:)=0;
        if ic<10
            lab(kl,1:3)=['qw',num2str(ic)];
        else
            lab(kl,1:4)=['qw',num2str(ic)];
        end
        %lab(kl+1,1:3)=['qb',num2str(ic)];
    end
end
% looks for 1s in fiCA, which contains layers LcC and LgC
fiCAV=find(fiCA==1)';
for kf=fiCAV
    lab(kf,:)=0;
    lab(kf,1:3)='Cav';
end

% -------------------------------------------------------------------------
% Torrelli Valerio, 03/11/2021
% I need to modify the following section in order to include multiple BTJs.
% Of course i would like to flag them like, TJ1, TJ2, TJ3, .... 


% Now, the following is Torrelli procedure
if flgStop==1
    'TO', keyboard
end
if isfield(ParDD,'iNEGF_multiple')==1 % if there are BTJs in the structure
%     
    nBTJ=ParDD.nBTJ; % total number of BTJs within the  VCSEL
    iNEGF_multiple=ParDD.iNEGF_multiple; % cell containing the position of the BTJs
    
    for indexBTJ=1:nBTJ
        nLayers=length(iNEGF_multiple{indexBTJ});
        
        if isfield(mireq,'n2') && ~isempty(mireq.n2.Nmesh)
            iNEGF_multiple{indexBTJ}(:)=iNEGF_multiple{indexBTJ}(:)+1-(Np_n2dd*NL_n2dd);
        elseif isfield(mireq,'p') && ~isempty(mireq.p.Nmesh)
            iNEGF_multiple{indexBTJ}(:)=iNEGF_multiple{indexBTJ}(:)+1-(Np_pdd*NL_pdd);
        end
        
        for indexLayer=1:nLayers
            LA=lab(iNEGF_multiple{indexBTJ}(indexLayer),:);
            
            if indexLayer==1 && indexBTJ==1
                lab(iNEGF_multiple{indexBTJ}(indexLayer),1:length(LA)+6)=[LA,'_TJ',num2str(indexBTJ),'_',num2str(indexLayer)];
            else
                LA(end-5:end)=['_TJ',num2str(indexBTJ),'_',num2str(indexLayer)];
                lab(iNEGF_multiple{indexBTJ}(indexLayer),1:length(LA))=LA;
            end           
            
        end
    end
    
    fprintf('Check the labels of the TJ layers!\n')%,keyboard
    
end

% -------------------------------------------------------------------------

% parte per plot

Ds=cumsum(D1);
ddd=[[0; Ds(1:end-1)] Ds];
xxd=[[0; X1(1:end-1)] X1];
nad=[[0; NA1(1:end-1)] NA1];
ndd=[[0; ND1(1:end-1)] ND1];



dd=reshape(ddd',prod(size(ddd)),1);
xx=reshape(xxd',prod(size(ddd)),1);
na=reshape(nad',prod(size(ddd)),1);
nd=reshape(ndd',prod(size(ddd)),1);

n=ParOpt.n(2:end-1,2);
%iraot=find(ParDD.rag(finq)>0 & n<2);


%'iraot', keyboard
%'iraot', keyboard
Prag=Rag;


%iraot=find(Prag>0 & X1_i(:,1)>1e-4);
%irarel=find(Prag>0 & X1_i(:,1)<1e-4);

iratot=find(Prag>0);



%iraot=finq(iraotd);
raot=unique(Prag(iratot)');


ratot=sort([raot ParDD.ragAdd]);
%raox=unique(Prag(iraot)');
%'rapt', keyboard
ieff=find(reprep==-1);

io=0;

if io==1
    
    lab(iraot,:)=0;
    lab(irarel,:)=0;
    for k=iraot'
        lab(k,1:2)='Ox';
    end
    for k=irarel'
        lab(k,1:3)='Rel';
    end
    
    lab(ieff,:)=0;
    for k=ieff'
        lab(k,1:2)='Eq';
    end
end

labo=lab;

%ratot=[ParVet(1)+[-ParDD.ragAdd 0 ParDD.ragAdd ] ParVet(2)];
kind=0;
for k=1:length(material)
    if findstr(material{k,2},'AlOx')==1
        kind=kind+1;
        Ks(kind)=k;
    end
end

iOxTutti=0;

if iOxTutti==0
    fidiv=find(diff(Ks)>1);
    RoxT=Rag(Ks,1);
    if length(fidiv)==0
        % Single ox case
        raotu=min(Rag(Ks,1));
    else
        % More than one OX case, maybe not correct for identical radi
        filoc=1:fidiv(1)+1;
        raotu=min(Rag(Ks(filoc),1));
        for kfi=2:length(fidiv)-1
            filoc=fidiv(kfi-1)+2:fidiv(kfi)+1;
            raotu=[raotu min(Rag(Ks(filoc),1))];
        end
        filoc=fidiv(end)+2:length(Ks);
        raotu=[raotu min(Rag(Ks(filoc),1))];
    end
    raver=[];
    for kot=1:length(raot)
        rai=raot(kot);
        if length(find(RoxT==rai))==0
            raver=[raver rai];
        end
    end
    % "unique" is the modification: it is very SILLY, maybe not compatible
    % with other structures, BE CAREFUL PLEASE!
    raver=sort(unique([raver raotu]));  
else
    raver=raot;
end


%'Ks', keyboard

if exist('Ks')
    fi=find(Rag(Ks,1)>0);
    Rox=min(Rag(Ks(fi),1));
else
    Rox=0;
end

R_contact=Rag(1,:);

if Str1D.quasi1D==1
    ParDD.ragAdd=[0 0];
end
raoADD=Rox+[-ParDD.ragAdd(1) 0 ParDD.ragAdd(2)];


fira=find(diff([0 raoADD])>0);
raoADD=raoADD(fira);

fiO=find(raver==Rox);

raIns=[raver(1:fiO-1) raoADD raver(fiO+1:end)];

diffRA=diff([0 raIns]);
mesh_r=[Contact.Nrmesh];
fiMze= find(diffRA==0);



%'qoiaer\', keyboard

if length(fiMze)==1
    mesh_r0=mesh_r(fiMze);
    fiS=fiMze;
    fiRao=find(raIns(fiS)~=raoADD);
    raoADD=raoADD(fiRao);
    fiMze= find(diffRA>0);
    raIns=raIns(fiMze);
    mesh_r=mesh_r(fiMze);
    if fiS==1
        mesh_r(fiS)=mesh_r(fiS)+mesh_r0;
    else
        mesh_r(fiS-1)=mesh_r(fiS-1)+mesh_r0;
    end
end


diffRA=diff([0 raIns]);

fiMze= find(diffRA<0);

if length(fiMze)==1
    mesh_r0=mesh_r(fiMze);
    fiS=fiMze;
    % 'qui', keyboard
    fiRao=find(raIns(fiS)~=raoADD);
    raoADD=raoADD(fiRao);
    fiMze= find(diffRA>0);
    raIns=raIns(fiMze);
    mesh_r=mesh_r(fiMze);
    if fiS==1
        mesh_r(fiS)=mesh_r(fiS)+mesh_r0;
    else
        mesh_r(fiS-1)=mesh_r(fiS-1)+mesh_r0;
    end
end

ra_dd=raIns;
sezionir=diff([0 ra_dd]);
if flgStop==1
     'qui sezioni', keyboard
end

%step_r=real(mesh_r(1:size(rag,2));
%punti_r=imag(mesh_r);

%if length(ra_dd)==1

%step_r=real(mesh_r(1:size(rag,2)+1));
%punti_r=imag(mesh_r(1:size(rag,2)+1));
%else
 step_r=real(mesh_r);
 punti_r=imag(mesh_r);
%end

if length(punti_r)~=length(ra_dd) & Str1D.quasi1D==0
 ' parametri radiali non allineati !!!'
 ra_dd
 mesh_r
 keyboard
end

if isfield(Str1D,'quasi1D')
if Str1D.quasi1D==1
    step_r(1)=0.5;
    punti_r(1)=1;
end
end
%ra_dd=sort([ratot Contact.radius_i Contact.radius_e]);
NPrad0=ceil(sezionir./step_r);
NPrad=NPrad0;
fival=find(NPrad-punti_r>0);
if length(fival)>0
    NPrad(fival)=punti_r(fival);
end
mesh_r=NPrad;

%'qui', keyboard
%   disp(['guardiamo raot']), keyboard

% mesh_r=[Contact.Nr_i Contact.Nr_m Contact.Nr_e];


%'raggi', keyboard


if length(mesh_r)~=length(ra_dd)
    'sezioni radiali e mesh radiali di diversa lunghezza: sistemare il file str'
    if Str1D.quasi1D==0
        keyboard
    end
    if length(mesh_r)>length(ra_dd)
        finox=find([0 ra_dd]~=Rox);
        mesh_save=mesh_r;
        mesh_r=mesh_save(finox);
    else
        'sezioni radiali e mesh radiali di diversa lunghezza: sistemare il file str', keyboard
    end
end

if length(mesh_r)~=length(ra_dd)
    if fiO==1
        me0=mesh_r(fiO);
        mesh_r=mesh_r(fiO+1:end);
        mesh_r(fiO)=mesh_r(fiO)+me0;
    end
end


if length(mesh_r)~=length(ra_dd) & Str1D.quasi1D==0
    'sezioni radiali e mesh radiali di diversa lunghezza: sistemare il file str',
	keyboard
end

if iplot==1
    figure, plot(dd,xx,'.-'),
    a=axis;
    a(3)=-.1;
    a(4)=1.1;
    axis(a)
    pausak
    
    figure, semilogy(dd,na,'.r-'),
    hold on
    semilogy(dd,nd,'.g-'),
    pausak
    
    figure, plot(dd,na,'.r-'),
    hold on
    plot(dd,nd,'.g-'),
    pausak
end

Dsc=cumsum(Dc);
ddd=[[0; Dsc(1:end-1)] Dsc];
xxd=[[0; Xc(1:end-1)] Xc];
nad=[[0; NAc(1:end-1)] NAc];
ndd=[[0; NDc(1:end-1)] NDc];

ddc=reshape(ddd',prod(size(ddd)),1);
xxc=reshape(xxd',prod(size(ddd)),1);
nac=reshape(nad',prod(size(ddd)),1);
ndc=reshape(ndd',prod(size(ddd)),1);

if iplot==1
    figure, plot(ddc,xxc),
    a=axis;
    a(3)=-.1;
    a(4)=1.1;
    axis(a)
    
    figure, semilogy(ddc,nac,'r'),
    hold on
    semilogy(ddc,ndc,'g'),
    
    figure, plot(ddc,nac,'r'),
    hold on
    plot(ddc,ndc,'g'),
    pausak
end

RagN=zeros(size(material));
for kr=1:length(ra_dd)
    for kc=1:size(Rag,2)
        fio1=find(abs(Rag(:,kc))==ra_dd(kr));
        RagN(fio1,kr)=ra_dd(kr);
    end
end
if Rox>0
    fio=find(ra_dd==raoADD(1));
    for Ksi=Ks
        RagN(Ksi,fio+[0:length(raoADD)-1])=raoADD;
    end
end

%' prima Caddd ', keyboard
nrig=length(material);
ncol=length(ra_dd);
materialO=material;
for kc=1:ncol
    for kr=1:nrig
        material{kr,kc}=materialO{kr,1};
    end
end

if flgStop==1
'qui ROTIAOP', keyboard
end

inew=1;
if inew==0
    Script_mate
else
    Script_mate1
end


fiM=find(D1_i>0);
NA_dd=flipud(NA1_i(fiM,:));
ND_dd=flipud(ND1_i(fiM,:));
d_dd=flipud(D1_i(fiM));
mesh_dd=flipud(Mesh1_i(fiM));
x_dd=flipud(X1_i(fiM,:));
lab=flipud(lab(fiM,:));
materialC=flipud(material(fiM,:));

if length(ParDD.iDD)==1
    material(fiM+1,:)=materialC;
    for kc=1:size(materialC,2)
        material{1,kc}='Ground';
    end
elseif length(ParDD.iDD)==0
    material(fiM+1,:)=materialC;
    for kc=1:size(materialC,2)
        material{1,kc}='Line';
    end
    for kc=1:size(material,2)
        material{size(material,1),kc}='vacuum';
    end
    fiC=find(ra_dd==Contact.radius_i);
    material{size(material,1),fiC+1}='Line';
else
    material=materialC;
end

if length(ParDD.iDD)==1
    indici=1:length(d_dd)-1;
    indiciR=length(d_dd):-1:2;
elseif length(ParDD.iDD)==2
    indici=2:length(d_dd)-1;
    indiciR=length(d_dd)-1:-1:2;
else
    indici=1:length(d_dd);
    indiciR=length(d_dd):-1:1;
end

StrDD.d_dd=d_dd(indici);
%StrDD.raggi=Rag(end:-1:1,:);
StrDD.mesh_dd=mesh_dd(indici);
StrDD.raggi=RagN(indiciR,:);
StrDD.x_dd=fliplr(x_dd(indici,:));
StrDD.NA_dd=fliplr(NA_dd(indici,:));
StrDD.ND_dd=fliplr(ND_dd(indici,:));


% Salvo gli indici e il numero delle BTJ ----------------------------------
if isfield(ParDD,'iNEGF_multiple')
    StrDD.nBTJ = nBTJ;
    StrDD.iNEGF_multiple = iNEGF_multiple;
end
% -------------------------------------------------------------------------



StrDD.lab=lab(indici,:);
StrDD.material=material;
StrDD.mesh_r=mesh_r;
StrDD.ra_dd=ra_dd;
StrDD.ParMore=ParMore;

if length(fiCAV)==0
    fiCAV=fix(length(d)/3)*[1 2];
    fiCav=fiCAV;
end
StrDD.fiPassiv1=fiCAV(end);

%StrDD.fiPassiv1=1;
% 'fine strDD', keyboard
if isfield(StrDD.ParMore,'depth_passiv')
    Depth=StrDD.ParMore.depth_passiv;
    Dtot=cumsum(d_dd);
    Dmes=Dtot(end)-Depth*1000;
    [du,imi]=min(abs(Dtot-Dmes));
    StrDD.fiPassiv2=imi;
    %figure, plot(cumsum(d_dd)/1000,'.'), hold on, plot(imi,Dmes/1000,'go')
    % 'fine strDD', keyboard
    
end


%      'fine strDD UNIFIED NEW', keyboard