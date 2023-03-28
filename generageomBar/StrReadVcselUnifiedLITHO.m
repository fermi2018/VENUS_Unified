flgStop=0;
if flgStop==1
'UNIFIED', keyboard
end
Buffer=ParMore.Buffer;
Mesa=ParMore.mesa;
Contact=ParMore.contact;

NPsottratto=1;  % con TJ litho bisogna sottrarre 1 paio per la rottura di periodicità della TJ


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
	
	mireq.n2.P_i=mireq.n2.P_i+NPsottratto;
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

clear fiS np
if isfield(mireq,'n2')
%     pun=[1:length(repe)]';   % indeces from 1 to number of layers in the str file
    
    fiP=find(Nd(:,1)>0);    % indeces for Nd>1
    fi=find(repe(fiP,1)>1);  % indeces for Nd>1 && #repetitions>1
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

fiR=find(repe(:,1)==1);
repe(fiR,2)=0;

rep=repe; % repe contains now the upgraded repe: equivalent layers are labeled with NEGATIVE number of reps
% Mesh contains the number of mesh points of each str layer;
% In this way the number of mesh points in the reps "excluded" by the eq
% medium (in case mireq.P_i>0) are assigned to a correct number of mesh
% points (therefore mesh_yd is longer than mesh). Similar procedure is
% applied to ragv, x, xmat, Nadd, Nddd

reTJ=rep(:,1);
fiT=find(reTJ==-500);
rep(fiT,1)=1;
repe(fiT,1)=1;

if flgStop==1
'fine rep', keyboard
end
%'fine rep', keyboard

mesh_yd=mesh(pun); 
ragv=rag(pun,:);

xold=x;
Naold=Na;
Ndold=Nd;
d=d(pun);
% Loop sulle colonne
for iCol=1:2

x=xold(pun,:,iCol);

xt=[[0 0]; x];

dsu=diff(sum(xt,2));
ddud=diff(xt,[],2);
ddu=ddud(2:end);
fiTJ=find(reTJ==-500);
fiug=find(dsu(1:fiTJ(end))==0 & ddu(1:fiTJ(end))~=0);
xo=x;


%for kg=fiug'
% xi=xo(kg-1,:);
% xu=xo(kg,:);
% xd=d(kg+[-1 0]);
% dx=diff(xi);
% xinterf=xi(1)+dx*xd(1)/sum(xd);
% x(kg-1,2)=xinterf;
% x(kg,1)=xinterf;
% 'qi', keyboard
%end

%MM=[[1:length(dsu)]' x dsu ddu d];


%'sono x;', keyboard
xmat=ParDD.xm(pun,:);
%firag=find(ragv.*xmat<0);
firag=find(ragv.*xmat(:,1:size(ragv,2))<0);
raot=ragv(firag);

xdd=x;

klv=1:length(rep);  % indeces of the updated repe

Nadd=Naold(pun,:,iCol);
Nddd=Ndold(pun,:,iCol);

% Old values are overwritten by the updated versions

Na=Nadd;
Nd=Nddd;



mesh=mesh_yd;

layerStacking_LITHO


DentroPazzia

StrDD.x_dd{iCol}=fliplr(x_dd(indici,:,:));
StrDD.NA_dd{iCol}=fliplr(NA_dd(indici,:,:));
StrDD.ND_dd{iCol}=fliplr(ND_dd(indici,:,:));
if iCol==1
    StrDD.lab=lab(indici,:);
end

iCol
% 'dentro pazzia',keyboard

end



%'dopo pazzia',keyboard
% figure(1),hold on,plot(X),keyboard

StrDD.flgBTJ_lithographic=Str1D.flgBTJ_lithographic;


% Salvo gli indici e il numero delle BTJ ----------------------------------
if isfield(ParDD,'iNEGF_multiple')
    StrDD.nBTJ = nBTJ;
    % indexes not related to StrDD.NA_dd, StrDD.ND_dd and StrDD.lab!
    StrDD.iNEGF_multiple = iNEGF_multiple; 
end
% -------------------------------------------------------------------------

if length(raoADD)==3
    for iVac=1:length(material)
        if strcmp(material{iVac,end},'vacuum')==1
            material{iVac,end}='AlGaAs';
            break
        end
    end
end

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
