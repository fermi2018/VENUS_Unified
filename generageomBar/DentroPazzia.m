if flgStop==1
'qui basta MATERIAL', keyboard
end
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

if flgStop==1
'FUBE ', keyboard
end

reprep=[];
k=1;
kc=1;
fiCA=[];
while k<=length(rep)
    
%    rei=rep(k,1);
    rei=reTJ(k,1);
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
        if rei==-500
            reprep=[reprep; rei];
        else
            reprep=[reprep; -1];
        end
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

if flgStop==1
    'reprep', keyboard
    'reprep', keyboard
end

% looks for imaginary reps, corresponding to QWs
lab=Dch;
fiq=find(reprep==1i);
finq=find(reprep>=1);
fiTJ=find(reprep==-500);
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


% The following was the way it was done when just 1 BTJ was included in the
% simulation

if flgStop==1
 'qui lab', keyboard
end
% if isfield(ParDD,'iNEGF')==1
%     if ~isempty(mireq.n2.Nmesh)
%         iNEGF=ParDD.iNEGF + 1 - (Np_n2dd*NL_n2dd);
%     else
%         iNEGF=ParDD.iNEGF;
%     end
%     for kn=1:length(iNEGF)
%         LA=lab(iNEGF(kn),:);
%         if kn==1
%             ib=ib+1;
%             lab(iNEGF(kn),1:length(LA)+4)=[LA,'_TJ',num2str(ib)];
%         else
%             ib=ib+1;
%             LA(end-3:end)=['_TJ',num2str(ib)];
%             lab(iNEGF(kn),:)=LA;
%         end
%     end
% end

% Now, the following is Torrelli procedure
if isfield(ParDD,'iNEGF_multiple')==1 && iCol==1 % if there are BTJs in the structure
    %
    clear iNEGF_multiple
    iNEGF_multiple=ParDD.iNEGF_multiple;
    nBTJ=ParDD.nBTJ; % total number of BTJs within the  VCSEL
    %    iNEGF_multiple=ParDD.iNEGF_multiple; % cell containing the position of the BTJs
    if nBTJ>1
        fprintf('da sistemare con nBTJ>1\n'),keyboard
    end
    
    for indexBTJ=1:nBTJ
        nLayers=length(ParDD.iNEGF_multiple{indexBTJ});
        %        nLayers=length(iNEGF_multiple{indexBTJ});
        
        if isfield(mireq,'n2') && ~isempty(mireq.n2.Nmesh)
            %            iNEGF_multiple{indexBTJ}(:)=iNEGF_multiple{indexBTJ}(:)+1-(Np_n2dd*NL_n2dd);
%             iNEGF_multiple{indexBTJ}(:)=fiTJ(1:nLayers)';
%            iNEGF_multiple{indexBTJ}(:)=fiTJ(1:nLayers-1)'-1;
            iNEGF_multiple{indexBTJ}(:)=fiTJ';
        elseif isfield(mireq,'p') && ~isempty(mireq.p.Nmesh)
			'TJ con specchio p, ATTENZIONE', keyboard
%            iNEGF_multiple{indexBTJ}(:)=iNEGF_multiple{indexBTJ}(:)+1-(Np_pdd*NL_pdd);
        end
        
%        for indexLayer=1:nLayers-1
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



dd=reshape(ddd',numel(ddd),1);
xx=reshape(xxd',numel(ddd),1);
na=reshape(nad',numel(ddd),1);
nd=reshape(ndd',numel(ddd),1);

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
if exist('Ks')
%if iOxTutti==0
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

if flgStop==1
'Ks', keyboard
end

if exist('Ks')
    fi=find(Rag(Ks,1)>0);
    Rox=min(Rag(Ks(fi),1));
else
    Rox=Rag(fiTJ(1));
end

R_contact=Rag(1,:);

if Str1D.quasi1D==1
    ParDD.ragAdd=[0 0];
end
raoADD=Rox+[ParDD.ragAdd(1) 0 ParDD.ragAdd(2)];

NUcolTJ=1+length(find(ParDD.ragAdd<0));

StrDD.NUcolTJ=NUcolTJ;


%'RAdd', keyboard


fira=find(diff([0 raoADD])>0);
raoADD=raoADD(fira);

fiO=find(raver==Rox);
fiOx=find(raoADD==Rox);

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

ddc=reshape(ddd',numel(ddd),1);
xxc=reshape(xxd',numel(ddd),1);
nac=reshape(nad',numel(ddd),1);
ndc=reshape(ndd',numel(ddd),1);

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
if Rox>0 & exist('Ks')
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
'qui ROTIAOP new', keyboard
end

inew=1;
if inew==0
    Script_mate
else
    Script_mate2Pazzia
end


fiM=find(D1_i>0);
NA_dd=flipud(NA1_i(fiM,:,:));
ND_dd=flipud(ND1_i(fiM,:,:));
d_dd=flipud(D1_i(fiM));
mesh_dd=flipud(Mesh1_i(fiM));
x_dd=flipud(X1_i(fiM,:,:));
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
