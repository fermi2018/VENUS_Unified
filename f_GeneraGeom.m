
function f_GeneraGeom(CellInfo,geom,StrTT,StrDD,ParVet,structure)
% clear
% close all
% clc
%%%%%%%%%%%%%%%%%%%%%%%%% Dati in ingresso %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Read data structure from 'filename.txt' build with scrividati8
%%%% Produce geom_filename.mat
%%%% Support graded molar fraction region
%%%% Support 2D structure

%%%%%%%%%%%%%%%%%%%%%%%%%% Dati in uscita %%%%%%%%%%%%%%%%
% Slayer.dopd=[0,0,0,1e18];
% Slayer.dopa=[0,0,1e18,0];
% Slayer.label={'line','ground','p-layer','n-layer'};
% Slayer.mater={'Au','Au','Si','Ge'};
% Slayer.cont=[1 1];
% Slayer.workf=[0 0];
% Slayer.X=[1 2 3 4 4 3 2 1; 10 11 12 13 13 12 11 10; 4 5 6 7 7 6 5 4; 7 8 9 10 10 9 8 7];
% Slayer.Y=[15 15 15 15 1 1 1 1; 15 15 15 15 1 1 1 1; 15 15 15 15 1 1 1 1; 15 15 15 15 1 1 1 1];
% Slayer.XP={'-t' '-t+x' '-x' '0' 'x' 't_p-x' 't_p' 't_p+x' 't_n+t_p-x' 't_n+t_p' 't_n+t_p+x' 't_n+t_p+t-x' 't_n+t_p+t'};
% Slayer.YP={'0' '' '' '' '' '' '' '' '' '' '' '' '' '' 'b'};
% Slayer.par={'t=0.01e-4;' 'b=0.05e-4;' 'x=1e-8;'};
% Slayer.divx=[1 1 1 1 100 1 1 100 1 1 1 1];
% Slayer.divy=[1];
% Slayer.ethj=[];
% Slayer.reg2=[1 2 0 0];
% Slayer.neum=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
% Slayer.grid=[13 15]
%
% Created by Marco Calciati, last revision 1-07-2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

%'in Genera;', keyboard
%%%% Set number of columns
ncs=size(CellInfo,2);
for indCol=1:ncs
    %%%% Set columns thickness
    Slayer.colthick(indCol)=CellInfo(1,indCol).thick_x;
    %%%% Set columns mesh
    Slayer.colmesh(indCol)=CellInfo(1,indCol).mesh_x;
end

%%%% Set number of layers
nls=numel(CellInfo);
nls2D=nls/ncs;

%%%% Set numero massimo punti [H]
nxp=1;
for i=1:nls2D
nxp=nxp+CellInfo(i,1).hpoint;
Slayer.hpoint(i)=CellInfo(i,1).hpoint;
Slayer.hthick(i)=CellInfo(i,1).hthick;
end

%%%% Set altezza colonna geometrica
hei=6;

%%%% set grid [H]
grid=[hei*ncs nxp];

Slayer.grid=grid;

%%%% Define draw paramteres
draw.x=[0 0]; draw.y=[0 0]; % X- and Y-coordinates of boundary segments
draw.X=[]; draw.Y=[]; % X- and Y-coordinates of patch vertices

%%%% Define geom
nx=Slayer.grid(1); ny=Slayer.grid(2); % set geom.grid
nLut=geom.GLUTm;
geom=struct('gd',[],'dl',[],'dgm',[],'label',cell(1,1),'material',cell(1,1),'color',[],'par',cell(1,1), ...
            'grid',[nx ny],'X',cell(1,1),'Y',cell(1,1),'bspec',[],'bspmc',[],'symflg',[], ...
            'electrode',[],'semiconductor',[],'reg2contact',[],'contact_type',[],'workfun',[], ...
            'ethj',[],'ethj_sd',[],'ethjflg',0, ...
            'div_x',[],'div_y',[],'mate',[],'function_dop',[], ...
            'edge',[],'nbs',[],'nd',[],'nc',[],'nm',[]);
geom.X(1:nx,1)=cellstr(''); geom.Y(1:ny,1)=cellstr('');
geom.GLUTm=nLut;

%%%%%%%%%%%%% DrawMode Set geom.gd %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Set X [H]
 t=0;
for i=1:nls2D %% number of layers
    for j=1:Slayer.hpoint(i)+1
        y(i,j)=t+j;
        xone(i,j)=1;
    end
    t=y(i,j)-1;
end

 y(all(y==0,2),:)=[]; %%% delete zeroes rows

 y=[y fliplr(y)]; %%%% mirror rows
 
 for i=1:nls2D
     yzero=find(y(i,:)==0);
     if(not(isempty(yzero)))
         indy=yzero(1);
         for j=1:(indy-1)
             y(i,indy) = y(i,indy+length(yzero));
             y(i,indy+length(yzero))=0;
             indy=indy+1;
         end
         celly{i}=y(i,1:indy-1);
     else
         celly{i}=y(i,:);
     end
 end

 for j=1:ncs-1
     for i=1:nls2D
     celly{i+nls2D*j}=celly{i};
     end
 end

Slayer.Y=celly;

%%%% set Y [H]
%es y=[1 1 1 1 11 11 11 11; 1 1 1 1 11 11 11 11; 1 1 1 1 11 11 11 11; 1 1 1 1 11 11 11 11; 11 11 11 11 20 20 20 20; 11 11 11 11 20 20 20 20; 11 11 11 11 20 20 20 20; 11 11 11 11 20 20 20 20];

uno=xone;
xt=xone*hei;

for j=1:ncs-1
  xone=[xone;hei*j*uno];
  xt=[xt;hei*(j+1)*uno];
end

xone=[xone,xt];

for i=1:nls2D*ncs
     yzero=find(xone(i,:)==0);
     if(not(isempty(yzero)))
         indy=yzero(1);
         for j=1:(indy-1)
             xone(i,indy) = xone(i,indy+length(yzero)/2);
             xone(i,indy+length(yzero)/2)=0;
             indy=indy+1;
         end
         cellx{i}=xone(i,1:indy-1);
     else
         cellx{i}=xone(i,:);
     end
end
Slayer.X=cellx;

for i=1:length(Slayer.X)
    X=Slayer.X{i}; Y=Slayer.Y{i};
    temp=[2 length(X) X Y]';
    geom.gd(1:size(temp,1),size(geom.gd,2)+1)=temp;
end
disp(' ')
disp('Generating geometry file');
disp(' ')
%keyboard

tic
[geom.dl,geom.bt]=decsg(geom.gd);
toc
nsd=max([geom.dl(6,:) geom.dl(7,:)]);
geom.label=cell(nsd,1);
geom.label(1:nsd)=cellstr('');
geom.material=cell(nsd,1);
geom.material(1:nsd)=cellstr('');
geom.mate=zeros(1,nsd);
geom.color=0.95*ones(nsd,3);
I1=ismember(geom.dl(6,:),0); I2=ismember(geom.dl(7,:),0);
geom.bspec=[]; geom.bspmc=find(I1|I2);

%%%%%%%%%%%%%%%%%%%%% Physical Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%da getPhMat Set mater label dopa dopd  %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[pp,ee,tt]=initmesh(geom.dl,'hmax',Inf,'MesherVersion','R2013a');
in1=tt(1,:); in2=tt(2,:); in3=tt(3,:);
x1=pp(1,in1) ;x2=pp(1,in2); x3=pp(1,in3); y1=pp(2,in1); y2=pp(2,in2); y3=pp(2,in3);
r12x=x2-x1; r23x=x3-x2; r31x=x1-x3; r12y=-y2+y1; r23y=-y3+y2; r31y=-y1+y3;
Ae=abs(r31x.*r23y-r31y.*r23x)./2; % triangle area
nbs=size(geom.dl,2);
nsd=max([geom.dl(6,:) geom.dl(7,:)]);
for ii=1:nsd,
    tri=find(tt(4,:)==ii);
    Amax=find(Ae(tri)==max(Ae(tri))); Amax=Amax(1);
    xm=mean(pp(1,tt(1:3,tri(Amax)))); ym=mean(pp(2,tt(1:3,tri(Amax))));
    xdata(ii)=xm;
    ydata(ii)=ym;
end

% [B,orderx]=sort(xdata); %%% ordina x
% ydata=ydata(orderx);    %%% ordina y
% 
% for i=1:ncs:nls2D*ncs
%     [ydata(i:i+ncs-1),Iy]=sort(ydata(i:i+ncs-1));
%     L=orderx(i:i+ncs-1);
%     orderx(i:i+ncs-1)=L(Iy);
% end
% 
% %% structure ordex
% %% 2 4 6 8
% %% 1 3 5 7
% %% structure data
% %% 1 2 3 4 5 6 7 8
% %% stucture output
% %% 1 3 5 7 2 4 6 8
% 
% %% structure ordex
% %% 3 6 9 12      3 6 9 12 15 18
% %% 2 5 8 11      2 5 8 11 14 17
% %% 1 4 7 10      1 4 7 10 13 16
% %% structure data
% %% 1 2 3 4 5 6 7 8 9 10 11 12
% %% stucture output
% %% 1 4 7 10 2 5 8 11 3 6 9 12
% 
% for j=1:ncs
%     for i=1:nls2D
%         orderx2(i+nls2D*(j-1))=orderx(ncs*i-(ncs-1)+(j-1));
%     end
% end
% orderx=orderx2;

[B,ordery]=sort(ydata); %%% ordina y
xdata=xdata(ordery);    %%% ordina x

for i=1:ncs:nls2D*ncs
    [xdata(i:i+ncs-1),Iy]=sort(xdata(i:i+ncs-1));
    L=ordery(i:i+ncs-1);
    ordery(i:i+ncs-1)=L(Iy);
end

for j=1:ncs
    for i=1:nls2D
        ordery2(i+nls2D*(j-1))=ordery(ncs*i-(ncs-1)+(j-1));
    end
end
ordery=ordery2;

indCell=0;
for indCol=1:ncs
    for indLay=1:nls2D
        indCell=indCell+1;
        %-- Sorted
        Slayer.label{ordery(indCell)}=CellInfo(indLay,indCol).label;
        Slayer.mater{ordery(indCell)}=CellInfo(indLay,indCol).material;
        Slayer.dtype(ordery(indCell))=CellInfo(indLay,indCol).dtype;
        gvet{ordery(indCell)}=CellInfo(indLay,indCol).gvet;
        dgvet{ordery(indCell)}=CellInfo(indLay,indCol).dgvet;
        dtype{ordery(indCell)}=CellInfo(indLay,indCol).dtype;
        %-- NOT Sorted
        Slayer.thickness(indCell,1)=CellInfo(indLay,indCol).thick_y; %% set layer thickness
        Slayer.tempmesh(indCell,1)=CellInfo(indLay,indCol).mesh_y;  %% set layer mesh
    end
end

geom.label=Slayer.label';         % set geom label
geom.material=Slayer.mater';      % set geom.material

%%%% find contacts position
% L'ordine dei materiali influenza:
% label, material, matercom, dopd, dopa, thickness, divx, divy, cont,
% workf, neum, reg2, par, XP
% L'ordine dei materiali non influenza:
% X, Y, YP, grid
lnp=1;
for i=1:nls
    if(strcmp(geom.label(i),'line'))
        linep(lnp)=i; %%% position 'line' contact
        lnp=lnp+1;
    end
end
lnp=lnp-1; %%% numero di layer 'line'

lng=1;
for i=1:nls
    if(strcmp(geom.label(i),'ground'))
        groundp(lng)=i; % position 'ground' contact
        lng=lng+1;
    end
end
lng=lng-1; %%% numero di layer 'ground'

lnel=1;
for i=1:nls
    if(strcmp(geom.label(i),'extraline'))
        extralinep(lnel)=i; %%%position 'ground' contact
        lnel=lnel+1;
    end
end
lnel=lnel-1; %%% numero di layer 'ground'

%%%% Set contacts type
Slayer.cont(1)=dgvet{linep(1)}(2);
Slayer.cont(2)=dgvet{groundp(1)}(2);
if(exist('extralinep') && not(isempty(extralinep))) % if extra-line has been used...
    Slayer.cont(3)=dgvet{extralinep(1)}(2);
end    

%%%% Set contacts work function
Slayer.workf(1)=dgvet{linep(1)}(3);
Slayer.workf(2)=dgvet{groundp(1)}(3);
if(exist('extralinep') && not(isempty(extralinep))) % if extra-line has been used...
    Slayer.workf(3)=dgvet{extralinep(1)}(3);
end    

%%%% Set doping for Contacts
for i=1:lnp
    dgvet{linep(i)}(2)=0;
    dgvet{linep(i)}(3)=0;
end
for i=1:lng
    dgvet{groundp(i)}(2)=0;
    dgvet{groundp(i)}(3)=0;
end
for i=1:lnel
    dgvet{extralinep(i)}(2)=0;
    dgvet{extralinep(i)}(3)=0;
end

%%%% Set geom color electrode semiconductor
% 'qui[', keyboard
for sdl=1:nls
    str=char(geom.material(sdl));
    [geom.semiconductor(sdl),geom.electrode(sdl),geom.color(sdl,:)]=f_GraphicLibrary(str);
end

%%%% set XP
% es: xp={'-t';'-t+x';'-x';'0';'x';'t1-x';'t1';'t1+x';'t2-x';'t2';'t2+x';'t2+t-x';'t2+t'};

xp={'t1'}; %%% primo punto
indt=2;
indx=1;
for i=1:nls2D
    if (Slayer.hpoint(i)==1)
        xp(indt)={['t' num2str(i+1)]};
        indt=indt+1;
    elseif(Slayer.hpoint(i)==2)
        if(Slayer.hpoint(i-1)<3)
            xp(indt)={['t' num2str(i+1) '-x' num2str(indx)]};
            xp(indt+1)={['t' num2str(i+1)]};
            indt=indt+2; indx=indx+1;
        else
            xp(indt)={['t' num2str(i) '+x' num2str(indx)]};
            xp(indt+1)={['t' num2str(i+1)]};
            indt=indt+2; indx=indx+1;
        end
    else
        xp(indt)={['t' num2str(i) '+x' num2str(indx)]};
        xp(indt+1)={['t' num2str(i+1) '-x' num2str(indx)]};
        xp(indt+2)={['t' num2str(i+1)]};
        indt=indt+3; indx=indx+1;
    end
end

geom.Y=xp';

%%%% set YP
%es: yp={'0' '' '' '' '' '' '' '' '' '' 'b1' '' '' '' '' '' '' '' '' 'b2'};
yp={'0'}; %%% primo punto
for i=2:hei*ncs
    yp=[yp {''}];
end
bisn=1;
for i=hei:hei:ncs*hei
    yp{i}=strcat('b', num2str(bisn));
    bisn=bisn+1;
end

geom.X=yp';

%%%%%%%%%%%%%%%%%%%%% Boundary Conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% set boundary condiction (da fare)
nneum=(nxp-1)*(ncs+1)+(nls2D+1)*ncs;
neum=zeros(1,nneum);
Slayer.neum=neum;

dl=geom.dl;
I1=ismember(dl(6,:),0); I2=ismember(dl(7,:),0); I=I1|I2;
nbs=size(geom.dl,2);
bcvet=ones(1,nbs);
for bs=1:nbs
    if(I(bs))
        value=Slayer.neum(bs);
        if(value), bcvet(bs)=2; % Neumann
        else, bcvet(bs)=3; % Dirichlet
        end
    end
end
geom.bspec=find(bcvet==2);
geom.bspmc=find(bcvet==3);

%%%%%%%%%%%% geom.par
%%%new [H]
par={['t1=' num2str(-Slayer.thickness(1)) ';' 'b1=' num2str(Slayer.colthick(1)) ';']};
for i=2:nls2D+1
    par=[par; {strcat('t', num2str(i), '=t', num2str(i-1), '+', num2str(Slayer.thickness(i-1),15), ';')}];
end
for i=2:ncs
    par=[par; {strcat('b', num2str(i), '=b' , num2str(i-1), '+', num2str(Slayer.colthick(i),15), ';')}];
end
indx=1;
for i=1:nls2D
    if (Slayer.hpoint(i)>1)
        par=[par; {strcat('x', num2str(indx), '=' , num2str(Slayer.hthick(i),15), ';')}];
        indx=indx+1;
    end
end

geom.par=par;
geom.gvet=gvet;
geom.dgvet=dgvet;
geom.dtype=dtype;
%%%%%%%%%%%%%%%%%% computedgm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dl=geom.dl;
eval(char(geom.par)')
for nn=1:size(dl,2),
    dl(2,nn)=eval(char(geom.X(round(dl(2,nn)))));
    dl(3,nn)=eval(char(geom.X(round(dl(3,nn)))));
    dl(4,nn)=eval(char(geom.Y(round(dl(4,nn)))));
    dl(5,nn)=eval(char(geom.Y(round(dl(5,nn)))));
end
geom.dgm=dl;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% GEOMODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

geom.nbs=size(geom.dl,2); % nbr. of boundary segments
geom.nd=length(geom.material); % nbr. of regions
if(isfield(geom,'reg2contact') && isempty(geom.reg2contact))
    geom.reg2contact=zeros(1,geom.nd);
    line=find(strcmp('line',geom.label)); ground=find(strcmp('ground',geom.label));
    geom.nm=length(line); geom.nc=geom.nm; % number of active electrodes
    geom.reg2contact(line)=1:geom.nm;
    if(not(isempty(ground))), geom.reg2contact(ground)=geom.nm+1; geom.nc=geom.nc+1; end
    geom.workfun=zeros(1,geom.nc); geom.contact_type=ones(1,geom.nc); end
%
electrode=find(geom.electrode); % all electrodes
semiconductor=find(geom.semiconductor); % find semiconductor regions
I5=ismember(geom.dgm(6,:),0); I6=ismember(geom.dgm(7,:),0);
I7=ismember(geom.dgm(6,:),electrode); I8=ismember(geom.dgm(7,:),electrode);
I9=ismember(geom.dgm(6,:),semiconductor); I10=ismember(geom.dgm(7,:),semiconductor);
% set Monte Carlo boundary conditions
if(isfield(geom,'edge') && isempty(geom.edge))
    geom.edge=zeros(1,geom.nbs);
    geom.edge(xor(I9,I10))=1;
    geom.edge(and(or(I7,I8),or(I9,I10)))=2;
end
if(isfield(geom,'div_x') && isempty(geom.div_x)),
    geom.div_x=ones(1,length(find(not(cellfun('isempty',geom.X))))-1);
end
if(isfield(geom,'div_y') && isempty(geom.div_y)),
    geom.div_y=ones(1,length(find(not(cellfun('isempty',geom.Y))))-1);
end
% set symflg
if(isfield(geom,'symflg') && isempty(geom.symflg)), geom.symflg=0; end

%%%%%%%%%%%%%%%%%%% Contacts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% set mesh divx
% ATTENZIONE: Si aspetta una struttura 'line n1 n2 nN ground' NOT ORDERED
indt=1;
for i=1:nls2D
    if (Slayer.hpoint(i)==1)
        mshx(indt)=Slayer.tempmesh(i);
        indt=indt+1;
    elseif(Slayer.hpoint(i)==2)
        if(Slayer.hpoint(i-1)<3)
            mshx(indt)=Slayer.tempmesh(i);
            mshx(indt+1)=1;
            indt=indt+2;
        else
            mshx(indt)=1;
            mshx(indt+1)=Slayer.tempmesh(i);
            indt=indt+2;
        end
    else
        mshx(indt)=1;
        mshx(indt+1)=Slayer.tempmesh(i);
        mshx(indt+2)=1;
        indt=indt+3;
    end
end

geom.div_y=mshx;

%%%% set mesh divy
geom.div_x=Slayer.colmesh;

%%%% set etheroj (da fare)
geom.ethj=[];

%%%% set region 2 contact
reg2=zeros(1,nls);
for i=1:lnp
    reg2(linep(i))=1;
end
for i=1:lng
    reg2(groundp(i))=2;
end
for i=1:lnel
    reg2(extralinep(i))=1; % setting extraline as line contacts...
end

geom.reg2contact=reg2;

%%%% Set contacts
geom.contact_type=Slayer.cont;

%%%% Set work function
geom.workfun=Slayer.workf;

lab=StrDD.lab;

x_ddC=StrDD.x_dd;
if iscell(x_ddC)
x_dd=x_ddC{1};
else
x_dd=x_ddC;
end



%%%% Set quantum well molar fraction
for indLayer=1:size(lab,1)
    strtmp=lab(indLayer,:);
    VerQW=strfind(strtmp,'qw'); % looking for QW
    if(not(isempty(VerQW)))
        indMQW=str2double(strtmp(3));
        geom.QWxmol{indMQW}=x_dd(indLayer,1);
    end
end

assignin('base','geom',geom)

toc

% clear geom StrTT temp Slayer CellInfo x y celly X Y cellx
%clear geom x y Slayer
