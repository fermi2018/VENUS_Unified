contStr=1;

igrigliaSTOP=IPLOT;
Gs.QWorientation='X'; % orient. quantum well (parallel to x (rho): VCSEL)
%%% open file
structure=structureName;

iplot=0; % 1 for plotting VELM str file parameters
iBulk=0; % if 0 the full quantum model is adopted
iHetero=0; % if 1 puts the "heteropoints"; if 0 doesn't
Str1D.flgBTJ_lithographic=mode.flgBTJ_lithographic;

if isfield(mode,'quasi1D')
    Str1D.quasi1D=mode.quasi1D;
else
    Str1D=[];
end

if isfield(mode,'stairs')
    Str1D.stairs=mode.stairs;
else
    Str1D.stairs=0;
end

global flgStop

% flgStop=1;
flgStop=0;

if flgStop==1
    'entro StrRead ', keyboard
end

if mode.flgBTJ_lithographic<=1
    [StrDD,StrTT,ParOpt]=StrRead(ParVet,structure,iplot,Str1D);
else    
    [StrDD,StrTT,ParOpt]=StrRead_Litho(ParVet,structure,iplot,Str1D);
    VelmOptions.PaTJ=ParOpt.radii.PaTJ;
end

StrDD.LITHO=mode.flgBTJ_lithographic;

if flgStop==1
    'dopo StrRead ', keyboard
    'dopo StrRead ', keyboard
    'dopo StrRead ', keyboard
    'dopo StrRead ', keyboard
end
%%%%% TIBALDI OGGI %%%%%

%StrDD.r_passiv=5;
%StrDD.mesh_passiv=10;
% StrDD.r_passiv=[];
% StrDD.mesh_passiv=[];

r_passiv=abs(StrDD.r_passiv);
mesh_passiv=StrDD.mesh_passiv;
%%%%%%%%%%%%%%%%%

ParMore=StrDD.ParMore;


% tContact=0.005e-4; % contact thickness, cm
tContact=0.2e-4; % contact thickness, cm
% tContact=0.0001e-4; % contact thickness, cm (RESISTOR)
heteromode=mode.heteromode; % 1 is no points, 2 is 1 point per interface, 3 is 2 points per interface
heterojunctionth=1e-9; % if heteromode is on, defines the distance of the interface points


iQW=0;
clear CellInfo QWInfo



if flgStop==1
    'inizio CellInfo;', keyboard
end

% Bottom contact

if isfield(StrDD,'flgBTJ_lithographic')
    if StrDD.flgBTJ_lithographic==2
        ScriptCell_Lito2
    else
        indLong=1;
        xini=0;
        PreCellOld
        ScriptCell_Lito1
    end
else
    indLong=1;
    xini=0;
    PreCellOx
%     'i material;', keyboard
%     material=material(:,1:end-1);
    ScriptCell_oxide
end


if contStr==1 && flgStop==1
    'Genera', keyboard
end

f_GeneraGeom(CellInfo,geom,StrTT,StrDD,ParVet,structure)

if isfield(StrDD,'iNEGF')
    geom.iNEGF=1;
end

fis= strfind(structure,'\');
strName=structure(fis(end)+1:end);
DirName=structure(1:fis(end));
%%%% save geom
save([DirName,'geom_' strName],'geom','StrTT','ParVet','StrDD')
