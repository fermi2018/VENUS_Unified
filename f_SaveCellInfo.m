
function [CellInfo]=f_SaveCellInfo(label,material,molarini,molarfin,dopini,dopfin,doptype,xini,thick_x,yini,thick_y,mesh_x,mesh_y,heteropoint,heterospace)
%-- CellInfo: per ciascun rettangolo
%-- Estensione grading 2d: molarini (punto in basso a sinistra nella
%cella), molarfin_x, molarfiny: così ho il triangolo e posso calcolare
%a,b,c !!!! =)
% heteropoint: se 1 non si mettono punti aggiuntivi all'interfaccia; se 2,
% si mette un punto solo in fondo al layer corrente (a destra), appena
% prima della PROSSIMA interfaccia. se 3, si mettono 2 punti (EVITARE)
% ha senso quindi mettere 1 se non ci sono eterointerfacce, 2 se il layer
% SUCCESSIVO ha un cambio di frazione molare o di materiale.

CellInfo.label=label;
CellInfo.material=material;
CellInfo.mesh_x=mesh_x;
CellInfo.mesh_y=mesh_y;
CellInfo.thick_x=thick_x;
CellInfo.thick_y=thick_y;
CellInfo.xini=xini;
CellInfo.xfin=xini+thick_x;
CellInfo.yini=yini;
CellInfo.yfin=yini+thick_y;
CellInfo.hpoint=heteropoint;
CellInfo.hthick=heterospace;

%%% molar fraction
a=0;
b=(molarfin-molarini)./thick_y;
c=molarini-(molarfin-molarini)./thick_y.*yini;
CellInfo.gvet=[a,b,c];

%%% doping
if strcmp(doptype,'C')
    d=0;
    e=dopini;
    f=dopfin;
else
    d=0;
    e=(dopfin-dopini)./thick_y;
    f=dopini-(dopfin-dopini)./thick_y.*yini;
end

CellInfo.dgvet=[d,e,f];
CellInfo.dtype=doptype;

return