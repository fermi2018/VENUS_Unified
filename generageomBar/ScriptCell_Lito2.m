% Gestisco le due colonne, dentro e fuori


% Dentro


NCin=1
NCend=StrDD.NUcolTJ;

x_dd=StrDD.x_dd{1};
Na_dd=StrDD.NA_dd{1}/mode.CarrierNorm;
Nd_dd=StrDD.ND_dd{1}/mode.CarrierNorm;
d_dd=StrDD.d_dd;
mesh_dd=StrDD.mesh_dd;
lab=StrDD.lab;
material=StrDD.material;
if isfield(StrDD,'nBTJ')
    mode.nBTJ=StrDD.nBTJ;
end
 
mesh_r=StrDD.mesh_r;
ra_dd=StrDD.ra_dd;
ra_long=StrDD.raggi;

 NUMcol=size(material,2);
if  isfield(mode,'quasi1D') 
 if mode.quasi1D==1 %|| mode.dimReduction==1
  NUMcol=1;
 end
end

PreCell



if strcmp(mode.ContactPosition,'left')==1   % switch the contact labels
    for indTrans=1:NUMcol
        if(strcmp(material{indLong,indTrans},'Ground')==1)
            material{indLong,indTrans}='Line';
        elseif(strcmp(material{indLong,indTrans},'Line')==1)
            material{indLong,indTrans}='Ground';
        end
    end
end

% 'QUI prima', keyboard
Sub_CellNEW

% 'QUI 1 dopo', keyboard
% FUori

NCin=NCend+1;
NCend=NUMcol;

x_dd=StrDD.x_dd{2};
Na_dd=StrDD.NA_dd{2}/mode.CarrierNorm;
Nd_dd=StrDD.ND_dd{2}/mode.CarrierNorm;

PreCell

% 'QUI 2 prima', keyboard
    
Sub_CellNEW


% 'QUI 2 dopo', keyboard
% 
% Tutte le colonne

% Top contact
indLong=size(material,1);
xini=0;

if strcmp(mode.ContactPosition,'left')==1   % switch the contact labels
    for indTrans=1:NUMcol
        if(strcmp(material{indLong,indTrans},'Ground')==1)
            material{indLong,indTrans}='Line';
        elseif(strcmp(material{indLong,indTrans},'Line')==1)
            material{indLong,indTrans}='Ground';
        end
    end
end

for indTrans=1:NUMcol;%size(material,2)
    if(strcmp(material{indLong,indTrans},'Ground')==1)
        mat='Au';
        CellInfo(indLong,indTrans)=f_SaveCellInfo('ground',mat,0,0,1,0,'C',xini,thickcolumn(indTrans),CellInfo(indLong-1,indTrans).yfin,tContact,meshcolumn(indTrans),1,1,heterojunctionth);
    elseif(strcmp(material{indLong,indTrans},'Line')==1)
        mat='Au';
        if(indTrans==size(material,2)-flagPassive)
            CellInfo(indLong,indTrans)=f_SaveCellInfo('line',mat,0,0,1,0,'C',xini,thickcolumn(indTrans),CellInfo(indLong-1,indTrans).yfin,tContact,meshcolumn(indTrans),1,1,heterojunctionth);
        else
            CellInfo(indLong,indTrans)=f_SaveCellInfo('extraline',mat,0,0,1,0,'C',xini,thickcolumn(indTrans),CellInfo(indLong-1,indTrans).yfin,tContact,meshcolumn(indTrans),1,1,heterojunctionth);
        end
    else
        mat=material{indLong,indTrans};
        CellInfo(indLong,indTrans)=f_SaveCellInfo('Layer_top',mat,0,0,1,0,'C',xini,thickcolumn(indTrans),CellInfo(indLong-1,indTrans).yfin,tContact,meshcolumn(indTrans),1,1,heterojunctionth);
    end
    xini = xini+thickcolumn(indTrans);
end


% Check: if cell's material is polyamide, DO NOT tag it as qw!!!
for indLong=2:size(CellInfo,1)-1
    for indTrans=1:NUMcol;%size(CellInfo,2)
        is_QW_Cond = not(isempty(strfind(CellInfo(indLong,indTrans).label,'qw')));
        is_Passiv_Cond = not(isempty(strfind(CellInfo(indLong,indTrans).material,'Polyamide')));
        if(is_QW_Cond & is_Passiv_Cond)
            indQW=str2num(CellInfo(indLong,indTrans).label(3));
            CellInfo(indLong,indTrans).label = ['PassivationQW',num2str(indQW)];
        end
    end
end

if NUMcol==1
    if strcmp(mode.ContactPosition,'left')==1   % switch the contact labels
        CellInfo(1).label='line';
    else
        CellInfo(end).label='line';
    end
end 

if mode.flgBTJ_lithographic==2     
    for iL=size(CellInfo,1):-1:1
        for iT=1:NUMcol
            if strcmp(CellInfo(iL,iT).material,'AlOx')==1
                if mode.etch.growth==1
                    CellInfo(iL,iT).material=CellInfo(iL+1,iT).material;
                    CellInfo(iL,iT).dtype=CellInfo(iL+1,iT).dtype;
                    CellInfo(iL,iT).dgvet(3)=CellInfo(iL+1,iT).dgvet(3);
                    CellInfo(iL,iT).gvet(3)=CellInfo(iL+1,iT).gvet(3);
                else
                    CellInfo(iL,iT).material=mode.etch.mat;
                    CellInfo(iL,iT).dtype=mode.etch.dtype;
                    CellInfo(iL,iT).dgvet(3)=mode.etch.dop;
                    CellInfo(iL,iT).gvet(3)=mode.etch.xm;
                end
            end
        end
    end
end

%'FINE CellNEW; ', keyboard