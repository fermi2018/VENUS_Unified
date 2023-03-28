
% ' inizio Cell Old', keyboard


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
        CellInfo(indLong,indTrans)=f_SaveCellInfo('ground',mat,0,0,1,0,'C',xini,thickcolumn(indTrans),-tContact,tContact,meshcolumn(indTrans),1,1,heterojunctionth);
    elseif(strcmp(material{indLong,indTrans},'Line')==1)
        mat='Au';
        if(indTrans==size(material,2)-flagPassive)
            CellInfo(indLong,indTrans)=f_SaveCellInfo('line',mat,0,0,1,0,'C',xini,thickcolumn(indTrans),-tContact,tContact,meshcolumn(indTrans),1,1,heterojunctionth);
        else
            CellInfo(indLong,indTrans)=f_SaveCellInfo('extraline',mat,0,0,1,0,'C',xini,thickcolumn(indTrans),-tContact,tContact,meshcolumn(indTrans),1,1,heterojunctionth);
        end
    else
        mat=material{indLong,indTrans};
        CellInfo(indLong,indTrans)=f_SaveCellInfo('Layer_bottom',mat,0,0,1,0,'C',xini,thickcolumn(indTrans),-tContact,tContact,meshcolumn(indTrans),1,1,heterojunctionth);
    end
    xini = xini+thickcolumn(indTrans);
end

% ' All the other layers', keyboard

for indLong=2:size(material,1)-1 % loop on longitudinal (y) direction
    xini=0; % initializing x of the transverse layer
    
    label=lab(indLong-1,:); % label of the layer
    xmol=x_dd(indLong-1,:); % initial and ending molar fraction
    
    dop_a=Na_dd(indLong-1,:); % initial and ending acceptor doping
    dop_d=Nd_dd(indLong-1,:); % initial and ending donor doping
    dopType='N'; % doping type
    dop=dop_d;
    if(max(dop_a)>max(dop_d)) % if we have more acceptor than donor
        dopType='P';
        dop=dop_a;
    end
    
    if(iBulk==0) % if "bulk mode" is deactivated, i.e., quantum model on
        QWcond=strfind(label,'qw');
        if(QWcond==1)
            iQW=iQW+1;
            QWInfo(iQW).xmol=xmol(1); % saving xmol for QW
            xmol=x_dd(indLong-2,:); % modifying QW cmol
        end
    end
    
    % treatment of heterojunctions, if necessary
    hetero=heteromode;
    if(iHetero==1)
        if(iBulk==0)
            if(indLong>2 & isempty(QWcond))
                if(abs(x_dd(indLong-2,2)-x_dd(indLong-1,1))>1e-8)
                    labforLastQWCheck=lab(indLong-2,:);
                    if(isempty(strfind(labforLastQWCheck,'qw')))
                        CellInfo(indLong-1,indTrans).hpoint=2;
                    end
                end
            end
        elseif(indLong>2 & iBulk==1)
            if(abs(x_dd(indLong-2,2)-x_dd(indLong-1,1))>1e-8)
                labforLastQWCheck=lab(indLong-2,:);
                if(isempty(strfind(labforLastQWCheck,'qw')))
                    CellInfo(indLong-1,indTrans).hpoint=2;
                end
            end
        end
    end
    
    if(iBulk==1)
        label=lab(indLong-1,:);
    elseif(iBulk==0)
        if(isempty(QWcond));
            label=lab(indLong-1,:);
        end
    end
    
%' All the other layers NU', keyboard	
	
    for indTrans=1:NUMcol;%size(material,2) % loop on transversal (x) direction
        mat=material{indLong,indTrans};
%         if indLong==size(material,1)-1
%             CellInfo(indLong,indTrans)=f_SaveCellInfo(label,mat,xmol(1),xmol(2),dop(1),dop(2),dopType,xini,thickcolumn(indTrans),...
%                 CellInfo(indLong-1,indTrans).yfin,d_dd(indLong-1)*1e-7,meshcolumn(indTrans),mesh_dd(indLong-1),1,heterojunctionth);
%         else
        
            CellInfo(indLong,indTrans)=f_SaveCellInfo(label,mat,xmol(1),xmol(2),dop(1),dop(2),dopType,xini,thickcolumn(indTrans),...
                CellInfo(indLong-1,indTrans).yfin,d_dd(indLong-1)*1e-7,meshcolumn(indTrans),mesh_dd(indLong-1),hetero,heterojunctionth);
%         end
        xini=xini+thickcolumn(indTrans);
    end % transversal loop
    
end % longitudinal loop

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