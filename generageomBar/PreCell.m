%

if isfield(mode,'xmolPiatto')
    if mode.xmolPiatto>=1
        [du,fiM]=max(d_dd(2:end));
        fiM=fiM+1;
        xM=x_dd(fiM,1);
        fiA=find(x_dd(:,1)>0);
        x_dd(fiA,1)=xM;
        x_dd(fiA,2)=xM;
        if flgStop==1
            'qui', keyboard
        end
    end
end

if isfield(mode,'xmolPiatto')
    if mode.xmolPiatto==2
%         [du,fiM]=max(d_dd(2:end));
%         fiM=fiM+1;
%         xM=x_dd(fiM,1);
        fiA=find(Na_dd(:,1)>0);
        Na=mean(Na_dd(fiA,1));
        Na_dd(fiA,1)=Na;
        Na_dd(fiA,2)=Na;
        fiA=find(Nd_dd(:,1)>0);
        Nd=mean(Nd_dd(fiA,1));
        Nd_dd(fiA,1)=Nd;
        Nd_dd(fiA,2)=Nd;
        if flgStop==1
            'qui dop', keyboard
        end
    end
end

if igrigliaSTOP==1
    mesh_r
    punti_X=sum(mesh_r)
    % 'Stop grigliato e material', keyboard
end
if(isempty(StrDD.r_passiv))
    thickcolumn=[ra_dd(1),diff(ra_dd)]*1e-4;
    meshcolumn=[mesh_r];
    if size(material,2)>length(mesh_r)
        meshcolumn=[mesh_r 1];
        thickcolumn=[ra_dd(1) .1]*1e-4;
    end
    flagPassive=0;
else
    thickcolumn=[ra_dd(1),diff(ra_dd),r_passiv]*1e-4;
    meshcolumn=[mesh_r,mesh_passiv];
    material(:,end+1) = material(:,end);
%'dopo', keyboard
    ms=material;
%    material(end-StrDD.fiPassiv1:end,end) = {'Polyamide'};
    material(StrDD.fiPassiv2:end,end) = {'Polyamide'};
    material(end,end) = {'Line'};
%    material(end,end) = {'vacuum'};
    flagPassive=1;
end
numcolumn=length(thickcolumn);
%'thick ', keyboard

if contStr==1 && flgStop==1
    'dopo StrRead', keyboard
end

