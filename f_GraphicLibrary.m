function [semiconductor,metal,color]=f_GraphicLibrary(material)
%
switch material
    case {'XY-LiNbO3','YX-LiNbO3'} % X-cut Y-propagating & Y-cut X-propagating LiNbO3
        semiconductor=0; metal=0; color=[0 1 0];
        
    case {'ZX-LiNbO3','ZY-LiNbO3'} % Z-cut X-propagating & Z-cut Y-propagating LiNbO3
        semiconductor=0; metal=0; color=[0 1 0];
        
    case {'XZ-LiNbO3','YZ-LiNbO3'} 
        semiconductor=0; metal=0; color=[0 1 0];
        
    case 'SiO2' % silicon-oxide
        semiconductor=0; metal=0; color=[0 1/2 1];
        
    case 'MgO' % magnesium-oxide
        semiconductor=0; metal=0; color=[1 0 0];
        
    case 'ZnO' % zinc-oxide
        semiconductor=1; metal=0; color=[1 0 0];
        
    case 'Si' % silicon
        semiconductor=1; metal=0; color=[0 1 0];
        
    case 'Ge' % germanium
        semiconductor=1; metal=0; color=[0.2 0.6 0.0];
        
    case 'SiGe' %
        semiconductor=1; metal=0; color=[0 .5 0.5];
        
    case 'i-Si' % semi-insulating silicon
        semiconductor=0; metal=0; color=[0.7 1 0];
        
    case 'p-Si' % semi-insulating silicon
        semiconductor=0; metal=0; color=[0.7 1 0];
        
    case 'epi-Si' % semi-insulating silicon
        semiconductor=0; metal=0; color=[0.7 0.8 0];
        
    case 'vacuum' % free-space
        semiconductor=0; metal=0; color=[0 1 1];
         color=[1 1 1];
    case 'Au' % gold
        semiconductor=0; metal=1; color=[1 1 0];
        color=[0.8 0.7 0.2];
     %   color=[0 0 0];
        
    case 'Al' % alluminium
        semiconductor=0; metal=1; color=[1/4 1/2 1/3];
        
    case 'Cu' % copper
        semiconductor=0; metal=1; color=[1 1 1/2];
        
    case 'bcb' % BCB
        semiconductor=0; metal=0; color=[0 1/2 0.8];
        
    case 'SiN' % silicon nitride
        semiconductor=0; metal=0; color=[0 0 0];
        
    case 'Si3N4' % silicon nitride
        semiconductor=0; metal=0; color=[0.2 0.7 0];
        
    case 'InP' % indium phosphide
        semiconductor=1; metal=0; color=[1 0.5 0.5];
        
    case 'i-InP' % indium phosphide
        semiconductor=0; metal=0; color=[0.8 0.5 0.5];
        
    case 'InGaAsP' %
        semiconductor=1; metal=0; color=[1 0.7 0.7];
        
    case 'GaAsP' %
        semiconductor=1; metal=0; color=[0.7 0.7 0.7];
        
    case 'AlGaAs' % aluminum gallium arsenide
%         semiconductor=1; metal=0; color=[0.8 0.7 0.2];
        semiconductor=1; metal=0; color=[1 1 1]*.9;

    case 'InGaAs' % indium gallium arsenide
        semiconductor=1; metal=0; color=[0.7 0.6 0.1];

    case 'AlOx' % oxide of AlAs (AlGaAs,xmol=0)
        semiconductor=0; metal=0; color=[1 0 0];

    case 'Polyamide' % oxide of AlAs (AlGaAs,xmol=0)
        semiconductor=0; metal=0; color=[1 0 1];

    case 'BaTiO3' % barium titanate
        semiconductor=0; metal=0; color=[0 .5 0];
        
    case 'BaTiO3-PLD31' % barium titanate
        semiconductor=0; metal=0; color=[0 0 0];
        
    case 'BaTiO3-corea' % barium titanate
        semiconductor=0; metal=0; color=[0 0 0];
        
    case 'BaTiO3-PLD33' % barium titanate
        semiconductor=0; metal=0; color=[0 0 0];
        
    case 'BaTiO3-PLD40' % barium titanate
        semiconductor=0; metal=0; color=[0 0 0];
        
    case 'sapphire' %
        semiconductor=0; metal=0; color=[0 0 0];
        
    case 'GaN' % gallium nitride
        semiconductor=1; metal=0; color=[0 1 1/2];        
        
    case 'AlGaN' % gallium nitride
        semiconductor=1; metal=0; color=[0 1 1/3];        
        
    case 'HgCdTe' % mercury cadmium telluride
        semiconductor=1; metal=0; color=[0 1 1/2];
    otherwise,
        error('Unknown material'),
%        ('Unknown material'),
 %      semiconductor=0; metal=0; color=[0 1 1];
 %        color=[1 1 1];
end
% *********************************************************************************************100
