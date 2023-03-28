
function [mesh] = f_EvalDopingMobility(mesh,mode)
%
% This function computes the zero-field mobilities introducing the effect
% of the scattering from doping impurities. As of this moment, scattering
% is assumed for all the introduced impurities, so not just from the
% ionized ones.
%
% Two possibilities are currently introduced:
%
% - 'none' : no model is introduced, and mobility is assumed to be equal to
%            the intrinsic one
%
% - 'Hilsum' : model from C. Hilsum, "Simple empirical relationship between
%              mobility and carrier concentration", Electronics Letters,
%              vol. 10, no. 13, June 1974.
%
% The input quantities are assumed to be defined on nodes, whereas the
% output ones will be defined on elements.
%
% Alberto Tibaldi, 21/10/2016
%

if isfield(mode,'NxCoe')
 Tt=pdeintrp(mesh.node,mesh.triangle(1:4,:),mesh.T.'); % T on triangles
 Tn=mesh.T;
 T300=mode.T300;
 NXt=mode.CoeffHilsum0+(Tt-T300)*mode.NxCoe*1e18; % Coefficient N_X for the Hilsum model
 NXn=mode.CoeffHilsum0+(Tn-T300)*mode.NxCoe*1e18; % Coefficient N_X for the Hilsum model
else
 NXt=mode.CoeffHilsum;
 NXn=mode.CoeffHilsum;
end

NXt=NXt/mode.CarrierNorm;
NXn=NXn/mode.CarrierNorm;
%NX=mode.CoeffHilsum;
%'Hilsum', keyboard

switch(mode.DopingMobilityModel)
    case{'none'}
        mobn0_n=mesh.mobnint_n;
        mobp0_n=mesh.mobpint_n;
        mobn0_t=mesh.mobnint_t;
        mobp0_t=mesh.mobpint_t;
        %
    case{'Hilsum'}
        NA_t=mesh.dop_a_t;
        ND_t=mesh.dop_d_t;
        NA_n=mesh.dop_a;
        ND_n=mesh.dop_d;
        Denominator_t=1+sqrt((NA_t+ND_t)./NXt);
        Denominator_n=1+sqrt((NA_n+ND_n)./NXn);
        mobn0_t=(mesh.mobnint_t)./Denominator_t;
        mobp0_t=(mesh.mobpint_t)./Denominator_t;
        mobn0_n=(mesh.mobnint_n)./Denominator_n;
        mobp0_n=(mesh.mobpint_n)./Denominator_n;
        %
    case{'Roland'}
        NA_t=mesh.dop_a_t;
        ND_t=mesh.dop_d_t;
        NA_n=mesh.dop_a;
        ND_n=mesh.dop_d;
        Denominator_tE=1+((NA_t+ND_t)./NXt).^0.35;
        Denominator_tH=1+((NA_t+ND_t)./(2*NXt)).^0.35;   
        Denominator_nE=1+((NA_n+ND_n)./NXn).^0.35;
        Denominator_nH=1+((NA_n+ND_n)./(2*NXn)).^0.35;           
        mobn0_t=(mesh.mobnint_t)./Denominator_tE;
        mobp0_t=(mesh.mobpint_t)./Denominator_tH;
        mobn0_n=(mesh.mobnint_n)./Denominator_nE;
        mobp0_n=(mesh.mobpint_n)./Denominator_nH;        
    case{'RolandFC'} % not using dopants, but free-carriers!
    %'FC', keyboard

      if isfield(mode,'elec')
%        'sono qui', keyboard      
        elec_t=pdeintrp(mesh.node,mesh.triangle(1:4,:),mode.elec');
        hole_t=pdeintrp(mesh.node,mesh.triangle(1:4,:),mode.hole');
        NA_t=hole_t;
        ND_t=elec_t;
        NA_n=mode.hole;
        ND_n=mode.elec;
%        'passo FULL', keyboard
       else 
        NA_t=mesh.dop_a_t;
        ND_t=mesh.dop_d_t;
        NA_n=mesh.dop_a;
        ND_n=mesh.dop_d;
        'passo Zero',        
       end         
        Denominator_tE=1+((NA_t+ND_t)./NXt).^0.35;
        Denominator_tH=1+((NA_t+ND_t)./(2*NXt)).^0.35;   
        Denominator_nE=1+((NA_n+ND_n)./NXn).^0.35;
        Denominator_nH=1+((NA_n+ND_n)./(2*NXn)).^0.35;         
        mobn0_t=(mesh.mobnint_t)./Denominator_tE;
        mobp0_t=(mesh.mobpint_t)./Denominator_tH;
        mobn0_n=(mesh.mobnint_n)./Denominator_nE;
        mobp0_n=(mesh.mobpint_n)./Denominator_nH;        
    case{'RolandIon'}
      if isfield(mode,'dop_am')
%        'sono qui', keyboard      
        dop_am_t=pdeintrp(mesh.node,mesh.triangle(1:4,:),mode.dop_am');
        dop_dp_t=pdeintrp(mesh.node,mesh.triangle(1:4,:),mode.dop_dp');
        NA_n=mode.dop_am;
        ND_n=mode.dop_dp;
        NA_t=dop_am_t;
        ND_t=dop_dp_t;

%        'passo FULL', keyboard
       else 
        NA_t=mesh.dop_a_t;
        ND_t=mesh.dop_d_t;
        NA_n=mesh.dop_a;
        ND_n=mesh.dop_d;

        'passo Zero',        
       end 
        Denominator_tE=1+((NA_t+ND_t)./NXt).^0.35;
        Denominator_tH=1+((NA_t+ND_t)./(2*NXt)).^0.35;   
        Denominator_nE=1+((NA_n+ND_n)./NXn).^0.35;
        Denominator_nH=1+((NA_n+ND_n)./(2*NXn)).^0.35;
        mobn0_t=(mesh.mobnint_t)./Denominator_tE;
        mobp0_t=(mesh.mobpint_t)./Denominator_tH;
        mobn0_n=(mesh.mobnint_n)./Denominator_nE;
        mobp0_n=(mesh.mobpint_n)./Denominator_nH;        
    otherwise
        error('P-please choose a doping-dependent mobility model >_>')
end

% plots doping e Free Carrier
% figure,plot(mesh.ygrid,mode.dop_dp(1:mesh.nny),mesh.ygrid,mesh.dop_d(1:mesh.nny),mesh.ygrid,mode.elec(1:mesh.nny)+10,'k'),set(gca,'yscale','log')
% figure,plot(mesh.ygrid,mode.dop_dp(1:mesh.nny),mesh.ygrid,mesh.dop_d(1:mesh.nny),mesh.ygrid,mode.elec(1:mesh.nny)+10,'k'),legend('ND+','ND','elec'),set(gca,'yscale','log')
% figure,plot(mesh.ygrid,mode.dop_dp(1:mesh.nny),mesh.ygrid,mesh.dop_d(1:mesh.nny),mesh.ygrid,mode.elec(1:mesh.nny)+10,'k'),set(gca,'yscale','log'),ylim([1e17,1e20])

mesh.mobn0_t=mobn0_t;
mesh.mobp0_t=mobp0_t;
mesh.mobn0_n=mobn0_n;
mesh.mobp0_n=mobp0_n;
     
return