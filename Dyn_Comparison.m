% close all

Z0=50;

%% VENUS results
if isfield(mode,'AM')
    
    % Parasitic Parameters
    Rm = 0;
    Cp = 0;
    Radd = 0;
    
    omega = mode.fvet*2*pi;
    
%     % Differential conductance and resistance coming from VENUS
%     G_ss=mode.Gss;
%     C_ss=mode.Css;
%     Yv=zeros(size(C_ss));
%     
%     for iF = 1:length(mode.fvet)
%         Yv(:,iF) = 1i*omega(iF)*C_ss(:,iF)+G_ss(:,iF);
%         Ys(:,iF) = 1./(1./(Yv(:,iF))+Rm);
%         Ypad(:,iF) = Ys(:,iF) + 1i*omega(iF)*Cp;
%     end
% 
%     Zpad = 1./Ypad;
%     Zin = Zpad + Radd;
%     Z_ss = Zin;
    
%     Z_ssB = 1./real(modePlot.YssB(:,1));
    Z_ss = 1./real(modePlot.Yss(:,1));
    Z_ss2 = 1./real(modePlot.Yss(:,2));
    Z_ss3 = 1./real(modePlot.Yss(:,3));
    Z_ss4 = 1./real(modePlot.Yss(:,4));
    
    % Reflection coefficient
    S11=(Z_ss-Z0)./(Z_ss+Z0);
    
    AM=mode.AM;
    fGH=mode.fvet/1e9;
    
    %% VENUS1 Quasi-static resistance
    
    figure(1),subplot(234),hold on
    plot(mode.CurDyn,Z_ss,[colo(kpar),'*'])
    plot(mode.CurDyn,Z_ss2,'*')
    plot(mode.CurDyn,Z_ss3,'*')
    plot(mode.CurDyn,Z_ss4,'*')
%     plot(mode.CurDyn,real(Z_ssB),'bo')
    
    figure(6)
    plot(mode.CurDyn,Z_ss,[colo(kpar),'*'])
    plot(mode.CurDyn,Z_ss2,'*')
    plot(mode.CurDyn,Z_ss3,'*')
    plot(mode.CurDyn,Z_ss4,'*')
   
    
    %% Extract experimental results from SmallSignalJul19 folder
    if input('Do you want to plot dynamic experimental results? 1:YES, Enter:NO')==1
    
    wafer = '01103_448C';
    device = 'C134R233';    % available: C134R233, C134R235, C134R237
    
    % Tvet = [20:10:60];      % available: 20, 30, 40, 50, 60
    % Ivet = [1];             % available: 1, 2, 3, 4, 5, 6
    
    Tvet = [20];            % available: 20, 30, 40, 50, 60
    Tvet = input('Select Temperature to be plotted (20:60) ');
%     Ivet = [1:6];           % available: 1, 2, 3, 4, 5, 6
    Ivet = [1 2 3 4 5 6];           % available: 1, 2, 3, 4, 5, 6
    
    colorMatrix = [[0 0.4470 0.7410],
          	[0.8500 0.3250 0.0980],
          	[0.9290 0.6940 0.1250],
          	[0.4940 0.1840 0.5560],
          	[0.4660 0.6740 0.1880],
          	[0.3010 0.7450 0.9330],
          	[0.6350 0.0780 0.1840]];
        
    for indT = 1:length(Tvet)
        for indI = 1:length(Ivet)
            fileName = ['MarkusNov18\SmallSignalJul19\',...
                wafer,'_',num2str(Tvet(indT)),'\',wafer,'_',device,'_I_',num2str(Ivet(indI)),'mA_',num2str(Tvet(indT)),'°C.s2p'];
            Sstr = sparameters(fileName);
            vFreq = Sstr.Frequencies;
            
            % Scattering parameters extraction
            S21_exp = squeeze(Sstr.Parameters(2,1,:));
            S11_exp = squeeze(Sstr.Parameters(1,1,:));
            
            Z_exp = Z0*(1+S11_exp)./(1-S11_exp);
            Z_ssEXP(indT,indI) = real(Z_exp(1));
            
        end
    end
    
    figure(1),subplot(234),hold on
            plot(Ivet,Z_ssEXP,[colo(kpar),'diamond'])
            
            figure(6),plot(Ivet,Z_ssEXP,[colo(kpar),'diamond'])
    end
    
    figure(1),subplot(234)
    legend('exp-static','VENUS1-static','VENUS1-dyn','VENUS1-dyn (Broy)','exp-dyn')
    
    
end