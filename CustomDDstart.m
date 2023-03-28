% Torrelli Valerio: modification 13/10/2022
% Model to invent the guess for the electrostatic potential.

% This was tested mainly for the multi-TJ case, sepcifically for the 2 and
% 3 TJ case.

if mode.DDLater == 1 && indv == 2
    
    indv = find(mode.v0_dd == mode.WantedV); % Setting the voltage index corresponding to the wanting voltage.
    PhiEquil = mode.phi(1:mesh.nny); % Saving the equilibrium electrostatic potential.
    
    % Determining the width and the position of the ARs
    
    dop = mesh.dop(1:mesh.nny); % doping profile at the central column of the VCSEL. I can extract the AR info looking at where the doping is zero!
    Ec = mode.ecb(1:mesh.nny); % I need to know where this variable is defined NaN in order to exclude the contacts.
    IndNan = find(isnan(Ec));
    dop(IndNan) = NaN; % Now I have excluded the contacts
    IndZeroDop = find(dop == 0); % Indices identifying the ARs, now I have to separate them.
    IndJumps = find(diff(IndZeroDop)>1); % Plot IndZeroDop and diff(IndZeroDop) to understand which "jumps" I am referring to.
    NAR = length(IndJumps) + 1; % #AR    
    zAR = zeros(1,NAR);    
    
    for iAR = 1:NAR % loop over the active regions
        if iAR == 1
            zAR(iAR) = mesh.ygrid(round(0.5*(IndZeroDop(1)+IndZeroDop(IndJumps(1)))));
            WidthAR = mesh.ygrid(IndZeroDop(IndJumps(1))) - mesh.ygrid(IndZeroDop(1));
        elseif iAR == NAR
            zAR(iAR) = mesh.ygrid(round(0.5*(IndZeroDop(IndJumps(end)+1)+IndZeroDop(end))));
        else
            zAR(iAR) = mesh.ygrid(round(0.5*(IndZeroDop(IndJumps(iAR-1)+1)+IndZeroDop(IndJumps(iAR)))));
        end
    end
    
    mode.zAR = zAR;
    mode.WidthAR = WidthAR*mode.AdjustWAR;
    
    % Starting the evaluation of the correction
    
    if mode.flgBTJ == 1
        
        % Determining the width and the position of the TJs
        
        zTJ = zeros(size(mode.zAR));
        for inTJ = 1:length(mode.zAR)
            zTJ(inTJ) = mesh.ygrid(mesh.iLBTJ{inTJ,2}(1));
            WidthTJ = abs(mesh.ygrid(mesh.iLBTJ{inTJ,3}(1)) - mesh.ygrid(mesh.iLBTJ{inTJ,1}(1)));
        end
        mode.zTJ = zTJ;
        mode.WidthTJ = WidthTJ*mode.AdjustWTJ;
        
        DeltaPhiFit = zeros(size(mesh.ygrid));
        
        for inTJ = 1:length(mode.zAR)
            Corr = zeros(size(mesh.ygrid));
            Corr = (mode.WantedV/length(mode.zAR))*mode.FracAR*(1+tanh((mesh.ygrid-mode.zAR(inTJ))/mode.WidthAR))/2;
            DeltaPhiFit = DeltaPhiFit + Corr;
            Corr = zeros(size(mesh.ygrid));
            Corr = (mode.WantedV/length(mode.zTJ))*mode.FracTJ*(1+tanh((mesh.ygrid-mode.zTJ(inTJ))/mode.WidthTJ))/2;
            DeltaPhiFit = DeltaPhiFit + Corr;
        end
        
    else
        mode.FracAR = 1;        
        DeltaPhiFit = zeros(size(mesh.ygrid));
        
        for inTJ = 1:length(mode.zAR)
            Corr = zeros(size(mesh.ygrid));
            Corr = (mode.WantedV/length(mode.zAR))*mode.FracAR*(1+tanh((mesh.ygrid-mode.zAR(inTJ))/mode.WidthAR))/2;
            DeltaPhiFit = DeltaPhiFit + Corr;
        end
        
    end
    
    PhiFit = PhiEquil + DeltaPhiFit;
    
    
    % Updating the initial guess
    
    for indCol = 1:mesh.nnx
        mode.phi((1:mesh.nny)+mesh.nny*(indCol-1)) = PhiFit;
    end
    uvet(1:nn)=mode.phi;
    
    
end

