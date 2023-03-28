nn=mesh.nn;                             % -> elec eqs.
pp=nn+mode.nflg*nn;                     % -> hole eqs.
tt=pp+mode.pflg*nn;                     % -> trap eqs.
qq=tt+mode.tflg*nl*nn;                  % -> circuit eqs.
rr=qq+nm;                               % -> current eqs.
vv=rr+nm;                               % -> n2D eqs.
ww=vv+length(mesh.xgrid)*NQW;           % -> p2D eqs.
ss=ww+length(mesh.xgrid)*NQW;           % -> optical rate eqs.

neq=nn+mode.nflg*nn+mode.pflg*nn+mode.tflg*nl*nn+2*nm+mode.oflg*not(mode.firstrun)*2*NQW*length(mesh.xgrid)+mode.oflg*mode.nmodes;
tvet = sparse(neq,1);
tvet(qq+1) = 1;

iomega = 1;
omega = 1; % static analysis
    
%% Matrices for Broyden's scheme
invJmatDC = inv(JmatDC);
yvet = rvet_new - rvet_old;
den = deltau.'*invJmatDC*yvet;
invJmatDC = invJmatDC + (1./den).*(deltau - invJmatDC*yvet)*deltau.'*invJmatDC;

Y_ssB = zeros(1,length(mode.fvet));
G_ssB = zeros(1,length(mode.fvet));
C_ssB = zeros(1,length(mode.fvet));
Pst_ssB = zeros(mode.nmodes,length(mode.fvet));

deltaBRO = invJmatDC*tvet;
    

Y_B=deltaBRO(rr+1); % small-signal admittance (Broyden's scheme), ohm/cm^2

Y_ssB(iomega)=Y_B;
G_ssB(iomega)=real(Y_B); % small-signal conductance (Broyden's scheme), ohm/cm^2
C_ssB(iomega)=imag(Y_B/omega);  % small-signal capacitance (Broyden's scheme), ohm/cm^2

Pst_ssB(:,iomega) = deltaBRO(ss+[1:mode.nmodes]);


%% Matrix for Newton's scheme
JmatDC = Kmat0 + Jmat0 + Jmat2;

Y_ss = zeros(1,length(mode.fvet));
G_ss = zeros(1,length(mode.fvet));
C_ss = zeros(1,length(mode.fvet));
Pst_ss = zeros(mode.nmodes,length(mode.fvet));


% for iomega = 1:length(mode.fvet)
% parfor iomega = 1:length(mode.fvet)
%     omega = mode.fvet(iomega);

%     omega = 2*pi*mode.fvet(iomega); % full fvet small-signal analysis
    
    Jmat = JmatDC + j.*omega.*Jmat1;
    
    delta = Jmat\tvet;
    
    
    Y=delta(rr+1); % small-signal admittance, ohm/cm^2
    
    Y_ss(iomega)=Y;
    G_ss(iomega)=real(Y); % small-signal conductance, ohm/cm^2
    C_ss(iomega)=imag(Y/omega);  % small-signal capacitance, ohm/cm^2

    Pst_ss(:,iomega) = delta(ss+[1:mode.nmodes]);
    
    disp(['Small-signal analysis ',num2str(iomega),' of ',num2str(length(mode.fvet)),' executed'])
    
% end

AM = abs(sum(Pst_ss)./Y_ss);
AM = AM./AM(1);

