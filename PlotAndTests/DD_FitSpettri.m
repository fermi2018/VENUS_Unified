
if exist('modePlot')==0
    modePlot=MODEplot{1};
    % % xQW=P.mesh.xgrid(1:P.mesh.nnxQW{1})*1e4;
    xQW=mesh.xgrid(1:mesh.nnxQW{1})*1e4;
end

if exist('CORRENTI')==0
    CORRENTI=[1 4 8];
end
CorLav=CORRENTI;
corrente=modePlot.ii_dd*1000;

pu=[];
for kC=1:length(CorLav)
    [~,puk]=min(abs(corrente-CorLav(kC)));
    pu=[pu puk];
end
%pu=1:10:indv;
%pu=[40:20:120];
%pu=[30:15:indv];

x=modePlot.x;

mode=modePlot;
Oxide=modePlot.Oxide;

Glut4Dinc(mode)

clear Gv Nv

%pu=1:indv;
kc=0;
for kv=pu
    kc=kc+1;
    Ver_Gain
    Gv(kc,:)=G0;
    Nv(kc,:)=DeltaN0;
end

hp=figure,

set(hp,'pos',[160         428        1390         525])
subplot(121)

E=squeeze(modePlot.E2(pu,1,:));
plot(xQW,modePlot.matgain(pu,:),'linewidth',2), hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(xQW,Gv,'.'),
ax.ColorOrderIndex = 1;
plot(x,E/max(E(:,1))*1000,'--','linewidth',2), hold on
% plot(xQW,E(:,1:mesh.nnxQW{1})/max(E(:,1)).*modePlot.matgain(pu,:),'--','linewidth',2), hold on
ylabel(' QW gain, 1/cm ')
xlabel(' \rho, \mum ')
title(['Correnti= ',num2str(CorLav)])
grid
subplot(122)
plot(xQW,modePlot.dn(pu,:),'linewidth',2), hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(xQW,Nv,'.'),
ylabel('\Deltan QW ')
xlabel(' \rho, \mum ')
grid
pausak

DLa=modePlot.Deltalam;
DeltaLam=[-20:5:20]+DLa;

clear GvS DelPlo

fiox=find(x<=Oxide*2);
%'qui gain', keyboard

%N2lin=modePlot.Elqw(kv,fiox);
%P2lin=modePlot.Hoqw(kv,fiox);
%T2=modePlot.Tqw(kv,fiox);
%L2=modePlot.lambda(1,kv);

kinc=0;
for kv=pu
    E=squeeze(modePlot.E2(kv,:,:));
    if size(E,2)==1
        E=E';
    end
    kinc=kinc+1;
    for kla=1:length(DeltaLam)
        DeltalamLocal=DeltaLam(kla);
        Spet_Gain
        N2v(kinc,:)=N2lin;
        P2v(kinc,:)=P2lin;
        T2v(kinc,:)=T2;
        L2v(kinc)=L2;
        dx=diff(x);
        Gamma=modePlot.Gamma;
        Gxdx=sum(Gamma)*1e-8*(G0.*x(fiox).*dx(fiox)*2*pi).';
        %Gxdx=1e-8*(x(fiox).*dx(fiox)*2*pi).';
        
        Gm=E(:,fiox)*Gxdx;
        
        
        %GvS(kinc,kla)=mean(G0);
        GvS(kinc,kla)=max(Gm);
    end
    DelPlo(kinc,:)=L2+DeltaLam;
end
% xme=x(fiox);
% hp=figure,
% set(hp,'pos',[160         428        1190         525])
% subplot(121)
% plot(xme,N2v,'linewidth',2), hold on,
% ax = gca;       ax.ColorOrderIndex = 1;
% plot(xme,P2v,'.'),
% ylabel(' Electron (--) /Hole (\cdot) density, 1/cm^3 ')
% xlabel(' \rho, \mum')
% title(['Correnti= ',num2str(CorLav)])

% An alternative is plotting directly the quantities stored in modePlot
xme=x(fiox);
hp=figure,
set(hp,'pos',[160         428        1190         525])
subplot(121)
hold on
for iI=1:length(pu)
    plot(xme,mode.nQW{pu(iI)}{2}(fiox),'linewidth',2) % indicates the central QW!
end
ax = gca;       ax.ColorOrderIndex = 1;
for iI=1:length(pu)
    plot(xme,mode.pQW{pu(iI)}{2}(fiox),'o')
end
ylabel(' Electron (--) /Hole (\cdot) density, 1/cm^3 ')
xlabel(' \rho, \mum')
title(['Correnti= ',num2str(CorLav)])


subplot(122)
plot(xme,T2v,'linewidth',2),
ylabel(' Temperature, K')
xlabel(' \rho, \mum')
pausak

Sym='o+ps';
coloriWhite
hg=figure,
set(hg,'pos',[410         428        1490         525])
subplot(121)
plot(DelPlo',GvS','linewidth',2)
lamPlo=modePlot.lambda(pu);
hold on
ax = gca;
ax.ColorOrderIndex = 1;
[du,idlaa]=min(abs(DLa-DeltaLam));
for k=1:length(pu)
    plot(DelPlo(k,idlaa),GvS(k,idlaa),'o'),
end
grid
NOMELUT=modePlot.LUT{modePlot.iLUT};
for kl=1:length(NOMELUT)
    if strcmp(NOMELUT(kl),'_')
        NOMELUT(kl)='-';
    end
end
fis= strfind(NOMELUT,'\');
NLU=NOMELUT(fis(end)+1:end);
%DirName=structure(1:fis(end));
title(['Correnti= ',num2str(CorLav),'      LUT= ',NLU])
legend(num2str(CorLav'),'location','Best')
xlabel(' Wavelength, nm')
ylabel(' Gain QW')
subplot(122)

%plot(DelPlo,GvS,'o-','linewidth',2)
plot(DelPlo,GvS,'linewidth',2)
hold on
for kCu=1:length(CorLav)
    plot(DelPlo(kCu,:),GvS(kCu,:),[Sym(kCu),'k'])
end
lamPlo=modePlot.lambda(pu);


[du,idlaa]=min(abs(DLa-DeltaLam));
plot(DelPlo(:,idlaa),GvS(:,idlaa),'k','linewidth',2),
DP=[min(min(DelPlo)) max(max(DelPlo)) ];
Los=min(modePlot.Lmod(:,end));
%'Los', keyboard
Gth=[Los Los];
hold on, plot(DP,Gth,'k', 'linewidth',2)
%legend(num2str(DeltaLam'),'location','best')
grid
xlabel(' Wavelength, nm')
ylabel(' Gain QW')
title(['\Delta\lambda= ',num2str(DeltaLam),' nm'])
pausak

%lamPlo=modePlot.lambda(pu);
%lamPlot=[1; 1]*lamPlo;
%yplo=[-1e4; 1e4]*ones(size(lamPlo));
%plot(lamPlot,yplo)

saltoVL=1;

if saltoVL==0
    ' Vario le LUT'
    
    
    LUT{1}='LUT4D_Feb_Lorentzian';
    LUT{2}='LUT4D_Feb_nMark_20';
    LUT{3}='LUT4D_Feb_nMark_40';
    
    
    h=figure;
    set(h,'pos',[207          21        1700         963])
    for kL=1:3
        subplot(3,2,kL)
        NOMELUT=LUT{kL};
        mode.GLUT=[NOMELUT,'_Der.mat'];
        Glut4Dinc(mode)
        
        kinc=0;
        for kv=pu
            kinc=kinc+1;
            for kla=1:length(DeltaLam)
                DeltalamLocal=DeltaLam(kla);
                Spet_Gain
                N2v(kinc,:)=N2lin;
                P2v(kinc,:)=P2lin;
                T2v(kinc,:)=T2;
                L2v(kinc)=L2;
                %'stop', keyboard
                
                Gxdx=Gamma*1e-8*(G0.*x(fiox).*dx(fiox)*2*pi).';
                %Gxdx=1e-8*(x(fiox).*dx(fiox)*2*pi).';
                
                Gm=E(:,fiox)*Gxdx;
                GvS(kinc,kla)=max(Gm);
            end
            DelPlo(kinc,:)=L2+DeltaLam;
        end
        
    end %saltoVL
    
    figure(h)
    subplot(3,2,2*kL-1)
    plot(DelPlo',GvS','linewidth',2)
    lamPlo=modePlot.lambda(pu);
    hold on
    ax = gca;
    ax.ColorOrderIndex = 1;
    [du,idlaa]=min(abs(DLa-DeltaLam));
    for k=1:length(pu)
        plot(DelPlo(k,idlaa),GvS(k,idlaa),'o'),
    end
    grid
    for kl=1:length(NOMELUT)
        if strcmp(NOMELUT(kl),'_')
            NOMELUT(kl)='-';
        end
    end
    title(['Correnti= ',num2str(CorLav),'      LUT= ',NOMELUT])
    figure(h)
    subplot(3,2,2*kL)
    
    plot(DelPlo,GvS,'linewidth',2)
    hold on
    for kCu=1:length(CorLav)
        plot(DelPlo(kCu,:),GvS(kCu,:),[Sym(kCu),'k'])
    end
    hold on, plot(DP,Gth,'k', 'linewidth',2)
    lamPlo=modePlot.lambda(pu);
    hold on
    ax = gca;
    ax.ColorOrderIndex = 1;
    [du,idlaa]=min(abs(DLa-DeltaLam));
    %for kCu=1:length(CorLav)
    plot(DelPlo(:,idlaa),GvS(:,idlaa),'k-','linewidth',2),
    %end
    title(['DeltaLam= ',num2str(DeltaLam)])
    %legend(num2str(DeltaLam'),'location','best')
    %legend(num2str(DeltaLam'),'location','best','orientation','horizontal')
    grid
end
