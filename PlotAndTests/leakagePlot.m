JN_X=mode.JXn;
JN_Y=(mode.JYn);
JP_X=mode.JXp;
JP_Y=(mode.JYp);

% modePlot=MODEplot{1};
%
% JN_X=modePlot.JXn;
% JN_Y=(modePlot.JYn);
% JP_X=modePlot.JXp;
% JP_Y=(modePlot.JYp);
%
% CORRENTI=[2 4];
% CorLav=CORRENTI;
% corrente=modePlot.ii_dd*1000;
% pu=[];
% for kC=1:length(CorLav)
% [du,puk]=min(abs(corrente-CorLav(kC)));
% pu=[pu puk];
% end
% fipst=pu;
%
% Vcor=pu;
% fipst=pu;
% %
% x=mesh.xgrid*1e4;
% z=mesh.ygrid*1e4;
%

z=y;
ycm=z*1e-4;
xcm=x*1e-4;

ipcor=0;
Fall=0;
Fall=Fall+2;
figure(Fall)
%set(hc,'pos',[80         421        1800         525 ])

set(Fall,'pos',[141         136        1275         829 ])
subplot(221)
xqw=1000*(modePlot.y-modePlot.yQW);
plot(xqw,modePlot.Ec(fipst,:),'linewidth',2)
%plot(1000*(modePlot.y-modePlot.yQW),modePlot.Ec(fipst,:)-modePlot.Ev(fipst,:),'linewidth',2)
ax = gca;
ax.ColorOrderIndex = 1;
hold on, plot(xqw,modePlot.Ev(fipst,:)+1.8,'--','linewidth',2)

% figure
% xqw=(modePlot.y)-350;
% plot(xqw,modePlot.Ec(fipst,:)+1.8,'linewidth',1)
%  ax = gca;
%  ax.ColorOrderIndex = 1;
% hold on, plot(xqw,modePlot.Ev(fipst,:)+1.8,'r','linewidth',1)
% xlim([-1 9])
% plot(xqw,modePlot.EFc(fipst,:)+1.8,'b','linewidth',2)
% plot(xqw,modePlot.EFv(fipst,:)+1.8,'r','linewidth',2)

%figure, plot(y,Rc-Rv)
dX=80;
xlim([-dX,dX])
title(' Ec (cont), Ev+1.8 (dashed) in Cavity')
xlabel(' z around QW (nm) ')
grid

%keyboard


for pcor=Vcor
    ipcor=ipcor+1;
    for indy=1:length(z)
        
        jn_y=squeeze(JN_Y(pcor,indy,:))';
        jp_y=squeeze(JP_Y(pcor,indy,:))';
        %    figure, plot(x,jn_y,x,jp_y), pausak
        curr_n(indy)=1000*trapz(xcm,2.*pi.*xcm.*jn_y);
        curr_p(indy)=1000*trapz(xcm,2.*pi.*xcm.*jp_y);
        
    end
    
    zcm=ycm;
    clear curx_n curx_p
    jnQ_x=sum(squeeze(mode.JnQW(pcor,:,:)),2)*WQW*1000;
    jpQ_x=sum(squeeze(mode.JpQW(pcor,:,:)),2)*WQW*1000;
    jQn=[jnQ_x' zeros(1,length(xcm)-length(xQW))];
    jQp=[jpQ_x' zeros(1,length(xcm)-length(xQW))];
    for indx=1:size(JN_X,3)
        
        jn_x=squeeze(JN_X(pcor,:,indx))';
        jp_x=squeeze(JP_X(pcor,:,indx))';
        %    figure, plot(zcm,jn_y,zcm,jp_y), pausak
        curx_n(indx)=1000*trapz(zcm,jn_x);
        curx_p(indx)=1000*trapz(zcm,jp_x);
        
    end
    
    inoplot=0
    if inoplot==0
        xcD=xcm*1e4;
        figure(99)
        set(99,'pos',[ 358         419        1500         525])
        
        % QW
        subplot(131)
        plot(xcD,2*pi*xcm.*curx_n,xcD,2*pi*xcm.*curx_p,'o-'),
        hold on, grid on
        ax = gca;
        ax.ColorOrderIndex = 1;
        plot(xcD,jQn,'--',xcD,jQp,'--','linewidth',2)
        legend('3D elec','3D hole','2D elec','2D hole')
        xlabel(' \rho, \mum')
        ylabel(' Current \rho, mA')
        
        subplot(132)
        plot(xQW(1:end-1),mode.efield_rho(pcor,:))
        grid on
        xlabel(' \rho, \mum')
        ylabel(' E_\rho, V/cm')
        
        subplot(133)
        yqw=y-modePlot.yQW;
        plot(yqw(1:end-1),mode.efield_z)
        grid on
        xlabel(' z, \mum')
        ylabel(' E_z, V/cm')
        pausak
        
        JTcil=2*pi*xcm.*curx_n+2*pi*xcm.*curx_p+jQn+jQp;
        
        figure(88), plot(xcm*1e4,JTcil),
        title(' Total Current normal to cilinder vs. cil. radius')
        xlabel('\rho, \mum')
        grid on
        pausak
        
    end
    
    fprintf('QW section! \n'),pausak
    
    ze=y(end);
    yo=ze-mode.zox;
    figure(Fall),subplot(222)
    yqw=y-modePlot.yQW;
    
    if isfield(mode,'El2D')
        semilogy(yqw,mode.El(pcor,:)+mode.El2D(pcor,:),yqw,mode.Ho(pcor,:)+mode.Ho2D(pcor,:))
        %  semilogy(yqw,mode.El2D(pcor,:),'--',yqw,mode.Ho2D(pcor,:),'--')
        %  hold off
    else
        semilogy(yqw,mode.El(pcor,:),yqw,mode.Ho(pcor,:))
    end
    %hold on
    ylim([1e17 1e20])
    %xlim(y(end)*[.98 1])
    xlim([-1 1]/4)
    ylabel(' N-P, 1/cm')
    xlabel(' z around QW, \mum ')
    grid
    
    
    
    NQW=mode.NMQW;
    yMQW=mode.yMQW;
    vWMQW=mode.vWMQW;
    
    Zl=[];
    Ind=[];
    %dX=50e-7;
    for kn=1:NQW
        [val,ind]=min(abs(ycm-(yMQW{kn}-vWMQW{kn}/2)));
        Zl=[Zl xqw(ind)];
        Ind=[Ind ind];
        [val,ind]=min(abs(ycm-(yMQW{kn}+vWMQW{kn}/2)));
        ind=ind+1;
        Zl=[Zl xqw(ind)];
        Ind=[Ind ind];
        %zIn=(yMQW{3}-dX)*1e4;
        %zFi=(yMQW{1}+dX)*1e4;
        %zLim=[zIn zFi];
    end
    fP=abs((curr_n)+(curr_p));
    figure(Fall)
    subplot(224)
    %semilogy(z,curr_n*1000,z,curr_p*1000), hold on,
    %semilogy(Zl,curr_n(Ind)*1000,'o',Zl,curr_p(Ind)*1000,'o'),
    %semilogy(z,fP,'g','linewidth',2),
    semilogy(xqw,abs(curr_n),xqw,abs(curr_p)), hold on,
    semilogy(Zl,abs(curr_n(Ind)),'o',Zl,abs(curr_p(Ind)),'o'),
    semilogy(xqw,fP,'g','linewidth',2),
    xlim([-dX,dX])
    hold off
    title(['Current = ',num2str(Cpl(Vcor(ipcor)),2), ' , mA'])
    grid
    %xlim(zLim)
    ylabel(' N-P currents, mA')
    xlabel(' z around QW, nm ')
    
    
    
    subplot(223)
    
    %     nQW=0; pQW=0;
    %     for indQW=1:mesh.NMQW
    %         nQW=nQW+modePlot.nQW{end}{indQW};
    %         pQW=pQW+modePlot.pQW{end}{indQW};
    %     end
    
    plot(xQW,modePlot.nQW{pcor}{2},'b.-',xQW,modePlot.pQW{pcor}{2},'r.-')
    %hold on
    %plot(x,modePlot.Cn(end,:),'c.',x,modePlot.Cn(end,:),'m.')
    %    K>> figure, plot(x,mode.Cp*8e-7*1e-9/9)
    %K>> figure, plot(x,mode.Cn*8e-7*1e-9/9)
    
    xlim([xQW(1),xQW(end)])
    xlabel('\rho, nm')
    ylabel('N-P 2D density 1e12/cm^2)')
    hold off
    grid
    %     pausak
    
    %% TJ
    if isfield(mesh,'yBTJ')==1
        fprintf('TJ section! \n'),pausak
        
        yTJ=y-mesh.yBTJ{2}*1e4;
        xTJ=1000*yTJ;

        figure(Fall+15),
        set(Fall+15,'pos',[141         136        1275         829 ])
        subplot(221)
        plot(xTJ,modePlot.Ec(fipst,:),'linewidth',2)
        ax = gca;
        ax.ColorOrderIndex = 1;
        hold on, plot(xTJ,modePlot.Ev(fipst,:)+1.8,'--','linewidth',2)
        
        
        %figure, plot(y,Rc-Rv)
        dX=80;
        xlim([-dX,dX])
        title(' Ec (cont), Ev+1.8 (dashed) in Cavity')
        xlabel(' z around TJ (nm) ')
        grid
        
        subplot(222)
        if isfield(mode,'El2D')
            semilogy(yTJ,mode.El(pcor,:)+mode.El2D(pcor,:),yTJ,mode.Ho(pcor,:)+mode.Ho2D(pcor,:))
        else
            semilogy(yTJ,mode.El(pcor,:),yTJ,mode.Ho(pcor,:))
        end
        %hold on
        ylim([1e17 1e20])
        %xlim(y(end)*[.98 1])
        xlim([-1 1]/4)
        ylabel(' N-P  (1/cm)')
        xlabel(' z around TJ (\mum) ')
        grid
        
        fprintf('Qui per plot corrente TJ!\n'),keyboard
        NBTJ=length(mesh.BTJcell);
        yBTJ=mesh.yBTJ;
        vBTJ=mesh.vBTJ;
        
        Zl=[];
        Ind=[];
        %dX=50e-7;
        for kn=1:NBTJ
            [val,ind]=min(abs(ycm-(yBTJ{kn}-vBTJ{kn}/2)));
            Zl=[Zl xTJ(ind)];
            Ind=[Ind ind];
            [val,ind]=min(abs(ycm-(yBTJ{kn}+vBTJ{kn}/2)));
            ind=ind+1;
            Zl=[Zl xTJ(ind)];
            Ind=[Ind ind];
            %zIn=(yMQW{3}-dX)*1e4;
            %zFi=(yMQW{1}+dX)*1e4;
            %zLim=[zIn zFi];
        end
        fP=abs((curr_n)+(curr_p));
        figure(Fall+15)
        subplot(224)
        %semilogy(z,curr_n*1000,z,curr_p*1000), hold on,
        %semilogy(Zl,curr_n(Ind)*1000,'o',Zl,curr_p(Ind)*1000,'o'),
        %semilogy(z,fP,'g','linewidth',2),
        semilogy(xTJ,abs(curr_n),xTJ,abs(curr_p)), hold on,
        semilogy(Zl,abs(curr_n(Ind)),'o',Zl,abs(curr_p(Ind)),'o'),
        semilogy(xTJ,fP,'g','linewidth',2),
        xlim([-dX,dX])
        hold off
        title(['Current = ',num2str(Cpl(Vcor(ipcor)),2), ' (mA)'])
        grid
        %xlim(zLim)
        ylabel(' N-P currents (mA)')
        xlabel(' z around TJ (nm) ')
        
        %     % Jz along a r>rTJ (xTJ=yTJ*1e3)
        %     subplot(223)
        %     semilogy(xTJ,abs(curr_n),xTJ,abs(curr_p)), hold on,
        %     semilogy(Zl,abs(curr_n(Ind)),'o',Zl,abs(curr_p(Ind)),'o'),
        %     semilogy(xTJ,fP,'g','linewidth',2),
        %     xlim([-dX,dX])
        %     hold off
        %     title(['Current = ',num2str(Cpl(Vcor(ipcor)),2), ' (mA)'])
        %     grid
        %     %xlim(zLim)
        %     ylabel(' N-P currents (mA)')
        %     xlabel(' z around TJ (nm) ')
        
        %             PlotRadialTJ
        
        clear irho
        irho=input('Want to see radial current details? [1: Yes - ENTER: No]\n');
        if irho==1
            plot_JrhoTJ
        end
    end
end

F101=1;
IncFig=222;

xo=1000*(modePlot.y-yo);
F101=F101+IncFig;
figure
plot(xo,modePlot.EcH(fipst,:),'linewidth',2)

ax = gca;
ax.ColorOrderIndex = 1;
hold on, plot(xo,modePlot.EvH(fipst,:)+1.8,'--','linewidth',2)

F100=F101+1;

dX=180;
xlim([-dX,dX])
title(' Ec (cont), Ev+1.8 (dashed) below Oxide')
xlabel(' z around Oxide (nm) ')
F100=F100+IncFig;
figure(F100)
set(F100,'pos',[ 188          84        1100         893])
subplot(221)
yqw=1000*(y-modePlot.yQW);
[du,iz]=min(abs(yqw));
p=diag(mode.EFc(Vcor,iz))*ones(size(mode.EFc(Vcor,:)));
plot(yqw,mode.EFc(Vcor,:)-p,'-','linewidth',2)
grid
hold on
ax = gca;
ax.ColorOrderIndex = 1;
% bande
%    p=diag(mode.EFv(Vcor,iz))*ones(size(mode.EFv(Vcor,:)));
plot(yqw,mode.Ec(Vcor,:)-p,'.','linewidth',2)
% valenza
%    p=diag(mode.EFv(Vcor,iz))*ones(size(mode.EFv(Vcor,:)));
%  plot(yqw,mode.EFv(Vcor,:)-p,'--','linewidth',2)
xlabel(' z, nm')
ylabel(' Cond. Band diagram related to zQW, eV')
title('Dots Ec, Cont. Ef, both at (r=0), eV')
xlim([-150 150])


subplot(222)
yqw=1000*(y-modePlot.yQW);
[du,iz]=min(abs(yqw));
p=diag(mode.EFv(Vcor,iz))*ones(size(mode.EFv(Vcor,:)));
plot(yqw,mode.EFv(Vcor,:)-p,'-','linewidth',2)
grid
hold on
ax = gca;
ax.ColorOrderIndex = 1;
% bande
%    p=diag(mode.EFv(Vcor,iz))*ones(size(mode.EFv(Vcor,:)));
plot(yqw,mode.Ev(Vcor,:)-p,'.','linewidth',2)
% valenza
%    p=diag(mode.EFv(Vcor,iz))*ones(size(mode.EFv(Vcor,:)));
%  plot(yqw,mode.EFv(Vcor,:)-p,'--','linewidth',2)
xlabel(' z, nm')
ylabel(' Val. Band diagram related to zQW, eV')
title('Dots Ev, Cont. Efv, both at (r=0), eV')
xlim([-150 150])

subplot(223)
semilogy(yqw,mode.El(Vcor,:),'linewidth',2)
hold on
ylim([1e16 2e18])
%xlim(y(end)*[.98 1])
xlim([-150 150])
ylabel(' N  (1/cm)')
xlabel(' z around QW (um) ')
grid

subplot(224)
semilogy(yqw,mode.Ho(Vcor,:),'linewidth',2)
hold on
ylim([1e16 2e18])
%xlim(y(end)*[.98 1])
xlim([-150 150])
ylabel(' P  (1/cm)')
xlabel(' z around QW (\mum) ')
grid

pausak

figure(F101)
set(F101,'pos',[ 188         70        1100         381])
subplot(121)
plot(xQW(1:end-1),mode.efield_rho(Vcor,:),'linewidth',2)
hold on
xlabel(' \rho, \mum')
ylabel(' E_\rho (z=QW), V/cm')
title(['Correnti = ',num2str(CORRENTI)])
grid
subplot(122)
%yqw=y-modePlot.yQW;
plot(yqw(1:end-1),mode.efield_z(Vcor,:),'linewidth',2)
hold on
grid
xlabel(' z, nm')
ylabel(' Ez (r=0), V/cm')
xlim([-150 150])

pausak
F99=F100+5;
F99=F99+IncFig;
figure(F99),
set(F99,'pos',[  118         171        1600         431])
subplot(131)
xcD=xcm*1e4;
plot(xcD,2*pi*xcm.*curx_n,xcD,2*pi*xcm.*curx_p),
hold on
%ax = gca;
%ax.ColorOrderIndex = 1;
plot(xcD,jQn,'--',xcD,jQp,'--','linewidth',2)
legend('3D elec','3D hole','2D elec','2D hole')
title(' @ last current')
xlabel(' \rho, \mum')
ylabel(' Integral over z of Current \rho, mA')
grid

subplot(132)
plot(yqw,squeeze(JN_Y(Vcor,:,1)),'linewidth',2)
hold on
%  ax = gca;
%  ax.ColorOrderIndex = 1;
plot(yqw,squeeze(JP_Y(Vcor,:,1)),'--','linewidth',2)
xlim([-50 +50])
xlabel(' zQW, \mum')
ylabel(' J_z, mA/cm^2')
grid

if isfield(mesh,'yBTJ')
    subplot(133)
    plot(yTJ*1e3,squeeze(JN_Y(Vcor,:,1)),'linewidth',2)
    hold on
    %  ax = gca;
    %  ax.ColorOrderIndex = 1;
    plot(yTJ*1e3,squeeze(JP_Y(Vcor,:,1)),'--','linewidth',2)
    xlim([-50 +50])
    xlabel(' zTJ, \mum')
    ylabel(' J_z, mA/cm^2')
    grid
end

pausak

F1002=F100+10;
XL=[-30 30];
F1002=F1002+IncFig;
vorticeUNO
pausak

F1003=F1002+10;
F1003=F1003+IncFig;
XL=[-300 300];
vorticeUNO1
pausak