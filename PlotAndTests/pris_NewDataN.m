% clear
clear global
close all
colordef white
dbstop if error

addpathVENUS

INEW=1;

FLO=1;

IMAX=15;
VMAX=2.6;
% CORRENTI=[1 5 10];
CORRENTI=[2 4];
iDeLa=0;
TTT000=273;



N{1}='Tdic2'; %
%N{1}='Tdic1'; %
%N{1}='TvenNots9'; %
%N{1}='Tdic2n'; %
%N{1}='Tdic2n'; %
N{1}='TdicP'; %
N{1}='TdicP1'; %
%N{1}='TdicP2'; %
%N{1}='Tlin'; %
%N{1}='TdicP3'; %


N{1}='Tvar'; %
N{1}='Tok'; %
N{1}='Tok1'; %
N{1}='Tok2'; %

N{1}='Tok4_110'; %
N{1}='Tok5_110'; %
N{1}='Tok5_20'; %
N{1}='Tok_og'; %
N{1}='Tok_og1'; %
N{1}='Tok_ogu'; %
N{1}='Tok_ogu1'; %
N{1}='Tok_oguv'; %
N{1}='Tok_not110'; %
N{1}='dueT'; %
N{1}='tau'; %
N{1}='nlg1'; %
N{1}='tau1'; %
%N{1}='perd'; %

%N{1}='TvarOld'; %
N{1}='plog1'; %
N{1}='fc'; %

N{1}='LD4'; % Los exp
N{1}='LD5'; % Auger val
N{1}='LD4a'; % Auger exp (2)
N{1}='LD5a'; % Auger val
N{1}='LD7'; % Coef Nperd
%N{1}='LD6'; % Fat Perd c10
%N{1}='LD6a'; % Fat Perd c2

N{1}='TpapN3'; %
%N{1}='DueT'; %


N{1}='nyT'; %
N{1}='nyT20d'; %
N{1}='nyT20e'; %
N{1}='nyT20e1'; %
N{1}='nyT110f'; %
N{1}='nyT110h'; %
N{1}='nyT110k'; %

%N{1}='nyT20g'; %

N{1}='nyT20kL'; %

N{1}='VenT20'; %
N{1}='VenT20Dif'; %
N{1}='VenT20Cter'; %

N{1}='VenT110per'; %
N{1}='VenT110ex'; %

N{1}='VenT110TERn2';
N{1}='VenT4';

N{1}='VenTnotte1';
N{1}='VenTnotte3';
N{1}='LunTnlg';
%N{1}='VenT110TERref'


IncFig=1;
F1002=1001;
F1003=2001;
F99=98;
F101=100;
F100=99;
Fcur=9;
Fall=19;

iLoad=input(' Vuoi ultima simulazione? [Enter: NO;  0: precedente] ')

% if length(iLoad)==0
%     load ultimonome
%     % save nomeLoopFinito nomeSave
%     fibar=strfind(nomeSave,'\');
%     nomeSave=nomeSave(fibar(end)+1:end);
%     %keyboard
%     eval(['load  ',nomeSave,'_pmat'])
%     fibar=strfind(nomeSave,'\');
%     nomeSave=nomeSave(fibar(end)+1:end);
% else
%     load penultimonome
%     fibar=strfind(nomeSW,'\');
%     nomeSW=nomeSW(fibar(end)+1:end);
%     % save nomeLoopFinito nomeSave
%     nomeSave=[nomeSW,N{1}];
%     eval(['load  ',nomeSave,'_pmat'])
%     fibar=strfind(nomeSave,'\');
%     nomeSave=nomeSave(fibar(end)+1:end);
% end


nomeSave
pausak
nome=nomeSave;

iuno=1;

%IPvet=9

for IPAR=IPvet
    
    TMa=[];
    
%     clear MODEplot
    
    if IPAR==4
        Fpar=1e-18;
    else
        Fpar=1;
    end
    % eval(['load All',num2str(IPAR),'.mat MODE MESH VInf VInp'])
%     eval(['exi=exist( ''',nome,num2str(IPAR),'.mat',''');'])
    
    if exi==0
        return
    end
%     eval(['load ',nome,num2str(IPAR)])
    
    
    % keyboard
    colo='rgbcmrb';
    %col{1}='k'; col{2}='r'; col{3}='g'; col{4}='b'; col{5}='c'; col{6}='m'; col{7}='r.'; col{8}='b.';
    col{1}='r'; col{2}='g'; col{3}='b'; col{4}='c'; col{5}='m'; col{6}='k'; col{7}='y';
    if IPAR>0
        ppar=PMAT{IPAR}
    else
        ppar=[];
    end
    %keyboard
    hh=figure;
    set(hh,'position',[263          89        1517         891])
    
    
    %  legend(TITE{1},'location','best')
    %hh1=figure;
    %set(hh1,'position',[ 400 50 1400 400])
    %'ver', keyboard
    for kpar=1:length(MODEplot)
        kpar
        mode=MODEplot{kpar};
        if isfield(mode,'FLos')
            FLos=mode.FLos;
        else
            FLos=FLO;
        end
        modePlot=MODEplot{kpar};
        if isfield(modePlot,'Isize')
            Isize=modePlot.Isize;
            T0=modePlot.T0-TTT000;
            %    'prima di load', keyboard
            %fibar=strfind(nomeSR,'\');
            %nomeSR=nomeSR(fibar(end)+1:end);
            if INEW==1
                eval(['load MarkusN_',num2str(Isize),'_T',num2str(T0),'.mat'])
            else
                eval(['load Markus_',num2str(Isize),'_T',num2str(T0),'.mat'])
            end
            PDB{kpar}=P_dB;
            CU{kpar}=Cur;
        else
            if irel==1
                eval(['load ',nomeSR,'Meas_ReliefVCSEL.mat'])
            else
                eval(['load ',nomeSR,'Meas_StandardVCSEL.mat'])
            end
        end
        CUR=Cur;
        Rmeas=diff(Vmeas)./diff(Imeas/1000);
        Rmeas=[Rmeas(1:end-2)];
        Imeas_res=Imeas(1:end-3);
        Imax=max(Imeas);
        Imax=IMAX;
        Pmax=max(Lmeas)*1.2;
        subplot(231)
        pa=1:5:length(Imeas);
        plot(Vmeas(pa),Imeas(pa),[colo(kpar),'.'],'markersize',10)
        hold on
        
        
        subplot(233)
        plot(Imeas(pa),Lmeas(pa),[colo(kpar),'.'],'markersize',10)
        hold on
        xlim([0 IMAX])
        %'tit', keyboard
        title(['settings-vari-',mode.settings(2:end)])
        
        if IPAR>0
            
            if kpar<length(ppar)
                tita=[tit{IPAR},'  Valore attuale= ',num2str(Fpar*ppar(kpar))];
            else
                tita=tit{IPAR};
            end
            % 'ver', keyboard
            % [fPES Gamma_z fPdif  E2 lambda TempVELM anti_gui]
            
            if IPAR==30
                titeff=['EFFETTI: ',TITE{PMAT{IPAR}(kpar)+2}]
            else
                titeff=[''];
            end
        else
            tita='Set 0';
            titeff='Set 0';
        end
        %mesh=MESH{kpar};
        I=mode.ii_dd*1000;
        %'pa', keyboard
        V=mode.vv_dd;
        P=sum(mode.Pst_dd,1)+mode.Psp_dd;
        Imax=max([IMAX max(I)]);
        Pmax=max([max(Lmeas) max(P)]) *1.1 ;
        Imod{kpar}=modePlot.ii_dd*1000;
        Pmod{kpar}=10*log10(modePlot.Pst_dd);
        
        figure(hh)
        
        subplot(231)
        hold on
        grid on
        box on
        axis([1.4 VMAX 0 IMAX])
        plot(V,I,[col{kpar}],'linewidth',2)
        ylabel('Current (mA)')
        xlabel('Voltage (V)')
        title(tita)
        % figure(11)
        subplot(233)
        hold on
        grid on
        box on
        if kpar==1
            axis([0 Imax 0 Pmax*1.2])
            Pold=Pmax;
        else
            if Pmax>Pold
                axis([0 Imax 0 Pmax*1.2])
            end
        end
        plot(I,P,[col{kpar},'.-'],'linewidth',2)
        
        xlabel('Current, mA')
        ylabel('Optical power, mW')
        
        %legend(num2str([1:kpar]'),'location','best')
        %  legend([titeff0, TITE{PMAT{IPAR}(1:kpar)+1}],'location','best')
        
        subplot(232)
        hold on
        grid on
        box on
        %            xlim([0 20])
        len=1:length(mode.ii_dd);
        %     plot(1000*mode.ii_dd,1e-12*mode.nMaxVet(len),[colo(kpar),'o-'],...
        %         1000*mode.ii_dd,1e-12*mode.pMaxVet(len),[colo(kpar),'+--'])
        plot(1000*mode.ii_dd,1e-12*mode.nMaxVet(len)/mode.NMQW,[colo(kpar),'o-'],...
            1000*mode.ii_dd,1e-12*mode.pMaxVet(len)/mode.NMQW,[colo(kpar),'+--'])
        hold on
        plot(1000*mode.ii_dd,1e-18*mode.n3MaxVet(len),[colo(kpar),'-'],...
            1000*mode.ii_dd,1e-18*mode.p3MaxVet(len),[colo(kpar),'--'],'linewidth',2)
        xlabel('Current (mA)')
        ylabel('N-P density (PER well) 1e12/cm^2)')
        legend('Electrons','Holes','location','best')
        if IPAR==30
            if kpar<=ppar(end)
                legend(TITE{[1 ppar(1:kpar)+2]},'location','best')
            else
                legend(TITE{[1 ppar(1):ppar(end)+2]},'location','best')
            end
        else
            titeff=nomeSave;
        end
        if length(titeff)>20
            title(titeff(15:end))
        else
            title(titeff)
        end
        % title(tit{IPAR})
        %figure(hh1)
        subplot(236)
        hold on
        grid on
        box on
        mI=mean(1000*mode.ii_dd);
        [du,iMii]=min(abs(1000*mode.ii_dd-mI));
        mI=1000*mode.ii_dd(iMii);
        [du,iMi]=min(abs(CUR-mI));
        LaC=mode.lambda(1,iMii);
        LaX=LAM(iMi,1);
        
        if iDeLa==1
            DeLam=LaX-LaC;
        else
            DeLam=0
            '  !!!!!!!!!!!!!!!!!!!!!! '
        end
        if isnan(LaX)
            DeLam=0;
        end
        TMa=[TMa mode.DeltaTmax(end)];
        fic=find(1000*mode.ii_dd>.1);
        plot(1000*mode.ii_dd(fic),mode.lambda(1,fic)+DeLam,col{kpar},'linewidth',2),
        title(['DeLam = ',num2str(DeLam),'    TMa ',num2str(TMa,'%0.f- ')])
        plot(CUR,LAM(:,1),[colo(kpar),'o'],'Markersize',5)
        xlabel('Current, mA')
        ylabel('Wavelength, nm')
        subplot(234)
        %         hold on
        %          grid on
        %          box on
        
        Resistence=diff(mode.vv0_dd)./diff(mode.ii_dd);
        Resistence=[Resistence(1),Resistence];
        % plot(mode.ii_dd*1000,Resistence,Imeasinterp,Resmeasinterp)
        %            plot(Imeasinterp,Resmeasinterp,'r-','LineWidth',2)
        pu=1:4:length(Rmeas);
        semilogy(Imeas_res(pu),Rmeas(pu),[colo(kpar),'o'],'Markersize',3)

        hold on
        lres=1:length(Resistence);
        %semilogy(mode.ii_dd(lres)*1000,Resistence(lres),[colo(kpar),'.-'],'Markersize',10)
        semilogy(mode.ii_dd(lres)*1e3,Resistence(lres),[colo(kpar),'.-'],'linewidth',2)
        
%         ylim([30,500])
        ylim([30,5e3])
        grid on
        xlim([0 Imax])
        %            xlim([0.5,1000*mode.ii_dd(end)+.1])
        xlabel('Current, mA')
        ylabel('Differential resistence, \Omega')
        %            keyboard
        subplot(235)
        hold on
        grid on
        box on
        if kpar==1
            CurSav=Cur;
        end
        if FLos==1
            %figure(111)
            
            Cur=CurSav;
            
            [du,icu]=min(abs(Imeas-Cur(fix(end/2))));
            Smis='o+s';
            k=kpar;
            if IPAR==11 || IPAR==41
                
                P_dB=PDB{k};
                Cur=CU{k};
                P_lin=10.^(P_dB/10);
                Pls=P_lin(:,1);
                for kp=2:size(P_dB,2)
                    fi=isnan(P_lin(:,kp))==0;
                    if length(fi)>0
                        Pls(fi)=Pls(fi)+P_lin(fi,kp);
                    end
                end
                FatPow=Lmeas(icu)/Pls(fix(end/2));
                scLo=10*log10(FatPow);
                plot(Cur,P_dB+scLo,[Smis(3),'--'])
%                 plot(Cur,P_dB(:,1)+scLo,[Smis(k),'--'])
                hold on
                ax = gca;
                ax.ColorOrderIndex = 1;
            else
                P_lin=10.^(P_dB/10);
                Pls=P_lin(:,1);
                for kp=2:size(P_dB,2)
                    fi=isnan(P_lin(:,kp))==0;
                    if length(fi)>0
                        Pls(fi)=Pls(fi)+P_lin(fi,kp);
                    end
                end
                FatPow=Lmeas(icu)/Pls(fix(end/2));
                scLo=10*log10(FatPow);
                
%                 plot(Cur,P_dB+scLo,'o','LineWidth',2)
                plot(Cur,P_dB(:,1)+scLo,'o','LineWidth',2)
                hold on
            end
            
            
            ax = gca;
            ax.ColorOrderIndex = 1;
            lines{1}='-'; lines{2}='--'; lines{3}='-.';  lines{4}='.';  lines{5}='+';  lines{6}='o';
            lines{7}='p'; lines{8}='s';
            plot(Imod{k},Pmod{k},lines{mod(k,length(lines))},'linewidth',2)
            hold on
            ax = gca;
            ax.ColorOrderIndex = 1;
            %axis([0 Imod{end}(end)+1 -40 5])
            axis([0 15 -7 5])
            xlabel('Current, mA')
            ylabel('Optical power, mW')
            
        else
            
            plot(mode.ii_dd*1000,mode.DeltaTmax,col{kpar},'LineWidth',2)
            xlabel('Current, mA')
            ylabel('Temperature rise, K')
        end
        pausak
        
%         PlotDyn
        Dyn_Comparison
        
        iDDPlot=input(' Vuoi DD_PLOT? [0, NO, Enter, Si] ')
        
        if length(iDDPlot)==1
            DD_PLOT_pris
            pausak
            % close(1229:1231)
            % close(1234:1238)
        end
    end
    pausak
    
    Cur=CurSav;
    
    [du,icu]=min(abs(Imeas-Cur(fix(end/2))));
    
    figure(111)
    Smis='o+s';
    if IPvet==11 || IPvet==41
        
        for k=1:length(MODEplot)
            P_dB=PDB{k};
            Cur=CU{k};
            P_lin=10.^(P_dB/10);
            Pls=P_lin(:,1);
            for kp=2:size(P_dB,2)
                fi=isnan(P_lin(:,kp))==0;
                if length(fi)>0
                    Pls(fi)=Pls(fi)+P_lin(fi,kp);
                end
            end
            FatPow=Lmeas(icu)/Pls(fix(end/2));
            scLo=10*log10(FatPow);
            plot(Cur,P_dB+scLo,[Smis(3),'--'])
%             plot(Cur,P_dB+scLo,[Smis(k),'--'])
            hold on
            ax = gca;
            ax.ColorOrderIndex = 1;
        end
    else
        P_lin=10.^(P_dB/10);
        Pls=P_lin(:,1);
        for kp=2:size(P_dB,2)
            fi=isnan(P_lin(:,kp))==0;
            if length(fi)>0
                Pls(fi)=Pls(fi)+P_lin(fi,kp);
            end
        end
        FatPow=Lmeas(icu)/Pls(fix(end/2));
        scLo=10*log10(FatPow);
        
        plot(Cur,P_dB+scLo,'o','LineWidth',2)
        hold on
    end
    
    
    ax = gca;
    ax.ColorOrderIndex = 1;
    lines{1}='-'; lines{2}='--'; lines{3}='-.';  lines{4}='.';  lines{5}='+';  lines{6}='o';
    for k=1:length(MODEplot)
        plot(Imod{k},Pmod{k},lines{mod(k,6)},'linewidth',2)
        hold on
        ax = gca;
        ax.ColorOrderIndex = 1;
    end
    %axis([0 Imod{end}(end)+1 -40 5])
    axis([0 15 -7 5])
    xlabel('Current, mA')
    ylabel('Optical power, mW')
    pausak
    % close all
end