% figure(50),hold on
figure(51),hold on
% figure(52),hold on

istart=16;
Tsoglia=10;

for iV=istart:length(mode.vv_dd)
    
    iVp=iV+[-5:0];
    
%     xstart=-5;
    polGrad=2;
    xstart=-2-polGrad;

%     xVoltP=mode.vv0_dd(iV+[xstart:-1]);
    xVoltI=mode.vv_dd(iV+[xstart:-1]);
%     yCurrent=mode.ii_dd(iV+[xstart:-1]);
%     yPst=mode.Pst_dd(iV+[xstart:-1]);
    yPdiss=mode.PDiss(iV+[xstart:-1]);
    
%     coeffI=polyfit(xVoltI,yCurrent,2);
%     coeffPst=polyfit(xVoltP,yPst,2);
    
    coeffPdiss=polyfit(xVoltI,log10(yPdiss),1);
    coeffPdissLIN=polyfit(xVoltI,yPdiss,1);

    Pdissfit=10.^polyval(coeffPdiss,mode.vv_dd(iV));
    PdissfitLIN=polyval(coeffPdissLIN,mode.vv_dd(iV));
    
%     mode.Pstfit(iV+1)=Pstfit;
    mode.Pdissfit(iV)=Pdissfit;
    mode.PdissfitLIN(iV)=PdissfitLIN;
    
%     xVolt1 = linspace(xVoltP(1),mode.vv0_dd(iV+1),101);
%     Ifit1 = polyval(coeffI,xVolt1);
%     Pstfit1 = polyval(coeffPst,xVolt1);
%     Pstfit1 = polyval(coeffPst,xVolt1);
        
%    	figure(50),chold
%     plot(-mode.vv0_dd(iVp),mode.Pst_dd(iVp),'r.-','linewidth',3)
%     plot(-xVoltP,yPst,'o','markersize',5)
%     plot(-mode.vv0_dd(19:iV+1),mode.Pstfit(19:iV+1),'*')
% %     plot(-xVolt1,Pstfit1,'k')
%     ylabel('Pst')
%     set(gcf,'Position',[1344         566         566         425])

%     figure(51),chold
%     plot(mode.vv0_dd(iVp),mode.ii_dd(iVp)*1e3,'r.-','linewidth',3)
%     plot(xVoltI,yCurrent*1e3,'o','markersize',5)
%     plot(mode.vv0_dd(23:iV+1),mode.Ifit(23:iV+1)*1e3,'*')
%     plot(xVolt1,Ifit1*1e3,'k')
%     ylabel('Current')
%     set(gcf,'Position',[689    88   560   420])
%     drawnow
    
    figure(51),chold
    plot(abs(mode.vv_dd(iVp)),mode.PDiss(iVp),'r.-','linewidth',3)
    plot(xVoltI,yPdiss,'o','markersize',5)
    plot(mode.vv_dd(istart:iV),mode.PdissfitLIN(istart:iV),'k*')
%     plot(-xVolt1,Pstfit1,'k')
    ylabel('PDiss')
    title('LIN')
    set(gcf,'Position',[1271          81         560         419])
    
    
%     figure(52),chold
%     plot(abs(mode.vv0_dd(iVp)),mode.PDiss(iVp),'r.-','linewidth',3)
%     plot(xVoltI,yPdiss,'o','markersize',5)
%     plot(mode.vv0_dd(istart:iV),mode.Pdissfit(istart:iV),'k*')
% %     plot(-xVolt1,Pstfit1,'k')
%     ylabel('PDiss')
%     title('LOG')
%     set(gcf,'Position',[689    88   560   420])

%     drawnow
%     keyboard

%    if mode.DeltaTmax(iV)>Tsoglia
%        FatTEMP=1./abs(mode.PDiss(2:iV)./mode.PDiss(1:iV-1));
%        
% %        xT=-5;
%        polGradT=1;
%        xT=-3-polGradT;
%        
%        xVT=mode.vv0_dd(iV+[xT:-1]);
%        yT=FatTEMP(iV+[xT:-1]-1);
%        
%        coeffT=polyfit(xVT,yT,1);
%        
%        FatTEMPfit=polyval(coeffT,mode.vv0_dd(iV));
%        
%        mode.FatTEMP(iV)=1/FatTEMPfit;
%         
%        PDissT(iV)=mode.PDiss(iV-1)/FatTEMPfit;
%    end
    pausak
end

figure,
plot(mode.ii_dd(istart:end)*1e3,abs(1-mode.PdissfitLIN(istart:end)./mode.PDiss(istart:end))*100,'.-')
hold on,grid on
plot(mode.ii_dd(istart:end)*1e3,abs(1-mode.Pdissfit(istart:end)./mode.PDiss(istart:end))*100,'.-')

% figure,
% plot(abs(mode.PDiss(istart:end)-PDissT(istart:end)),'.-')
% hold on,grid on
% 
% figure
% hold on,grid on
% plot(1./FatTEMP(10:end))
% plot(mode.FatTEMP(10:end))