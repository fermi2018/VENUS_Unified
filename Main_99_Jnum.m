
% clear
clear geom mesh mode uvet
close all
clc

% load('WorkSpace_Jacobian.mat')
% load('WorkSpace_Jacobian_BTJ.mat')
load(['WorkSpace_Jacobian_',strSave,'_',nomeSav,'.mat'])

nn = mesh.nn;
uvet0 = uvet;

if mode.mp==0
    delta=1e-5;
    [Kmat0,Jmat0,Jmat2,~,rvet0,mode0]=assem_GBT(geom,mesh,mode,uvet,mode.Vnum);    % VENUS
else
    delta=mp('1e10');
    [Kmat0,Jmat0,Jmat2,~,rvet0,mode0]=assem_mp(geom,mesh,mode,uvet,mode.Vnum);    % VENUS (mp)
end

Jmat_anal = Jmat0+Kmat0+Jmat2;

Jmat_num=sparse(length(uvet0),length(uvet0));

% In case of TJ indexes: left index in the {} brackets is the TJ index; 
% the index in the () brackets identifies the column index (radial)
iTJp=mesh.iLBTJ{1,end}(1);  % right BTJ node, first column
iTJn=mesh.iRBTJ{1,1}(1); % left BTJ node, first column
% iTJ=iTJp:iTJn;
iTJ=iTJn;
iTJ=[iTJp iTJn];

i_num=[iTJ iTJ+nn iTJ+2*nn];

% 
% iQWL=min(mesh.inMQW{3});  % right QW node, first column
% iQWR=min(mesh.inMQW{1});  % left QW node, first column
% iQW=iQWL-2:iQWR+2;
% 
% i_num=[iQW iQW+nn iQW+2*nn];

for ind = i_num
%     for ind = 1:length(uvet0)
    
    uvet = uvet0;
    uvet(ind) = uvet(ind)*(1+delta);
    
    if mode.mp==0
        [~,~,~,~,rvet_num,mode_num]=assem_GBT(geom,mesh,mode,uvet,mode.Vnum);    % VENUS
    else
        [~,~,~,~,rvet_num,mode_num]=assem_mp(geom,mesh,mode,uvet,mode.Vnum);    % VENUS (mp)
    end
    
    rvet_inc = rvet_num;
    Jmat_num(:,ind) = (rvet_inc-rvet0)/(delta*uvet0(ind));
    
    figure(19)
    plot(abs(Jmat_num(:,ind)),'bo')
    hold on
    plot(abs(Jmat_anal(:,ind)),'r.')
    title(['Abs - ',num2str(ind),' of ',num2str(3*nn)]),set(gca,'yscale','log')
    set(gcf,'Position',[0 511 640 485])
    hold off
    legend('numerico','analitico')
    
    figure(20)
    plot(abs(1-abs(Jmat_num(:,ind))./abs(Jmat_anal(:,ind))),'gd')
    title(['Relative error - ',num2str(ind),' of ',num2str(3*nn)]),set(gca,'yscale','log')
    set(gcf,'Position',[640 511 640 485])
    hold off
    
    drawnow
    pausak
end


% return

% figure,plot((Jmat_num(:,ind)),'b*'),hold on,plot((Jmat_anal(:,ind)),'ro'),title('Segno +'),set(gca,'yscale','log')
% figure,plot(-(Jmat_num(:,ind)),'b*'),hold on,plot(-(Jmat_anal(:,ind)),'ro'),title('Segno -'),set(gca,'yscale','log')
% 
% figure,plot(abs(Jmat_num(:,ind)),'b*'),hold on,plot(abs(Jmat_anal(:,ind)),'ro'),title('Abs'),set(gca,'yscale','log')
% 
% figure,plot(abs(Jmat_num),'b*'),hold on,plot(abs(Jmat_anal(:,i_num)),'ro'),title('Abs'),set(gca,'yscale','log')


figure,spy(Jmat_num,'bo')
hold on
spy(Jmat_anal,'r.')

