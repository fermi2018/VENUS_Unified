if indv>1
    JN_X=reshape(mode.Jn_x,mesh.nny,mesh.nnx);
    JN_Y=-reshape(mode.Jn_y,mesh.nny,mesh.nnx);
    JP_X=reshape(mode.Jp_x,mesh.nny,mesh.nnx);
    JP_Y=-reshape(mode.Jp_y,mesh.nny,mesh.nnx);
    
    x=mesh.xgrid;
    z=mesh.ygrid*1e7;
    
    for indy=1:mesh.nny
        
        jn_y=JN_Y(indy,:);
        jp_y=JP_Y(indy,:);
        curr_n(indy)=trapz(x,2.*pi.*x.*jn_y);
        curr_p(indy)=trapz(x,2.*pi.*x.*jp_y);
        
    end
    
    iploleakage=0;
    
    Zl=[];
    Ind=[];  % left and rightmost nodes of each QW; (mesh.inMQW: central QW node)
    for kn=1:NQW
        [val,ind]=min(abs(mesh.ygrid-(mesh.yMQW{kn}-mesh.vWMQW{kn}/2)));
        Zl=[Zl z(ind)];
        Ind=[Ind ind];
        [val,ind]=min(abs(mesh.ygrid-(mesh.yMQW{kn}+mesh.vWMQW{kn}/2)));
        ind=ind+1;
        Zl=[Zl z(ind)];
        Ind=[Ind ind];
        zIn=(mesh.yMQW{3}-50e-7)*1e7;
        zFi=(mesh.yMQW{1}+50e-7)*1e7;
        zLim=[zIn zFi];
    end
    if iploleakage==1
        figure(222)
        plot(z,curr_n*1000,z,curr_p*1000)
        xlim(zLim)
        
        figure(222)
        semilogy(z,curr_n*1000,z,curr_p*1000)
        hold on
        semilogy(Zl,curr_n(Ind)*1000,'o',Zl,curr_p(Ind)*1000,'o')
        xlim(zLim)
        pausak
    end
    
    QWlea=length(Ind)-1;
    Fleak1=1/(curr_n(Ind(QWlea)-1)/curr_p(Ind(QWlea)-1));
    QWlea=2;
    Fleak2=(curr_n(Ind(QWlea)+1)/curr_p(Ind(QWlea)+1));
    Fleak=min([Fleak1 Fleak2]);
    Fleak=sum([Fleak1 Fleak2]);
    mode.Fleak(indv)=Fleak;
else
    mode.Fleak(indv)=0;
end