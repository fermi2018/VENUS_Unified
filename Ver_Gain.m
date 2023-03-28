global N2lin P2lin N2i P2i T2 L2 indGan
global G0 dgE0 dgH0 rsp0 drspE0 drspH0 Rsp0 dRspE0 dRspH0

%'vedo NUOVA ver', keyboard
%tic
%[N2lin,P2lin,N2i,P2i,T2,L2,indGan]=vet_mat(mode,mesh);

N2lin=modePlot.Elqw(kv,:);
P2lin=modePlot.Hoqw(kv,:);
T2=modePlot.Tqw(kv,:);
L2=modePlot.lambda(1,kv);

% 15/12/2017
% Prima di oggi, facevamo tutto vettoriale (avevamo che L2 era un
% vettore/matrice Nmodi x Nelementi) perché non avevamo capito benissimo
% come usare interpn. Usando una/due lunghezze d'onda, conviene fare come
% ora. Vedi NUOVA_Gain_OLD_BACKUP per come si faceva prima.
% In questo modo è molto più veloce.
% Il codice quindi è un po' raffazzonato perché è stato modificato meno
% possibile, ma si può mettere più elegante (evitare di mettere L2 come
% vettore lungo / matrice). Allo stesso modo, non ha più senso che N2lin,
% P2lin, T2 siano matrici: è sufficiente che siano sempre vettori di
% lunghezza Nelementi.
if mean(N2lin)>1e9

G0=zeros(size(N2lin));
dgE0=G0;
dgH0=G0;
rsp0=G0;
drspE0=G0;
drspH0=G0;
DeltaN0=G0;
%        'ver dent', keyboard
    [vG0,vDeltaN0] =...
        f_InterpGain4Dver(N2lin,P2lin,modePlot.Deltalam+L2,T2);

        indlam=1;
%     toc
    for ind_diag=1:length(N2lin)
        G0(indlam,ind_diag)=vG0(ind_diag,ind_diag,1,ind_diag);
%        dgE0(indlam,ind_diag)=vdgE0(ind_diag,ind_diag,1,ind_diag);
%        dgH0(indlam,ind_diag)=vdgH0(ind_diag,ind_diag,1,ind_diag);
%        rsp0(indlam,ind_diag)=vrsp0(ind_diag,ind_diag,1,ind_diag);
%        drspE0(indlam,ind_diag)=vdrspE0(ind_diag,ind_diag,1,ind_diag);
%        drspH0(indlam,ind_diag)=vdrspH0(ind_diag,ind_diag,1,ind_diag);
        DeltaN0(indlam,ind_diag)=vDeltaN0(ind_diag,ind_diag,1,ind_diag);
    end


else

G0=zeros(size(N2lin));
DeltaN0=zeros(size(N2lin));
end