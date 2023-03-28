function [ex,ey,ez,hx,hy,hz]=f_EvalFloquetModes_z(x,k0,csi,ky,d,nindex)

%
% Autore: Renato Orta
% Torino, 02 aprile 2010
% ------------------------------------------------------
% Riga di comando:
%
% [hFloqAn, eFloqAn] = FloqAnComp(x1, d1, n1, x2, d2, n2,Nmodi, polariz);
%
% ------------------------------------------------------
% Funzione per il calcolo dell'andamento dei modi di Floquet(sia
% campo elettrico che magnetico) TE e TM rispetto a z!!!!!!!!!!!!!!!
% ------------------------------------------------------
% Variabili in ingresso:
%	x1 = vettore contenente i valori di x tra 0 e d1;
%	x2 = vettore contenente i valori di x tra d1 e d;
%	d1 = larghezza dello strato 1
%	d2 = larghezza dello strato 2
%	n1 = costante dielettrica strato 1
%	n2 = costante dielettrica strato 2
%	polariz = polarizzazione dell'onda piana incidente
%				= 0 --> TE
%				= 1 --> TM
%
% Variabili in uscita:
%	hFloqAn = matrice contenente 
%		l'andamento delle autofunzioni h:
%  ind. riga --> ordine modo
%  ind. colonna --> valori di x
%	eFloqAn= matrice contenente 
%		l'andamento delle autofunzioni e
%  ind. riga --> ordine modo
%  ind. colonna --> valori di x
%-------------------------------------------------------

Nmodi=length(csi);
kt=f_psqrt(csi.^2+ky^2);
Nmodi=length(csi);
kxsukt=(csi./kt);
kysukt=(ky./kt);

ind0=find(kt==0);
kxsukt(ind0)=1;
kysukt(ind0)=0;

kz=f_psqrt(k0^2*nindex^2-kt.^2);
for ind=1:Nmodi
    Floq=exp(-j*csi(ind)*x)/sqrt(d);
    % indice di riga : modo di Floquet
    % indice di colonna: punto x
    
    exTE(ind,:)=j*kysukt(ind)*Floq;
    eyTE(ind,:)=-j*kxsukt(ind)*Floq;
    ezTE(ind,:)=zeros(1,length(x));
    hxTE(ind,:)=j*kxsukt(ind)*Floq;
    hyTE(ind,:)=j*kysukt(ind)*Floq;
    hzTE(ind,:)=-j*kt(ind)/kz(ind)*Floq;
    exTM(ind,:)=j*kxsukt(ind)*Floq;
    eyTM(ind,:)=j*kysukt(ind)*Floq;
    ezTM(ind,:)=-j*kt(ind)/kz(ind)*Floq;
    hxTM(ind,:)=-j*kysukt(ind)*Floq;
    hyTM(ind,:)=j*kxsukt(ind)*Floq;
    hzTM(ind,:)=zeros(1,length(x));
end
    
ex=[exTE;exTM];
ey=[eyTE;eyTM];
ez=[ezTE;ezTM];
hx=[hxTE;hxTM];
hy=[hyTE;hyTM];
hz=[hzTE;hzTM];



    
