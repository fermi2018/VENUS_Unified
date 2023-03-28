
%
% EQ P generation
%decP=Na.*repe; % p doping of each layer is multiplied by the #repetitions
%'qui basta impilamento Nuomo', keyboard

decP=Na.*repmat(repe(:,1),1,size(Na,2));
fieq=find(decP(:,1)<0); % Finds the p-side reps by looking for NEGATIVE values
if length(fieq)>0
    % Repetition layers are NOT removed, but simply are made all equal!
    Na(fieq,:)=ParMore.mireq.p.Dop_eq;
    x(fieq,:)=ParMore.mireq.p.x_eq;
    mesh(fieq,:)=-ParMore.mireq.p.Nmesh;
end
%
% Similar procedure for n-side repetitions
% EQ N - BOTTOM DBR
%decN=Nd.*repe;  % identifies the eq n by looking for negative values 
decN=Nd.*repmat(repe(:,1),1,size(Nd,2));
fieq=find(decN(:,1)<0);
% Similar procedure already used to identify the BOTTOM n DBR 
difi=diff(fieq);          % understand the position of "holes" in various repetitions
fdis=find(difi>1);      % understand the position of "holes" in various repetitions
if length(fieq)>0
    if isfield(mireq,'n2')
        pip=0;
        for kk=1:length(fdis)
            fiS{kk}=fieq(pip+1:fdis(kk));
            pip=fdis(kk);
            np(kk)=repe(fiS{kk}(1));
        end
        kk=kk+1;
        fiS{kk}=fieq(pip+1:end); % cells containing the #repetitions of each "block"
        np(kk)=repe(fiS{kk}(1));     % # of repetition of each block identified in fiS
        [~,npMa]=max(abs(np));  % idenntifies the max number of repetitions and its position in np
        fieq=fiS{npMa};       % indeces of the str layers that repeat the most among the p-doped ones
    end
    Nd(fieq,:)=ParMore.mireq.n.Dop_eq;
    x(fieq,:)=ParMore.mireq.n.x_eq;
    mesh(fieq,:)=-ParMore.mireq.n.Nmesh;
end

if isfield(mireq,'n2')
    % Similar procedure for n-side repetitions
    % EQ N - TOP DBR
    decN=Nd.*repmat(repe(:,1),1,size(Nd,2));	
%    decN=Nd.*repe;
    fieq=find(decN(:,1)<0);

    difi=diff(fieq);        % understand the position of "holes" in various repetitions
    fdis=find(difi>1);      % understand the position of "holes" in various repetitions
    if length(fieq)>0
        pip=0;
        for kk=1:length(fdis)
            fiS{kk}=fieq(pip+1:fdis(kk));
            pip=fdis(kk);
            np(kk)=repe(fiS{kk}(1));
        end
        kk=kk+1;
        fiS{kk}=fieq(pip+1:end); % cells containing the #repetitions of each "block"
        np(kk)=repe(fiS{kk}(1));     % # of repetition of each block identified in fiS
        [~,npMa]=min(abs(np));  % idenntifies the max number of repetitions and its position in np
        fieq=fiS{npMa};       % indeces of the str layers that repeat the most among the p-doped ones

        Nd(fieq,:)=ParMore.mireq.n2.Dop_eq;
        x(fieq,:)=ParMore.mireq.n2.x_eq;
        mesh(fieq,:)=-ParMore.mireq.n2.Nmesh;
    end
end

%'VER', keyboard
if flgStop==1
'qui dopo', keyboard
end

%fiCavet=d<0;    % identifies the cavity (i.e., the LcC or LgC layers)
%fiCav=find(d<0);
%d=abs(d);       % brings the values of d to positive once that the cavity has been identified        

if iCol==1
fiCavet=d<0;
fiCav=find(d<0);
end
% fiCav=fiCav([1 end]);
% fiCavet=fiCavet([1 end]);
%' qui Pmat ', keyboard

d=abs(d);       % brings the values of d to positive once that the cavity has been identified        

%' qui ', keyboard

% Generation of the StrDD from the data coming from each layer in the str
X=[];
Ra=[];
D=[];
NA=[];
ND=[];

X_i=[];
Ra_i=[];
D_i=[];
NA_i=[];
ND_i=[];
Mesh_i=[];


Xi=[];
Di=[];
NAi=[];
NDi=[];

Xc=[];
Dc=[];
NAc=[];
NDc=[];

fiCA=[];

kink=1;

ks=1;
NPA=1042;
repip=abs(rep(ks,1));
while ks<=length(d)
    repi=abs(rep(ks,:));
	if ks>NPA
	repi
    ks,  pausak
	end
    Xl=[];
    Ral=[];
    Dl=[];
    NAl=[];
    NDl=[];
    Xli=[];
    Rali=[];
    Mli=[];
    Dli=[];
    NAli=[];
    NDli=[];
    
    if repi(1)>1    % check wheter the layer has repetitions or not
        
        for kk=1:repi(2)
		    k=ks+kk-1;		
%			'entro rep', keyboard
%            repi=abs(rep(k,:)); % vector containing the #reps and #layers involved
            di=d(k);    % length of k-th layer
            Srep=sign(rep(k,1));    % sign of repetetion, NEGATIVE in case of eq layer
            ragi=ragv(k,:)*Srep;
            na=Na(k,:); % doping levels of the k-th layer (already imposed at mireq value)
            nd=Nd(k,:);
            xi=x(k,:);  % molar fraction of the k-th layer
            mi=(mesh(k));   % number of mesh points, sign included
            % Assign these values to the matrices
            Mli=[Mli; mi];
            Dli=[Dli; di];
            Xli=[Xli; xi];
            Rali=[Rali; ragi];
            NAli=[NAli; na];
            NDli=[NDli; nd];
            
            ndisc=abs(mesh(k));
            ndiscS=sign(mesh(k));   % takes the sign of the mesh points: NEGATIVE for eq layers
            % Assign values to the xmol matrices checking for graded layers
            if xi(1)==xi(2)
                Xl=[Xl; ndiscS*ones(ndisc,1)*xi(1)];
            else
                Xl=[Xl; ndiscS*linspace(xi(1),xi(2),ndisc)'];
            end
            Dl=[Dl; ones(ndisc,1)*di/ndisc];    % length of each mesh element
            % Assign values to the doping matrices checking for graded layers
            if na(1)==na(2)
                NAl=[NAl; ones(ndisc,1)*na(1)]; % differently from the mesh, always POSITIVE
            else
                NAl=[NAl; linspace(na(1),na(2),ndisc)'];
            end
            if na(1)==na(2)
                NDl=[NDl; ones(ndisc,1)*nd(1)];
            else
                NDl=[NDl; linspace(nd(1),nd(2),ndisc)'];
            end
            
            % Check if the repetition belong or not to the eq medium
            if rep(k,2)>0   % Positive: NO eq
                %     for iext=1:repi(2)
                %      for iint=1:repi(1)
                Pdu=ParDD.mat{pun(k)};
%                Pdu=ParDD.mat{pun(k)}{1};
                Pmat{kink,1}=Pdu;
%                Pmat{kink,1}{1}=Pdu;
                kink=kink+1;
                %      end
                %     end
            else            % Negative: repetitions BELONG to equivalent medium
                if isop==1  % enters here if an actual mireq.p exists
				   for kkk=1:length(Pmat{pun(k)})
                    Pmat{kink,1}{kkk}=ParMore.mireq.p.Mateq;
				   end 	
                    kink=kink+1;
                    isop=0;
                elseif isop==-1
				   for kkk=1:length(Pmat{pun(k)})
                    Pmat{kink,1}{kkk}=ParMore.mireq.n.Mateq;
				   end 				
                    kink=kink+1;
                    isop=0;
                end
            end
%			Pmat{end}
%							'Pmat 1', keyboard
            %       'ver0', keyboard
%            k=k+1;  % k is increased since kk loops over layers of the repetition
            repip=repi(1);
            
        end  % end of kk loop (loop over the number of layers in the repetition)
        kinkS=kink-1;
        ks=ks+repi(2);       
        if rep(k-1,2)>0 % Check again if the last kk layer (k-1) belongs to an eq medium
            for iext=1:repi(1)-1    % looop over the number of layers in the repetition (repi(1))
%							'Pmat 2', keyboard
                for iint=kinkS-repi(2)+[1:repi(2)]
                    Pmat{kink,1}=Pmat{iint,1};  % Assign material to the layers between kink and k (kink=k-kk)
                    kink=kink+1;
                end
            end
%			Pmat{end}
%							'Pmat 2', keyboard
        else
            isop=-1;
        end
        %            'ver0', keyboard
        
        Dc=[Dc; Dl];
        Xc=[Xc; Xl];
        NDc=[NDc; NDl];
        NAc=[NAc; NAl];
        
        Dlt=repmat(Dl,repi(1),1);
        Xlt=repmat(Xl,repi(1),1);
        NDt=repmat(NDl,repi(1),1);
        NAt=repmat(NAl,repi(1),1);
        
        Dlti=repmat(Dli,repi(1),1);
        Mlti=repmat(Mli,repi(1),1);
        Xlti=repmat(Xli,repi(1),1);
        Rati=repmat(Rali,repi(1),1);
        NDti=repmat(NDli,repi(1),1);
        NAti=repmat(NAli,repi(1),1);
        
        D=[D; Dlt];
        X=[X; Xlt];
        ND=[ND; NDt];
        NA=[NA; NAt];
        
        
%        'ver NUOVO', keyboard
        
        Mesh_i=[Mesh_i; Mlti];
        D_i=[D_i; Dlti];
        X_i=[X_i; Xlti];
        Ra_i=[Ra_i; Rati];
        ND_i=[ND_i; NDti];
        NA_i=[NA_i; NAti];
        
        %   'qui', keyboard

%           'qui periodico NUOVO', keyboard        
		
    else % in case of no repetition of the layer k
	    k=ks;	
		if ks>NPA
		 'unostrato', keyboard
		end
        di=d(k);    % length of layer k
        Srep=sign(rep(k,1));    % sign of the repetition: if NEGATIVE, it is an eq. layer!
        ragi=ragv(k,:)*Srep;    % radii of the layer
        na=Na(k,:);             % doping of the k-th layer
        nd=Nd(k,:);
        xi=x(k,:)*Srep;         % molar fraction of the k-th layer
%		ks
%		xi
%		keyboard
        ndisc=abs(mesh(k));     % number of mesh points in the k-th layer
        mi=mesh(k);             % includes also the sign to understand if it belongs to an eq layer!

        Mli=[Mli; mi];
        Dli=[Dl; di];
        Xli=[Xli; xi];
        NAli=[NAli; na];
        NDli=[NDli; nd];
        % Filled with the same quanities as previous ones
        Mesh_i=[Mesh_i; mi];
        D_i=[D_i; di];
        X_i=[X_i; xi];
        Ra_i=[Ra_i; ragi];
        ND_i=[ND_i; nd];
        NA_i=[NA_i; na];
        
        % The material of the kink-th (k-th) layer is assigned from
        % ParDD.mat, which similarly to original Na,Nd,x,... has its
        % original length: however pun indeces are needed to "double" take
        % the values of the eq layers
        Pmat{kink,1}=ParDD.mat{pun(k)};
%						'Pmat 3', %keyboard
        kink=kink+1;
        
        % xmol is assigned, based on the number of mesh points (ndisc)
%		'X', keyboard
%		'X', keyboard
        if xi(1)==xi(2) % region with flat xmol
            X=[X; ones(ndisc,1)*xi(1)];
            Xc=[Xc; ones(ndisc,1)*xi(1)];
        else            % region with graded xmol
            X=[X; linspace(xi(1),xi(2),ndisc)'];
            Xc=[Xc; linspace(xi(1),xi(2),ndisc)'];
        end
		
%		X, 		keyboard
        % thickness of each mesh element is assigned, based on the number of mesh points (ndisc)
        D=[D; ones(ndisc,1)*di/ndisc];
        Dc=[Dc; ones(ndisc,1)*di/ndisc];
        if na(1)==na(2) % region with flat Na
            NA=[NA; ones(ndisc,1)*na(1)];
            NAc=[NAc; ones(ndisc,1)*na(1)];
        else            % region with graded Na
            NA=[NA; linspace(na(1),na(2),ndisc)'];
            NAc=[NAc; linspace(na(1),na(2),ndisc)'];
        end
        if nd(1)==nd(2) % region with flat Nd
            ND=[ND; ones(ndisc,1)*nd(1)];
            NDc=[NDc; ones(ndisc,1)*nd(1)];
        else            % region with flat Nd    
            ND=[ND; linspace(nd(1),nd(2),ndisc)'];
            NDc=[NDc; linspace(nd(1),nd(2),ndisc)'];
        end
        ks=ks+1;
        %   k,   pausak
    end
    
    %pausak
    
    
end

if flgStop==1
 'dopp', keyboard
 'dopp', keyboard
end

% riduco per struttura a PUNTI
Ds=cumsum(D);
veq=X.*(NA+ND); % multiplies doping (always positive) and xmol (negative in eq medium)
fineq=find(veq>=0); % looks for region w/out eq medium

% eq p (mesh)
veq=X.*NA;  % multiplies p-doping and xmol (negative in eq medium)
fieqA=find(veq<0);  % look for eq p

if length(fieqA)>0
    if fieqA(1)>1
        din_A=Ds(fieqA(1)-1);   % distance from origin of eq p
    else
        din_A=Ds(fieqA(1));
    end
    % multiplies Nadd (always positive) to repe, whose entries are negative in case of eq p 
    % veq now is not depending on the mesh, but on the overall #layers
    veq=repe.*Nadd; 
    fiA=find(veq(:,1)<0);   % find eq p LAYERS indeces
    LA=sum(d(fiA).*abs(repe(fiA,1)));   % Total LENGTH of the eq p
    DLa=LA/ParMore.mireq.p.Nmesh;   % length of each element of eq p mesh
    dLa_A=din_A+[DLa:DLa:LA]';  % creates the new mesh for eq p!
else
    dLa_A=[];
    DLa=0;
end


% Very similar procedure of eq p!
% eq n (mesh) - BOTTOM

if flgStop==1
'Unified LITO', keyboard	
'VEQ 0', keyboard
end




veq=X.*ND;
fieqD=find(veq<0);  % look for eq n

% Similar procedure already used to identify the BOTTOM n DBR 
difi=diff(fieqD);          % understand the position of "holes" in various repetitions
fdis=find(difi>1);      % understand the position of "holes" in various repetitions
%fdis=fdis(2);

%'dLa_D   NUOVO   ********', keyboard

%CCCCCCCCC
%if length(fieqD)>0

if length(fieqD)>0
    if isfield(mireq,'n2')
        pip=0;
        for kk=1:length(fdis)
            fiS{kk}=fieqD(pip+1:fdis(kk));
            pip=fdis(kk);
            np(kk)=veq(fiS{kk}(1));
        end
        kk=kk+1;
        fiS{kk}=fieqD(pip+1:end); % cells containing the #repetitions of each "block"
        np(kk)=veq(fiS{kk}(1));     % # of repetition of each block identified in fiS
%         [~,npMa]=max(abs(np));  % idenntifies the max number of repetitions and its position in np
        npMa=length(np); % idenntifies the max number of repetitions and its position in np
        fieqD=fiS{npMa};       % indeces of the str layers that repeat the most among the p-doped ones
        %
    end
    din_D=Ds(fieqD(1)-1);   % distance from origin of eq n
    
    if flgStop==1
	'KK',keyboard
    end
    %
%    veq=repe.*Nddd;
    veq=Nddd.*repmat(repe(:,1),1,size(Nddd,2));		
    if isfield(mireq,'n2')
        fieqD=find(veq(:,1)<0);  % look for eq n
        
        difi=diff(fieqD);          % understand the position of "holes" in various repetitions
        fdis=find(difi>1);      % understand the position of "holes" in various repetitions
        if length(fieqD)>0
            pip=0;
            for kk=1:length(fdis)
                fiS{kk}=fieqD(pip+1:fdis(kk));
                pip=fdis(kk);
                np(kk)=repe(fiS{kk}(1));
            end
            kk=kk+1;
            fiS{kk}=fieqD(pip+1:end); % cells containing the #repetitions of each "block"
            np(kk)=repe(fiS{kk}(1));     % # of repetition of each block identified in fiS
            [~,npMa]=max(abs(np));  % idenntifies the max number of repetitions and its position in np
            fiD=fiS{npMa};       % indeces of the str layers that repeat the most among the n-doped ones
        end
    else
        fiD=find(veq(:,1)<0);
    end
    %
    LD=sum(d(fiD).*abs(repe(fiD,1)));   % Total LENGTH of the eq n
    DLd=LD/ParMore.mireq.n.Nmesh;   % length of each element of eq n mesh
    dLa_D=din_D+[DLd:DLd:LD]';  % creates the new mesh for eq n!
else
    dLa_D=[];
    DLd=0;
end

if flgStop==1
'find VEQ0', keyboard
'VEQ', keyboard
end

% eq n - TOP
veq=X.*ND;
fieqD=find(veq<0);  % look for eq n

% Similar procedure already used to identify the TOP n DBR 
difi=diff(fieqD);          % understand the position of "holes" in various repetitions
fdis=find(difi>1);      % understand the position of "holes" in various repetitions
if isfield(mireq,'n2')
    if length(fieqD)>0
        pip=0;
        for kk=1:length(fdis)
            fiS{kk}=fieqD(pip+1:fdis(kk));
            pip=fdis(kk);
            np(kk)=veq(fiS{kk}(1));
        end
        kk=kk+1;
        fiS{kk}=fieqD(pip+1:end); % cells containing the #repetitions of each "block"
        np(kk)=veq(fiS{kk}(1));     % # of repetition of each block identified in fiS
        %     [~,npMa]=max(abs(np));  % idenntifies the max number of repetitions and its position in np
        npMa=1; % idenntifies the max number of repetitions and its position in np
        fieqD=fiS{npMa};       % indeces of the str layers that repeat the most among the p-doped ones
        %
        din_D2=Ds(fieqD(1)-1);   % distance from origin of eq n
        %
        REPE=repmat(repe(:,1),1,size(Nd,2))	;
        veq=REPE.*Nddd;				
%        veq=repe.*Nddd;
        fieqD=find(veq(:,1)<0);  % look for eq n
        
        difi=diff(fieqD);          % understand the position of "holes" in various repetitions
        fdis=find(difi>1);      % understand the position of "holes" in various repetitions
        if length(fieqD)>0
            pip=0;
            for kk=1:length(fdis)
                fiS{kk}=fieqD(pip+1:fdis(kk));
                pip=fdis(kk);
                np(kk)=repe(fiS{kk}(1));
            end
            kk=kk+1;
            fiS{kk}=fieqD(pip+1:end); % cells containing the #repetitions of each "block"
            np(kk)=repe(fiS{kk}(1));     % # of repetition of each block identified in fiS
            [~,npMa]=min(abs(np));  % idenntifies the max number of repetitions and its position in np
            fiD=fiS{npMa};       % indeces of the str layers that repeat the most among the n-doped ones
            
        end
        
        %
        LD2=sum(d(fiD).*abs(repe(fiD,1)));   % Total LENGTH of the eq n
        DLd2=LD2/ParMore.mireq.n2.Nmesh;   % length of each element of eq n mesh
        dLa_D2=din_D2+[DLd2:DLd2:LD2]';  % creates the new mesh for eq n!
    else
        dLa_D2=[];
        DLd2=0;
        
    end
else
    dLa_D2=[];
    DLd2=0;
end

if flgStop==1
'VEQ fine isoooooo', keyboard
end

% Orders the elements length and extract the needed indeces (iso)
[~,iso]=sort([Ds(fineq); dLa_A; dLa_D2; dLa_D]);

% length of each element ordered relying on sort
du=[D(fineq); ones(size(dLa_A))*DLa; ones(size(dLa_D2))*DLd2; ones(size(dLa_D))*DLd;];
D1=du(iso); 

% xmol of each element ordered relying on sort
du=[X(fineq); ones(size(dLa_A))*ParMore.mireq.p.x_eq; ones(size(dLa_D2))*ParMore.mireq.n2.x_eq; ones(size(dLa_D))*ParMore.mireq.n.x_eq];
X1=du(iso); 

% NA of each element ordered relying on sort
du=[NA(fineq); ones(size(dLa_A))*ParMore.mireq.p.Dop_eq; zeros(size(dLa_D2)); zeros(size(dLa_D))];
NA1=du(iso);

% ND of each element ordered relying on sort
du=[ND(fineq); zeros(size(dLa_A)); ones(size(dLa_D2))*ParMore.mireq.n2.Dop_eq; ones(size(dLa_D))*ParMore.mireq.n.Dop_eq];
ND1=du(iso);


% riduco per struttura a LAYERS DD
Ds=cumsum(D_i); % length of each layer progressively added
veq=Mesh_i(:,1).*(NA_i(:,1)+ND_i(:,1)); % multiplies doping and MESH POINTS
fineq=find(veq>=0); % looks for region w/out eq medium

% eq p
%
veq=Mesh_i.*NA_i(:,1);  % multiplies p-doping and mesh points (negative in eq medium)
fieqA=find(veq<0);  % look for eq p
if length(fieqA)>0
    if fieqA(1)>1
        din_A=Ds(fieqA(1)); % distance from origin of eq p
    else
        din_A=Ds(fieqA(1)+1);
    end
    % multiplies Nadd (always positive) to mesh points, whose entries are
    % negative in case of eq p layers
    veq=Mesh_i.*NA_i(:,1);
    fiA=find(veq<0);    % find eq p LAYERS indeces
    DLA=sum(D_i(fiA));  % total length of eq p
    DLa=LA/ParMore.mireq.p.Nmesh;   % length of each element of eq p mesh (same as before)
    %  dLa_A=din_A+DLA;
    dLa_A=din_A;
else
    DLa=0;
    DLA=[];
    dLa_A=[];
end



% eq n - BOTTOM
% very similar to eq p 
veq=Mesh_i.*ND_i(:,1);
fieqD=find(veq<0);

difi=diff(fieqD);          % understand the position of "holes" in various repetitions
fdis=find(difi>1);      % understand the position of "holes" in various repetitions
if length(fieqD)>0
    if isfield(mireq,'n2')
    pip=0;
    for kk=1:length(fdis)
        fiS{kk}=fieqD(pip+1:fdis(kk));
        pip=fdis(kk);
        np(kk)=veq(fiS{kk}(1));
    end
    kk=kk+1;
    fiS{kk}=fieqD(pip+1:end); % cells containing the #repetitions of each "block"
    np(kk)=veq(fiS{kk}(1));     % # of repetition of each block identified in fiS
    [~,npMa]=max(abs(np));  % identifies the min number of pmesh points and its position in np
    fiD=fiS{npMa};       % indeces of the str layers that repeat the most among the n-doped ones
    
    din_D=Ds(fiD(1));
    else
    din_D=Ds(fieqD(1));
    veq=Mesh_i.*ND_i(:,1);
    fiD=find(veq<0);
    end
    DLD=sum(D_i(fiD));
    %    dLa_D=din_D+DLD;
    dLa_D=din_D;
else
    DLa=0;
    DLD=[];
    dLa_D=[];
end

% eq n - TOP
% very similar to eq p 
veq=Mesh_i.*ND_i(:,1);
fieqD=find(veq<0);

difi=diff(fieqD);          % understand the position of "holes" in various repetitions
fdis=find(difi>1);      % understand the position of "holes" in various repetitions
if isfield(mireq,'n2')
    if length(fieqD)>0
        pip=0;
        for kk=1:length(fdis)
            fiS{kk}=fieqD(pip+1:fdis(kk));
            pip=fdis(kk);
            np(kk)=veq(fiS{kk}(1));
        end
        kk=kk+1;
        fiS{kk}=fieqD(pip+1:end); % cells containing the #repetitions of each "block"
        np(kk)=veq(fiS{kk}(1));     % # of repetition of each block identified in fiS
        [~,npMa]=min(abs(np));  % identifies the min number of pmesh points and its position in np
        fiD=fiS{npMa};       % indeces of the str layers that repeat the most among the n-doped ones
        
        din_D2=Ds(fiD(1));
        DLD2=sum(D_i(fiD));
        %    dLa_D=din_D+DLD;
        dLa_D2=din_D2;
    else
        DLa2=0;
        DLD2=[];
        dLa_D2=[];
    end
else
    DLa2=0;
    DLD2=[];
    dLa_D2=[];
end

%'iso ;;;', keyboard

% Insert in Ds the distance at which eq p and eq n begin
[du,iso]=sort([Ds(fineq); dLa_A; dLa_D2; dLa_D]);
Ds1=du;

% Insert to the length of each layer (D_i), the length of eq p and eq n
du=[D_i(fineq); DLA; DLD2; DLD];
D1_i=du(iso);

OO=ones(size(X_i(1,:,:)));
% Check if n and p eq medium are described or not
if length(ParMore.mireq.n2.x_eq)>0 && length(ParMore.mireq.n.x_eq)>0
    dux=[X_i(fineq,:,:); ParMore.mireq.n2.x_eq*OO; ParMore.mireq.n.x_eq*OO];
    dura=[Ra_i(fineq,:); zeros(1,size(Ra_i,2)); zeros(1,size(Ra_i,2))];
    duNA=[NA_i(fineq,:,:); 0*OO; 0*OO];
    duND=[ND_i(fineq,:,:); ParMore.mireq.n2.Dop_eq*OO; ParMore.mireq.n.Dop_eq*OO];
elseif length(ParMore.mireq.p.x_eq)>0 && length(ParMore.mireq.n.x_eq)>0
    dux=[X_i(fineq,:,:); ParMore.mireq.p.x_eq*OO; ParMore.mireq.n.x_eq*OO];
    dura=[Ra_i(fineq,:); zeros(1,size(Ra_i,2)); zeros(1,size(Ra_i,2))];
    duNA=[NA_i(fineq,:,:); ParMore.mireq.p.Dop_eq*OO; 0*OO];
    duND=[ND_i(fineq,:,:); 0*OO; ParMore.mireq.n.Dop_eq*OO];
elseif length(ParMore.mireq.n2.x_eq)>0
    dux=[X_i(fineq,:,:); ParMore.mireq.p.x_eq*OO];
    dura=[Ra_i(fineq,:);  zeros(1,size(Ra_i,2))];
    duNA=[NA_i(fineq,:,:); ParMore.mireq.p.Dop_eq*OO];
    duND=[ND_i(fineq,:,:); ParMore.mireq.n.Dop_eq*OO];
elseif length(ParMore.mireq.n.x_eq)>0
    dux=[X_i(fineq,:,:); ParMore.mireq.n.x_eq*OO];
    dura=[Ra_i(fineq,:); zeros(1,size(Ra_i,2))];
    duNA=[NA_i(fineq,:,:); ParMore.mireq.p.Dop_eq*OO];
    duND=[ND_i(fineq,:,:); ParMore.mireq.n.Dop_eq*OO];
else
    dux=[X_i(fineq,:,:)];
    dura=[Ra_i(fineq,:)];
    duNA=[NA_i(fineq,:,:)];
    duND=[ND_i(fineq,:,:)];
end

% molar fraction (dux) sorted between fineq layers and eq ones 
X1_i=dux(iso,:,:);
% Radius (dura)
Rag=dura(iso,:);
% Acceptor Doping
NA1_i=duNA(iso,:,:);
% Donor Doping
ND1_i=duND(iso,:,:);
% Mesh points 
du=[Mesh_i(fineq); ParMore.mireq.p.Nmesh; ParMore.mireq.n2.Nmesh; ParMore.mireq.n.Nmesh];
Mesh1_i=du(iso);