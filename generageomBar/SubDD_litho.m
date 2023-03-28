function [x_out,d_out,puRag]=SubDD_litho(x_in,L_in,puTJ,pu_Below,ifp)

ifi=0;
if ifp==-10
ifi=1;
end

%L_i=abs(L_in);
L_i=(L_in);

if ifi==1
Cs=cumsum(L_i);
figure, plot(Cs,x_in,'.-')


hold on

D=sum(L_i(1:puTJ(1)-1));
%plot(cumsum(L_i(puTJ))+D,x_in(puTJ),'r','linewidth',2)
plot(Cs(puTJ),x_in(puTJ),'or','linewidth',2)
pausak
end


Detch=sum(L_i(puTJ));

x_lateral=-.1;
pui=[1:puTJ(1)-1];



Le=[Detch; L_i(pui)];
%'nex', keyboard

nex=[x_lateral*ones(size(x_in(1,:))); x_in(pui,:)];
Ltop=flipud(Le);
Les=cumsum(Ltop);
ntopE=flipud(nex);

LeInter=[[0;Les(1:end-1)+1e-2] Les];

neInter=[ntopE ntopE];
Lint_ext=reshape(LeInter',2*length(Les),1);
nint_ext=reshape(neInter.',2*length(Les),1);


puiTJ=[1:puTJ(end)];
LtopTJ=flipud(L_i(puiTJ));
ntopTJ=flipud(x_in(puiTJ));
Les=cumsum(LtopTJ);
ntopI=ntopTJ;



LeInter=[[0;Les(1:end-1)+1e-3] Les];
LeInterF=[[0;Les(1:end-1)] Les-1e-3];
%LeInter=[[0;Les(1:end-1)+1e-5] Les];
%LeInter=[[0;Les(1:end-1)] Les];
neInter=[ntopI ntopI];



Lint_int=reshape(LeInter',2*length(Les),1);
Lint_intF=reshape(LeInterF',2*length(Les),1);
nint_int=reshape(neInter.',2*length(Les),1);



Dtot=[0; cumsum(LtopTJ); cumsum(Ltop)];
[Ds,iso]=sort(Dtot);

dL=diff([-10; Ds]);
fiVeri=find(dL>0);


Dfinal=Ds(fiVeri);
Lfinal=diff(Dfinal);

n_intI=interp1(Lint_int,nint_int,Dfinal);
n_extI=interp1(Lint_ext,nint_ext,Dfinal);

if ifp==-10
'VEDI FIT', keyboard
end


%de=-1
%n_intIp=interp1(Lint_int,nint_int,1564+de)
%n_extIp=interp1(Lint_ext,nint_ext,1564+de)

%  figure, plot(Dfinal,real(n_intI))
%  hold on, plot(Dfinal,real(n_extI),'r')

LeInter=[[Dfinal(1:end-1)] Dfinal(2:end)];
neInter=[n_intI(2:end) n_intI(2:end) ];
Lplot=reshape(LeInter',2*length(LeInter),1);
np_int=reshape(neInter.',2*length(LeInter),1);

neInter=[n_extI(2:end) n_extI(2:end) ];
np_ext=reshape(neInter.',2*length(LeInter),1);

%'qui; SubDD_litho', keyboard

x_outd=[n_intI(2:end) n_extI(2:end)];
x_out=[flipud(x_outd); repmat(x_in(pu_Below,:),1,2)];


d_out=[flipud(Lfinal); L_i(pu_Below)];

fi=find(abs(d_out)>1e-4);
d_out=d_out(fi);
x_out=x_out(fi,:);
%'qui;', keyboard

if ifp==-10


figure, plot(Lint_ext,nint_ext,'b','linewidth',2)
hold on, plot(Lint_int,nint_int,'r','linewidth',2)
title(' RED interno, BLU esterno')

%figure, 
plot(Lplot,real(np_ext),'bo')
hold on, plot(Lplot,np_int,'ro')

%plot(Lplot,real(np_ext),'g')
%hold on, plot(Lplot,np_int,'m')
'prima di uscire'
pausak
end

