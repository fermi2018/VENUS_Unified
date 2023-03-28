function [N2,P2,N2i,P2i,T2,L2,indGan]=vet_mat_Ver(modePlot,ind)

n2D_vet=[];
p2D_vet=[];
TQW_vet=[];
n2Di_vet=[];
p2Di_vet=[];

lxgr=length(mesh.xgrid);
indGain=1:lxgr;
NQW=mesh.NMQW;
for kqw=1:NQW
 n2D_vet=[n2D_vet mode.n2D{kqw}];
 p2D_vet=[p2D_vet mode.p2D{kqw}];
 n2Di_vet=[n2Di_vet mode.n2Di{kqw}];
 p2Di_vet=[p2Di_vet mode.p2Di{kqw}]; 
 inQW = mesh.inMQW{kqw};
 TQW=mesh.T(inQW);  
 TQW_vet=[TQW_vet TQW];
 indGan{kqw}=indGain+(kqw-1)*lxgr;
end
OT=ones(size(TQW_vet));
OL=ones(size(mode.vlambda));

%keyboard
L2=mode.vlambda*OT;

N2=OL*n2D_vet;
P2=OL*p2D_vet;
T2=OL*TQW_vet;
N2i=OL*n2Di_vet;
P2i=OL*p2Di_vet;





