function [g,dgE,dgH,rsp,drspE,drspH]=f_InterpGain_lin(vportE,vportH,indQW,indLam)


global N2lin P2lin indGan 
global G0 dgE0 dgH0 rsp0 drspE0 drspH0


% 'in gain', keyboard

ind=indGan{indQW};
dE=vportE-N2lin(indLam,ind);
dH=vportH-P2lin(indLam,ind);

dgE=dgE0(indLam,ind);
dgH=dgH0(indLam,ind);
g=G0(indLam,ind)+dgE.*dE+dgH.*dH;

drspE=drspE0(indLam,ind);
drspH=drspH0(indLam,ind);
rsp=rsp0(indLam,ind)+drspE.*dE+drspH.*dH; % from Rtot, (1/s)






