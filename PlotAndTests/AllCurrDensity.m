%clear all
close all
%load du1f
%load du2
%load prova
%load du1

 figure, plot(reshape(mode.Jn_y,mesh.nny,mesh.nnx))
 hold on
 plot(JN_Y,'.'), pausak

[Jn_xR,Jn_yR,Jp_xR,Jp_yR] = f_EvalCurrentDensityPost(geom,mesh,mode,iga,Tar);
% [Jn_xR,Jn_yR,Jp_xR,Jp_yR] = f_EvalCurrentDensity(geom,mesh,mode);
J_XNv=reshape(Jn_xR,mesh.nny,mesh.nnx);
J_YNv=reshape(Jn_yR,mesh.nny,mesh.nnx);
J_XPv=reshape(Jp_xR,mesh.nny,mesh.nnx);
J_YPv=reshape(Jp_yR,mesh.nny,mesh.nnx);
figure, plot(J_YNv)