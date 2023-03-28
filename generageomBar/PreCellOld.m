
x_dd=StrDD.x_dd;
Na_dd=StrDD.NA_dd/mode.CarrierNorm;
Nd_dd=StrDD.ND_dd/mode.CarrierNorm;
d_dd=StrDD.d_dd;
mesh_dd=StrDD.mesh_dd;
lab=StrDD.lab;
material=StrDD.material;
if isfield(StrDD,'nBTJ')
    mode.nBTJ=StrDD.nBTJ;
end
 
mesh_r=StrDD.mesh_r;
ra_dd=StrDD.ra_dd;
ra_long=StrDD.raggi;

 NUMcol=size(material,2);
if  isfield(mode,'quasi1D') 
 if mode.quasi1D==1 %|| mode.dimReduction==1
  NUMcol=1;
 end
end


PreCell