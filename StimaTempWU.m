if mode.Pmin_Pfit<100
    Fat_STIMA_Temp=1;
elseif mode.Pmin_Pfit==100
    Fat_STIMA_Temp=1.05;
end
return

DTr=reshape(DeltaT,mesh.nny,mesh.nnx);
%'entro Stima', keyboard
mode.Tmed(indv)=mean(mean(DTr(20:fix(mesh.nny*2/3),1:fix(mesh.nnx/3))));
%Tm=mode1.Tmed;

%Tmed=mean(mean(DeltaTold));

if indv>8 & mode.Tmed(end)>30
   Tmd=mode.Tmed(2:end);
   fim=find(Tmd>.2);
%'qui', keyboard 
   Tm=Tmd(fim); 
   Vm=mode.vv_dd(fim);
   dT=diff(Tm)./diff(Vm)./Tm(2:end);
   Vf=Vm(2:end);
   coT=polyfit(Vf,dT,1);
   Vmf=[v0_dd];
   dT_est=polyval(coT,Vmf);  
   fiV=find(dd==v0_dd);
   dV=diff(dd(fiV+[-1 0]));
   
 %Vmp=Vm(end-5:end);
 %Tmp=Tm(end-5:end);
 Fat_STIMA_Temp=1+dT_est*dV
 figure, plot(Vf,dT), pausak 

 FSmax=1.1;

 if Fat_STIMA_Temp>FSmax
  Fat_STIMA_Temp=FSmax;
 end
 
 mode.Fat_STIMA(indv)=Fat_STIMA_Temp;
 mode.T_STIMA(indv)=Fat_STIMA_Temp*Tmd(end);
 if Fat_STIMA_Temp<1
   Fat_STIMA_Temp=1
 end
 %keyboard

  if isfield(mode,'TempEst')
   Tadd=[mode.TempEst T_est];
   mode.TempEst=Tadd;
   Vadd=[mode.VTempEst v0_dd];
   mode.VTempEst=Vadd;   
  else
   mode.TempEst=T_est;
   mode.VTempEst=v0_dd;
  end
else
 Fat_STIMA_Temp=1;
end 

if isfield(mode,'STIMA_Temp')
 if mode.STIMA_Temp==0
  Fat_STIMA_Temp=1;
 end
end