%      'entro coe_zin', keyboard

if dos>0

%      'entro coe_zin', keyboard
if dos>1 
%      'entro coe_zin', keyboard
end
Liz=Lizi;
Lizi_i=Lizi;
 if shailoc(istr0)==6 & igraef_new==2
  Lizi_i=Li(istr0);
 end
% if shailoc(istr0)==6 & igraef_new==2
%  Lizi_i=Li(istr0);
% end 
ioxi=0;
ncompari=0;
     Oos=Oo;
     Dos=dos;
     Amods=Amod;     
%     'inizio coe_zi: Dos=, istr=',
%     [Dos istr],
     istr
     iAmod=0;
     if istr>=9 & dovesono==10
          iAmod=1;

          'inizio coe_zi: Dos=, istr=',
           [Dos istr],
           
          keyboard
     end     
%     pausak
%     Dos
%     if istr==32
%      'istr coe_zi;  ', pausak
%      keyboard
%     end

     zip=zi(icfiez);
      le=length(find(imag(nitr)<-1));
      if le>0 & Dos/Liz>1
       Lizi_i=Dos;
       Liz=Lizi_i;
%'Liz grosso', keyboard
      end     
     if Dos>Liz & ifiez==1

%      ncompari=nitr(1)+10*aitr(1);
      if length(aitr)+1==length(nitr)
       ncompari=sum(nitr(1:end-1).*(aitr+1))+nitr(end);
      else
       ncompari=sum(nitr.*(aitr+1));      
      end

%      if ncompari>10
%       'ossido', keyboard
%       ioxi=1;
%      end
      finc=find(abs(ncompari-ncompar)<1e-6);
%     'valuto Oo  1', keyboard
      if length(finc)==0
%       'valuto Oo '
%       pausak
       ncompar(incdiv)=ncompari;
%       ' prima di emme', keyboard
       eval(emme)
%       ' dopo emme', keyboard       
       Treus{incdiv}=Oo;
       incdiv=incdiv+1;
      else
       Oo=Treus{finc};
      end
      if icfiez>47
%       'icfiez', keyboard
      end
      icfiez=icfiez+1;

      zi(icfiez)=zi(icfiez-1)+Lizi_i;

      Nz(icfiez)=nitr(1);
%      ' uno '
%      pausak
      if iAmod==1
       'Oo Amod', pausak
      end
      Amod=Oo*Amod;
      Acoz(:,icfiez)=Amod;
%      ' uno ', keyboard
%%%%%%%%%%%      Lis=2*Lizi_i;
      Lis=Lizi_i;
      ifirst=0;
      izemp=0;
%      if Dos>10
%       ' DOs'
%       Dos=100;
%%       keyboard
%      end
%      tic
%  Lis
%  Dos
% 'Lis ', pausak
%      while Lis<Dos

      decis=Lis-Dos;
      if decis>=0
%decis,      keyboard
      end

      le=length(find(imag(nitr)<-10));
%      while abs(Lis-Dos)>1e-6
      while decis+Liz<0
%         Lis-Dos
%       ' dentro while Lis', pausak
%      while Lis<Dos-Liz
%%       if Lis<10*Dos & ifirst==1
%%        izemp=1;
%%        pu1=1:length(KK);
%%        pu2=pu1+length(KK);
%%        Od0=diag(Oo,0);
%%        Odu=diag(Oo,length(KK));
%%        Odd=diag(Oo,-length(KK));
%%        O11=Od0(pu1);
%%        O12=Odu;
%%        O22=Od0(pu2);
%%        O21=Odd;;
%%        ifirst=0;
%%        'qui', keyboard
%%       end
if icfiez>48
 %'icfiez dentro', keyboard
end
       icfiez=icfiez+1;
%       'altro'
%       'sue ', keyboard
       zi(icfiez)=zi(icfiez-1)+Lizi_i;
       Nz(icfiez)=nitr(1);
       if le==0
%       pausak
%%       if izemp==0
      if iAmod==1
       'Oo Amod DENTRO', pausak
      end

         Amod=Oo*Amod;
%         istr
%         decis
%         pausak
       else
      if iAmod==1
       'da vedere Oo Amod', pausak
      end       
       
         Liz=Lis;
         eval(emme)
         Amod=Oo*Amods;
%         'passo nuovo', pausak
         Liz=Lizi_i;

       end  
%%       else
%%        Amod(pu1)=O11.*Amod(pu1)+O12.*Amod(pu2);
%%        Amod(pu2)=O21.*Amod(pu1)+O22.*Amod(pu2);
%        ' proviamo ', keyboard
%        'vel'
%%       end
       Acoz(:,icfiez)=Amod;
       Lis=Lis+Lizi_i;
       decis= Lis-Dos;
       if abs(Lis-Dos)<=1e-6 | decis>0
        decis=0;
       end
%       ' decis '
%       decis
%       pausak
      end
%      toc
%      keyboard
     end
      icfiez=icfiez+1;
%     'altro 1',     'tre', keyboard

      if iAmod==1
       ' FINE Oo Amod', pausak
      end
      Nz(icfiez)=nitr(1);
      zi(icfiez)=zip+Dos;
%     pausak
      Amod=Oos*Amods;
      Acoz(:,icfiez)=Amod;

%      zi(end-2:end)
%     'fine coe_zi', pausak
%if dos<1e-5
% ' controlla coezi ', keyboard
%end
end
%if istr==48
% ' controlla coezi ', keyboard
%end
Liz=Lizi;

if dovesono==2
%istr
% ' controlla coezi ', pausak
end 