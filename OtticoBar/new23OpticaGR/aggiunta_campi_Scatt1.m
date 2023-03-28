   
   Refz=[];
   Lz=[];
   LzRatio=[];
   Layz=[];
   
   ik=0;
   irep=1;
   while ik<ficrit(end)
    ik=ik+irep
    ilay=fmlsdu(ik,1);
    pulay=ik;
    if ilay>1
     irep=fmlsdu(ik,2);
     pulay=pulay+[0:ilay-1];
    else
     irep=1;
    end 

%'qui ver'
     reps=fix(Litn_sav(pulay)/Lizi);
     Lhv=[];
     Thv=[];
     ThR=[];
     Rhv=[];
     for kit=1:length(pulay)
       Lhd=repmat(pulay(kit),reps(kit),1);
       Thdr=repmat(Lizi/Litn_sav(pulay(kit)),reps(kit),1);
       Thd=repmat(Lizi,reps(kit),1);
       Rhd=repmat(nitn(1,pulay(kit)),reps(kit),1);
       Thur=(Litn_sav(pulay(kit))-Lizi*reps(kit))/Litn_sav(pulay(kit));
       Thu=(Litn_sav(pulay(kit))-Lizi*reps(kit));
       if Thu>0
        Lhd=[Lhd; pulay(kit)];
        Thd=[Thd; Thu];
        Thdr=[Thdr; Thur];
        Rhd=[Rhd; nitn(1,pulay(kit))];
       end       
       Lhv=[Lhv; Lhd];       
       Thv=[Thv; Thd];
       ThR=[ThR; Thdr];
       Rhv=[Rhv; Rhd];
     end
    pulay
    
     Layz=[Layz; repmat(Lhv,irep,1)];
     Lz=[Lz; repmat(Thv,irep,1)];
     LzRatio=[LzRatio; repmat(ThR,irep,1)];
     Refz=[Refz; repmat(Rhv,irep,1)];

%'vedo', keyboard
   end