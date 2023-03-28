 if Effetti>0
  effetti(Effetti)=1 
  if Effetti==8
   FAT_Diff_E=2;
   FAT_Diff_H=2;  
  elseif Effetti==9
   iStruttura=6;
  elseif Effetti==10
   iLUT=6;
  end


  
 mode.effetti=effetti;
  %'effetti', keyboard
  
   structureName=STR{iStruttura};
   NOMELUT=LUT{iLUT};
   mode.GLUT=[NOMELUT,'_Der.mat'];
   geom.GLUTm=[NOMELUT,'_more.mat'];
   mode.FAT_idiffusioneQW_E=FAT_Diff_E;
   mode.FAT_idiffusioneQW_H=FAT_Diff_H;

%    save Contributi effetti IPAR PMAT Svar kpar Effetti LUT STR -append
    eval(['save ',nomeSW,'Contributi_',nomeSav,' ',' effetti IPAR PMAT Svar kpar Effetti LUT STR -append'])

 end     
