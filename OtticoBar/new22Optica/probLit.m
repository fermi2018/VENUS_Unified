d0=ParDD0.d;
figure, plot(cumsum(flipud(abs(d))),'o')
 hold on
 plot(cumsum(flipud(abs(L_i))),'r.')
 plot(cumsum(flipud(abs(d0))),'g+')

 Cfd=cumsum(flipud(abs(d)));
 CfL=cumsum(flipud(abs(L_i)));
  Cfd(end)-CfL(end)
 
 fd=flipud(d);
 fL=flipud(L_i);
 P=[fd(44:104) fL(44+[0:60])]