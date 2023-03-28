load PL.txt -ascii

figure, plot(PL(:,1),PL(:,2)-linspace(0,1000,length(PL(:,1)))')