
function Sub_Fig_Divider
%Version 1: 8 August 2018
%Make sure that the figure number of the subplot you want to divide is
%figure(1) or change the value of i (Line 11 of the code), also close all other figure which you dont need.
%HOW TO USE THIS
%Keep the figure(1) open (figure(i) if you changed the number)
%Next run this file (Sub_Fig_Divider).
%Copyright (c) 2018, Akshay Kumar 
%All rights reserved.
clear set;
i=1;%%i= The figure number (CHANGE IF REQUIRED)
FIG=figure(i);
g=figure(2);
plot_new=copyobj(get(FIG,'Children'),g);
[R,~]=size(plot_new);
ii=R;
k=0;
while ii~=0
    
    TH1=plot_new(ii);
    TH2=plot_new(ii-1);
    if isempty(TH1.Tag)
        pause(0.00001)
        set(plot_new(ii),'Position', get(0, 'DefaultAxesPosition'))        
        if isempty(TH2.Tag)        
           PLO=copyobj(g.Children(ii),figure(i+2+k));
           ii=ii-1;
        else
            pause(0.00001);
            PLO=copyobj([g.Children(ii),g.Children(ii-1)],figure(i+2+k));
           ii=ii-2;           
        end
    end
    k=k+1;
end
close 2;
end
