%load campo

%'qui ZOOM', keyboard
clear nma
if ~exist('ISTO')
 ISTO=  Ps.ifpstop;
 Ps.ifpstop=0;
 ifp=-4;
end 

isav_Az0=isav_Az;



if isfield(Ps,'izoom')==1
 izo=Ps.izoom;
else
 izo=0;
end


if izo>=1 & imod>=0

iORTA=0;
IPass=1;
isav_Az=0;
df=diff(fre_camp);

fre_camp_old=fre_camp;

%fre_camp=ze+([-1 0 1])*df(1)/50;
fre_camp=ze+([-1  1])*df(1)/20;

%fre_camp=ze+([-1 0 1])*df(1)/100;
%fre_camp=ze+([-1 0 1])*df(1)/30;
%fre_camp=fre_camp(3:5);
%'primo zoom',  keyboard
%kcav_at=kcav;
istopemme=0;
Ps.isav_Az=0;
imod=0;
nmasce=1
sub_zoom
%sub_old

%'dopo zoom', keyboard
clear nma
Ps.isav_Az=isav_Az0;
isav_Az=isav_Az0;

end

if izo>=2 & imod>=0
IPass=IPass+1;
isav_Az=0;
df=diff(fre_camp);

fre_camp_old=fre_camp;

IC=0;
imod=-1;
%while imod>=0 & IC<4

%IC=IC+1

%fre_camp=ze+([-1  1])*df(1)/500*IC;
fre_camp=ze+([-1  1])*df(1)/50;

'prima zoom', % keyboard
%kcav_at=kcav;
istopemme=0;
Ps.isav_Az=0;
imod=0;
nmasce=1
sub_zoom
%sub_old

%end %while

'dopo zoom 2', % keyboard
clear nma
Ps.isav_Az=isav_Az0;
isav_Az=isav_Az0;

end

if izo>=3 & imod>=0

while IPass<izo 
%while IPass<izo & imod>=0
IPass=IPass+1

isav_Az=0;
df=diff(fre_camp);

fre_camp_old=fre_camp;

IC=0;
imod=-1;
%while imod<0 & IC<4

%IC=IC+1

fre_camp=ze+([-1 1])*df(1)/50;
%fre_camp=ze+([-1 0 1])*df(1)/30;
%fre_camp=fre_camp(3:5);
'prima zoom', % keyboard
%kcav_at=kcav;
istopemme=0;
Ps.isav_Az=0;
imod=0;
nmasce=1
sub_zoom
%sub_old

%end %while IC

IPass
'dopo zoom i', %pausak
clear nma
Ps.isav_Az=isav_Az0;
isav_Az=isav_Az0;
end %while

end

if imod==-1
[du,imr]=min(abs(alvet(:,1)));
[du,imc]=min(abs(alvet(imr,:)));
 fris=fre_camp(imc);
end

'dopo zoom i', %keyboard
    sub_Ez_old

    
if ISTO==1
     figure, plot(zi,abs(Acoz(10,:))), 
    keyboard 
end     
imod=1;

clear iORTA