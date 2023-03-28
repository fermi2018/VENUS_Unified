puindic=diff(pimu([1 end]))+1;
 A=zeros(1,puindic,puindic);
 B=A;
 Bpp=A;
 Bmm=A;
 Amp=A;
 Apm=A;
 
 for imu=pimu
  jmu=imu-meun;
  mu=mbv(imu);
  for inu=pimu
   jnu=inu-meun;
   nu=mbv(inu);
%   if (nu+mu)/2-fix((nu+mu)/2)==0
    dmn0=mu-nu;
    dpn0=mu+nu;
%   [mu,nu], pausak    
%    if sdmn==0
%     keyboard
%    end

% A
    dmn=dmn0;
    fim=find(dmn==muv);
    if length(fim)==1
     A(:,jmu,jnu)=AB(:,fim);
     if ipolar==0
      sdmn=sign(mu-nu);
      C(:,jmu,jnu)=sdmn*CD(:,fim);
     end
    end

% B
    if length(fim)==1
     if fim==1 & mu==0
      B(:,jmu,jnu)=AB(:,fim);
     end 
    end    
    
% Amp    
    dmn=dmn0-2;
    fim=find(dmn==muv);
    if length(fim)==1
     Amp(:,jmu,jnu)=AB(:,fim);
     if ipolar==0
      D(:,jmu,jnu)=CD(:,fim);
     end
    end  %if   
   
% Apm    
    dmn=dmn0+2;
    fim=find(dmn==muv);
    if length(fim)==1
     Apm(:,jmu,jnu)=AB(:,fim);
     if ipolar==0
      D(:,jmu,jnu)=CD(:,fim);
     end
    end  %if 

% Bmm    
    dmn=dpn0-2;
    fim=find(dmn==muv);
    if length(fim)==1
     Bmm(:,jmu,jnu)=AB(:,fim);
     if ipolar==0
      D(:,jmu,jnu)=CD(:,fim);
     end
    end  %if   
   

  end
 end
 

%'fineMPNB', keyboard