function [G] = mpGji(rj,ri,threshold)
G=mp('1'); %-- Excluding calculation of this function
igi=mp('0');
if igi==mp('1')
 G = (rj-ri)./log(rj./ri).*2./(rj+ri);
 indDiag = abs(rj-ri)<threshold; % diagonal (origin excluded)
 G(indDiag) = 1;
 indAxis_rj = abs(ri)<threshold; % rj axis (ri=0)
 G(indAxis_rj) = threshold;
 indAxis_ri = abs(rj)<threshold; % ri axis (rj=0)
 G(indAxis_ri) = threshold;
%-- Origin: no lim exists!
 indZer = abs(rj)<threshold & abs(ri)<threshold; 
 G(indZer) = threshold; % arbitrary choice!
 'passo Gji', keyboard
end 