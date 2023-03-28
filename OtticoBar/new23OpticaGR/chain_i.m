%iold_1=0;
global iold_1


if length(iold_1)==0
 iold_1=1;
end

%'CHAI', keyboard

if iold_1==0
 chain_i_mod
else
 chain_i_vecchio
end