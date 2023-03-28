function [v_rid]=VectorUnpacker(v,PUN)

v_rid(PUN.purid)=v(1:end-2);
v_rid(PUN.purid1)=v(1:end-2);
v_rid(PUN.pulast)=v(end-1:end);
