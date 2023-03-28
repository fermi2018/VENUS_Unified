function [M_rid]=MatrixCompacter(M,PUN)

if isvector(M)==0
    A=M(PUN.purid,PUN.purid)+M(PUN.purid1,PUN.purid1);
    B=M(PUN.purid,PUN.pulast)+M(PUN.purid1,PUN.pulast);
    C=M(PUN.pulast,PUN.purid)+M(PUN.pulast,PUN.purid1);
    D=M(PUN.pulast,PUN.pulast);
    
    M_rid=[A./2 B./2;C D];
else
    M_rid=M(PUN.purid);
    M_rid(end+1:end+2)=M(PUN.pulast);
end

