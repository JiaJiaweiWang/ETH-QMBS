function [D]=sdistance(A,B,p)
C=(A/snorm(A,p))-(B/snorm(B,p));
D=snorm(C,p);
end