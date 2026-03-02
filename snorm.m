function [sn] = snorm(A,p)
Enn=eig(A);
len=size(Enn,1);
ss=0;
for i=1:len
    ss=ss+abs(Enn(i,1))^p;
end
sn=ss^(1/p);
end