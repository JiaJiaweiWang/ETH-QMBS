function [beta]=BetaC(E0,E)
num=size(E,1);
beta=fzero(@EcEnn,0);
function [dE]=EcEnn(beta1)
    fu=0;
    fd=0;
    for l=1:num
        fu=fu+exp(-beta1*E(l))*E(l);
        fd=fd+exp(-beta1*E(l));
    end
    Ec=fu/fd;
    dE=Ec-E0;
end
end