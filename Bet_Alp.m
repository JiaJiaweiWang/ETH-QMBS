function bet_alp=Bet_Alp(E0,N0,E,N)
num=size(E,1);
x0=[0,0];
bet_alp=fsolve(@disEN,x0);
    function dEdN=disEN(bet_alp)
        Efu=0;
        Nfu=0;
        fd=0;
        for i=1:num
            expon=bet_alp(1)*E(i)-bet_alp(2)*N(i);
            Efu=Efu+exp(-expon)*E(i);
            Nfu=Nfu+exp(-expon)*N(i);
            fd=fd+exp(-expon);
        end
        Eg=Efu/fd;
        Ng=Nfu/fd;
        dEdN(1)=real(Eg-E0);
        dEdN(2)=real(Ng-N0);
    end
end