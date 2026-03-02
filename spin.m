function [spinstate]=spin(k)
    if k==-1
        spinstate=sparse([0;0;1]);
    elseif k==0
        spinstate=sparse([0;1;0]);
    elseif k==1
        spinstate=sparse([1;0;0]);
    end
end