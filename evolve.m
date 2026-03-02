function [O_exp]=evolve(E,V,psi0,O,T)
    % T=T(:);
    % M = length(T);
    % O_exp = zeros(M, 1);
    % O_eig = V' * O * V;
    % c = V' * psi0;
    % c = c / sqrt(c'*c);
    % for k=1:M
    %     ct=c .* exp(-1j * E * T(k));
    % 
    %     O_exp(k)=ct'*O_eig*ct;
    % end
    % O_exp = real(O_exp);
    T = T(:);

    c = V' * psi0;
    c = c / sqrt(c'*c);

    c_time = c .* exp(-1j * E * T');

    O_eig = V' * O * V;

    O_exp = sum(conj(c_time) .* (O_eig * c_time), 1).';
    O_exp = real(O_exp);
end