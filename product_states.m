function productstates=product_states(perm,N)
    productstates=spin(perm(1));
    for i=2:N
        productstates=kron(productstates,spin(perm(i)));
    end
end