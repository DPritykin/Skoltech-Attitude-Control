function q_conj = quatConjugate(q)
    q_conj = zeros(4, 1);
    q_conj(1) = q(1); 
    q_conj(2:4) = -q(2:4);
end