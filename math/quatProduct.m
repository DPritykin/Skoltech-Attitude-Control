function q = quatProduct(q1, q2)

    q = zeros(4, 1);
    q(1) = q1(1) * q2(1) - q1(2) * q2(2) - q1(3) * q2(3) - q1(4) * q2(4);
    q(2:4) = q1(1) * q2(2:4) + q2(1) * q1(2:4) + ...
             crossProduct(q1(2:4), q2(2:4));

end