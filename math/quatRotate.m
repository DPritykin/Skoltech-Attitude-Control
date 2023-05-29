function r = quatRotate(q, vec)

    q = quatConjugate(q);

    qxvec = crossProduct(q(2:4), vec);

    r = dotProduct(q(2:4), vec) * q(2:4) + ...
        q(1)^2 * vec + 2 * q(1) * qxvec + ...
        crossProduct(q(2:4), qxvec);
end