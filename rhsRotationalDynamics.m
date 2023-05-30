function res = rhsRotationalDynamics(t, x, sat, orb, envB, mCtrl, trqDist)

    q = x(1:4) / norm(x(1:4));
    omega = x(5:7);

    %% Gravity-gradient torque
    ez_b = quatRotate(q, [0; 0; 1]);
    trqGrav = 3 * orb.meanMotion^2 * crossProduct(ez_b, sat.J * ez_b);

    %% Magnetic torque
    if ~exist('mCtrl', 'var')
        mCtrl = [0; 0; 0];
    end    
    trqMagn = crossProduct(mCtrl, envB);

    %% Disturbance torque
    if ~exist('trqDist', 'var')
        trqDist = [0; 0; 0];
    end

    %% Right-hand side equations
    Omega = omega - quatRotate(q, [0; orb.meanMotion; 0]);

    dq = 0.5 * quatProduct(q, [0; Omega]);

    dOmega = sat.invJ * (- crossProduct(omega, (sat.J) * omega) + ...
                           trqMagn + trqGrav + trqDist);

    res = [dq; dOmega];
end

