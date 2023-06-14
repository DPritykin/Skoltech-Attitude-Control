function res = rhsRotationalDynamics(t, x, sat, orb, trqCtrl, trqDist)

    q = x(1:4) / norm(x(1:4));
    omega = x(5:7);

    rwControl = length(x) == 10;
    if rwControl
        h = x(8:10);
    else
        h = [0; 0; 0];
    end

    %% Gravity-gradient torque
    ez_b = quatRotate(q, [0; 0; 1]);
    trqGrav = 3 * orb.meanMotion^2 * crossProduct(ez_b, sat.J * ez_b);

    %% Magnetic torque
    if ~exist('trqCtrl', 'var')
        trqCtrl = [0; 0; 0];
    end    

    %% Disturbance torque
    if ~exist('trqDist', 'var')
        trqDist = [0; 0; 0];
    end

    %% Right-hand side equations
    Omega = omega - quatRotate(q, [0; orb.meanMotion; 0]);

    dq = 0.5 * quatProduct(q, [0; Omega]);

    dOmega = sat.invJ * (- crossProduct(omega, (sat.J) * omega + h) + ...
                           trqCtrl + trqGrav + trqDist);

    if rwControl
        dh = - trqCtrl;
        res = [dq; dOmega; dh];
    else
        res = [dq; dOmega];
    end
end