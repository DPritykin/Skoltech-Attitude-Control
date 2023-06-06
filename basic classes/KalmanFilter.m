classdef KalmanFilter < handle

    properties(SetAccess = protected, GetAccess = public)
        odeOptions = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

        sat % satellite object
        orb % orbit object

        P % error covariance matrix
        R % sensor noise covariance matrix
        Q % process noise covariance matrix
    end

    methods
        function this = KalmanFilter(options)
            arguments
                options.sat
                options.orb
                options.env
                options.sigmaQ0 {mustBeNumeric} = []
                options.sigmaOmega0 {mustBeNumeric} = []
            end

            this.sat = options.sat;
            this.orb = options.orb;
            this.initProcessCovariance(options.env.distTorqueSigma);
            this.initErrorCovariance(options.sigmaQ0, options.sigmaOmega0);
            this.initMeasurementsCovariance();
        end

        function estimatedX = estimate(this, t0, x0, mCtrl, bModel0, bmodelT, bSensor,omegaSensor)

            [predictedX, predictedP] = this.prediction(t0, x0, bModel0, mCtrl);

            estimatedX = this.correction(predictedX, predictedP, bmodelT, bSensor,omegaSensor);
        end

        function [predictedX, predictedP] = prediction(this, t0, x0, bModel0, mCtrl)
            bModel = quatRotate(x0(1:4), bModel0);

            % magnetorquers on
            timeInterval = [t0, t0 + this.sat.controlParams.tCtrl];
            [ ~, stateVec ] = ode45(@(t, x) rhsRotationalDynamics(t, x, this.sat, this.orb, bModel, mCtrl), ...
                                    timeInterval, x0, this.odeOptions);

            % magnetorquers off
            x0 = stateVec(end, 1:7)';
            timeInterval = [t0 + this.sat.controlParams.tCtrl, t0 + this.sat.controlParams.tLoop];
            [ ~, stateVec ] = ode45(@(t, x) rhsRotationalDynamics(t, x, this.sat, this.orb, bModel), ...
                                    timeInterval, x0, this.odeOptions);

            predictedX = stateVec(end, 1:7)';
            predictedX(1:4) = predictedX(1:4) / vecnorm(predictedX(1:4));

            Phi = this.calcEvolutionMatrix(x0, bModel, mCtrl);
            predictedP = Phi * this.P * Phi' + this.Q;
        end

        function estimatedX = correction(this, predictedX, predictedP, bModelT, bSensor,omegaSensor)
            bModelBody = quatRotate(predictedX(1:4), bModelT);
            bModelNorm = vecnorm(bModelBody);
            omegaNorm = vecnorm(omegaSensor);

            b_meas = bSensor / bModelNorm;
            w_meas = omegaSensor / omegaNorm;
            z = [b_meas ; w_meas];
             
            Hx_B = [bModelBody / bModelNorm];
            Hx_w = [omegaSensor / omegaNorm];
            Hx = [Hx_B ; Hx_w];

            H = this.calcObservationMatrix (Hx_B,Hx_w);
            K = this.calcKalmanGain(predictedP, H, bModelNorm,omegaNorm);

            correctedX = K * (z - Hx);
            qCor = vec2unitQuat(correctedX(1:3));
            
            estimatedX = zeros(6, 1);
            estimatedX(1:4) = quatProduct(predictedX(1:4), qCor);
            estimatedX(5:7) = predictedX(5:7) + correctedX(4:6);

            this.P = (eye(6) - K * H) * predictedP;
        end
    end

    methods (Access = private)

        function initErrorCovariance(this, sigmaQ0, sigmaOmega0)
            this.P = blkdiag(eye(3) * sigmaQ0^2, eye(3) * sigmaOmega0^2);
        end

        function initProcessCovariance(this, distTorqueSigma)
            G = [zeros(3); this.sat.invJ];
            D = eye(3) * distTorqueSigma^2;

            this.Q = G * D * G' * this.sat.controlParams.tCtrl;
        end

        function initMeasurementsCovariance(this)
            this.R = [diag(this.sat.mtm.noiseSigma .* this.sat.mtm.noiseSigma), zeros(3);
                      zeros(3) diag(this.sat.gyro.noiseSigma .* this.sat.gyro.noiseSigma)] ;
        end

        function Phi = calcEvolutionMatrix(this, state, bModel, mCtrl)
            q = state(1:4);
            omega = state(5:7);
            omegaRel = omega - quatRotate(q, [0; this.orb.meanMotion; 0]);

            e3 = quatRotate(q, [0, 0, 1]);

            Fgrav = 6 * this.orb.meanMotion^2 * (skewSymm(e3)* this.sat.J * skewSymm(e3) - ...
                                                 skewSymm(this.sat.J * e3) * skewSymm(e3));

            Fgyr = skewSymm(this.sat.J * omega) - skewSymm(omega) * this.sat.J;

            Fmagn = 2 * skewSymm(mCtrl) * skewSymm(bModel);

            F1 = [-skewSymm(omegaRel), 0.5 * eye(3)];
            F2 = [this.sat.invJ * (Fgrav + Fmagn), this.sat.invJ * Fgyr];

            F = [F1; F2];

            Phi =  eye(6) + F * this.sat.controlParams.tLoop;
        end

        function H = calcObservationMatrix(this, bModel,omegaSensor)
            H = [2 * skewSymm(bModel),zeros(3);zeros(3)  skewSymm(omegaSensor)];
        end

        function K = calcKalmanGain(this, P, H, bModelNorm,omegaNorm)
             S = H * P * H' + this.R / [(bModelNorm^2)*eye(3) , zeros(3) ; zeros(3) , (omegaNorm^2)*eye(3)];

            K =  P * H' / S;
        end

    end
end
