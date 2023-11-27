classdef MagneticKalmanFilter < AbstractKalmanFilter
% magnetometer-based EKF for magnetic attitude control

    methods

        function this = MagneticKalmanFilter(options)
            arguments
                options.sat {mustBeA(options.sat, 'Satellite')}
                options.orb {mustBeA(options.orb, 'CircularOrbit')}
                options.env {mustBeA(options.env, 'Environment')}
                options.sigmaQ0 {mustBeNumeric} = 1e-10
                options.sigmaOmega0 {mustBeNumeric} = 1e-10
            end

            this.sat = options.sat;
            this.orb = options.orb;

            this.initProcessCovariance(options.env.distTorqueSigma);
            this.initErrorCovariance(options.sigmaQ0, options.sigmaOmega0);
            this.initMeasurementsCovariance();
        end

        function estimatedX = estimate(this, t0, x0, bModel0, bmodelT, bSensor, mCtrl)
            arguments
                this
                t0 {mustBeNumeric}
                x0 {mustBeNumeric}
                bModel0(3, 1) {mustBeNumeric}
                bmodelT(3, 1) {mustBeNumeric}
                bSensor(3, 1) {mustBeNumeric}
                mCtrl(3, 1) {mustBeNumeric} = [0; 0; 0];                
            end

            [predictedX, predictedP] = this.prediction(t0, x0, bModel0, mCtrl);

            estimatedX = this.correction(predictedX, predictedP, bmodelT, bSensor);
        end

    end

    methods (Access = protected)
        
        % EKF predictor
        function [predictedX, predictedP] = prediction(this, t0, x0, bModel0, mCtrl)
            bModel = quatRotate(x0(1:4), bModel0);

            % magnetorquers on
            timeInterval = [t0, t0 + this.sat.controlParams.tCtrl];
            ctrlTorque = crossProduct(mCtrl, bModel);
            [ ~, stateVec ] = ode45(@(t, x) rhsRotationalDynamics(t, x, this.sat, this.orb, ctrlTorque), ...
                                    timeInterval, x0, this.odeOptions);

            % magnetorquers off
            x1 = stateVec(end, 1:7)';
            timeInterval = [t0 + this.sat.controlParams.tCtrl, t0 + this.sat.controlParams.tLoop];
            [ ~, stateVec ] = ode45(@(t, x) rhsRotationalDynamics(t, x, this.sat, this.orb), ...
                                    timeInterval, x1, this.odeOptions);

            predictedX = stateVec(end, 1:7)';
            predictedX(1:4) = predictedX(1:4) / vecnorm(predictedX(1:4));

            Phi = this.calcEvolutionMatrix(x0, bModel, mCtrl);
            predictedP = Phi * this.P * Phi' + this.Q;
        end

        % EKF corrector
        function estimatedX = correction(this, predictedX, predictedP, bModelT, bSensor)
            bModelBody = quatRotate(predictedX(1:4), bModelT);
            bModelNorm = vecnorm(bModelBody);

            z = bSensor / bModelNorm;
            Hx = bModelBody / bModelNorm;

            H = MagneticKalmanFilter.calcObservationMatrix(Hx);
            K = this.calcKalmanGain(predictedP, H, bModelNorm);

            correctedX = K * (z - Hx);
            qCor = vec2unitQuat(correctedX(1:3));
            
            estimatedX = zeros(6, 1);
            estimatedX(1:4) = quatProduct(predictedX(1:4), qCor);
            estimatedX(5:7) = predictedX(5:7) + correctedX(4:6);

            this.P = (eye(6) - K * H) * predictedP;
        end

%% initialization methods
        function initErrorCovariance(this, sigmaQ0, sigmaOmega0)
            this.P = blkdiag(eye(3) * sigmaQ0^2, eye(3) * sigmaOmega0^2);
        end

        function initProcessCovariance(this, distTorqueSigma)
            G = [zeros(3); this.sat.invJ];
            D = eye(3) * distTorqueSigma^2;

            this.Q = G * D * G' * this.sat.controlParams.tLoop;
        end

        function initMeasurementsCovariance(this)
            this.R = diag(this.sat.mtm.noiseSigma .* this.sat.mtm.noiseSigma);
        end

%% state-dependent EKF matrices
        function Phi = calcEvolutionMatrix(this, state, bModel, mCtrl)
            q = state(1:4);
            omega = state(5:7);

            e3 = quatRotate(q, [0; 0; 1]);

            Fgrav = 6 * this.orb.meanMotion^2 * (skewSymm(e3)* this.sat.J * skewSymm(e3) - ...
                                                 skewSymm(this.sat.J * e3) * skewSymm(e3));

            Fgyr = skewSymm(this.sat.J * omega) - skewSymm(omega) * this.sat.J;

            Fmagn = 2 * skewSymm(mCtrl) * skewSymm(bModel);

            F1 = [-skewSymm(omega), 0.5 * eye(3)];
            F2 = [this.sat.invJ * (Fgrav + Fmagn), this.sat.invJ * Fgyr];

            F = [F1; F2];

            Phi =  eye(6) + F * this.sat.controlParams.tLoop;
        end

        function K = calcKalmanGain(this, P, H, bModelNorm)
            S = H * P * H' + this.R / bModelNorm^2;

            K =  P * H' / S;
        end

    end

    methods (Static)

        function H = calcObservationMatrix(bModel)
            H = [2 * skewSymm(bModel), zeros(3)];
        end

    end

end
