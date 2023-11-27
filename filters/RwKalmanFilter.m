classdef RwKalmanFilter < AbstractKalmanFilter
% EKF for rw attitude control

    methods

        function this = RwKalmanFilter(options)
            arguments
                options.sat Satellite
                options.orb CircularOrbit
                options.env Environment
                options.sigmaQ0 {mustBeNumeric} = 1e-10
                options.sigmaOmega0 {mustBeNumeric} = 1e-10
            end

            this.sat = options.sat;
            this.orb = options.orb;

            this.initProcessCovariance(options.env.distTorqueSigma);
            this.initErrorCovariance(options.sigmaQ0, options.sigmaOmega0);
            this.initMeasurementsCovariance();
        end

        function estimatedX = estimate(this, t0, x0, bModelOrb, bSensor, sunModelOrb, sunSensor, rwCtrl)
            arguments
                this
                t0 {mustBeNumeric}
                x0 {mustBeNumeric}
                bModelOrb(3, 1) {mustBeNumeric}
                bSensor(3, 1) {mustBeNumeric}
                sunModelOrb(3, 1) {mustBeNumeric}
                sunSensor(3, 1) {mustBeNumeric}
                rwCtrl(3, 1) {mustBeNumeric} = [0; 0; 0];
            end

            [predictedX, predictedP] = this.prediction(t0, x0, rwCtrl);

            estimatedX = this.correction(predictedX, predictedP, bModelOrb, bSensor, sunModelOrb, sunSensor);
        end

    end

    methods (Access = protected)
        
        % EKF predictor
        function [predictedX, predictedP] = prediction(this, t0, x0, ctrlTorque)

            timeInterval = [t0, t0 + this.sat.controlParams.tLoop];
            [ ~, stateVec ] = ode45(@(t, x) rhsRotationalDynamics(t, x, this.sat, this.orb, ctrlTorque), ...
                                    timeInterval, x0, this.odeOptions);

            predictedX = stateVec(end, 1:10)';
            predictedX(1:4) = predictedX(1:4) / vecnorm(predictedX(1:4));

            Phi = this.calcEvolutionMatrix(x0, ctrlTorque);
            predictedP = Phi * this.P * Phi' + this.Q;
        end

        % EKF corrector
        function estimatedX = correction(this, predictedX, predictedP, bModelOrb, bSensor, sunModelOrb, sunSensor)
            bModelBody = quatRotate(predictedX(1:4), bModelOrb);
            sunModelBody = quatRotate(predictedX(1:4), sunModelOrb);
            bModelNorm = vecnorm(bModelBody);

            z = [bSensor / bModelNorm; sunSensor];
            Hx = [bModelBody / bModelNorm; sunModelBody];

            H = RwKalmanFilter.calcObservationMatrix(bModelBody / bModelNorm, sunModelBody);
            K = this.calcKalmanGain(predictedP, H, bModelNorm);

            correctedX = K * (z - Hx);
            qCor = vec2unitQuat(correctedX(1:3));
            
            estimatedX = zeros(6, 1);
            estimatedX(1:4) = quatProduct(predictedX(1:4), qCor);
            estimatedX(5:7) = predictedX(5:7) + correctedX(4:6);
            estimatedX(8:10) = predictedX(8:10);

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
            this.R = diag([this.sat.mtm.noiseSigma .* this.sat.mtm.noiseSigma; ...
                           this.sat.sunSensors.sensors(1).noiseSigma .* this.sat.sunSensors.sensors(1).noiseSigma]);
        end

%% state-dependent EKF matrices
        function Phi = calcEvolutionMatrix(this, state, rwCtrl)
            q = state(1:4);
            omega = state(5:7);

            F1 = [-skewSymm(omega), 0.5 * eye(3)];

            if isequal(rwCtrl, [0; 0; 0])
                e3 = quatRotate(q, [0; 0; 1]);

                Fgrav = 6 * this.orb.meanMotion^2 * (skewSymm(e3)* this.sat.J * skewSymm(e3) - ...
                                                     skewSymm(this.sat.J * e3) * skewSymm(e3));

                Fgyr = skewSymm(this.sat.J * omega) - skewSymm(omega) * this.sat.J;

                F2 = [this.sat.invJ * Fgrav, this.sat.invJ * Fgyr];
            else
                F2 = [-this.sat.controlParams.kQ  * eye(3), - this.sat.controlParams.kW  * eye(3)];
            end

            F = [F1; F2];

            Phi =  eye(6) + F * this.sat.controlParams.tLoop;
        end

        function K = calcKalmanGain(this, P, H, bModelNorm)
            measurementsCov = this.R;
            measurementsCov(1:3, 1:3) = this.R(1:3, 1:3) / bModelNorm^2;

            S = H * P * H' + measurementsCov;

            K =  P * H' / S;
        end

    end

    methods (Static)

        function H = calcObservationMatrix(bModel, sunModel)
            H = [2 * skewSymm(bModel), zeros(3);
                 2 * skewSymm(sunModel), zeros(3)];
        end

    end

end