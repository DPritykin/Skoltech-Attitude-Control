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
                options.sigmaQ0 {mustBeNumeric} = 1e-12
                options.sigmaOmega0 {mustBeNumeric} = 1e-12
                options.sigmaResDipole0 {mustBeNumeric} = 1e-12
                options.sigmaMtmBias0 {mustBeNumeric} = 1e-12
                options.sigmaGyroBias0 {mustBeNumeric} = 1e-12
            end

            this.sat = options.sat;
            this.orb = options.orb;
            this.initProcessCovariance(options.env.distTorqueSigma);
            this.initErrorCovariance(options.sigmaQ0, options.sigmaOmega0, options.sigmaResDipole0, options.sigmaMtmBias0, options.sigmaGyroBias0);
            this.initMeasurementsCovariance();
        end

        function estimatedX = estimate(this, t0, x0, mCtrl, bModel0, bmodelT, bSensor, omegaSensor)

            [predictedX, predictedP] = this.prediction(t0, x0, bModel0, mCtrl);

            estimatedX = this.correction(predictedX, predictedP, bmodelT, bSensor, omegaSensor);
        end

        function [predictedX, predictedP] = prediction(this, t0, x0, bModel0, mCtrl)
            bModel = quatRotate(x0(1:4), bModel0);

            mResEst = x0(8:10);
            predictedMtmBias = x0(11:13);
            predictedGyroBias = x0(14:16);
            
            % magnetorquers on
            timeInterval = [t0, t0 + this.sat.controlParams.tCtrl];
            [ ~, stateVec ] = ode45(@(t, x) rhsRotationalDynamics(t, x, this.sat, this.orb, bModel, mCtrl, mResEst), ...
                                    timeInterval, x0(1:7), this.odeOptions);

            % magnetorquers off
            x1 = stateVec(end, 1:7)';
            x1(1:4) = x1(1:4) / vecnorm(x1(1:4));
            timeInterval = [t0 + this.sat.controlParams.tCtrl, t0 + this.sat.controlParams.tLoop];
            [ ~, stateVec ] = ode45(@(t, x) rhsRotationalDynamics(t, x, this.sat, this.orb, bModel, [0; 0; 0], mResEst), ...
                                    timeInterval, x1(1:7), this.odeOptions);

            predictedX = [stateVec(end, 1:7)'; mResEst; predictedMtmBias; predictedGyroBias];
            predictedX(1:4) = predictedX(1:4) / vecnorm(predictedX(1:4));

            Phi = this.calcEvolutionMatrix(x0, bModel, mCtrl);
            predictedP = Phi * this.P * Phi' + this.Q;
        end

        function estimatedX = correction(this, predictedX, predictedP, bModelT, bSensor, omegaSensor)
            bModelBody = quatRotate(predictedX(1:4), bModelT);
            bModelNorm = vecnorm(bModelBody);

            mtmBias = predictedX(11:13);
            gyroBias = predictedX(14:16);

            if ~isempty(this.sat.gyro)
                z = [bSensor / bModelNorm; omegaSensor];
                Hx = [(bModelBody + mtmBias) / bModelNorm; predictedX(5:7) + gyroBias];   % magnetometer + gyro
            elseif isempty(this.sat.gyro)
                z = bSensor / bModelNorm;
                Hx = (bModelBody + mtmBias) / bModelNorm ;  % magnetometer only
            end

            H = this.calcObservationMatrix(bModelBody / bModelNorm);
            K = this.calcKalmanGain(predictedP, H, bModelNorm);

            correctedX = K * (z - Hx);
            qCor = vec2unitQuat(correctedX(1:3));

            estimatedX = zeros(16, 1);
            estimatedX(1:4) = quatProduct(predictedX(1:4), qCor);      % quaternion
            estimatedX(5:16) = predictedX(5:16) + correctedX(4:15);

            this.P = (eye(15) - K * H) * predictedP;
        end
    end

    methods (Access = private)

        function initErrorCovariance(this, sigmaQ0, sigmaOmega0,  sigmaResDipole0, sigmaMtmBias0, sigmaGyroBias0)
            this.P = blkdiag(eye(3) * sigmaQ0^2, ...
                             eye(3) * sigmaOmega0^2, ...
                             eye(3) * sigmaResDipole0^2, ...
                             eye(3) * sigmaMtmBias0^2, ...
                             eye(3) * sigmaGyroBias0^2);
        end

        function initProcessCovariance(this, distTorqueSigma)
            G = [zeros(3); this.sat.invJ; zeros(9, 3)];
            D = eye(3) * distTorqueSigma^2;

            this.Q = G * D * G' * this.sat.controlParams.tLoop;
        end

        function initMeasurementsCovariance(this)
            if ~isempty(this.sat.gyro)  
                this.R = [diag(this.sat.mtm.noiseSigma .* this.sat.mtm.noiseSigma), zeros(3);
                          zeros(3) diag(this.sat.gyro.noiseSigma .* this.sat.gyro.noiseSigma)];
            elseif isempty(this.sat.gyro) 
                this.R = diag(this.sat.mtm.noiseSigma .* this.sat.mtm.noiseSigma);
            end
        end

        function Phi = calcEvolutionMatrix(this, state, bModel, mCtrl)
            q = state(1:4);
            omega = state(5:7);

            e3 = quatRotate(q, [0, 0, 1]);

            Fgrav = 6 * this.orb.meanMotion^2 * (skewSymm(e3)* this.sat.J * skewSymm(e3) - ...
                                                 skewSymm(this.sat.J * e3) * skewSymm(e3));

            Fgyr = skewSymm(this.sat.J * omega) - skewSymm(omega) * this.sat.J;
            Fmagn = 2 * skewSymm(mCtrl) * skewSymm(bModel);
            Fmres = -this.sat.invJ * skewSymm(bModel);

            F = zeros(15);
            F(1:3, 1:6) = [-skewSymm(omega), 0.5 * eye(3)];
            F(4:6, 1:9) = [this.sat.invJ * (Fgrav + Fmagn), this.sat.invJ * Fgyr, Fmres];

            Phi =  eye(15) + F * this.sat.controlParams.tLoop;
        end

        function H = calcObservationMatrix(this, bModel)
            if ~isempty(this.sat.gyro)
                H = [2 * skewSymm(bModel), zeros(3, 6), eye(3), zeros(3);
                    zeros(3), eye(3), zeros(3, 6), eye(3)];                % magnetometer + gyroscope
            elseif isempty(this.sat.gyro)
                H = [2 * skewSymm(bModel), zeros(3, 6), eye(3)];           % only magnetometer
            end
        end

        function K = calcKalmanGain(this, P, H, bModelNorm)
            measCovariance = this.R;
            measCovariance(1:3, 1:3) = this.R(1:3, 1:3) / bModelNorm^2;

            if ~isempty(this.sat.gyro)
                S = H * P * H' + measCovariance;
            elseif isempty(this.sat.gyro)
                S = H * P * H' + measCovariance;
            end

            K =  P * H' / S;
        end

    end
end
