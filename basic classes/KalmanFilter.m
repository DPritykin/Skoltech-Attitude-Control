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
                options.sigmaResDipole0 {mustBeNumeric} = []
                options.sigmaBias0_mtm {mustBeNumeric} = []
                options.sigmaBias0_gyro {mustBeNumeric} = []
            end

            this.sat = options.sat;
            this.orb = options.orb;
            this.initProcessCovariance(options.env.distTorqueSigma);
            this.initErrorCovariance(options.sigmaQ0, options.sigmaOmega0,options.sigmaResDipole0,options.sigmaBias0_mtm,options.sigmaBias0_gyro);
            this.initMeasurementsCovariance();
        end

        function estimatedX = estimate(this, t0, x0, mCtrl, bModel0, bmodelT, bSensor,omegaSensor)

            [predictedX, predictedP] = this.prediction(t0, x0, bModel0, mCtrl);

            estimatedX = this.correction(predictedX, predictedP, bmodelT, bSensor,omegaSensor);
        end

        function [predictedX, predictedP] = prediction(this, t0, x0, bModel0, mCtrl)
            bModel = quatRotate(x0(1:4), bModel0);
            bInduced = this.sat.calcResidualMagnFieldAtPosition(this.sat.residualDipolePos);

            m_resEst = x0(8:10);
            PredictedBias_mtm = x0(11:13);
            PredictedBias_gyro = x0(14:16);
            
            % magnetorquers on
            timeInterval = [t0, t0 + this.sat.controlParams.tCtrl];
            [ ~, stateVec ] = ode45(@(t, x) rhsRotationalDynamics(t, x, this.sat, this.orb, bModel, mCtrl, m_resEst), ...
                                    timeInterval, x0(1:7), this.odeOptions);

            % magnetorquers off
            x0 = stateVec(end, 1:7)';
            timeInterval = [t0 + this.sat.controlParams.tCtrl, t0 + this.sat.controlParams.tLoop];
            [ ~, stateVec ] = ode45(@(t, x) rhsRotationalDynamics(t, x, this.sat, this.orb, bModel), ...
                                    timeInterval, x0(1:7), this.odeOptions);

            predictedX = [stateVec(end, 1:7)' ; m_resEst ;PredictedBias_mtm; PredictedBias_gyro];
            predictedX(1:4) = predictedX(1:4) / vecnorm(predictedX(1:4));

            Phi = this.calcEvolutionMatrix(x0, bModel, mCtrl, m_restEst);
            predictedP = Phi * this.P * Phi' + this.Q;
        end

        function estimatedX = correction(this, predictedX, predictedP, bModelT, bSensor,omegaSensor)
            bModelBody = quatRotate(predictedX(1:4), bModelT);
            bModelNorm = vecnorm(bModelBody);

            mtm_bias = predictedX(11:13);
            gyro_bias = predictedX(14:16);
           
            if ~isempty(this.sat.gyro) 
                z = [bSensor / bModelNorm ; omegaSensor]; 
                Hx = [(bModelBody + mtm_bias) / bModelNorm  ; predictedX(5:7) + gyro_bias];   % magnetometer + gyro          
            
            elseif isempty(this.sat.gyro) 
                 z = bSensor / bModelNorm;
                 Hx = (bModelBody + mtm_bias) / bModelNorm ;  % magnetometer only 
                 
            end   

            H = this.calcObservationMatrix (Hx);
            K = this.calcKalmanGain(predictedP, H, bModelNorm);

            correctedX = K * (z - Hx);
            qCor = vec2unitQuat(correctedX(1:3));
            
            estimatedX = zeros(15, 1);
            estimatedX(1:4) = quatProduct(predictedX(1:4), qCor); % quaternion
            estimatedX(5:7) = predictedX(5:7) + correctedX(4:6);  % angular velocity
            estimatedX(8:10) = predictedX(8:10) + correctedX(7:9); %  residual dipole
            estimatedX(11:13) = predictedX(11:13) + correctedX(10:12); % magnetometer bias
            estimatedX(14:16) = predictedX(14:16) + correctedX(13:15); % magnetometer bias

            this.P = (eye(15) - K * H) * predictedP;
        end
    end

    methods (Access = private)

        function initErrorCovariance(this, sigmaQ0, sigmaOmega0,  sigmaResDipole0, sigmaBias0_mtm, sigmaBias0_gyro)
            this.P = blkdiag(eye(3) * sigmaQ0^2, eye(3) * sigmaOmega0^2, eye(3)*sigmaResDipole0, eye(3)*sigmaBias0_mtm, eye(3)*sigmaBias0_gyro);
        end

        function initProcessCovariance(this, distTorqueSigma)
            G = [zeros(3); this.sat.invJ; zeros(3); zeros(3); zeros(3)];
            D = eye(3) * distTorqueSigma^2;

            this.Q = G * D * G' * this.sat.controlParams.tCtrl;
        end

        function initMeasurementsCovariance(this)
            if ~isempty(this.sat.gyro)  
                this.R = [diag(this.sat.mtm.noiseSigma .* this.sat.mtm.noiseSigma), zeros(3);
                      zeros(3) diag(this.sat.gyro.noiseSigma .* this.sat.gyro.noiseSigma)];

            elseif isempty(this.sat.gyro) 
                this.R = diag(this.sat.mtm.noiseSigma .* this.sat.mtm.noiseSigma);
            end
        end

        function Phi = calcEvolutionMatrix(this, state, bModel, mCtrl, mRes)
            q = state(1:4);
            omega = state(5:7);
            omegaRel = omega - quatRotate(q, [0; this.orb.meanMotion; 0]);

            e3 = quatRotate(q, [0, 0, 1]);

            Fgrav = 6 * this.orb.meanMotion^2 * (skewSymm(e3)* this.sat.J * skewSymm(e3) - ...
                                                 skewSymm(this.sat.J * e3) * skewSymm(e3));

            Fgyr = skewSymm(this.sat.J * omega) - skewSymm(omega) * this.sat.J;

            Fmagn = 2 * skewSymm(mCtrl + mRes) * skewSymm(bModel);

            F1 = [-skewSymm(omegaRel), 0.5 * eye(3), zeros(3), zeros(3), zeros(3)];
            F2 = [this.sat.invJ * (Fgrav + Fmagn), this.sat.invJ * Fgyr, zeros(3), zeros(3), zeros(3)];
            F3 = [zeros(3), zeros(3), zeros(3), zeros(3), zeros(3)]; % residual magnetization model linearization
            F4 = [zeros(3), zeros(3), zeros(3), zeros(3), zeros(3)]; % * bias eqn linearization = 0?
            F5 = [zeros(3), zeros(3), zeros(3), zeros(3), zeros(3)];
            F = [F1; F2; F3; F4; F5];

            Phi =  eye(15) + F * this.sat.controlParams.tLoop;
        end

        function H = calcObservationMatrix(this, bModel)
            if ~isempty(this.sat.gyro)
                H = [2 * skewSymm(bModel) zeros(3,6) eye(3) zeros(3);zeros(3) eye(3) zeros(3,6) eye(3)]; % magnetometer+gyroscope

            elseif isempty(this.sat.gyro)
                 H = [2*skewSymm(bModel),zeros(3,6),eye(3)];  % only magnetometer             
            end
        end

        function K = calcKalmanGain(this, P, H, bModelNorm)
             if ~isempty(this.sat.gyro)
                S = H * P * H' + this.R / [(bModelNorm^2)*eye(3) , zeros(3) ; zeros(3) , eye(3)]; 
                 
            elseif isempty(this.sat.gyro) 
                S = H * P * H' + this.R / bModelNorm^2;
            end

            K =  P * H' / S;
        end

    end
end
