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

        function estimatedX = estimate(this, t0, x0, mCtrl, bModel0, bmodelT, bSensor,SS_Vec_Sensor,SS_Vec_Model0,SS_Vec_ModelT)

            [predictedX, predictedP] = this.prediction(t0, x0, bModel0, mCtrl);

            estimatedX = this.correction(predictedX, predictedP, bmodelT, bSensor,SS_Vec_ModelT,SS_Vec_Sensor);
        end

        function [predictedX, predictedP] = prediction(this, t0, x0, bModel0, mCtrl)
            bModel = quatRotate(x0(1:4), bModel0);

            % magnetorquers on
            timeInterval = [t0, t0 + this.sat.controlParams.tCtrl];
            ctrlTorque = crossProduct(mCtrl, bModel);
            [ ~, stateVec ] = ode45(@(t, x) rhsRotationalDynamics(t, x, this.sat, this.orb, ctrlTorque), ...
                                    timeInterval, x0, this.odeOptions);

            % magnetorquers off
            x0 = stateVec(end, 1:7)';
            timeInterval = [t0 + this.sat.controlParams.tCtrl, t0 + this.sat.controlParams.tLoop];
            [ ~, stateVec ] = ode45(@(t, x) rhsRotationalDynamics(t, x, this.sat, this.orb), ...
                                    timeInterval, x0, this.odeOptions);

            predictedX = stateVec(end, 1:7)';
            predictedX(1:4) = predictedX(1:4) / vecnorm(predictedX(1:4));

            Phi = this.calcEvolutionMatrix(x0, bModel, mCtrl);
            predictedP = Phi * this.P * Phi' + this.Q;
        end

        function estimatedX = correction(this, predictedX, predictedP, bModelT, bSensor,SS_Vec_ModelT,SS_Vec_Sensor)

            bModelBody = quatRotate(predictedX(1:4), bModelT);
            SS_Vec_ModelT_body = quatRotate(predictedX(1:4),SS_Vec_ModelT); 

            bModelNorm = vecnorm(bModelBody);
            SunModelNorm = vecnorm(SS_Vec_ModelT_body);
            
            if ~isempty(this.sat.ss) 
                z = [bSensor / bModelNorm ; SS_Vec_Sensor/norm(SS_Vec_Sensor)]; 
                Hx = [bModelBody / bModelNorm ; SS_Vec_ModelT_body/SunModelNorm];        
            
            elseif isempty(this.sat.ss)
                z = bSensor / bModelNorm;
                Hx = bModelBody / bModelNorm;
            end 

            H = this.calcObservationMatrix(bModelBody/bModelNorm,SS_Vec_ModelT_body/SunModelNorm);
            K = this.calcKalmanGain(predictedP, H, bModelNorm,SunModelNorm);

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

            this.Q = G * D * G' * this.sat.controlParams.tLoop;
        end

        function initMeasurementsCovariance(this)
            if ~isempty(this.sat.ss)  
                this.R = [diag(this.sat.mtm.noiseSigma .* this.sat.mtm.noiseSigma), zeros(3);
                          zeros(3), diag(this.sat.ss.noiseSigma .* this.sat.ss.noiseSigma)];

            elseif isempty(this.sat.ss) 
                this.R = diag(this.sat.mtm.noiseSigma .* this.sat.mtm.noiseSigma);
            end 
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

        function H = calcObservationMatrix(this, bModel,Sun_Model)
            if ~isempty(this.sat.ss)
                H = [2 * skewSymm(bModel) zeros(3); 2*skewSymm(Sun_Model) zeros(3)];  

            elseif isempty(this.sat.ss)
                H = [2 * skewSymm(bModel), zeros(3)];
            end 
        end

        function K = calcKalmanGain(this, P, H, bModelNorm, SunModelNorm)
            S = H * P * H' + [this.R(1:3,:)/bModelNorm^2; this.R(4:6,:)];

            K =  P * H' / S;
        end

    end
end