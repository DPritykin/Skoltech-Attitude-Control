classdef Satellite < handle

    properties(SetAccess = protected, GetAccess = public)
        J             % inertia tensor
        invJ          % inversed inertia tensor

        controlParams % control parameters structure

        mtm           % magnetometer
        gyro          % gyroscope

        mtq           % magnetorquer array
        rw            % reaction wheels array
    end

    methods
        function this = Satellite(inertiaTensor)
            this.J = inertiaTensor;
            this.invJ = inv(inertiaTensor);
        end

        function setControlParams(this, parameters)
            arguments
                this
                % [s] sampling time step
                parameters.tMeas {mustBeNumeric} = 0;
                % [s] control time step
                parameters.tCtrl {mustBeNumeric} = 0;
                % [-] required orbital orientation
                parameters.qReq(4, 1) {mustBeNumeric} = [1; 0; 0; 0];
                % [-] required angular velocity (absolute)
                parameters.omegaReq(3, 1) {mustBeNumeric} = [0; 0; 0];                
                % [N * m * s / T^2] stabilizing parameter for PID-regulator
                parameters.kW {mustBeNumeric} = 0;
                % [N * m / T^2] orientating parameter for PID-regulator
                parameters.kQ {mustBeNumeric} = 0;
            end

            this.controlParams.tCtrl = parameters.tCtrl;
            this.controlParams.tMeas = parameters.tMeas;
            this.controlParams.kW = parameters.kW;
            this.controlParams.kQ = parameters.kQ;

            % [s] control loop duration
            this.controlParams.tLoop = parameters.tCtrl + parameters.tMeas;

            % conjugate quaternion
            this.controlParams.qReqCnj = quatConjugate(parameters.qReq);
            this.controlParams.omegaReq = parameters.omegaReq;
        end

        function setMagnetometer(this, mtm)
            if isa(mtm, 'Magnetometer')
                this.mtm = mtm;
            else
                error('Satellite:InvalidMtm', 'Invalid Magnetometer object!');
            end
        end

        function setGyroscope(this, gyro)
            if isa(gyro, 'Gyroscope')
                this.gyro = gyro;
            else
                error('Satellite:InvalidGyro', 'Invalid Gyroscope object!');
            end
        end

        function setMtqArray(this, mtqArray)
            if isa(mtqArray, 'MtqArray')
                this.mtq = mtqArray;
            else
                error('Satellite:InvalidMtq', 'Invalid Magnetorquer Array object!');
            end
        end

        function setRwArray(this, rwArray)
            if isa(rwArray, 'RwArray')
                this.rw = rwArray;
            else
                error('Satellite:InvalidRw', 'Invalid Reaction Wheels Array object!');
            end
        end

        function m = calcControlMagneticMoment(this, q, omegaRel, B)

            ctrl = this.controlParams;
            dq = quatProduct(ctrl.qReqCnj, q);

            if (dq(1) == 0)
                multiplier = 1;
            else
                multiplier = sign(dq(1));
            end

            S = 4 * multiplier * dq(2:4);

            m = (-ctrl.kW * crossProduct(B, omegaRel) - ctrl.kQ * crossProduct(B, S)) * ctrl.tLoop / ctrl.tCtrl;

            if any(abs(m) > this.mtq.maxMagneticMoment)
                maxRatio = max(abs(m) ./ this.mtq.maxMagneticMoment);
                m = m / maxRatio;
            end
        end

        function m = calcBdotMagneticMoment(this, magnField, angVelocity)
            m = cross(angVelocity, magnField);

            maxRatio = max(abs(m) ./ this.mtq.maxMagneticMoment);
            m = m / maxRatio;
        end

        function trq = calcRwControl(this, q, omega, externalTorqueToCompensate)
            quatErr = quatProduct(this.controlParams.qReqCnj, q);
            omegaErr = omega - this.controlParams.omegaReq;

            trqToActuate = -externalTorqueToCompensate ...
                           -this.controlParams.kW * this.J * omegaErr ...
                           -this.controlParams.kQ * this.J * quatErr(2:4);

            trq = this.rw.actuateCommand(trqToActuate, this.controlParams.tLoop);
        end
    end
end