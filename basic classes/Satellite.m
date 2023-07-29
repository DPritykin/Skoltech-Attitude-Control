classdef Satellite < handle

    properties(SetAccess = protected, GetAccess = public)
        env           % environment referenced
        J             % inertia tensor
        invJ          % inversed inertia tensor

        controlParams % control parameters structure

        mtm           % magnetometer
        gyro          % gyroscope

        mtq           % magnetorquer array

        residualDipole   % residual magnetization - a dipole that can be set externally
        residualDipolePos  % position of residual dipole wrt sat axes
    end

    methods
        function this = Satellite(inertiaTensor)
            this.J = inertiaTensor;
            this.invJ = inv(inertiaTensor);
        end

        function setResidualDipole(this, mRes, pos)
            arguments
                this
                mRes(3,1) {mustBeNumeric} = [0; 0; 0]
                pos(3,1) {mustBeNumeric} = [0; 0; 0]
            end
            this.residualDipole = mRes;
            this.residualDipolePos = pos;
        end
        
        function setControlParams(this, parameters)
            arguments
                this
                % [s] sampling time step
                parameters.tMeas {mustBeNumeric} = 0.5;
                % [s] control time step
                parameters.tCtrl {mustBeNumeric} = 2;
                % [-] required orbital orientation
                parameters.qReq {mustBeNumeric} = [1; 0; 0; 0];
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
        end

        function setEnvironment(this, env)
            if isa(env, 'Environment')
                this.env = env;
            else
                error('Satellite:InvalidEnvironment', 'Invalid Environment object!')
            end
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

        function m = calcControlMagneticMoment(this, q, omegaRel, B, mRes)

            ctrl = this.controlParams;
            dq = quatProduct(ctrl.qReqCnj, q);

            if (dq(1) == 0)
                multiplier = 1;
            else
                multiplier = sign(dq(1));
            end

            S = 4 * multiplier * dq(2:4);

            m = ((-ctrl.kW * crossProduct(B, omegaRel) - ctrl.kQ * crossProduct(B, S)) - mRes) * ctrl.tLoop / ctrl.tCtrl;

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

       function m = calcResidualDipoleMoment(this)
            m = this.residualDipole; % Consider changing for a more complicated and true-to-life model
        end

        function b = calcResidualMagnFieldAtPosition(this,pos)
            mRes= this.calcResidualDipoleMoment();
            pos3 = vecnorm(pos)^3;
            ePos = pos / vecnorm(pos);

            b = (this.env.mu0 / 4*pi) * (3 * (dot(mRes,ePos) * ePos - mRes)) / pos3; % dipole magnetic field formula
        end
    end
end
