classdef Satellite < handle

    properties(SetAccess = protected, GetAccess = public)
        J             % inertia tensor
        invJ          % inversed inertia tensor

        controlParams % control parameters structure

        mtm           % magnetometer
        gyro          % gyroscope

        mtq           % magnetorquer array

        maxMagneticMoment = 0.2 % [A * m^2] maximal dipole moment
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

            this.controlParams.tCtrl = parameters.tMeas;
            this.controlParams.tMeas = parameters.tCtrl;
            this.controlParams.kW = parameters.kW;
            this.controlParams.kQ = parameters.kQ;

            % [s] control loop duration
            this.controlParams.tLoop = parameters.tCtrl + parameters.tMeas;

            quat = parameters.qReq;
            % conjugate quaternion
            this.controlParams.qReqCnj = [quat(1); -quat(2:4)];
        end

        function setMagnetometer(this, mtm)
            if isa(mtm, 'Magnetometer')
                this.mtm = mtm;
            else
                error('Satellite:InvalidMtm', 'Invalid Magnetometer object!');
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

            if abs(max(m)) > this.maxMagneticMoment
                m = m * this.maxMagneticMoment / abs(max(m));
            end
        end        
    end
end