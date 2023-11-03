classdef ReactionWheel < AbstractActuator

    properties
        maxTorque       % [Nm], maximum torque the reaction wheel can produce
        maxAngMomentum  % [Nms], maximum angular momentum the reaction wheel is allowed to attain

        angularMomentum = 0
    end

    methods
        function this = ReactionWheel(parameters)
            arguments
                parameters.name {mustBeText} = '';
                parameters.state {mustBeNumericOrLogical} = true;
                parameters.axis {mustBeNumeric} = [0; 0; 1];
                parameters.dcm {mustBeNumeric} = eye(3);
                parameters.noiseSigma {mustBeNumeric} = 1e-10;
                parameters.maxTorque {mustBeNumeric} = 0.97e-3;
                parameters.maxAngMomentum {mustBeNumeric} = 1e-2;
            end

            this.name       = parameters.name;
            this.state      = parameters.state;
            this.axis       = parameters.axis;
            this.dcm        = parameters.dcm;
            this.axisBf     = this.dcm * this.axis;

            this.noiseSigma = parameters.noiseSigma;

            this.maxTorque       = parameters.maxTorque;
            this.maxAngMomentum  = parameters.maxAngMomentum;
        end

        function setMaxTorque(this, maxTorque)
            this.maxTorque = maxTorque;
        end

        function setMaxAngMomentum(this, maxAngMomentum)
            this.maxAngMomentum = maxAngMomentum;
        end        

        % simulates the control torque output along the output axis
        % given the required torque and the current angular momentum of RW
        function outputTorque = actuateControlTorque(this, requiredTorque, duration)
            if abs(requiredTorque) > this.maxTorque
                requiredTorque = sign(requiredTorque) * this.maxTorque;
            end

            if abs(this.angularMomentum + requiredTorque * duration) > this.maxAngMomentum
                requiredTorque = 0;
            end

            outputTorque = this.actuateRequiredAction(requiredTorque);
            this.angularMomentum = this.angularMomentum + sign(requiredTorque) * vecnorm(outputTorque) * duration;
            
        end
    end
end