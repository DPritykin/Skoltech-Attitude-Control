classdef AbstractActuator < handle

    properties(SetAccess = protected, GetAccess = public)
        name                % name of the actuator (e.g. mtqX or rwA)
        state = true        % whether or not the actuator is operable or faulted
        axis = [0; 0; 1]    % direction, in which the actuated output is produced (own axes)
        dcm                 % conversion to body frame dcm * axis yields direction in the sat body frame
        axisBf

        noiseSigma          % normal noise sigma
    end

    methods
        function this = AbstractActuator(parameters)
            arguments
                parameters.name {mustBeText} = '';
                parameters.state {mustBeNumericOrLogical} = true;
                parameters.axis {mustBeNumeric} = [0; 0; 1];
                parameters.dcm {mustBeNumeric} = eye(3);
                parameters.noiseSigma {mustBeNumeric} = 0;
            end

            this.name  = parameters.name;
            this.state  = parameters.state;
            this.axis   = parameters.axis;
            this.dcm    = parameters.dcm;
            this.axisBf = this.dcm * this.axis;

            this.noiseSigma = parameters.noiseSigma;
        end

        function setDcm(this, dcm)
            this.dcm = dcm;
            this.axisBf = this.dcm * this.axis;
        end

        function setAxis(this, axis)
            this.axis = axis;
            this.axisBf = this.dcm * this.axis;
        end

        function setName(this, name)
            this.name = name;
        end

        function actuatedAction = actuateRequiredAction(this, requiredAction)
            actuatedAction = this.state * this.axisBf * (requiredAction + normrnd(0, this.noiseSigma));
        end
    end
end