classdef AbstractSensor < handle

    properties(SetAccess = protected, GetAccess = public)
        bias = [0; 0; 0]      % normal noise mean
        noiseSigma            % normal noise sigma

        position = [0; 0; 0]  % position (in the body-frame) 
                              % at which measurements are taken
        dcm = eye(3)          % dcm from sensor's axes to the
                              % host satellite body frame

        fov = 120;            % field of view (FOV) in degrees, specific to sun sensor
    end

    methods (Abstract)
        [sensorOutput,intensity] = getSensorReadings(this, trueValue)
    end

    methods
        function this = AbstractSensor(parameters)
            arguments
                parameters.bias {mustBeNumeric} = [0; 0; 0];
                parameters.sigma {mustBeNumeric} = 1e-12;
                parameters.position {mustBeNumeric} = [0; 0; 0];
                parameters.dcm {mustBeNumeric} = eye(3);
                parameters.fov {mustBeNumeric} = 120;
            end

            this.bias = parameters.bias;

            if length(parameters.sigma) == 3
                this.noiseSigma = parameters.sigma;
            else
                this.noiseSigma = repmat(parameters.sigma(1), 3, 1);
            end
            this.position = parameters.position;
            this.dcm = parameters.dcm;

            if isfield(parameters, 'fov')
                this.fov = parameters.fov;
            end

        end

        function setBias(this, bias)
            arguments
                this
                bias {mustBeNumeric} = [0; 0; 0];
            end

            this.bias = bias;
        end
    end
end
