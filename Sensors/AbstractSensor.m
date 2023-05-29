classdef AbstractSensor < handle

    properties(SetAccess = protected, GetAccess = public)
        bias = [0; 0; 0]      % normal noise mean
        noiseSigma            % normal noise sigma

        position = [0; 0; 0]  %
        dcm = eye(3)          % dcm from sensor's axes to the
        % host satellite body frame
    end

    methods (Abstract)
        sensorOutput = getSensorReadings(this, trueValue)
    end

    methods
        function this = AbstractSensor(options)
            arguments
                options.bias {mustBeNumeric} = [0; 0; 0];
                options.sigma {mustBeNumeric} = 1e-12;
                options.position {mustBeNumeric} = [0; 0; 0];
                options.dcm {mustBeNumeric} = eye(3);
            end

            this.bias = options.bias;

            if length(options.sigma) == 3
                this.noiseSigma = options.sigma;
            else
                this.noiseSigma = repmat(options.sigma(1), 3, 1);
            end
            this.position = options.position;
            this.dcm = options.dcm;
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