classdef Magnetometer < AbstractSensor

    methods
        function val = getSensorReadings(this, trueValue, InducedField)
            if nargin < 3
                InducedField = [0; 0; 0];
            end

            updatedBias = this.getBias(InducedField);

            % trueValue - magnetic field in the satellite body frame
            generatedNoise = normrnd(updatedBias, this.noiseSigma, [3, 1]);

            val = trueValue + updatedBias + this.dcm * generatedNoise;
        end


        function mtmBias = getBias(this, InducedField)
            if nargin < 3
                InducedField = [0; 0; 0];
            end

            mtmBias = this.bias + InducedField;
        end
    end
end
