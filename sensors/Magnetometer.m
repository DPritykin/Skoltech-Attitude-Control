classdef Magnetometer < AbstractSensor

    methods
        function val = getSensorReadings(this, trueValue, inducedField)
            if nargin < 3
                inducedField = [0; 0; 0]; % TODO be careful what ref-frame the inducedField is in
            end

            updatedBias = this.getBias(inducedField);

            % trueValue - magnetic field in the satellite body frame
            generatedNoise = normrnd(updatedBias, this.noiseSigma, [3, 1]);

            val = trueValue + updatedBias + this.dcm * generatedNoise;
        end


        function mtmBias = getBias(this, inducedField)
            if nargin < 3
                inducedField = [0; 0; 0];
            end

            mtmBias = this.bias + inducedField;
        end
    end
end
