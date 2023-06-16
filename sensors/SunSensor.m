classdef SunSensor < AbstractSensor

    methods
        function val = getSensorReadings(this, trueValue)
            % trueValue - sun vector in the body frame 
            generatedNoise = normrnd(this.bias, this.noiseSigma, [3, 1]);

            val = trueValue + this.dcm * generatedNoise;
        end
    end
end