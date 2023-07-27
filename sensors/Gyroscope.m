classdef Gyroscope < AbstractSensor

    methods
        function [val,intG] = getSensorReadings(this, trueValue)
            % trueValue - angular velocity in the satellite body frame
            generatedNoise = normrnd(this.bias, this.noiseSigma, [3, 1]);

            val = trueValue + this.dcm * generatedNoise;
            intG = 0; 
        end
    end
end