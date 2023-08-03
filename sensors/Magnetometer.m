classdef Magnetometer < AbstractSensor

    methods
        function [val,intM] = getSensorReadings(this, trueValue)
            % trueValue - magnetic field in the satellite body frame
            generatedNoise = normrnd(this.bias, this.noiseSigma, [3, 1]);

            val = trueValue + this.dcm * generatedNoise; 
            intM = 0; 
        end
    end
end