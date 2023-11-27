classdef Gyroscope < AbstractSensor

    methods

        function val = getSensorReadings(this, trueValue)
            % trueValue - angular velocity in the satellite body frame
            generatedNoise = normrnd(this.bias, this.noiseSigma, [3, 1]);

            val = trueValue + this.dcm * generatedNoise;
        end

    end
end