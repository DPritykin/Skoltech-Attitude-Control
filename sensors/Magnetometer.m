classdef Magnetometer < AbstractSensor

    methods
       function val = getMagnetometerBias(this, InducedField)
            val = InducedField; % Induced Bias of the magnetometer
        end
        
        function val = getSensorReadings(this, trueValue)
            % trueValue - magnetic field in the satellite body frame
            generatedNoise = normrnd(this.bias, this.noiseSigma, [3, 1]);

            val = trueValue + this.dcm * generatedNoise;
        end
    end
end
