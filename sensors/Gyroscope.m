classdef Gyroscope < AbstractSensor

     methods
         function val = getSensorReadings(this, trueValue)
            gyroBias = this.getBias();

            % trueValue - angular velocity in the satellite body frame
            generatedNoise = normrnd(gyroBias, this.noiseSigma, [3, 1]);

            val = trueValue + this.dcm * generatedNoise;
         end

         function gyroBias = getBias(this)
            gyroBias = this.bias;  %TODO: Add random walk or other modifications
        end
    end  
          
end
