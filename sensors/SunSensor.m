classdef SunSensor < AbstractSensor

    methods
        function [val,intensity] = getSensorReadings(this, trueValue)
            % trueValue - sun vector in the body frame 

            randomVector = rand(3, 1);
            eAxis = cross(trueValue, randomVector);
            eAxis = eAxis / norm(eAxis);
            
            rotationAngle = normrnd(0, this.noiseSigma);
            
            % As all the three axes are assumed to have the same sigma 
            rotationQuaternion = [cos(rotationAngle(1)/2); eAxis*sin(rotationAngle(1)/2)];

            rotatedVector = quatRotate(rotationQuaternion, trueValue);
            val = rotatedVector + this.bias;

            sunSensorNormal = [0; 0; -1];
            cos_incident_angle = dot(rotatedVector, sunSensorNormal);

            % Calculating intensity (proportional to the incident angle cosine when in fov)
            if cos_incident_angle >= cos(this.fov)
                intensity = cos_incident_angle;
            else
                intensity = 0;
            end
        end
    end
end


          