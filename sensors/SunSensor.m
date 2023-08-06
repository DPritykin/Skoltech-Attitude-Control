classdef SunSensor < AbstractSensor
    properties
        fov                       % [deg] sun sensor field of view constraint 
        boresightDirection        % center of sun sensor's FOV  
        boresightDirectionSatRF   % center of sun sensor's FOV in Sat Ref Frame 
    end

    methods
        function this = SunSensor(parameters)
            arguments
                parameters.bias {mustBeNumeric} = [0; 0; 0];
                parameters.sigma {mustBeNumeric} = deg2rad(0.2);
                parameters.position {mustBeNumeric} = [2; 2; -3] * 1e-2;
                parameters.dcm {mustBeNumeric} = eye(3);
                parameters.fov {mustBeNumeric} = deg2rad(120);
                parameters.boresightDirection {mustBeNumeric} = [0; 0; -1];
            end

            this.bias                    = parameters.bias;

            if length(parameters.sigma) == 3
                this.noiseSigma = parameters.sigma;
            else
                this.noiseSigma = repmat(parameters.sigma(1), 3, 1);
            end

            this.position                = parameters.position;
            this.dcm                     = parameters.dcm;

            this.fov                     = parameters.fov;
            this.boresightDirection      = parameters.boresightDirection; 
            
            this.boresightDirectionSatRF = this.dcm * this.boresightDirection; 
        end

        function [val,intensity] = getSensorReadings(this, trueValue)
            % trueValue - sun vector in the body frame 

            % sunSensorNormal = this.boresightDirection;

                cos_incident_angle = dot(trueValue, this.boresightDirectionSatRF);

                randomVector = rand(3, 1);
                eAxis = cross(trueValue, randomVector);
                eAxis = eAxis / norm(eAxis);
            
                rotationAngle = normrnd(0, this.noiseSigma);
            
                % As all the three axes are assumed to have the same sigma 
                rotationQuaternion = [cos(rotationAngle(1)/2); eAxis*sin(rotationAngle(1)/2)];

                rotatedVector = quatRotate(rotationQuaternion, trueValue);
                
             if acos(cos_incident_angle) <= this.fov/2
                intensity = cos_incident_angle;   
                val = rotatedVector + this.bias;
             else
                val = [0;0;0];
                intensity = 0;  
             end
        end
    end
end


          