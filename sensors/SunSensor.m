classdef SunSensor < AbstractSensor

    methods
        function [val,intensity] = getSensorReadings(this, trueValue)
            % trueValue - sun vector in the body frame 
            generatedNoise = normrnd(this.bias, this.noiseSigma, [3, 1]);

            random_orthogonal = rand(3, 1); % Generating a random orthogonal vector for rotation axis
            random_orthogonal = random_orthogonal/norm(random_orthogonal);

            rotation_axis = cross(trueValue, random_orthogonal);
            rotation_axis = rotation_axis / norm(rotation_axis);

            cos_angles = cos(generatedNoise);
            sin_angles = sin(generatedNoise);

            % Extract individual angles for x, y, and z components
            cos_angle_x = cos_angles(1);
            cos_angle_y = cos_angles(2);
            cos_angle_z = cos_angles(3);
            
            sin_angle_x = sin_angles(1);
            sin_angle_y = sin_angles(2);
            sin_angle_z = sin_angles(3);
            
            % Create the 3x3 skew-symmetric matrices for the rotation axes
            Sx = [0, -rotation_axis(3), rotation_axis(2);
                  rotation_axis(3), 0, -rotation_axis(1);
                  -rotation_axis(2), rotation_axis(1), 0];
            
            Sy = [0, -rotation_axis(3), rotation_axis(2);
                  rotation_axis(3), 0, -rotation_axis(1);
                  -rotation_axis(2), rotation_axis(1), 0];
            
            Sz = [0, -rotation_axis(3), rotation_axis(2);
                  rotation_axis(3), 0, -rotation_axis(1);
                  -rotation_axis(2), rotation_axis(1), 0];
            
            % Compute the rotation matrices for x, y, and z axes using Rodrigues' rotation formula
            Rx = eye(3) + sin_angle_x * Sx + (1 - cos_angle_x) * (Sx * Sx);
            Ry = eye(3) + sin_angle_y * Sy + (1 - cos_angle_y) * (Sy * Sy);
            Rz = eye(3) + sin_angle_z * Sz + (1 - cos_angle_z) * (Sz * Sz);
            
            % Compute the overall rotation matrix R_noise by combining the rotations about x, y, and z axes
            R_noise = Rx * Ry * Rz; 
            val = R_noise * trueValue;

            %INCIDENT ANGLE CALCULATION
            sunSensorNormal = this.position/1e-2;
            cos_incident_angle = dot(val, sunSensorNormal);
            incident_angle_deg = acosd(cos_incident_angle);

            if incident_angle_deg < this.fov
                intensity = 1;
            else
                intensity = 0;
            end 
        end
    end
end