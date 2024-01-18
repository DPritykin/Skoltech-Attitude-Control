classdef SunSensor < AbstractSensor

    properties (SetAccess = protected)
        fov                       % [deg] sun sensor field of view constraint 
        boresightDirection        % center of sun sensor's FOV  
        boresightDirectionSatRF   % center of sun sensor's FOV in Sat Ref Frame 
    end

    methods

        function this = SunSensor(parameters)
            arguments
                parameters.bias(3, 1) {mustBeNumeric}  = [0; 0; 0];
                parameters.sigma {mustBeNumeric} = deg2rad(0.2);
                parameters.dcm(3, 3) {mustBeNumeric} = eye(3);
                parameters.fovDeg {mustBeNumeric} = 120;
                parameters.boresightDirection(3, 1) {mustBeNumeric} = [0; 0; 1];
            end

            this.bias = parameters.bias;
            this.noiseSigma = parameters.sigma;
            this.fov = deg2rad(parameters.fovDeg);
            this.dcm = parameters.dcm;
            this.boresightDirection  = parameters.boresightDirection;
        end

        function [val, intensity] = getSensorReadings(this, trueSunDirection, sunEclipse)
            arguments
                this
                trueSunDirection(3, 1) {mustBeNumeric}
                sunEclipse {mustBeNumericOrLogical} = 0
            end

            sunIncidentAngleCos = dot(trueSunDirection, this.boresightDirectionSatRF);

            randomVector = rand(3, 1);
            eAxis = cross(trueSunDirection, randomVector);
            eAxis = eAxis / norm(eAxis);

            rotationAngle = normrnd(0, this.noiseSigma(1));

            % As all the three axes are assumed to have the same sigma
            rotationQuaternion = [cos(rotationAngle / 2); eAxis * sin(rotationAngle / 2)];

            rotatedVector = quatRotate(rotationQuaternion, trueSunDirection);

            if acos(sunIncidentAngleCos) <= this.fov/2 && sunEclipse == 0
                intensity = sunIncidentAngleCos;
                val = rotatedVector + this.bias;
                val = val / vecnorm(val);
            else
                val = NaN;
                intensity = 0;
            end
        end

        function setDcm(this, dcm)
            arguments
                this
                dcm(3, 3) {mustBeNumeric}
            end

            if abs(det(dcm) - 1) > 1e-10
                error('SunSensor:NonOrthogonalDcm', 'Matrix is not orthogonal!');
            end

            this.dcm = dcm;
            this.updateBoresightDirectionRF();
        end

        function set.fov(this, fov)
            this.fov = fov;
        end

        function set.boresightDirection(this, vec)
            arguments
                this
                vec(3, 1) {mustBeNumeric} = [0; 0; 1];
            end

            this.boresightDirection =  vec / vecnorm(vec);
            this.updateBoresightDirectionRF();
        end

    end

    methods (Access = protected)

        function updateBoresightDirectionRF(this)
            this.boresightDirectionSatRF = this.dcm * this.boresightDirection;
        end

    end
end