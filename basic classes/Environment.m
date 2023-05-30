classdef Environment < handle

    properties(SetAccess = protected, GetAccess = public)
        muG = 3.986e+14;         % [m^3 / s^2] standard gravitational parameter of the Earth
        mu0 = 1.257e-6;          % [N / A^2] vacuum permeability
        muE = 7.94e+22;          % [A * m^2] magnetic dipole moment of the Earth
        earthRadius = 6371e+3;   % [m] radius of the Earth

        distTorqueSigma = 0;     % [N * m] disturbance of the torque
        magnFieldSigma = 0;      % [T] noise to create actual magnetic field
    end

    methods
        function this = Environment(options)
            arguments
                options.dstSigma {mustBeNumeric} = [0; 0; 0];
                options.magnSigma {mustBeNumeric} = [0; 0; 0];
            end

            this.distTorqueSigma = options.dstSigma;
            this.magnFieldSigma = options.magnSigma;
        end

        function B = directDipoleOrbital(this, argLat, inclination, orbitRadius)
            % [T] magnitude of the magnetic field on the orbit
            B0 = this.muE * this.mu0 / (4 * pi * orbitRadius^3);   

            B = B0 * [cos(argLat) * sin(inclination); 
                      cos(inclination); 
                      -2 * sin(argLat) * sin(inclination)];
        end

        function val = getMagneticField(this, onBoardModelValue)
            % onBoardModelValue - magnetic field in the satellite body frame
            % as known to the onboard model
            val = onBoardModelValue + normrnd(0, this.magnFieldSigma, [3, 1]);
        end

        function val = getDisturbanceTorque(this)
            val = normrnd(0, this.distTorqueSigma, [3, 1]);
        end        
    end
end