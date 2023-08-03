classdef Environment < handle

    properties(SetAccess = protected, GetAccess = public)

        muG = 3.986e+14;         % [m^3 / s^2] standard gravitational parameter of the Earth
        mu0 = 1.257e-6;          % [N / A^2] vacuum permeability
        muE = 7.94e+22;          % [A * m^2] magnetic dipole moment of the Earth
        earthRadius = 6371e+3;   % [m] radius of the Earth
        km2m = 1e3;

        distTorqueSigma = 0;     % [N * m] disturbance of the torque
        magnFieldSigma = 0;      % [T] noise to create actual magnetic field
    end

    methods
        function this = Environment(parameters)
            arguments
                parameters.distTorqueSigma {mustBeNumeric} = [0; 0; 0];
                parameters.magnFieldSigma {mustBeNumeric} = [0; 0; 0];
            end

            this.distTorqueSigma = parameters.distTorqueSigma;
            this.magnFieldSigma = parameters.magnFieldSigma;
        end

        function B = directDipoleOrbital(this, argLat, inclination, orbitRadius)
            % [T] magnitude of the magnetic field on the orbit
            B0 = this.muE * this.mu0 / (4 * pi * orbitRadius^3);   

            B = B0 * [cos(argLat) * sin(inclination); 
                      cos(inclination); 
                      -2 * sin(argLat) * sin(inclination)];
        end

        function S_Vec = SunVecCalc(this,time,startTime)

            baseTime = startTime;
            durations = seconds(time);
            timeValues = baseTime + durations;

            utc = datetime(timeValues, 'TimeZone', 'UTC');
            Mjd =  Mjday(year(utc), month(utc), day(utc), hour(utc), minute(utc), second(utc));
            Mjd_TDB = Mjday_TDB(Mjd);

            [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
            r_Neptune,r_Pluto,r_Moon,r_Sun] = JPL_Eph_DE430(Mjd_TDB);
           
            S_Vec = r_Sun/norm(r_Sun);

        end 

        % picewise exponential atmospeheric model
        function rhoAtmo = getAtmosphericDensity(this, altitude)
            altitudeKm = altitude / this.km2m;

            baseAltitudes = [1000:-100:500, 450:-50:200, 180, 150:-10:30, 25, 0];
            scaleAltitudes = [268, 181.05, 124.64, 88.667, 71.835, 63.822, 60.828, 58.515, 53.298, 53.628, 45.546, 37.105, 29.740, ...
                22.523, 16.149, 12.636, 9.473, 7.263, 5.877, 5.382, 5.799, 6.549, 7.714, 8.382, 7.554, 6.682, 6.349, 7.249];
            nominalDensities = [3.019e-15, 5.245e-15, 1.170e-14, 3.614e-14, 1.454e-13, 6.967e-13, 1.585e-12, 3.725e-12, 9.518e-12, ...
                2.418e-11, 7.248e-11, 2.789e-10, 5.464e-10, 2.07e-9, 3.845e-9, 8.484e-9, 2.438e-8, 9.661e-8, 5.297e-7, 3.396e-6, ...
                1.905e-5, 8.77e-5, 3.206e-4, 1.057e-3, 3.972e-3, 1.774e-2, 3.899e-2, 1.225];

            intervalNumber = sum(altitudeKm <= baseAltitudes, 2) + 1;
            baseAltitude = baseAltitudes(intervalNumber);
            scaleAltitude = scaleAltitudes(intervalNumber);
            nominalDensity = nominalDensities(intervalNumber);

            rhoAtmo = nominalDensity * exp((baseAltitude - altitudeKm) / scaleAltitude);
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