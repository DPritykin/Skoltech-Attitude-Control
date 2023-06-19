classdef CircularOrbit < handle

    properties(SetAccess = protected, GetAccess = public)
        altitude
        inclination
        raan = 0

        orbitRadius
        meanMotion
    end

    methods
        function this = CircularOrbit(env, altitudeKm, inclinationDeg, raanDeg)
            arguments
                env
                altitudeKm {mustBeNumeric}
                inclinationDeg  {mustBeNumeric}
                raanDeg  {mustBeNumeric} = 0;
            end

            this.altitude =  altitudeKm * 1e3;
            this.inclination =  deg2rad(inclinationDeg);
            this.raan =  deg2rad(raanDeg);

            this.orbitRadius = env.earthRadius + this.altitude;
            this.meanMotion = sqrt(env.muG / this.orbitRadius^3);
        end

        function [ePos, eVel] = calcUnitPosVel(this, argLat)
            raanMatrix = [cos(this.raan), -sin(this.raan), 0
                          sin(this.raan), cos(this.raan), 0
                          0, 0, 1];
            inclMatrix = [1, 0, 0
                          0, cos(this.inclination), -sin(this.inclination)
                          0, sin(this.inclination), cos(this.inclination)];
            argLatMatrix =  [cos(argLat), -sin(argLat), 0
                             sin(argLat), cos(argLat), 0
                             0, 0, 1];

            ePos = raanMatrix * inclMatrix * argLatMatrix * [1; 0; 0];
            eVel = raanMatrix * inclMatrix * argLatMatrix * [0; 1; 0];
        end

        function orb2eciMatrix = orb2eci(this, argLat)
            [ePos, eVel] = this.calcUnitPosVel(argLat);

            orb2eciMatrix = zeros(3);
            orb2eciMatrix(:, 1) = eVel;
            orb2eciMatrix(:, 2) = cross(ePos, eVel);
            orb2eciMatrix(:, 3) = ePos;
        end

        function eci2orbMatrix = eci2orb(this, argLat)
            eci2orbMatrix = this.orb2eci(argLat)';
        end
    end
end