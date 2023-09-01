classdef SSArray < handle

    properties(SetAccess = protected, GetAccess = public)
        sunsensors    % array of sun sensors
    end

    methods
        function this = SSArray(parameters)
            arguments
                parameters.baselineSS SunSensor
                parameters.dcm {mustBeNumericOrLogical} 
                parameters.sensorCount {mustBeNumeric} 
            end
            
            SSArray(1, parameters.sensorCount) = SunSensor();

            for sensorIdx = 1:parameters.sensorCount
                SSArray(sensorIdx).setbias(parameters.baselineSS.bias);
                SSArray(sensorIdx).setnoiseSigma(parameters.baselineSS.noiseSigma);
                SSArray(sensorIdx).setDcm(parameters.dcm(-2+(3*sensorIdx):0+(3*sensorIdx),1:3));
                SSArray(sensorIdx).setfov(parameters.baselineSS.fov);
                SSArray(sensorIdx).setboresightDirection(parameters.baselineSS.boresightDirection);
            end

                this.sunsensors = SSArray;
        end

        function [direction, intensity] = getSensorReadingsArray(this, sunDirection, sunEclipse)
            nSensors = numel(this.sunsensors);
            intensity = zeros(nSensors, 1);
            direction = zeros(nSensors, 3);

            for sensorIdx = 1:nSensors
                [direction(sensorIdx, :), intensity(sensorIdx)] = this.sunsensors(sensorIdx).getSensorReadings(sunDirection, sunEclipse);
            end
        end
    end
end
