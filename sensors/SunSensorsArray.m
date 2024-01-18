classdef SunSensorsArray < handle

    properties(SetAccess = protected, GetAccess = public)
        sensors    % array of sun sensors
    end

    methods

        function this = SunSensorsArray(parameters)
            arguments
                parameters.baselineSensor SunSensor
                parameters.dcm(:, 3) {mustBeNumeric} = eye(3)
            end

            sensorsCount = size(parameters.dcm, 1) / 3;

            sunSensors(1, sensorsCount) = SunSensor();

            for sensorIdx = 1:sensorsCount
                idx = 3 * sensorIdx;
                sunSensors(sensorIdx) = SunSensor(bias = parameters.baselineSensor.bias, ...
                                                  sigma = parameters.baselineSensor.noiseSigma, ...
                                                  dcm = parameters.dcm((idx - 2):(idx), 1:3), ...
                                                  fovDeg = rad2deg(parameters.baselineSensor.fov), ...
                                                  boresightDirection = parameters.baselineSensor.boresightDirection);
            end

            this.sensors = sunSensors;
        end

        function [maxDirection, maxIntensity] = getSensorsReadings(this, sunDirection, sunEclipse)
            arguments
                this
                sunDirection(3, 1) {mustBeNumeric}
                sunEclipse {mustBeNumericOrLogical} = false
            end

            sensorsCount = length(this.sensors);
            signalIntensity = zeros(sensorsCount, 1);
            direction = zeros(3, sensorsCount);

            for sensorIdx = 1:sensorsCount
                [direction(:, sensorIdx), signalIntensity(sensorIdx)] = ...
                    this.sensors(sensorIdx).getSensorReadings(sunDirection, sunEclipse);
            end

            [~, maxIntencityIdx] = max(signalIntensity);
    
            % Extracting the direction and intensity of the sensor with the highest intensity
            maxDirection = direction(:, maxIntencityIdx);
            maxIntensity = signalIntensity(maxIntencityIdx);
        end

    end
end