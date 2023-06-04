classdef Magnetorquer < AbstractActuator

    properties
        area            % [m^2], average area of the current loops in the mtq winding
        nTurns          % [-], number of turns in the winding

        resistance      % [Ohm], wire resistance
        alpha = 0       % [-], resistance temperature coefficient
    end

    methods
        function this = Magnetorquer(parameters)
            arguments
                parameters.name {mustBeText} = '';
                parameters.state {mustBeNumericOrLogical} = true;
                parameters.axis {mustBeNumeric} = [0; 0; 1];
                parameters.dcm {mustBeNumeric} = eye(3);
                parameters.noiseSigma {mustBeNumeric} = 1e-10;
                parameters.area {mustBeNumeric} = 0.05^2;
                parameters.nTurns {mustBeNumeric} = 200;
                parameters.resistance {mustBeNumeric} = 25 
                parameters.alpha {mustBeNumeric} = 0;
            end

            this.name       = parameters.name;
            this.state      = parameters.state;
            this.axis       = parameters.axis;
            this.dcm        = parameters.dcm;
            this.axisBf     = this.dcm * this.axis;

            this.noiseSigma = parameters.noiseSigma;

            this.area       = parameters.area;
            this.nTurns     = parameters.nTurns;
            this.resistance = parameters.resistance;
            this.alpha      = parameters.alpha;
        end

        % simulates the magnetic moment output along the output axis
        % given the input voltage (PWM is not simulated)
        function magneticMoment = actuateMagneticMoment(this, voltage)
            res = this.resistance * (1 + this.alpha);
            pureMagneticMoment = voltage / res * (this.area * this.nTurns);

            magneticMoment = this.actuateRequiredAction(pureMagneticMoment);
        end

        % expected max magnetic moment given the max voltage
        function maxMagneticMoment = calcMaxMagneticMoment(this, maxVoltage)
            maxMagneticMoment = maxVoltage / this.resistance * (this.area * this.nTurns) * this.axisBf;
        end
    end
end