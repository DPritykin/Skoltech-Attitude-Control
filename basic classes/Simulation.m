classdef Simulation < handle

    properties
        startDate

        simulationTime

        odeOptions = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

        sat     % Satellite

        orb     % Orbit

        env     % environment

        ekf     % Kalman filter
    end

    methods
        function this = Simulation(simTime, startDate, odeOptions)
            arguments
                simTime {mustBeNumeric}
                startDate(1, 6) {mustBeNumeric} = [2000 1 1 12 0 0];
                odeOptions = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
            end

            this.startDate = datetime(startDate);
            this.simulationTime = simTime;

            if nargin > 1
                this.odeOptions = odeOptions;
            end
        end

        % adds Environment to simulation
        function setEnvironment(this, env)
            if isa(env, 'Environment')
                this.env = env;
            else
                error('Simulation:InvalidEnvironment', 'Invalid Environment object!')
            end
        end

        % adds Orbit to simulation
        function setOrbit(this, orb)
            if isa(orb, 'CircularOrbit')
                this.orb = orb;
            else
                error('Simulation:InvalidOrbit', 'Invalid Orbit object!')
            end
        end

        % adds Satellite to simulation
        function setSatellite(this, sat)
            if isa(sat, 'Satellite')
                this.sat = sat;
            else
                error('Simulation:InvalidSatelite', 'Invalid Satellite object!')
            end
        end

        % adds EKF to simulation
        function setFilter(this, ekf)
            if isa(ekf, 'AbstractKalmanFilter')
                this.ekf = ekf;
            else
                error('Simulation:InvalidFilter', 'Invalid EKF object!')
            end
        end

        % runs simulation
        function simResults = run(this, simulationType, omega0, q0)
            arguments
                this
                simulationType SimulationType
                omega0(3, 1) {mustBeNumeric} = 10 * this.orb.meanMotion * normrnd(0, 1, [3, 1])
                q0(4, 1) {mustBeNumeric} = rand(4, 1)
            end
            % initial conditions and integration settings
            q0 = q0 / vecnorm(q0);

            switch simulationType
                case SimulationType.noControl
                    simResults = this.simulateFreeRotation(q0, omega0);
                case SimulationType.fullMagneticControl
                    simResults = this.simulateThreeAxialControl(q0, omega0);
                case SimulationType.bDotControl
                    simResults = this.simulateBdotControl(q0, omega0);
                case SimulationType.rwControl
                    simResults = this.simulateRwControl(q0, omega0);
                otherwise
                    error('Simulation:InvalidSimulationType', 'Invalid Simulation type instruction!')
            end
        end
    end

    methods(Access = protected)
        
        function simResults = simulateFreeRotation(this, qDyn, omegaDyn)
            arguments
                this
                qDyn(4, 1) {mustBeNumeric}
                omegaDyn(3, 1) {mustBeNumeric}
            end

            tSpan = 0:1:this.simulationTime;
            [T, X] = ode45(@(t, x) rhsRotationalDynamics(t, x, this.sat, this.orb), ...
                tSpan, [qDyn; omegaDyn], this.odeOptions);

            simResults = [T, X];
            sunModelBody = zeros(size(simResults, 1), 3);

            for idx = 1:length(T)
                t = T(idx);
                argLat = this.orb.meanMotion * t;
                sunModelEci = Environment.getSunDirectionEci(datevec(this.startDate + t / 86400));
                sunModelOrb = this.orb.eci2orb(argLat) * sunModelEci;
                sunModelBody(idx, :) = quatRotate(X(idx, 1:4), sunModelOrb);
            end

            simResults = [simResults, sunModelBody]';
        end

        % rotational dynamics with three-axial magnetic controller
        function simResults = simulateThreeAxialControl(this, qDyn, omegaDyn)
            arguments
                this
                qDyn(4, 1) {mustBeNumeric}
                omegaDyn(3, 1) {mustBeNumeric}
            end

            t = 0;
            mCtrl = [0; 0; 0];
            stateEst = [1; 0; 0; 0; 0; 0; 0;];
            simResults = zeros(8, ceil(this.simulationTime / this.sat.controlParams.tLoop));

            for iterIdx = 1:size(simResults, 2)

                %% controlled dynamics (magnetorquers on)
                simTime = t + this.sat.controlParams.tCtrl;

                stateVec = this.integrate([t simTime], [qDyn; omegaDyn], mCtrl, SimulationType.fullMagneticControl);

                t = simTime;
                qDyn = stateVec(end, 1:4)' / vecnorm(stateVec(end, 1:4));
                omegaDyn = stateVec(end, 5:7)';

                %% measurements (magnetorquers off)
                simTime = t + this.sat.controlParams.tMeas;

                stateVec = this.integrate([t simTime], [qDyn; omegaDyn]);

                t = simTime;
                qDyn = stateVec(end, 1:4)' / norm(stateVec(end, 1:4));
                omegaDyn = stateVec(end, 5:7)';
                mtmMeasuredField = this.calcSensorMagnField(t, qDyn);

                %% state estimation (Extended Kalman filter)
                t0 = t - this.sat.controlParams.tLoop;
                bModel0 = this.env.directDipoleOrbital(this.orb.meanMotion * t0, this.orb.inclination, this.orb.orbitRadius);
                bmodelT = this.env.directDipoleOrbital(this.orb.meanMotion * t, this.orb.inclination, this.orb.orbitRadius);

                stateEst = this.ekf.estimate(t0, stateEst, bModel0, bmodelT, mtmMeasuredField, mCtrl);
                qEst = stateEst(1:4);
                omegaEst = stateEst(5:7);

                %% control moment for the next control loop (based on the Kalman estimate of the state)                
                omegaRel = omegaEst - quatRotate(qEst, [0; this.orb.meanMotion; 0]);
                mCtrl = this.sat.calcControlMagneticMoment(qEst, omegaRel, mtmMeasuredField);

                simResults(:, iterIdx) = [t; qDyn; omegaDyn];
            end
        end

        % rotational dynamics with b-dot controller
        function simResults = simulateBdotControl(this, qDyn, omegaDyn)
            arguments
                this
                qDyn(4, 1) {mustBeNumeric}
                omegaDyn(3, 1) {mustBeNumeric}
            end

            t = 0;
            mCtrl = [0; 0; 0];
            simResults = zeros(8, ceil(this.simulationTime / this.sat.controlParams.tLoop));

            for iterIdx = 1:size(simResults, 2)

                %% controlled dynamics (magnetorquers on)
                simTime = t + this.sat.controlParams.tCtrl;
                % TODO: need to make sure the control period is not enough to
                % start accelerating the sat

                stateVec = this.integrate([t simTime], [qDyn; omegaDyn], mCtrl, SimulationType.bDotControl);

                t = simTime;
                qDyn = stateVec(end, 1:4)' / vecnorm(stateVec(end, 1:4));
                omegaDyn = stateVec(end, 5:7)';

                %% measurements (magnetorquers off)
                simTime = t + this.sat.controlParams.tMeas;

                stateVec = this.integrate([t simTime], [qDyn; omegaDyn]);

                t = simTime;
                qDyn = stateVec(end, 1:4)' / norm(stateVec(end, 1:4));
                omegaDyn = stateVec(end, 5:7)';

                %% measurements
                bmodelT = this.env.directDipoleOrbital(this.orb.meanMotion * t, this.orb.inclination, this.orb.orbitRadius);
                bSensor = this.sat.mtm.getSensorReadings(quatRotate(qDyn, bmodelT));
                omegaSensor = this.sat.gyro.getSensorReadings(omegaDyn);

                %% control moment for the next control loop (based on the measurements)                
                mCtrl = this.sat.calcBdotMagneticMoment(bSensor, omegaSensor);

                simResults(:, iterIdx) = [t; qDyn; omegaDyn];
            end
        end

        % rotational dynamics with reaction wheels controller
        function simResults = simulateRwControl(this, qDyn, omegaDyn, rwAngMomentumDyn)
            arguments
                this
                qDyn(4, 1) {mustBeNumeric}
                omegaDyn(3, 1) {mustBeNumeric}
                rwAngMomentumDyn(3, 1) {mustBeNumeric} = zeros(3, 1)
            end

            t = 0;
            rwCtrl = zeros(3, 1);
            stateEst = zeros(10, 1);
            simResults = zeros(11, ceil(this.simulationTime / this.sat.controlParams.tLoop));

            for iterIdx = 1:size(simResults, 2)

                %% controlled dynamics
                simTime = t + this.sat.controlParams.tLoop;

                stateVec = this.integrate([t simTime], ...
                                          [qDyn; omegaDyn; rwAngMomentumDyn], ...
                                          rwCtrl, ...
                                          SimulationType.rwControl);

                t = simTime;
                qDyn = stateVec(end, 1:4)' / vecnorm(stateVec(end, 1:4));
                omegaDyn = stateVec(end, 5:7)';
                rwAngMomentumDyn = stateVec(end, 8:10)';

                %% control torque for the next control loop
                                
                if isempty(this.ekf) % no Kalman filter involved
                    qEst = qDyn;
                    omegaEst = omegaDyn;

                    calcControlFlag = true;
                else                 % run the Kalman filter and get the state estimate
                    mtmMeasuredField = this.calcSensorMagnField(t, qDyn);
                    argLat = this.orb.meanMotion * t;
                    bModelOrb = this.env.directDipoleOrbital(argLat, this.orb.inclination, this.orb.orbitRadius);
                    sunModelEci = Environment.getSunDirectionEci(datevec(this.startDate + t / 86400));
                    sunModelOrb = this.orb.eci2orb(argLat) * sunModelEci;
                    measuredSunDirection = this.sat.sunSensors.getSensorsReadings(quatRotate(qDyn, sunModelOrb)); % implement shadow

                    t0 = t - this.sat.controlParams.tLoop;
                    if t0 < 1e-4 % TODO find a reasonable threshold
                        stateEst(1:4) = triad(sunModelOrb, measuredSunDirection, bModelOrb, mtmMeasuredField);
                    end

                    stateEst = this.ekf.estimate(t0, stateEst, bModelOrb, mtmMeasuredField, sunModelOrb, measuredSunDirection, rwCtrl);
                    qEst = stateEst(1:4);
                    omegaEst = stateEst(5:7);

                    calcControlFlag = (this.ekf.traceP < 1e-4);
                end

                if calcControlFlag
                    ez_b = quatRotate(qEst, [0; 0; 1]);
                    trqGrav = 3 * this.orb.meanMotion^2 * crossProduct(ez_b, this.sat.J * ez_b);
                    externalTorqueToCompensate = trqGrav - ... 
                                                 crossProduct(omegaEst, (this.sat.J) * omegaEst + rwAngMomentumDyn);

                    rwCtrl = this.sat.calcRwControl(qEst, omegaEst, externalTorqueToCompensate);
                else
                    rwCtrl = [0; 0; 0];
                end

                simResults(:, iterIdx) = [t; qDyn; omegaDyn; rwAngMomentumDyn];
            end
        end

        % propagate the rotational dynamics of the satellite
        function stateVec = integrate(this, timeInterval, initialConditions, ctrlAction, simulationType)
            arguments
                this 
                timeInterval(1, 2) {mustBeNumeric}
                initialConditions(:, 1) {mustBeNumeric}
                ctrlAction(3, 1) {mustBeNumeric} = [0; 0; 0];
                simulationType SimulationType = SimulationType.noControl;
            end

            distTorque = this.env.getDisturbanceTorque();

            switch simulationType
                case SimulationType.noControl
                    ctrlTorque = [0; 0; 0];
                case SimulationType.rwControl
                    ctrlTorque = ctrlAction;
                case {SimulationType.bDotControl, SimulationType.fullMagneticControl}
                    envB = this.calcEnvMagnField(timeInterval(1), initialConditions(1:4));
                    ctrlTorque = crossProduct(ctrlAction, envB);
            end

            [ ~, stateVec ] = ode45(@(t, x) rhsRotationalDynamics(t, x, this.sat, this.orb, ctrlTorque, distTorque), ...
                timeInterval, initialConditions, this.odeOptions);
        end

        % environmental magnetic field (different from the onboard direct dipole model)
        function envMagnField = calcEnvMagnField(this, t, q)
            directDipole = this.env.directDipoleOrbital(this.orb.meanMotion * t, ...
                                                        this.orb.inclination, ...
                                                        this.orb.orbitRadius);
            onBoardModelField = quatRotate(q, directDipole);
            envMagnField = this.env.getMagneticField(onBoardModelField);
        end

        % magnetic field as measured by the sensor
        function sensedMagnField = calcSensorMagnField(this, t, q)
            envMagnField = this.calcEnvMagnField(t, q);
            sensedMagnField = this.sat.mtm.getSensorReadings(envMagnField);
        end
    end
end
