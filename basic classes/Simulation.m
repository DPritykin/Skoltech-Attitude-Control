classdef Simulation < handle

    properties
        simulationTime

        odeOptions = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

        sat     % Satellite

        orb     % Orbit

        env     % environment

        ekf     % Kalman filter
    end

    methods
        function this = Simulation(simTime, odeOptions)
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
            if isa(ekf, 'KalmanFilter')
                this.ekf = ekf;
            else
                error('Simulation:InvalidFilter', 'Invalid EKF object!')
            end
        end

        % runs simulation
        function simResults = run(this, simulationType, omega0, q0)
            % initial conditions and integration settings
            if ~exist('q0', 'var')
                q0 = rand(4, 1);
                q0 = q0 / vecnorm(q0);
            end

            if ~exist('omega0', 'var')
                omega0 = 10 * this.orb.meanMotion * normrnd(0, 1, [3, 1]);
            end

            switch simulationType
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

        function simResults = simulateThreeAxialControl(this, q0, omega0)
            t = 0;
            mCtrl = [0; 0; 0];
            stateEst = [1; 0; 0; 0; 0; 0; 0;];
            simResults = zeros(9, ceil(this.simulationTime / this.sat.controlParams.tLoop));
            startTime = datetime('2021-03-14 01:00:00', 'TimeZone', 'UTC');
            for iterIdx = 1:size(simResults, 2)

                %% controlled dynamics (magnetorquers on)
                simTime = t + this.sat.controlParams.tCtrl;

                stateVec = this.integrate([t simTime], [q0; omega0], mCtrl, SimulationType.fullMagneticControl);

                t = simTime;
                q0 = stateVec(end, 1:4)' / vecnorm(stateVec(end, 1:4));
                omega0 = stateVec(end, 5:7)';

                %% measurements (magnetorquers off)
                simTime = t + this.sat.controlParams.tMeas;

                stateVec = this.integrate([t simTime], [q0; omega0]);

                t = simTime;
                q0 = stateVec(end, 1:4)' / norm(stateVec(end, 1:4));
                omega0 = stateVec(end, 5:7)';

                SS_VecT_icrs = this.env.SunVecCalc(t, startTime);
                SS_VecT_icrs = SS_VecT_icrs/norm(SS_VecT_icrs);

                TransformedT = this.orb.eci2orb(this.orb.meanMotion * t);
                SS_Vec_ModelT = TransformedT * SS_VecT_icrs;

                [ePos, ~] = this.orb.calcUnitPosVel(this.orb.meanMotion * t);
                ePosAbs = ePos .* [this.orb.orbitRadius; 0; 0];

                eclipse = this.env.sunEclipse(TransformedT, ePosAbs, t, startTime);

                [mtmMeasuredField] = this.calcSensorMagnField(t, q0);
                [ssMeasuredVector, intSS] = this.calcSSVector(t, q0, startTime, eclipse);

                %% state estimation (Extended Kalman filter)
                t0 = t - this.sat.controlParams.tLoop;

                bModel0 = this.env.directDipoleOrbital(this.orb.meanMotion * t0, this.orb.inclination, this.orb.orbitRadius);
                bmodelT = this.env.directDipoleOrbital(this.orb.meanMotion * t, this.orb.inclination, this.orb.orbitRadius);

                stateEst = this.ekf.estimate(t0, stateEst, mCtrl, bModel0, bmodelT, mtmMeasuredField,ssMeasuredVector,SS_Vec_ModelT);
                qEst = stateEst(1:4);
                omegaEst = stateEst(5:7);

                %% control moment for the next control loop (based on the Kalman estimate of the state)                
                omegaRel = omegaEst - quatRotate(qEst, [0; this.orb.meanMotion; 0]);
                mCtrl = this.sat.calcControlMagneticMoment(qEst, omegaRel, mtmMeasuredField);

                simResults(:, iterIdx) = [t; q0; omega0; intSS];
            end
        end

        function simResults = simulateBdotControl(this, q0, omega0)
            t = 0;
            mCtrl = [0; 0; 0];
            simResults = zeros(8, ceil(this.simulationTime / this.sat.controlParams.tLoop));

            for iterIdx = 1:size(simResults, 2)

                %% controlled dynamics (magnetorquers on)
                simTime = t + this.sat.controlParams.tCtrl;
                % TODO: need to make sure the control period is not enough to
                % start accelerating the sat

                stateVec = this.integrate([t simTime], [q0; omega0], mCtrl, SimulationType.bDotControl);

                t = simTime;
                q0 = stateVec(end, 1:4)' / vecnorm(stateVec(end, 1:4));
                omega0 = stateVec(end, 5:7)';

                %% measurements (magnetorquers off)
                simTime = t + this.sat.controlParams.tMeas;

                stateVec = this.integrate([t simTime], [q0; omega0]);

                t = simTime;
                q0 = stateVec(end, 1:4)' / norm(stateVec(end, 1:4));
                omega0 = stateVec(end, 5:7)';

                %% measurements
                bmodelT = this.env.directDipoleOrbital(this.orb.meanMotion * t, this.orb.inclination, this.orb.orbitRadius);

                %% control moment for the next control loop (based on the measurements)                
                mCtrl = this.sat.calcBdotMagneticMoment(bSensor, omegaSensor);

                simResults(:, iterIdx) = [t; q0; omega0];
            end
        end

        function simResults = simulateRwControl(this, q0, omega0, rwAngMomentum0)

            if ~exist('rwAngMomentum0', 'var')
                rwAngMomentum0 = zeros(3, 1);
            end

            t = 0;
            rwCtrl = [0; 0; 0];
            stateEst = [1; 0; 0; 0; 0; 0; 0];
            simResults = zeros(12, ceil(this.simulationTime / this.sat.controlParams.tLoop));
            startTime = datetime('2021-03-14 01:00:00', 'TimeZone', 'UTC');

            for iterIdx = 1:size(simResults, 2)

                %% controlled dynamics
                simTime = t + this.sat.controlParams.tLoop;

                stateVec = this.integrate([t simTime], ...
                                          [q0; omega0; rwAngMomentum0], ...
                                          rwCtrl, ...
                                          SimulationType.rwControl);

                t = simTime;
                q0 = stateVec(end, 1:4)' / vecnorm(stateVec(end, 1:4));
                omega0 = stateVec(end, 5:7)';
                rwAngMomentum0 = stateVec(end, 8:10)';

                %% state estimation (Extended Kalman filter)
                t0 = t - this.sat.controlParams.tLoop;

                bModel0 = this.env.directDipoleOrbital(this.orb.meanMotion * t0, this.orb.inclination, this.orb.orbitRadius);
                bmodelT = this.env.directDipoleOrbital(this.orb.meanMotion * t, this.orb.inclination, this.orb.orbitRadius);

                SS_VecT_icrs = this.env.SunVecCalc(t, startTime);
                SS_VecT_icrs = SS_VecT_icrs / vecnorm(SS_VecT_icrs);

                TransformedT = this.orb.eci2orb(this.orb.meanMotion * t);
                SS_Vec_ModelT = TransformedT * SS_VecT_icrs;

                [ePos, ~] = this.orb.calcUnitPosVel(this.orb.meanMotion * t);
                ePosAbs = ePos .* [this.orb.orbitRadius; 0; 0];

                eclipse = this.env.sunEclipse(TransformedT, ePosAbs, t, startTime);

                [mtmMeasuredField] = this.calcSensorMagnField(t, q0);
                [ssMeasuredVector, intSS] = this.calcSSVector(t, q0, startTime, eclipse);

                stateEst = this.ekf.estimate(t0, stateEst, rwCtrl, bModel0, bmodelT, mtmMeasuredField, ssMeasuredVector, SS_Vec_ModelT);
                qEst = stateEst(1:4);
                omegaEst = stateEst(5:7);

                %% control torque for the next control loop (based on the Kalman estimate of the state)  
                ez_b = quatRotate(qEst, [0; 0; 1]);
                trqGrav = 3 * this.orb.meanMotion^2 * crossProduct(ez_b, this.sat.J * ez_b);
                externalTorqueToCompensate = trqGrav - - crossProduct(omegaEst, (this.sat.J) * omegaEst + rwAngMomentum0);
                rwCtrl = this.sat.calcRwControl(qEst, omegaEst, rwAngMomentum0, externalTorqueToCompensate);

                simResults(:, iterIdx) = [t; q0; omega0; rwAngMomentum0; intSS];
            end
        end


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

        function envMagnField = calcEnvMagnField(this, t, q)
            directDipole = this.env.directDipoleOrbital(this.orb.meanMotion * t, ...
                                                        this.orb.inclination, ...
                                                        this.orb.orbitRadius);
            onBoardModelField = quatRotate(q, directDipole);
            envMagnField = this.env.getMagneticField(onBoardModelField);
        end

        function [sensedMagnField] = calcSensorMagnField(this, t, q)
            envMagnField = this.calcEnvMagnField(t, q);
            [sensedMagnField,~] = this.sat.mtm.getSensorReadings(envMagnField);
        end

        function [sensedSunVec, intSS] = calcSSVector(this, t, q, startTime, eclipse)

            TransformedT = this.orb.eci2orb(this.orb.meanMotion * t);
            SunVecICRS = this.env.SunVecCalc(t,startTime);
            SunVecICRS = SunVecICRS / vecnorm(SunVecICRS);
            SunVec = TransformedT * SunVecICRS;
            [sensedSunVec,intSS] = this.sat.ss.getSensorReadings(quatRotate(q, SunVec), eclipse);

        end

    end
end
    
           