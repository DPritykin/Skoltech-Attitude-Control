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
        function [simResults,ekfResults] = run(this, simulationType, omega0, q0)
            % initial conditions and integration settings
            if ~exist('q0', 'var')
                q0 = rand(4, 1);
                q0 = q0 / vecnorm(q0);
            end

            if ~exist('omega0', 'var')
                omega0 = 10 * this.orb.meanMotion * normrnd(0, 1, [3, 1]);
            end

            switch lower(simulationType)
                case 'fullmagneticcontrol'
                    [simResults,ekfResults] = this.simulateThreeAxialControl(q0, omega0);
                case 'bdotcontrol'
                    simResults = this.simulateBdotControl(q0, omega0);
                otherwise
                    error('Simulation:InvalidSimulationType', 'Invalid Simulation type instruction!')
            end
        end
    end

    methods(Access = protected)

        function [simResults, ekfResults] = simulateThreeAxialControl(this, q0, omega0)
            t = 0;
            mCtrl = [0; 0; 0];
            
            stateEst = [1; 0; 0; 0; 0; 0; 0;]; % Initialization [q; w]
            if ~isempty(this.sat.gyro)
                stateEst(1:4) = q0;
                stateEst(5:7) =  this.sat.gyro.getSensorReadings(omega0);
            end
            mRes = zeros(3, 1);               % Initializing residual dipole
            mtmBias = zeros(3, 1);            % Initializing mtm bias 
            gyroBias = zeros(3, 1);           % Initializing gyro bias
            
            ekfStateEst = [stateEst; mRes; mtmBias; gyroBias]; % Initialization
            
            simResults = zeros(8, ceil(this.simulationTime / this.sat.controlParams.tLoop));
            ekfResults = zeros(17, ceil(this.simulationTime / this.sat.controlParams.tLoop));
             
            for iterIdx = 1:size(simResults, 2)

                %% controlled dynamics (magnetorquers on)
                simTime = t + this.sat.controlParams.tCtrl;

                stateVec = this.integrate([t simTime], [q0; omega0], mCtrl);

                t = simTime;
                q0 = stateVec(end, 1:4)' / vecnorm(stateVec(end, 1:4));
                omega0 = stateVec(end, 5:7)';

                %% measurements (magnetorquers off)
                simTime = t + this.sat.controlParams.tMeas;

                stateVec = this.integrate([t simTime], [q0; omega0], [0; 0; 0]);

                t = simTime;
                q0 = stateVec(end, 1:4)' / norm(stateVec(end, 1:4));
                omega0 = stateVec(end, 5:7)';

                % magnetometer
                mtmReadings = this.calcSensorMagnField(t, q0);

                % gyroscope (TODO: take several points, use sliding mean)
                if ~isempty(this.sat.gyro)
                    gyroReadings = this.sat.gyro.getSensorReadings(omega0);
                elseif isempty(this.sat.gyro)
                    gyroReadings = [];
                end

                %% state estimation (Extended Kalman filter)
                t0 = t - this.sat.controlParams.tLoop;
                bModel0 = this.env.directDipoleOrbital(this.orb.meanMotion * t0, this.orb.inclination, this.orb.orbitRadius);
                bmodelT = this.env.directDipoleOrbital(this.orb.meanMotion * t, this.orb.inclination, this.orb.orbitRadius);
                
                ekfStateEst = this.ekf.estimate(t0, ekfStateEst, mCtrl, bModel0, bmodelT, mtmReadings, gyroReadings);
                qEst = ekfStateEst(1:4);
                omegaEst = ekfStateEst(5:7);
                mResEst = ekfStateEst(8:10);
                mtmBiasEst = ekfStateEst(11:13);
                gyroBiasEst = ekfStateEst(14:16);
                
                ekfResults(:, iterIdx) = [t; qEst; omegaEst; mResEst; mtmBiasEst; gyroBiasEst];

                %% control moment for the next control loop (based on the Kalman estimate of the state)
                omegaRel = omegaEst - quatRotate(qEst, [0; this.orb.meanMotion; 0]);
                mCtrl = this.sat.calcControlMagneticMoment(qEst, omegaRel, mtmReadings, mResEst);

                simResults(:, iterIdx) = [t; q0; omega0];
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

                stateVec = this.integrate([t simTime], [q0; omega0], mCtrl);

                t = simTime;
                q0 = stateVec(end, 1:4)' / vecnorm(stateVec(end, 1:4));
                omega0 = stateVec(end, 5:7)';

                %% measurements (magnetorquers off)
                simTime = t + this.sat.controlParams.tMeas;

                stateVec = this.integrate([t simTime], [q0; omega0], [0; 0; 0]);

                t = simTime;
                q0 = stateVec(end, 1:4)' / norm(stateVec(end, 1:4));
                omega0 = stateVec(end, 5:7)';

                %% measurements
                bmodelT = this.env.directDipoleOrbital(this.orb.meanMotion * t, this.orb.inclination, this.orb.orbitRadius);
                bSensor = this.sat.mtm.getSensorReadings(quatRotate(q0, bmodelT));
                omegaSensor = this.sat.gyro.getSensorReadings(omega0);

                %% control moment for the next control loop (based on the measurements)                
                mCtrl = this.sat.calcBdotMagneticMoment(bSensor, omegaSensor);

                simResults(:, iterIdx) = [t; q0; omega0];
            end
        end

        function stateVec = integrate(this, timeInterval, initialConditions, mCtrl)
            envB = this.calcEnvMagnField(timeInterval(1), initialConditions(1:4));
            distTorque = this.env.getDisturbanceTorque();
            mRes = this.sat.calcResidualDipoleMoment();
            
            [ ~, stateVec ] = ode45(@(t, x) rhsRotationalDynamics(t, x, this.sat, this.orb, envB, mCtrl, mRes, distTorque), ...
                                    timeInterval, initialConditions, this.odeOptions);
        end

        function envMagnField = calcEnvMagnField(this, t, q)
            directDipole = this.env.directDipoleOrbital(this.orb.meanMotion * t, ...
                                                        this.orb.inclination, ...
                                                        this.orb.orbitRadius);
            onBoardModelField = quatRotate(q, directDipole);
            envMagnField = this.env.getMagneticField(onBoardModelField);
        end

        function sensedMagnField = calcSensorMagnField(this, t, q)
            envMagnField = this.calcEnvMagnField(t, q);
            inducedField = this.sat.calcResidualMagnFieldAtPosition(this.sat.mtm.position); 
            sensedMagnField = this.sat.mtm.getSensorReadings(envMagnField, inducedField);
        end
    end
end
