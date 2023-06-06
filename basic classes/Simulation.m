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

            switch lower(simulationType)
                case 'fullmagneticcontrol'
                    simResults = this.simulateThreeAxialControl(q0, omega0);
                case 'bdotcontrol'
                    simResults = this.simulateBdotControl(q0, omega0);
                otherwise
                    error('Simulation:InvalidSimulationType', 'Invalid Simulation type instruction!')
            end
        end

        function plotResults(this, simData)
            meanMotion = this.orb.meanMotion;
            [pitch, roll, yaw] = quat2angle(simData(2:5, 1:end)', 'YXZ');

            pitch = rad2deg(pitch);
            roll = rad2deg(roll);
            yaw = rad2deg(yaw);

            len = length(pitch);
            startIndex = round(len * 3 / 4);

            ref0 = zeros(len, 1);

            rmseRoll = sqrt(immse(roll(startIndex:end), ref0(startIndex:end)));
            rmsePitch = sqrt(immse(pitch(startIndex:end), ref0(startIndex:end)));
            rmseYaw = sqrt(immse(yaw(startIndex:end), ref0(startIndex:end)));
            rmseW1 = sqrt(immse(simData(6, startIndex:end)', ref0(startIndex:end)));
            rmseW2 = sqrt(immse(simData(7, startIndex:end)' - meanMotion, ref0(startIndex:end)));
            rmseW3 = sqrt(immse(simData(8, startIndex:end)', ref0(startIndex:end)));

            red = [203/255, 37/255, 37/255];
            green = [138/255, 181/255, 73/255];
            blue = [33/255, 144/255, 209/255];

            figure
            subplot(2, 3, 1)
            plot(simData(1, 1:end) / 3600, (pitch), 'Color', red, 'LineWidth', 2)
            hold on
            plot(simData(1, 1:end) / 3600, (roll), 'Color', green, 'LineWidth', 2)
            hold on
            plot(simData(1, 1:end) / 3600, (yaw), 'Color', blue, 'LineWidth', 2)
            grid on
            xlabel('Time in hours')
            ylabel('Euler Angles, [deg]')
            legend('pitch','roll','yaw');

            subplot(2, 3, 4)
            plot(simData(1, startIndex:end) / 3600, (pitch(startIndex:end)), 'Color', red, 'LineWidth', 2)
            hold on
            plot(simData(1, startIndex:end) / 3600, (roll(startIndex:end)), 'Color', green, 'LineWidth', 2)
            hold on
            plot(simData(1, startIndex:end) / 3600, (yaw(startIndex:end)), 'Color', blue, 'LineWidth', 2)
            grid on
            xlabel('Time in hours')
            ylabel('Euler Angles Errors, [deg]')
            legend_p = ['RMSE_{pitch} = ', num2str(rmsePitch, '%10.2e\n')];
            legend_r = ['RMSE_{roll} = ', num2str(rmseRoll, '%10.2e\n')];
            legend_y = ['RMSE_{yaw} = ', num2str(rmseYaw, '%10.2e\n')];
            legend(legend_p, legend_r, legend_y);

            subplot(2, 3, 2) % angular velocity
            plot(simData(1, 1:end) / 3600, simData(6, 1:end), 'Color', red, 'LineWidth', 2)
            hold on
            plot(simData(1, 1:end) / 3600, simData(7, 1:end), 'Color', green, 'LineWidth', 2)
            hold on
            plot(simData(1, 1:end) / 3600, simData(8, 1:end), 'Color', blue, 'LineWidth', 2)
            grid on
            xlabel('Time in hours')
            ylabel('Angular Velocity Components, [rad/sec]')
            legend('\omega_1','\omega_2','\omega_3');

            subplot(2, 3, 5) % angular velocity
            plot(simData(1, startIndex:end) / 3600, simData(6, startIndex:end), 'Color', red, 'LineWidth', 2)
            hold on
            plot(simData(1, startIndex:end) / 3600, (simData(7, startIndex:end) - meanMotion), 'Color', green, 'LineWidth', 2)
            hold on
            plot(simData(1, startIndex:end) / 3600, simData(8, startIndex:end), 'Color', blue, 'LineWidth', 2)
            grid on
            xlabel('Time in hours')
            ylabel('Angular Velocity Errors, [rad/sec]')
            legend_1 = ['RMSE_{\omega_1} = ', num2str(rmseW1, '%10.2e\n')];
            legend_2 = ['RMSE_{omega_2} = ', num2str(rmseW2, '%10.2e\n')];
            legend_3 = ['RMSE_{\omega_3} = ', num2str(rmseW3, '%10.2e\n')];
            legend(legend_1, legend_2, legend_3);

            subplot(2, 3, 3) % quaternion
            plot(simData(1, 1:end) / 3600, simData(2, 1:end), 'k', 'LineWidth', 2)
            hold on
            plot(simData(1, 1:end) / 3600, simData(3, 1:end), 'Color', red, 'LineWidth', 2)
            hold on
            plot(simData(1, 1:end) / 3600, simData(4, 1:end), 'Color', green, 'LineWidth', 2)
            hold on
            plot(simData(1, 1:end) / 3600, simData(5, 1:end), 'Color', blue, 'LineWidth', 2)
            grid on
            xlabel('Time in hours')
            ylabel('Quaternion Components')
            legend('q0','q1','q2','q3');
        end
    end

    methods(Access = protected)

        function simResults = simulateThreeAxialControl(this, q0, omega0)
            t = 0;
            mCtrl = [0; 0; 0];
            stateEst = [1; 0; 0; 0; 0; 0; 0;];
            simResults = zeros(8, ceil(this.simulationTime / this.sat.controlParams.tLoop));

            for iterIdx = 1:size(simResults, 2)

                %% controlled dynamics (magnetorquers on)
                simTime = t + this.sat.controlParams.tCtrl;

                stateVec = this.integrate([t simTime], ...
                                          [q0; omega0], ...
                                          mCtrl);

                t = simTime;
                q0 = stateVec(end, 1:4)' / vecnorm(stateVec(end, 1:4));
                omega0 = stateVec(end, 5:7)';

                %% measurements (magnetorquers off)
                simTime = t + this.sat.controlParams.tMeas;

                stateVec = this.integrate([t simTime], ...
                                          [q0; omega0], ...
                                          [0; 0; 0]);

                t = simTime;
                q0 = stateVec(end, 1:4)' / norm(stateVec(end, 1:4));
                omega0 = stateVec(end, 5:7)';
                mtmMeasuredField = this.calcSensorMagnField(t, q0);

                %% state estimation (Extended Kalman filter)
                t0 = t - this.sat.controlParams.tLoop;
                bModel0 = this.env.directDipoleOrbital(this.orb.meanMotion * t0, this.orb.inclination, this.orb.orbitRadius);
                bmodelT = this.env.directDipoleOrbital(this.orb.meanMotion * t, this.orb.inclination, this.orb.orbitRadius);
                omegaSensor = this.sat.gyro.getSensorReadings(omega0);
                
                stateEst = this.ekf.estimate(t0, stateEst, mCtrl, bModel0, bmodelT, mtmMeasuredField,omegaSensor);
                qEst = stateEst(1:4);
                omegaEst = stateEst(5:7);

                %% control moment for the next control loop (based on the Kalman estimate of the state)                
                omegaRel = omegaEst - quatRotate(qEst, [0; this.orb.meanMotion; 0]);
                mCtrl = this.sat.calcControlMagneticMoment(qEst, omegaRel, mtmMeasuredField);

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

                stateVec = this.integrate([t simTime], ...
                                          [q0; omega0], ...
                                          mCtrl);

                t = simTime;
                q0 = stateVec(end, 1:4)' / vecnorm(stateVec(end, 1:4));
                omega0 = stateVec(end, 5:7)';

                %% measurements (magnetorquers off)
                simTime = t + this.sat.controlParams.tMeas;

                stateVec = this.integrate([t simTime], ...
                                          [q0; omega0], ...
                                          [0; 0; 0]);

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

            [ ~, stateVec ] = ode45(@(t, x) rhsRotationalDynamics(t, x, this.sat, this.orb, envB, mCtrl, distTorque), ...
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
            sensedMagnField = this.sat.mtm.getSensorReadings(envMagnField);
        end
    end
end
