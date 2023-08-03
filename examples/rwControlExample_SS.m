% An example of a satellite (3U CubeSat) in a circular orbit
% required attitude - orbital
% adcs sensor - none
% adcs actuators - 3 reaction wheels
% adcs state determination - none

clc
clear

%% environment settings 

env = Environment(distTorqueSigma = 3e-9, ...  % [N * m] disturbance of the torque
                  magnFieldSigma = 2e-7);      % [T] noise to create actual magnetic field

%% orbit settings

orb = CircularOrbit(env, ... % Environment object
                    470, ... % [km] circular orbit altitude
                    98);     % [deg] orbit inclination

%% satellite settings

 sat = Satellite([0.05466 -0.04e-3 -0.06e-3; -0.04e-3 0.05531 0.29e-3; -0.06e-3 0.29e-3 0.01201]); % [kg * m^2] inertia tensor for satellite
%  sat = Satellite(diag([0.015 0.014 0.007])); % [kg * m^2] inertia tensor for satellite

% defining control parameters
sat.setControlParams(tMeas = 0.02, ...              % [s] sampling time step
                     tCtrl = 0.1, ...               % [s] control time step
                     qReq = [1; 0; 0; 0], ...       % [-] required orbital orientation
                     kQ = 3, ...                  % [-] P gain for PID-regulator
                     kW = 3.3);                     % [-] D gain for PID-regulator

% adding a magnetometer
mtm = Magnetometer(bias = [0; 0; 0;], ...           % [T] magnetometer bias
                   sigma = 1e-7, ...                % [T] magnetometer measurement deviation
                   position = [2; 2; -3] * 1e-2);   % [m] magnetometer position in the body-frame

sat.setMagnetometer(mtm);

% adding a sun sensor 
ss = SunSensor(bias = [0; 0; 0;], ...             % [rad] sun sensor bias
               sigma = deg2rad(0.2), ...          % [rad] sun sensor measurement deviation (CubeSense Gen-1 SS)
               position = [2; 2; -3] * 1e-2, ...  % [m] sun sensor position in the body-frame
               fov = 120);                        % [deg] sun sensor field of view constraint 

sat.setSunSensor(ss);

% defining a reaction wheel
rw = ReactionWheel(maxTorque = 1e-3, ...        % [Nm] max torque RW can produce
                   maxAngMomentum = 1e-2);      % [Nms] max angular momentum RW can have

% setting up an array of 3 identical RWs onboard of the sat
standardRwArray = RwArray(baselineRw = rw, ...         % a ReactionWheel object
                          doStandardXyzArray = true);  % 3 RWs along the X, y, and Z axes of the satellite body-frame

sat.setRwArray(standardRwArray);

%% EKF settings

ekf = KalmanFilter(sat = sat, ...      % Satellite object
                   orb = orb, ...      % CircularOrbit object 
                   env = env, ...      % Environment object
                   sigmaQ0 = 1, ...    % variance to initialize the error covariance matrix (quaternion part)
                   sigmaOmega0 = 0.1); % variance to initialize the error covariance matrix (omega part)

%% simulation settings

simulationTime = 100;
sim = Simulation(simulationTime);

sim.setEnvironment(env);
sim.setOrbit(orb);
sim.setSatellite(sat);
sim.setFilter(ekf);

%% simulation loop
simResults = sim.run(SimulationType.rwControl);

plotResults(simResults);

%% plotter
function plotResults(simData)
    [pitch, roll, yaw] = quat2angle(simData(2:5, 1:end)', 'YXZ');
    
    pitch = rad2deg(pitch);
    roll = rad2deg(roll);
    yaw = rad2deg(yaw);

    red = [203/255, 37/255, 37/255];
    green = [138/255, 181/255, 73/255];
    blue = [33/255, 144/255, 209/255];

    timeInSeconds = simData(1, :);
    
    figure
    subplot(1, 3, 1)
    plot(timeInSeconds, pitch, 'Color', red, 'LineWidth', 2)
    hold on
    plot(timeInSeconds, roll, 'Color', green, 'LineWidth', 2)
    hold on
    plot(timeInSeconds, yaw, 'Color', blue, 'LineWidth', 2)
    grid on
    xlabel('Time in seconds')
    ylabel('Euler Angles, [deg]')
    legend('pitch','roll','yaw');

    subplot(1, 3, 2) % angular velocity
    plot(timeInSeconds, simData(6, 1:end), 'Color', red, 'LineWidth', 2)
    hold on
    plot(timeInSeconds, simData(7, 1:end), 'Color', green, 'LineWidth', 2)
    hold on
    plot(timeInSeconds, simData(8, 1:end), 'Color', blue, 'LineWidth', 2)
    grid on
    xlabel('Time in seconds')
    ylabel('Angular Velocity Components, [rad/sec]')
    legend('\omega_1','\omega_2','\omega_3');
       
    subplot(1, 3, 3) % quaternion
    plot(timeInSeconds, simData(2, 1:end), 'k', 'LineWidth', 2)
    hold on
    plot(timeInSeconds, simData(3, 1:end), 'Color', red, 'LineWidth', 2)
    hold on
    plot(timeInSeconds, simData(4, 1:end), 'Color', green, 'LineWidth', 2)
    hold on
    plot(timeInSeconds, simData(5, 1:end), 'Color', blue, 'LineWidth', 2)
    grid on
    xlabel('Time in seconds')
    ylabel('Quaternion Components')
    legend('q0','q1','q2','q3');
end