% An example of a satellite (3U CubeSat) in a circular orbit
% required attitude - orbital
% adcs sensor - magnetometer
% adcs actuators - 3 magnetorquers
% adcs state determination - Extended Kalman Filter

clc
clear
% addpath('C:\Users\HOOR BANO\Desktop\Skoltech\Thesis\MATLAB Files\Skoltech-Attitude-Control-main\math\Sun_Vector');

%% environment settings 

env = Environment(distTorqueSigma = 3e-9, ...  % [N * m] disturbance of the torque
                  magnFieldSigma = 2e-7);      % [T] noise to create actual magnetic field

%% orbit settings

orb = CircularOrbit(env, ... % Environment object
                    470, ... % [km] circular orbit altitude
                    52);     % [deg] orbit inclination

%% satellite settings

sat = Satellite(diag([0.015 0.014 0.007])); % [kg * m^2] inertia tensor for satellite

% defining control parameters
sat.setControlParams(tMeas = 0.5, ...           % [s] sampling time step
                     tCtrl = 3, ...             % [s] control time step
                     qReq = [1; 0; 0; 0], ...   % [-] required orbital orientation
                     kQ = 12, ...               % [N * m / T^2] orientating parameter for PID-regulator
                     kW = 60 / orb.meanMotion)  % [N * m * s / T^2] stabilizing parameter for PID-regulator

% adding a magnetometer
mtm = Magnetometer(bias = [0; 0; 0;], ...         % [T] magnetometer bias
                   sigma = 1e-7, ...              % [T] magnetometer measurement deviation
                   position = [2; 2; -3] * 1e-2); % [m] magnetometer position in the body-frame

sat.setMagnetometer(mtm);

% adding a sun sensor 
ss = SunSensor(bias = [0; 0; 0;], ...             % [rad] sun sensor bias
               sigma = deg2rad(0.2), ...          % [rad] sun sensor measurement deviation (CubeSense Gen-1 SS)
               position = [2; 2; -3] * 1e-2);     % [m] sun sensor position in the body-frame

sat.setSunSensor(ss);


% defining a magnetorquer
mtq = Magnetorquer(area = 0.05^2, ...             % [m^2] area of a coil
                   nTurns = 200, ...              % [-] number of turns in the winding  
                   resistance = 25);              % [Ohm] wire resistance

% setting up an array of 3 identical magnetorquers onboard of the sat
standardMtqArray = MtqArray(baselineMtq = mtq, ...      % a Magnrtorquer object
                            maxInputVoltage = 5, ...    % [V] maximum innput voltage available in EPS
                            doStandardXyzArray = true); % 3 mtqs along the X, y, and Z axes of the satellite body-frame

sat.setMtqArray(standardMtqArray);

%% EKF settings

ekf = KalmanFilter(sat = sat, ...      % Satellite object
                   orb = orb, ...      % CircularOrbit object 
                   env = env, ...      % Environment object
                   sigmaQ0 = 1, ...    % variance to initialize the error covariance matrix (quaternion part)
                   sigmaOmega0 = 0.1); % variance to initialize the error covariance matrix (omega part)

%% simulation settings

simulationTime = 1 * 3600;
sim = Simulation(simulationTime);

sim.setEnvironment(env);
sim.setOrbit(orb);
sim.setSatellite(sat);
sim.setFilter(ekf);

%% simulation loop
simResults = sim.run('fullMagneticControl');

plotResults(simResults, sim.orb.meanMotion);

%% plotter
function plotResults(simData, meanMotion)
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

    timeInHours = simData(1, :) / 3600;
    
    figure
    subplot(2, 3, 1)
    plot(timeInHours, pitch, 'Color', red, 'LineWidth', 2)
    hold on
    plot(timeInHours, roll, 'Color', green, 'LineWidth', 2)
    hold on
    plot(timeInHours, yaw, 'Color', blue, 'LineWidth', 2)
    grid on
    xlabel('Time in hours')
    ylabel('Euler Angles, [deg]')
    legend('pitch','roll','yaw');
    
    subplot(2, 3, 4)
    plot(timeInHours(startIndex:end), pitch(startIndex:end), 'Color', red, 'LineWidth', 2)
    hold on
    plot(timeInHours(startIndex:end), roll(startIndex:end), 'Color', green, 'LineWidth', 2)
    hold on
    plot(timeInHours(startIndex:end), yaw(startIndex:end), 'Color', blue, 'LineWidth', 2)
    grid on
    xlabel('Time in hours')
    ylabel('Euler Angles Errors, [deg]')
    legend_p = ['RMSE_{pitch} = ', num2str(rmsePitch, '%10.2e\n')];
    legend_r = ['RMSE_{roll} = ', num2str(rmseRoll, '%10.2e\n')];
    legend_y = ['RMSE_{yaw} = ', num2str(rmseYaw, '%10.2e\n')];
    legend(legend_p, legend_r, legend_y);
    
    subplot(2, 3, 2) % angular velocity
    plot(timeInHours, simData(6, 1:end), 'Color', red, 'LineWidth', 2)
    hold on
    plot(timeInHours, simData(7, 1:end), 'Color', green, 'LineWidth', 2)
    hold on
    plot(timeInHours, simData(8, 1:end), 'Color', blue, 'LineWidth', 2)
    grid on
    xlabel('Time in hours')
    ylabel('Angular Velocity Components, [rad/sec]')
    legend('\omega_1','\omega_2','\omega_3');
    
    subplot(2, 3, 5) % angular velocity
    plot(timeInHours(startIndex:end), simData(6, startIndex:end), 'Color', red, 'LineWidth', 2)
    hold on
    plot(timeInHours(startIndex:end), (simData(7, startIndex:end) - meanMotion), 'Color', green, 'LineWidth', 2)
    hold on
    plot(timeInHours(startIndex:end), simData(8, startIndex:end), 'Color', blue, 'LineWidth', 2)
    grid on
    xlabel('Time in hours')
    ylabel('Angular Velocity Errors, [rad/sec]')
    legend_1 = ['RMSE_{\omega_1} = ', num2str(rmseW1, '%10.2e\n')];
    legend_2 = ['RMSE_{omega_2} = ', num2str(rmseW2, '%10.2e\n')];
    legend_3 = ['RMSE_{\omega_3} = ', num2str(rmseW3, '%10.2e\n')];
    legend(legend_1, legend_2, legend_3);
    
    subplot(2, 3, 3) % quaternion
    plot(timeInHours, simData(2, 1:end), 'k', 'LineWidth', 2)
    hold on
    plot(timeInHours, simData(3, 1:end), 'Color', red, 'LineWidth', 2)
    hold on
    plot(timeInHours, simData(4, 1:end), 'Color', green, 'LineWidth', 2)
    hold on
    plot(timeInHours, simData(5, 1:end), 'Color', blue, 'LineWidth', 2)
    grid on
    xlabel('Time in hours')
    ylabel('Quaternion Components')
    legend('q0','q1','q2','q3');
end
