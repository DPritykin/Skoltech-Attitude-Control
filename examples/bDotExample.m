clc
clear

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
                     tCtrl = 1);                % [s] control time step

% adding a magnetometer
mtm = Magnetometer(bias = [0; 0; 0;], ...         % [T] magnetometer bias
                   sigma = 1e-7, ...              % [T] magnetometer measurement deviation
                   position = [2; 2; -3] * 1e-2); % [m] magnetometer position in the body-frame

sat.setMagnetometer(mtm);

% adding a gyroscope
gyro = Gyroscope(bias = [0; 0; 0;], ...           % [s^-1] gyroscope bias
                 sigma = 1e-5);                   % [s^-1] gyroscope measurement deviation

sat.setGyroscope(gyro);

% defining a magnetorquer
mtq = Magnetorquer(area = 0.05^2, ...             % [m^2] area of a coil
                   nTurns = 200, ...              % [-] number of turns in the winding  
                   resistance = 25);              % [Ohm] wire resistance

% setting up an array of 3 identical magnetorquers onboard of the sat
standardMtqArray = MtqArray(baselineMtq = mtq, ...      % a Magnetorquer object
                            maxInputVoltage = 5, ...    % [V] maximum innput voltage available in EPS
                            doStandardXyzArray = true); % 3 mtqs along the X, y, and Z axes of the satellite body-frame

sat.setMtqArray(standardMtqArray);



%% simulation settings

simulationTime = 2 * 3600;
sim = Simulation(simulationTime);

sim.setEnvironment(env);
sim.setOrbit(orb);
sim.setSatellite(sat);

%% simulation loop
omega0 = normrnd(0, 0.5, [3, 1]);
simResults = sim.run(SimulationType.bDotControl, omega0);

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
    
    timeInHours = simData(1, :) / 3600;

    figure
    subplot(1, 2, 1)
    plot(timeInHours, pitch, 'Color', red, 'LineWidth', 2)
    hold on
    plot(timeInHours, roll, 'Color', green, 'LineWidth', 2)
    hold on
    plot(timeInHours, yaw, 'Color', blue, 'LineWidth', 2)
    grid on
    xlabel('Time in hours')
    ylabel('Euler Angles, [deg]')
    legend('pitch','roll','yaw');
       
    subplot(1, 2, 2) % angular velocity
    plot(timeInHours, simData(6, :), 'Color', red, 'LineWidth', 2)
    hold on
    plot(timeInHours, simData(7, :), 'Color', green, 'LineWidth', 2)
    hold on
    plot(timeInHours, simData(8, :), 'Color', blue, 'LineWidth', 2)
    grid on
    xlabel('Time in hours')
    ylabel('Angular Velocity Components, [rad/sec]')
    legend('\omega_1','\omega_2','\omega_3');    
end
