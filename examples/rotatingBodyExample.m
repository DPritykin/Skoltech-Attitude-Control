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

%% simulation settings

simulationTime = 1000;
sim = Simulation(simulationTime);

sim.setEnvironment(env);
sim.setOrbit(orb);
sim.setSatellite(sat);

%% simulation loop
omega0 = normrnd(0, 0.01, [3, 1]);
simResults = sim.run(SimulationType.noControl, omega0);

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
