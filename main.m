clc
clear

%% environment settings 

% dstSigma [N * m] disturbance of the torque
% magnSigma [T] noise to create actual magnetic field
env = Environment(dstSigma = 3e-9, magnSigma = 2e-7);

%% orbit settings

% environment, altitude [km], and inclination [deg]
orb = CircularOrbit(env, 470, 52);

%% satellite parameters structure

% [kg * m^2] inertia tensor for satellite
sat = Satellite(diag([0.015 0.015 0.007]));

sat.setControlParams(tMeas = 0.5, ...           % [s] sampling time step
                     tCtrl = 3, ...             % [s] control time step
                     qReq = [1; 0; 0; 0], ...   % [-] required orbital orientation
                     kQ = 12, ...               % [N * m / T^2] orientating parameter for PID-regulator
                     kW = 60 / orb.meanMotion)  % [N * m * s / T^2] stabilizing parameter for PID-regulator

mtm = Magnetometer(bias = [0; 0; 0;], ...         % [T] magnetometer bias
                   sigma = 1e-7, ...              % [T] magnetometer measurement deviation
                   position = [2; 2; -3] * 1e-2); % [T] magnetometer position in the body-frame

sat.setMagnetometer(mtm);


%% simulation settings

simulationTime = 10 * 3600;
sim = Simulation(simulationTime);

sim.setEnvironment(env);
sim.setOrbit(orb);
sim.setSatellite(sat);

%% simulation loop
simRunsCount = 1;
RMSE = zeros(6, simRunsCount);

for i = 1:simRunsCount
    simResults = sim.run();    
    sim.plotResults(simResults);
end