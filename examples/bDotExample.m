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
standardMtqArray = MtqArray(baselineMtq = mtq, ...      % a Magnrtorquer object
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
simResults = sim.run('bDotControl', omega0);

sim.plotResults(simResults);
