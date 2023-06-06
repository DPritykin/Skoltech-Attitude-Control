% An example of a satellite (3U CubeSat) in a circular orbit
% required attitude - orbital
% adcs sensor - magnetometer / gyroscope
% adcs actuators - 3 magnetorquers
% adcs state determination - Extended Kalman Filter

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
                     tCtrl = 3, ...             % [s] control time step
                     qReq = [1; 0; 0; 0], ...   % [-] required orbital orientation
                     kQ = 12, ...               % [N * m / T^2] orientating parameter for PID-regulator
                     kW = 60 / orb.meanMotion)  % [N * m * s / T^2] stabilizing parameter for PID-regulator

% adding a gyroscope
gyro = Gyroscope(bias = [0; 0; 0;], ...         % [s^-1] gyroscope bias
                 sigma = 1e-5);                   % [s^-1] gyroscope measurement deviatio
                  

sat.setGyroscope(gyro);

sat.setMagnetometer(mtm);

% adding a gyroscope
gyro = Gyroscope(bias = [0; 0; 0;], ...         % [T] magnetometer bias
                   sigma = 1.45e-6, ...              % [rad/s] magnetometer measurement deviation (0.3 deg/hr : STIM 377H Gyro)
                   position = [1; 1; -3] * 1e-2); % [m] magnetometer position in the body-frame

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

%% EKF settings

ekf = KalmanFilter(sat = sat, ...      % Satellite object
                   orb = orb, ...      % CircularOrbit object 
                   env = env, ...      % Environment object
                   sigmaQ0 = 1, ...    % variance to initialize the error covariance matrix (quaternion part)
                   sigmaOmega0 = 0.1); % variance to initialize the error covariance matrix (omega part)

%% simulation settings

simulationTime = 8 * 3600;
sim = Simulation(simulationTime);

sim.setEnvironment(env);
sim.setOrbit(orb);
sim.setSatellite(sat);
sim.setFilter(ekf);

%% simulation loop
simResults = sim.run('fullmagneticcontrol');

sim.plotResults(simResults);
