
format long g
global PC    % Planetary Coefficients
global const % Astronomical Constants
load DE436Coeff.mat
PC = DE436Coeff;
%% Defining UTC times 
% Define the start and end times in UTC format
start_utc = datetime('2021-03-14 01:00:00', 'TimeZone', 'UTC');
end_utc = datetime('2021-03-14 02:34:37', 'TimeZone', 'UTC');
                                                  
% Convert the start and end times to numeric values
start_num = datenum(start_utc);
end_num = datenum(end_utc);

% Initialize the output arrays
num_seconds = round((end_num - start_num) * 86400);
date_second = zeros(num_seconds, 1);
sun_icrs = zeros(num_seconds, 3);
sun_eci = zeros(num_seconds, 3);

%% Iterations 
% Loop over each second between the start and end times
for i = 1:num_seconds
    
    % Calculate the datetime for the current second
    date_second(i,1) = start_num + (i-1)/86400;
    utc(i,1) = datetime(date_second(i,1), 'ConvertFrom', 'datenum', 'TimeZone', 'UTC');
    
    Mjd =  Mjday(year(utc(i,1)), month(utc(i,1)), day(utc(i,1)), hour(utc(i,1)), minute(utc(i,1)), second(utc(i,1)));
    Mjd_TDB = Mjday_TDB(Mjd);
    [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, r_Neptune,r_Pluto,r_Moon,r_Sun] = JPL_Eph_DE430(Mjd_TDB);

    sun_icrs(i,1:3) = r_Sun; 

    [RA, Dec, r] = cart2sph(sun_icrs(i,1), sun_icrs(i,2), sun_icrs(i,3));
    GAST = mean_sidereal_time(Mjd); % Calculate Greenwich Apparent Sidereal Time
    hour_angle = GAST - RA;
    [x_eq, y_eq, z_eq] = sph2cart(hour_angle, Dec, r);

    theta = -GAST; % Rotation angle in radians
    R = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
    sun_eci(i,:) = R * [x_eq; y_eq; z_eq];

    Sun_Vb(i,:) = reshape(Rq(i,:,:), [], 3) * sun_eci(i,:).';
    Sun_Vb(i,:) = Sun_Vb(i,:)/norm(Sun_Vb(i,:));
end 
Sun_Vb = reshape(Sun_Vb, num_seconds, 3, []);

%% Transformations & Noise (S1)
Sun_Dist = [0.00001 0 0 ; 0 0.00001 0 ; 0 0 0.00001];
v = sqrtm(Sun_Dist)*(randn(3,num_seconds));
Sun1_Vt = Sun_Vb + v';

%% Plotting
figure(1)
plot(Sun1_Vt(:,1),'Linewidth',1)
hold on
plot(Sun1_Vt(:,2),'Linewidth',1)
hold on
plot(Sun1_Vt(:,3),'Linewidth',1)
grid on
xlabel('Time (sec)')
ylabel('Sun Vectors')
title('Sun Sensor S1')
legend('X-axis','Y-axis','Z-axis')