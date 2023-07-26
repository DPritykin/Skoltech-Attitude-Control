function GST = mean_sidereal_time(JD)
% This function calculates the mean Greenwich Apparent Sidereal Time
% (GST) for a given Julian Date (JD).

% Convert Julian Date to centuries since J2000.0
T = (JD - 2451545.0) / 36525;

% Calculate the mean Greenwich Sidereal Time (GMST) for the same time
GMST = 280.46061837 + 360.98564736629 * (JD - 2451545.0) ...
    + 0.000387933 * T^2 - T^3 / 38710000;

% Make sure GMST is within 0-360 degrees
GMST = mod(GMST, 360);

% Convert GMST to GST by adding the equation of the equinoxes
% The equation of the equinoxes takes into account the precession and nutation of the Earth's axis
% This implementation uses a simplified formula that assumes no change in precession or nutation
epsilon = 23.439291111; % Obliquity of the ecliptic for J2000.0, in degrees
EE = GMST + 0.00273790931 * sind(125.04 - 1934.136 * T); % Equation of the equinoxes, in degrees
GST = mod(EE + epsilon, 360);

end
