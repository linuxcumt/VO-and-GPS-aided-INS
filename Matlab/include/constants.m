%% Constants.m
% 
% Defines global constants for various calculations

global GRAVCONST MASSEARTH MUEARTH A_EARTH B_EARTH OMEGA_EARTH ECEF2LATLON_MAXITER ECEF2LATLON_PRECISION DEGRAD

% Earth and gravity parameters
GRAVCONST = 6.674e-11; % [m^3/kg/s^2] Gravitational Constant for Earth
MASSEARTH = 5.972e24; % [kg] Earth mass
MUEARTH = GRAVCONST*MASSEARTH; % [m^3/s^2] Earth Gravity constant
OMEGA_EARTH = 7.292115e-5; % [rad/s]

% Latitude and Longitude calculations
A_EARTH = 6378137.00000;   % [m]
B_EARTH = 6356752.31425;   % [m]
ECEF2LATLON_MAXITER = 100;
ECEF2LATLON_PRECISION = eps;
DEGRAD = pi/180;