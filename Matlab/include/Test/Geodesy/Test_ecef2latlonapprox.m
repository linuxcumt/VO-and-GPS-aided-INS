%% Test_ecef2latlon
%
% Tests the conversion functions to go from Latitude and Longitude to 
% ECEF coordinates.
%
% @author: Matt Marti
% @date: 2018-08-15

clear;


%% Test 1
ecef = [1061985.2311866253;-4837817.3795761978;4005245.7246151697];
latlonalt = ecef2latlon(ecef);
[ latlonapproxp, dlatlonapp_dr ] = ecef2latlonapprox( ecef, 1 );
latlonapproxp = latlonapproxp * 180/pi;

assert(abs(latlonapproxp(1) - latlonalt(1)) < 1e-4, 'Bad latitude')
assert(abs(latlonapproxp(2) - latlonalt(2)) < 1e-4, 'Bad longitude')


% Numerical derivative
dx = 1e1;
D_num = zeros(2,3);
for i = 1:3
    xp = ecef;
    xm = ecef;
    xp(i) = xp(i) + dx;
    xm(i) = xm(i) - dx;
    [ latlonapproxp ] = ecef2latlonapprox( xp, 0 );
    [ latlonapproxm ] = ecef2latlonapprox( xm, 0 );
    D_num(:,i) = (latlonapproxp - latlonapproxm) / (2*dx);
end

% Assertions
p = 2e-6;
D_diff = dlatlonapp_dr - D_num;
Dcheck = abs(D_diff) <= p;
assert(prod(prod(Dcheck)) == 1, 'Bad Derivative');


%% Output
fprintf('PASSED: Test_ecef2latlon\n');