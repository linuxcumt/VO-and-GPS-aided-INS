%% Test_latlon2ecef
%
% Tests the conversion functions to go from the ECEF coordinates to 
% Latitude and Longitude.
%
% @author: Matt Marti
% @date: 2018-08-05

clear;


%% Test 1 - nominal test
latlonalt = [39.1489974;-77.618957;138.280];
ecef = latlon2ecef(latlonalt);

assert(size(ecef,1) == 3, 'Incorrect size of output');
assert(size(ecef,2) == 1, 'Incorrect size of output');
assert(ecef(1) == 1061985.2311866253, 'Incorrect ECEF value');
assert(ecef(2) == -4837817.3795761978, 'Incorrect ECEF value');
assert(ecef(3) == 4005245.7246151697, 'Incorrect ECEF value');


%% Test 2 - transposed input
latlonalt = [39.1489974;-77.618957;138.280];
ecef = latlon2ecef(latlonalt');

assert(size(ecef,1) == 1, 'Incorrect size of output');
assert(size(ecef,2) == 3, 'Incorrect size of output');
assert(ecef(1) == 1061985.2311866253, 'Incorrect ECEF value');
assert(ecef(2) == -4837817.3795761978, 'Incorrect ECEF value');
assert(ecef(3) == 4005245.7246151697, 'Incorrect ECEF value');


%% Test 3 - multiple data
latlonalt = [39.1489974, 27;...
    -77.618957, 89;...
    138.280, 600];
ecef = latlon2ecef(latlonalt);

assert(size(ecef,1) == 3, 'Incorrect size of output');
assert(size(ecef,2) == 2, 'Incorrect size of output');
assert(ecef(1,1) == 1061985.2311866253, 'Incorrect ECEF value');
assert(ecef(2,1) == -4837817.3795761978, 'Incorrect ECEF value');
assert(ecef(3,1) == 4005245.7246151697, 'Incorrect ECEF value');
assert(ecef(1,2) == 99259.181292558031, 'Incorrect ECEF value');
assert(ecef(2,2) == 5686554.6877512438, 'Incorrect ECEF value');
assert(ecef(3,2) == 2878487.9708811892, 'Incorrect ECEF value');


%% Test 4 - data from school
latlonalt = [-80.42; 37.23; 630];
ecef = latlon2ecef(latlonalt);

assert(abs(ecef(1,1) -  0.8480108977792e6) < 1e-7, 'Incorrect ECEF value');
assert(abs(ecef(2,1) - 0.64437541898631e6) < 1e-7, 'Incorrect ECEF value');
assert(abs(ecef(3,1) -  -6.26813853517580e6) < 1e-7, 'Incorrect ECEF value');


%% Test 4 - bad input
latlonalt = [39.1489974;-77.618957;138.280; 0];
try
    latlon2ecef(latlonalt);
    error('No exception thrown');
catch err
    str = 'Incorrect size of input';
    assert(strcmp(err.message,str), 'Incorrect Exception thrown');
end


%% Test 5 - test that function's inverse works correctly
latlonalt = [39.1489974;-77.618957;138.280];
ecef = latlon2ecef(latlonalt);
lla = ecef2latlon(ecef);

assert(abs(latlonalt(1) - lla(1)) < 1e-9, 'Incorrect value');
assert(abs(latlonalt(2) - lla(2)) < 1e-9, 'Incorrect value');
assert(abs(latlonalt(3) - lla(3)) < 1e-9, 'Incorrect value');


%% Output
fprintf('PASSED: Test_latlon2ecef\n');