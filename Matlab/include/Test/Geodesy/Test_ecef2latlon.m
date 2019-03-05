%% Test_ecef2latlon
%
% Tests the conversion functions to go from Latitude and Longitude to 
% ECEF coordinates.
%
% @author: Matt Marti
% @date: 2018-08-15

clear;


%% Test 1 - nominal test
ecef = [1061985.2311866253;-4837817.3795761978;4005245.7246151697];
latlonalt = ecef2latlon(ecef);

assert(size(latlonalt,1) == 3, 'Incorrect size of output');
assert(size(latlonalt,2) == 1, 'Incorrect size of output');
assert(abs(latlonalt(1) - 39.1489974) < 1e-12, 'Incorrect value');
assert(abs(latlonalt(2) - -77.618957) < 1e-12, 'Incorrect value');
assert(abs(latlonalt(3) - 138.280) < 1e-8, 'Incorrect value');


%% Test 2 - transposed input
ecef = [1061985.2311866253;-4837817.3795761978;4005245.7246151697];
latlonalt = ecef2latlon(ecef');

assert(size(latlonalt,1) == 1, 'Incorrect size of output');
assert(size(latlonalt,2) == 3, 'Incorrect size of output');
assert(abs(latlonalt(1) - 39.1489974) < 1e-12, 'Incorrect value');
assert(abs(latlonalt(2) - -77.618957) < 1e-12, 'Incorrect value');
assert(abs(latlonalt(3) - 138.280) < 1e-8, 'Incorrect value');


%% Test 3 - multiple data
ecef = [1061985.2311866253, 99259.181292558031;...
             -4837817.3795761978, 5686554.6877512438;...
             4005245.7246151697, 2878487.9708811892];
latlonalt = ecef2latlon(ecef);

assert(size(latlonalt,1) == 3, 'Incorrect size of output');
assert(size(latlonalt,2) == 2, 'Incorrect size of output');
assert(abs(latlonalt(1,1) - 39.1489974) < 1e-12, 'Incorrect value');
assert(abs(latlonalt(2,1) - -77.618957) < 1e-12, 'Incorrect value');
assert(abs(latlonalt(3,1) - 138.280) < 1e-8, 'Incorrect value');
assert(abs(latlonalt(1,2) - 27) < 1e-12, 'Incorrect value');
assert(abs(latlonalt(2,2) - 89) < 1e-12, 'Incorrect value');
assert(abs(latlonalt(3,2) - 600) < 1e-8, 'Incorrect value');


%% Test 4 - data from school
ecef = [8.4801089777923e5; 6.4437541898631e5; -6.26813853517580e6];
latlonalt = ecef2latlon(ecef);

assert(abs(latlonalt(1,1) - -80.4200) < 1e-12, 'Incorrect value');
assert(abs(latlonalt(2,1) - 37.2300) < 1e-12, 'Incorrect value');
assert(abs(latlonalt(3,1) - 630) < 1e-8, 'Incorrect value');


%% Test 4 - bad input
latlonalt = [39.1489974;-77.618957;138.280; 0];
try
    ecef2latlon(latlonalt);
    error('No exception thrown');
catch err
    str = 'Incorrect size of input';
    assert(strcmp(err.message,str), 'Incorrect Exception thrown');
end


%% Test 5 - test that function's inverse works correctly
ecef = [1061985.2311866253;-4837817.3795761978;4005245.7246151697];
latlonalt = ecef2latlon(ecef);
ecef2 = latlon2ecef(latlonalt);

assert(abs(ecef(1) - ecef2(1)) < 1e-9, 'Incorrect value');
assert(abs(ecef(2) - ecef2(2)) < 1e-9, 'Incorrect value');
assert(abs(ecef(3) - ecef2(3)) < 1e-9, 'Incorrect value');


%% Output
fprintf('PASSED: Test_ecef2latlon\n');