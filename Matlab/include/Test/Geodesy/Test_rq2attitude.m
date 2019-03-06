%% Test_rq2attitude
%
% Tests the computation of attitude and latitude and longitude
%
% @author: Matt Marti
% @date: 2019-03-06

clear;


%% Test 1

% Data and function call
r = 6400e3*[ 1; 0; 0 ];
q = [ 1; 0; 0; 0 ];
q = q/norm(q);
[ latlonalt, attitude ] = rq2attitude( r, q );

% Don't need to check latitude and longitude too much, there's already a
% test
assert(latlonalt(1) == 0, 'Bad latitude')
assert(latlonalt(2) == 0, 'Bad longitude')

% Attitude
assert(attitude(1) == 90, 'Heading');
assert(attitude(2) == 0, 'Pitch');
assert(attitude(3) == -90, 'Roll');


%% Test 2

% Data and function call
r = 6400e3*[ 1; 0; 1 ];
q = [ 1; 0; 0; 0 ];
q = q/norm(q);
[ lla, attitude ] = rq2attitude( r, q ); %#ok

% Attitude
assert(attitude(1) == 90, 'Heading'); % Heading is undefined
assert(attitude(2) == 0, 'Pitch');
assert(abs(attitude(3) - -45) < 3e-1, 'Roll');


%% Test 3

% Data and function call
r = [ -10000e3; 0; 0 ];
q = [ 1; 0; 0; 0 ];
q = q/norm(q);
[ lla, attitude ] = rq2attitude( r, q );

% Attitude
assert(attitude(1) == 270, 'Heading');
assert(attitude(2) == 0, 'Pitch');
assert(attitude(3) == 90, 'Roll');


%% Output
fprintf('PASSED: Test_rq2attitude\n');