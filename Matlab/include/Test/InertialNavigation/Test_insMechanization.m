%% Test_imuMechanization
% 
% Ensures that imuMechanization.m works the way I expect it to. The first
% couple tests are in lieu of actual INS data

clear, clear global


%% Set global constants to 0 to simulate inertial reference frame

% Set gravity and Earth rottaion to 0
constants;
global GRAVCONST MASSEARTH MUEARTH OMEGA_EARTH
GRAVCONST = 0;
MASSEARTH = 0;
MUEARTH = 0;
OMEGA_EARTH = 0;

% Double check gravity and coriolis effect are now non-existent
r = 6400e3*[ 0.728; 0.729; .3];
g = gravitymodel(r);
ac = centrifugalmodel(r);
assert(norm(g) == 0, 'Gravity force not 0');
assert(norm(ac) == 0, 'Centrifugal force not 0');


%% Test 1: Stationary INS doesn't move or rotate

% Parameters
deltaflag = 0;
dt = 1e-2;
R_b_m = eye(3);

% Initialize state
x0 = zeros(25,1);
x0([1,4,7],1) = [1; 2; 3]; % Arbitrary position vector
x0([10,11,12,13],1) = [1; 0; 0; 0]; % Identity rotation matrix

% Initialize measurements
N = 100;
zhist = zeros(6,N);

% Run simulation
xhist = 1324*ones(25,N);
xk = x0;
for i = 1:N
    zk = zhist(:,i);
    [xkp1, latlonalt, attitude] ...
        = imuMechanization(dt, xk, zk, R_b_m, deltaflag);
    xhist(:,i) = xkp1;
    xk = xkp1;
end

% Assert state does not change for stationary massless Earth
p = 1e-9;
for i = 1:N
    for j = 1:25
        assert(abs(xhist(j,i) - x0(j)) < p, 'Bad state');
    end
end


%% Test 2: Constant Acceleration on x3 axis
% Acceleration on x3 axis when orientation is straight up and down, and
% position is x1 = 1; Basically, IMU goes up

% Parameters
deltaflag = 0;
dt = 1e-2; % [s]
R_b_m = eye(3);

% Initialize state
x0 = zeros(25,1);
x0([1,4,7],1) = [100; 0; 0]; % Arbitrary position vector
x0([2,5,8],1) = [2; 0; 0]; % Arbitrary velocity vector
q = [1; 0; 0; 0]; % Set heading, pitch, roll to 0
q = q/norm(q); % Normalize quaternion
x0([10,11,12,13],1) = q; % Attitude

% Acceleartion
T = 5; % [s]
acc = 3;
deltav = x0(3) + T * acc;
deltax = x0(2)*T + 0.5*acc*T^2;
xf = x0(1) + deltax;
vf = x0(2) + deltav;

% Initialize measurements
N = T / dt;
zhist = zeros(6,N);
zhist(1,:) = acc; % IMU accelerates up

% Run simulation
xhist = 1324*ones(25,N);
xk = x0;
for i = 1:N
    zk = zhist(:,i);
    [xkp1, latlonalt, attitude] ...
        = imuMechanization(dt, xk, zk, R_b_m, deltaflag);
    xhist(:,i) = xkp1;
    xk = xkp1;
end

% Assert state does not change for all but x axis parameters
p = 1e-9;
for i = 1:N
    for j = 4:25
        assert(abs(xhist(j,i) - x0(j)) < p, 'Bad state');
    end
end

% Position and velocity check
assert(abs(xhist(1,end) - xf) < 2e-1, 'Bad final position');
assert(abs(xhist(2,end) - vf) < 1e-1, 'Bad final velocity');

% Attitude in this case should be [ 90; 0; 90 ]
assert(abs(attitude(1) - 90) < 1e0, 'Bad attitude');
assert(abs(attitude(2) - 0) < 1e-3, 'Bad attitude');
assert(abs(attitude(3) - -90) < 1e-3, 'Bad attitude');


%% Test 2: Constant Rotation on theta_z
% Rotate about the body-frame y axis for T seconds to get roll the device
% to zero roll

% Parameters
deltaflag = 0;
dt = 1e-3;
R_b_m = eye(3);
T = 1.5;

% Initialize state
x0 = zeros(25,1);
x0([1,4,7],1) = 6400e3*[1; 0; 0]; % Arbitrary position vector
q = [1; 0; 0; 0]; % Set heading, pitch, roll to 0
q = q/norm(q); % Normalize quaternion
x0([10,11,12,13],1) = q; % Attitude

% Rotation Rate:
d1 = -90;
d2 = 0;
deg = d2 - d1;
dtheta = deg * pi/180;
dtheta_dt = dtheta / T;
z = dtheta_dt;

% Initialize measurements
N = T/dt;
zhist = zeros(6,N);
zhist(5,:) = z; % Apply rotation to y axis gyro

% Run simulation
xhist = 1324*ones(25,N);
attitudehist = zeros(3,N);
xk = x0;
for i = 1:N
    zk = zhist(:,i);
    [xkp1, latlonalt, attitude] ...
        = imuMechanization(dt, xk, zk, R_b_m, deltaflag);
    xhist(:,i) = xkp1;
    attitudehist(:,i) = attitude;
    xk = xkp1;
end
x0_90_0_0 = xhist(:,end);

% Assert state does not change for stationary massless Earth
p = 1e-9;
for i = 1:N
    for j = [1:9, 14:25]
        assert(abs(xhist(j,i) - x0(j)) < p, 'Bad state');
    end
end

% Check attitude
assert(abs(attitude(1) - 90) < 1e0, 'Bad heading');
assert(abs(attitude(2) - 0) < 1e-3, 'Bad pitch');
assert(abs(attitude(3) - 0) < 1e-1, 'Bad roll');

% Roll another 90 degrees to +90
for i = 1:N
    zk = zhist(:,i);
    [xkp1, latlonalt, attitude] ...
        = imuMechanization(dt, xk, zk, R_b_m, deltaflag);
    xhist(:,i) = xkp1;
    attitudehist(:,i) = attitude;
    xk = xkp1;
end
assert(abs(attitude(1) - 90) < 1e0, 'Bad heading');
assert(abs(attitude(2) - 0) < 1e-3, 'Bad pitch');
assert(abs(attitude(3) - 90) < 1e-1, 'Bad roll');

% Roll another 90 degrees to +180
for i = 1:N
    zk = zhist(:,i);
    [xkp1, latlonalt, attitude] ...
        = imuMechanization(dt, xk, zk, R_b_m, deltaflag);
    xhist(:,i) = xkp1;
    attitudehist(:,i) = attitude;
    xk = xkp1;
end
assert(abs(attitude(1) - 90) < 1e0, 'Bad heading');
assert(abs(attitude(2) - 0) < 1e-3, 'Bad pitch');
assert(abs(attitude(3) - 180) < 1e-1, 'Bad roll');


%% Test 4
% Rotate about x axis to achieve pitch angle of 45 degrees after last
% manuever

% Parameters
deltaflag = 0;
dt = 1e-3;
R_b_m = eye(3);
T = 0.5;

% Initialize state
x0 = x0_90_0_0;

% Rotation Rate:
d1 = 0;
d2 = 45;
deg = d2 - d1;
dtheta = deg * pi/180;
dtheta_dt = dtheta / T;
z = dtheta_dt;

% Initialize measurements
N = T/dt;
zhist = zeros(6,N);
zhist(4,:) = z; % Rotate about x axis

% Run simulation
xhist = 1324*ones(25,N);
attitudehist = zeros(3,N);
xk = x0;
for i = 1:N
    zk = zhist(:,i);
    [xkp1, latlonalt, attitude] ...
        = imuMechanization(dt, xk, zk, R_b_m, deltaflag);
    xhist(:,i) = xkp1;
    attitudehist(:,i) = attitude;
    xk = xkp1;
end

% Assert state does not change for stationary massless Earth
p = 1e-9;
for i = 1:N
    for j = [1:9, 14:25]
        assert(abs(xhist(j,i) - x0(j)) < p, 'Bad state');
    end
end

% Check attitude
assert(abs(attitude(1) - 90) < 1e0, 'Bad heading');
assert(abs(attitude(2) - d2) < 1e-1, 'Bad pitch');
assert(abs(attitude(3) - 0) < 1e-1, 'Bad roll');


%% Test 5
% Rotate about z axis to achieve heading angle of 45 degrees

% Parameters
deltaflag = 0;
dt = 1e-4;
R_b_m = eye(3);
T = 0.5;

% Initialize state
x0 = x0_90_0_0;

% Rotation Rate:
d1 = 0; % Corrosponds to 90 degrees heading
d2 = 90; % Corrosponds to 0 heading
deg = - (d2 - d1); % Heading is backwards, want ccw motion
dtheta = deg * pi/180;
dtheta_dt = dtheta / T;
z = dtheta_dt;

% Initialize measurements
N = T/dt;
zhist = zeros(6,N);
zhist(6,:) = z; % Rotate on z axis

% Run simulation
xhist = 1324*ones(25,N);
attitudehist = zeros(3,N);
xk = x0;
for i = 1:N
    zk = zhist(:,i);
    [xkp1, latlonalt, attitude] ...
        = imuMechanization(dt, xk, zk, R_b_m, deltaflag);
    xhist(:,i) = xkp1;
    attitudehist(:,i) = attitude;
    xk = xkp1;
end

% Assert state does not change for stationary massless Earth
p = 1e-9;
for i = 1:N
    for j = [1:9, 14:25]
        assert(abs(xhist(j,i) - x0(j)) < p, 'Bad state');
    end
end

% Check attitude
assert((abs(attitude(1) - 0) < 1e-1) || (abs(attitude(1) - 360) < 1e-1), ...
    'Bad heading'); % Annoying wraparound
assert(abs(attitude(2) - 0) < 1e-1, 'Bad pitch');
assert(abs(attitude(3) - 0) < 1e-1, 'Bad roll');

% Rotate to 270
for i = 1:N
    zk = zhist(:,i);
    [xkp1, latlonalt, attitude] ...
        = imuMechanization(dt, xk, zk, R_b_m, deltaflag);
    xhist(:,i) = xkp1;
    attitudehist(:,i) = attitude;
    xk = xkp1;
end
assert(abs(attitude(1) - 270) < 1e-1, 'Bad heading');
assert(abs(attitude(2) - 0) < 1e-1, 'Bad pitch');
assert(abs(attitude(3) - 0) < 1e-1, 'Bad roll');

% Rotate to 180
for i = 1:N
    zk = zhist(:,i);
    [xkp1, latlonalt, attitude] ...
        = imuMechanization(dt, xk, zk, R_b_m, deltaflag);
    xhist(:,i) = xkp1;
    attitudehist(:,i) = attitude;
    xk = xkp1;
end
assert(abs(attitude(1) - 180) < 1e-1, 'Bad heading');
assert(abs(attitude(2) - 0) < 1e-1, 'Bad pitch');
assert(abs(attitude(3) - 0) < 1e-1, 'Bad roll');


%% Output
clear global
fprintf('PASSED: Test_imuMechanization\n');
