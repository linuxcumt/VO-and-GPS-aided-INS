%% Filter_Test.m
% 
% 
% This script tests the dynamics model and measurement model of the INS
% algorithm I am working on for the Machine Learning project.
% 
%   x(k+1) = F*x(k) + Gamma*v(k)        - State Transition
%   z(k+1) = H*x(k+1) + w(k+1)          - Measurement
% 
% In this example, a vehicle will be driving on a road in two dimensions.
% The vehicle state will be modeled as the x1 and x2 position, velocity,
% and acceleration components.
% 
% The point of this script is to show what it means that the Process Noise
% Covariance (Q) needs to be tuned. It's not like the Measurement Noise
% Covariance (R) that can be measured by taking the standard deviation of
% sensor errors. There are rules of thumb to guess Q, but beyond that
% there's no theoretical way to do it. I think that it would probably be a
% good problem for machine learning.
% 
% DEPENDENCIES
% zhist.mat
% xhist.mat
% 
% @author: Matt Marti
% @date: 2018-12-16

clear, clc, clear global
global dt


%% Governing System Parameters
% Refer to the notebook for the models

% Measurement Delta Time
dt = .01;
tof = 100;
N = round(tof / dt);

% Initial state
x0 = [ 2; .3; -0.005; -1; -.2; .002; -3; .3; 1; zeros(6,1); ...
    0.001; 0.002; -0.001; 0.0001; 0.0005; -0.0002];

% Parameters
nx = 21;
nv = 12;
nz = 6;

% F - Dynamics Propagation function
f = @(k, x, u, v) dynamics_model(k, x, v);

% H - Measurement Model
h = @(k, x, w) measurement_model(k, x, w);

% R - Measurement Noise Covariance
R = 1e-3*eye(6);

% Q - Process Noise Covariance
Q = zeros(nv,nv);
Q(1,1) = .01; % Accel
Q(2,2) = .01; % Accel
Q(3,3) = .01; % Accel
Q(4,4) = .02; % Gyro
Q(5,5) = .02; % Gyro
Q(6,6) = .02; % Gyro
Q(7,7) = 0.00001; % Accel bias
Q(8,8) = 0.00001; % Accel bias
Q(9,9) = 0.00001; % Accel bias
Q(10,10) = 0.00001; % Gyro bias
Q(11,11) = 0.00001; % Gyro bias
Q(12,12) = 0.00001; % Gyro bias


%% Generate Data

% Generate True Measurment Noise
Sr = chol(R)';
whist = Sr*randn(nz,N);

% Generate True Process Noise
Sq = chol(Q)';
v0 = Sq*randn(nv,1);
vhist = Sq*randn(nv,N);

% Preallocate data arrays
xhist = zeros(nx,N);
zhist = zeros(nz,N);

% Generate data
xk = x0;
vk = v0;
for kp1 = 1:N
    k = kp1 - 1;
    
    % Current state data
    wkp1 = whist(:,kp1);
    
    % Compute process
    xkp1 = f(k, xk, 0, vk);
    zkp1 = h(kp1, xkp1, wkp1);
    
    % Save data
    xhist(:,kp1) = xkp1;
    zhist(:,kp1) = zkp1;
    
    % Iterate state data
    xk = xkp1;
    vk = vhist(:,kp1);
end


%% Run Filter on System

% Load measurement time history
N = size(zhist,2);

% Initial state estimate - assume we are starting from the first GPS
% measurement
xhat0 = x0;

% Initial Covariance estimate - This can be set pretty big and it will
% converge fast for a linear system. Don't want it to be too small if you
% think you're not accurate
P0 = 5*eye(nx);

% Preallocate data arrays
xhathist = zeros(nx,N);
Phist = zeros(nx,nx,N);
epsilonhist = zeros(1,N);

% Initialize loop variables
xhatk = xhat0;
Pk = P0;
Qk = Q;
Rkp1 = R;

% Run Filter
for kp1 = 1:N % Index by k+1
    k = kp1 - 1;
    
    % Obtain Measurement Data
    zkp1 = zhist(:,kp1);
    
    % State propagation
    [ xbarkp1, Fk, Gammak ] = f(k, xhatk, 0, zeros(nv,1));
    [ zbarkp1, Hkp1 ] = h(kp1, xbarkp1, 0);
    
    % Dynamic propagation of covariance
    Pbarkp1 = Fk*Pk*(Fk') + Gammak*Qk*(Gammak');

    % Kalman Gain calculation
    nukp1 = zkp1 - zbarkp1;
    Skp1 = Hkp1*Pbarkp1*(Hkp1') + Rkp1;
    invSkp1 = inv(Skp1);
    Wkp1 = Pbarkp1*(Hkp1')*invSkp1; %#ok

    % Innovation Statistic
    epsilonkp1 = (nukp1')*invSkp1*nukp1; %#ok

    % Measurement update of state covariance
    xhatkp1 = xbarkp1 + Wkp1*nukp1;
    Pkp1 = Pbarkp1 - Wkp1*Skp1*(Wkp1');

    % Save data
    xhathist(:,kp1) = xhatkp1;
    Phist(:,:,kp1) = Pkp1;
    epsilonhist(kp1) = epsilonkp1;
    
    % Iterate
    xhatk = xhatkp1;
    Pk = Pkp1;
end


%% Chi-Squared Distribution Test for Filter Consistency
% This is what's known as a consistency test. It makes sure that the Filter
% Innovation is sampled from a Chi-Squared distribution of degree equal to 
% the number of elements of the measurement vector z. Essentially, this is
% what you use to tell that the Filter is working correctly (In the absence
% of truth data).

% "1% of points may lie outside these bounds"
alpha = 0.01;

% Chi-Squared Distribution Bounds for Innovation Statistic
% These are displayed as red lines on the Innovation Statistic Mean Time
% History. A certain percentage of points must lie within these bounds.
nz = 4; % Number of Measurements
r1nu = chi2inv(alpha/2, nz);
r2nu = chi2inv(1-alpha/2, nz);

% Chi-Squared Distribution Bounds for Innovation Statistic Mean
% These are displayed as magenta lines on the Consistency Test Time
% Hisotry. The mean value of the Innovation Statistic must lie within these
% bounds.
r1nu_mean = chi2inv(alpha/2, N*nz)/N;
r2nu_mean = chi2inv(1-alpha/2, N*nz)/N;

% Chi-squared distribution test
passcount = zeros(N,1);
for k = 1:N
    passcount(k) = (r1nu <= epsilonhist(k)) && (epsilonhist(k) <= r2nu);
end
passrate = 100*sum(passcount)/length(passcount);
pass = passrate >= 100*(1-alpha);

% Filter consistency can also be measured by running Monte-Carlo
% simulations. However, for this example, we are assuming we don't have the
% true state time history.

% Display whether filter passes consistency test
if pass
    fprintf('Filter passed consistency test\n');
else
    fprintf('Filter failed consistency test\n');
end
fprintf('Filter pass rate: %.2f [%%]\n', passrate);
fprintf('Mean Innovation Statistic falls within bounds: %d\n', ...
    r1nu_mean <= mean(epsilonhist) && mean(epsilonhist) <= r2nu_mean);


%% Plot State Estimates

% Ground Track Plot
figure(1)
hold off
plot(xhist(1,:), xhist(4,:),'k.-', 'linewidth', 1.25, 'markersize', 10);
hold on
plot(zhist(1,:), zhist(3,:),'ro', 'markersize', 5);
plot(xhathist(1,:), xhathist(4,:),'b.-', 'linewidth', 1.25, 'markersize', 10);
title('Object Ground Track');
xlabel('x1');
ylabel('y1');
legend({'True', 'GPS','Kalman'});
grid on, grid minor

% Innovation Statistic Plot
figure(2);
hold off;
semilogy(epsilonhist', 'linewidth', 1.5);
hold on;
semilogy(r1nu*ones(size(epsilonhist)), 'r--', 'linewidth', 1.75);
semilogy(r2nu*ones(size(epsilonhist)), 'r--', 'linewidth', 1.75);
semilogy(r1nu_mean*ones(size(epsilonhist)), 'm--', 'linewidth', 1.75);
semilogy(r2nu_mean*ones(size(epsilonhist)), 'm--', 'linewidth', 1.75);
semilogy(mean(epsilonhist)*ones(size(epsilonhist)), 'b-.', ...
    'linewidth', 1);
semilogy(nz*ones(size(epsilonhist)), 'k-.', 'linewidth', 1);
hold off;
title('Innovation Statistic Consistency Test Time History');
ylabel('Innovation Statistic');
xlabel('Index of Innovation Statistic k');
grid on, grid minor;

% Acceleration Time History
figure(3);
hold off
plot(xhist([3,6,9],:)', 'b');
hold on
plot(xhathist([3,6,9],:)','r');
title('Acceleration Time History');
legend({'True', '', '', 'Kalman', '', ''});
grid on, grid minor

% Position Error
figure(4);
hold off
plot((xhathist([1,4,7],:) - xhist([1,4,7],:))');
title('Position Error Time History');
grid on, grid minor

% Velocity Error
figure(5);
hold off
plot((xhathist([2,5,8],:) - xhist([2,5,8],:))');
title('Velocity Error Time History');
grid on, grid minor

% Acceleration Error
figure(6);
hold off
plot((xhathist([3,6,9],:) - xhist([3,6,9],:))');
title('Acceleration Error Time History');
grid on, grid minor