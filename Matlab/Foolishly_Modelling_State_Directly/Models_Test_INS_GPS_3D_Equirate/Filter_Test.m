%% Filter_Test.m
% 
% This script tests the dynamics model and measurement model of the INS
% algorithm I am working on for the Machine Learning project. The INS will
% implement an extended Kalman Filter that estimates the state goverened
% by:
% 
%   x(k+1) = f( k, x(k), v(k) )                - State Transition
%   z(k+1) = h( k+1, x(k+1), w(k+1) )          - Measurement
% 
% Note that the results of these tests are kind of bad because of the
% initial guess of the state covariance matrix. I didn't put any effort
% into making an educated guess, so it's a poor guess of the covariance.
% 
% This script is for funzies. I want to see how well the INS works if the
% GPS happens at the same rate and if I fudge the Loosely coupled filter.
% I'm doing this to procrastinate implementing the one from the guy's
% thesis.
% 
% @dependencies
% dynamics_model.m    v 2019-03-01
% measurement_model.m v 2019-03-01
% 
% @author: Matt Marti
% @date: 2019-03-01

clear, clc, clear global


%% Governing System Parameters
% Refer to the notebook for the models

% Measurement Delta Time
dt = 1;
tof = 1000;
N = round(tof / dt);

% Initial state
x0 = [ 2; .3; -0.005; -1; -.2; .002; 6700e3; .3; 1; zeros(6,1); ...
    0.001; 0.002; -0.001; 0.0001; 0.0005; -0.0002];

% Parameters
nx = 21;
nv = 12;
nz_INS = 6;
nz_GPS = 6;
nz = nz_INS + nz_GPS;

% F - Dynamics Propagation function
f = @(k, x, u, v) dynamics_model( dt, x, v );

% H - Measurement Model
h = @(k, x, w) measurement_model( x, w );

% R - Measurement Noise Covariance
R_INS = 1e-3*eye(6);
R_GPS = diag([5;.1;5;.1;7;.2]);
R = [R_GPS, zeros(nz_GPS,nz_INS); zeros(nz_INS,nz_GPS), R_INS];

% Q - Process Noise Covariance
Q = zeros(nv,nv);
Q(1,1) = .01; % Accel
Q(2,2) = .01; % Accel
Q(3,3) = .01; % Accel
Q(4,4) = .002; % Gyro
Q(5,5) = .002; % Gyro
Q(6,6) = .002; % Gyro
Q(7,7) = 0.0000001; % Accel bias
Q(8,8) = 0.0000001; % Accel bias
Q(9,9) = 0.0000001; % Accel bias
Q(10,10) = 0.0000001; % Gyro bias
Q(11,11) = 0.0000001; % Gyro bias
Q(12,12) = 0.0000001; % Gyro bias


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
P0 = 1e-3*eye(nx);

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
plot(xhathist(1,:), xhathist(4,:),'b.-', 'linewidth', 1.25, 'markersize', 10);
plot(zhist(1,:), zhist(3,:),'ro', 'markersize', 5);
title('Object Ground Track');
xlabel('x1');
ylabel('y1');
legend({'True', 'Kalman', 'GPS'});
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

% Acceleration Bias Error
figure(4);
hold off
plot(xhathist(16:18,:)', 'b');
hold on
plot(xhist(16:18,:)', 'r');
title('Acceleration Bias Time History');
legend({'True', '', '', 'Kalman', '', ''});
grid on, grid minor

% Gyro Bias Error
figure(5);
hold off
plot(xhathist(19:21,:)', 'b');
hold on
plot(xhist(19:21,:)', 'r');
title('Gyro Bias Time History');
legend({'True', '', '', 'Kalman', '', ''});
grid on, grid minor

% Position Error
figure(6);
hold off
plot((xhathist([1,4,7],:) - xhist([1,4,7],:))');
title('Position Error Time History');
grid on, grid minor

% Velocity Error
figure(7);
hold off
plot((xhathist([2,5,8],:) - xhist([2,5,8],:))');
title('Velocity Error Time History');
grid on, grid minor

% Acceleration Error
figure(8);
hold off
plot((xhathist([3,6,9],:) - xhist([3,6,9],:))');
title('Acceleration Error Time History');
grid on, grid minor

% Acceleration Bias Error
figure(9);
hold off
plot(xhathist(16:18,:)' - xhist(16:18,:)');
title('Acceleration Bias Error Time History');
grid on, grid minor

% Gyro Bias Error
figure(10);
hold off
plot(xhathist(19:21,:)' - xhist(19:21,:)');
title('Gyro Bias Error Time History');
grid on, grid minor