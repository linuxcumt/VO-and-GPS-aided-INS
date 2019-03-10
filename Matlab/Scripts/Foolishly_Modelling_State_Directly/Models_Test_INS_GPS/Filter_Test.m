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
% This script implements a loosely coupled GPS aided INS. The GPS sample
% rate is lower than that of the INS, so using the observations to estimate
% the state is a little tricky. GPS provides position measurements.
% 
% Note that the results of these tests are kind of bad because of the
% initial guess of the state covariance matrix. I didn't put any effort
% into making an educated guess, so it's a poor guess of the covariance.
% 
% @dependencies
% dynamics_model.m        v 2019-03-01
% measurement_model_INS.m v 2019-03-01
% measurement_model_GPS.m v 2019-03-02
% 
% @author: Matt Marti
% @date: 2019-03-08

clear, clc, clear global


%% Governing System Parameters
% Refer to the notebook for the models

% Measurement Delta Time
dt_INS = .04;
dt_GPS = 1;
tof = 200;
Nx = round(tof / dt_INS);
Nz_INS = Nx;
Nz_GPS = round(tof / dt_GPS);

% Initial state
x0 = [ 2; .3; -0.005; -1; -.2; .002; 6700e3; .3; 1; zeros(6,1); ...
    .7; 0.2; -0.3; 0.1; 0.5; -0.2];

% Parameters
nx = 21;
nv = 12;
nz_INS = 6;
nz_GPS = 6;

% F - Dynamics Propagation function
f = @(k, x, u, v) dynamics_model( dt_INS, x, v );

% H - Measurement Model
h_INS = @(k, x, w) measurement_model_INS( x, w );
h_GPS = @(k, x, w) measurement_model_GPS( x, w );

% R - Measurement Noise Covariance
R_INS = 1e-3*eye(6);
R_GPS = diag([5;.1;5;.1;7;.2]);

% Q - Process Noise Covariance
Q = zeros(nv,nv);
Q(1,1) = 1e-4; % Accel
Q(2,2) = 1e-4; % Accel
Q(3,3) = 1e-4; % Accel
Q(4,4) = 1e-4; % Gyro
Q(5,5) = 1e-4; % Gyro
Q(6,6) = 1e-4; % Gyro
Q(7,7) = 1e-10/dt_INS; % Accel bias
Q(8,8) = 1e-10/dt_INS; % Accel bias
Q(9,9) = 1e-10/dt_INS; % Accel bias
Q(10,10) = 1e-10/dt_INS; % Gyro bias
Q(11,11) = 1e-10/dt_INS; % Gyro bias
Q(12,12) = 1e-10/dt_INS; % Gyro bias


%% Generate Data

% Generate INS True Measurment Noise
Sr_INS = chol(R_INS)';
whist_INS = Sr_INS*randn(nz_INS,Nz_INS);

% Generate INS True Measurment Noise
Sr_GPS = chol(R_GPS)';
whist_GPS = Sr_GPS*randn(nz_GPS,Nz_GPS);

% Generate True Process Noise
Sq = chol(Q)';
v0 = Sq*randn(nv,1);
vhist = Sq*randn(nv,Nx);

% Preallocate data arrays
xhist = zeros(nx,Nx);
zhist_INS = zeros(nz_INS,Nz_INS);
zhist_GPS = zeros(nz_GPS,Nz_GPS);

% Generate data
xk = x0;
vk = v0;
kp1gps = 1;
for kp1 = 1:Nx
    k = kp1 - 1;
    
    % State Propagation
    xkp1 = f(k, xk, 0, vk);
    xhist(:,kp1) = xkp1;
    
    % INS Observation
    wkp1_INS = whist_INS(:,kp1);
    zhist_INS(:,kp1) = h_INS(kp1, xkp1, wkp1_INS);
    
    % GPS Observation
    if ~mod(kp1-1,round(dt_GPS/dt_INS))
        wkp1_GPS = whist_GPS(:,kp1gps);
        zhist_GPS(:,kp1gps) = h_GPS(kp1gps, xkp1, wkp1_GPS);
        kp1gps = kp1gps + 1;
    end
    
    % Iterate state data
    xk = xkp1;
    vk = vhist(:,kp1);
end


%% Run Filter on System

% Initial state estimate - assume we are starting from the first GPS
% measurement
xhatk = x0;

% Initial Covariance estimate - This can be set pretty big and it will
% converge fast for a linear system. Don't want it to be too small if you
% think you're not accurate
Pk = 1e0*eye(nx);

% Preallocate data arrays
xhathist = zeros(nx,Nx);
Phist = zeros(nx,nx,Nx);
epsilonhist = zeros(1,Nx);

% Initialize loop variables
Qk = Q;

% Run Filter
kp1gps = 1;
progressbar = waitbar(0, 'Progress');
for kp1 = 1:Nx % Index by k+1
    k = kp1 - 1;
    gpsflag = ~mod(kp1-1,round(dt_GPS/dt_INS));
    
    % State propagation
    [ xbarkp1, Fk, Gammak ] = f(k, xhatk, 0, zeros(nv,1));
    
    % Inertial measurments
    zkp1 = zhist_INS(:,kp1);
    [ zbarkp1, Hkp1 ] = h_INS(kp1, xbarkp1, 0);
    Rkp1 = R_INS;
    
    % Obtain Measurement Data
    if gpsflag
        
        % GPS measurments
        zkp1_GPS = zhist_GPS(:,kp1gps);
        [ zbarkp1_GPS, Hkp1_GPS ] = h_GPS(kp1gps, xbarkp1, 0);
        kp1gps = kp1gps + 1;
        
        % Combine measurements
        zkp1 = [ zkp1_GPS; zkp1 ]; %#ok
        zbarkp1 = [ zbarkp1_GPS; zbarkp1 ]; %#ok
        Hkp1 = [ Hkp1_GPS; Hkp1 ]; %#ok
        Rkp1 = [ R_GPS,                zeros(nz_GPS,nz_INS) ; ...
                 zeros(nz_INS,nz_GPS), Rkp1                 ]; %#ok
    end
    
    % Dynamic propagation of covariance
    Pbarkp1_INS = Fk*Pk*(Fk') + Gammak*Qk*(Gammak');

    % Kalman Gain calculation
    nukp1 = zkp1 - zbarkp1;
    Skp1 = Hkp1*Pbarkp1_INS*(Hkp1') + Rkp1;
    invSkp1 = inv(Skp1);
    Wkp1 = Pbarkp1_INS*(Hkp1')*invSkp1; %#ok

    % Innovation Statistic
    epsilonkp1 = (nukp1')*invSkp1*nukp1; %#ok

    % Measurement update of state covariance
    xhatkp1 = xbarkp1 + Wkp1*nukp1;
    Pkp1 = Pbarkp1_INS - Wkp1*Skp1*(Wkp1');
    
    % Save data
    xhathist(:,kp1) = xhatkp1;
    Phist(:,:,kp1) = Pkp1;
    epsilonhist(kp1) = epsilonkp1;
    
    % Iterate
    xhatk = xhatkp1;
    Pk = Pkp1;
    if ~mod(kp1/Nx*100-1,1)
        waitbar(kp1/Nx, progressbar);
    end
end
close(progressbar);


%% Chi-Squared Distribution Test for Filter Consistency
% This is what's known as a consistency test. It makes sure that the Filter
% Innovation is sampled from a Chi-Squared distribution of degree equal to 
% the number of elements of the measurement vector z. Essentially, this is
% what you use to tell that the Filter is working correctly (In the absence
% of truth data).

% "1% of points may lie outside these bounds"
alpha = 0.01;

% Weighted average order of chi-squared distribution
nz = (nz_INS/dt_INS+(nz_GPS+nz_INS)/dt_GPS)/(1/dt_INS+1/dt_GPS);

% Chi-Squared Distribution Bounds for Innovation Statistic
% These are displayed as red lines on the Innovation Statistic Mean Time
% History. A certain percentage of points must lie within these bounds.
r1nu = chi2inv(alpha/2, nz);
r2nu = chi2inv(1-alpha/2, nz);

% Chi-Squared Distribution Bounds for Innovation Statistic Mean
% These are displayed as magenta lines on the Consistency Test Time
% Hisotry. The mean value of the Innovation Statistic must lie within these
% bounds.
r1nu_mean = chi2inv(alpha/2, Nx*nz)/Nx;
r2nu_mean = chi2inv(1-alpha/2, Nx*nz)/Nx;

% Chi-squared distribution test
passcount = zeros(Nx,1);
for k = 1:Nx
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
plot(zhist_GPS(1,:), zhist_GPS(3,:),'ro', 'markersize', 5);
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

% Acceleration Bias
figure(4);
hold off
plot(xhist(16:18,:)', 'r');
hold on
plot(xhathist(16:18,:)', 'b');
title('Acceleration Bias Time History');
legend({'True', '', '', 'Kalman', '', ''});
grid on, grid minor

% Gyro Bias
figure(5);
hold off
plot(xhist(19:21,:)', 'r');
hold on
plot(xhathist(19:21,:)', 'b');
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

% Covariance condition number
condhist = zeros(size(Phist,3),1);
for i = 1:length(condhist)
    condhist(i) = cond(squeeze(Phist(:,:,i)));
end
figure(11);
hold off
semilogy(condhist');
title('Covariance Condition Number Time History');
grid on, grid minor