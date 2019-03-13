%% INS_simulation_NoRotation
% 
% Simulation of estimating a trajectory of an object using a GPS aided
% Inertial Navigation System. Based on the thesis by Eduardo Infante
% "Development and Assessment of Loosely-Coupled INS Using Smartphone
% Sensors"
% 
% Rotation and misalignment are not considered in this model so as to help
% me develop the filter.
% 
% @author: Matt Marti
% @date: 2019-03-12

clear, clc, clear global
% constants;


%% Parameters
% Refer to the notebook for the models

% Parameters
nx = 15; % Dimension of error state
ny = 12; % Dimension of IMU state
nv = 9; % Number of noise terms in error state
nvy = 3; % Number of noise terms for manuevering noise
nz_IMU = 3; % Number of measurements for IMU
nz_GPS = 6; % Number of measurements for GPS
nInt = 10; % Number of integration steps

% Measurement Times
dt_IMU = .01;.04; % IMU Sampling time period
dt_GPS = 1; % GPS Sampling time period
tof = 50; % Time Of Flight
Nx = round(tof / dt_IMU); % Number of IMU Samples
Nz_IMU = Nx; % Number of IMU observations
Nz_GPS = round(tof / dt_GPS); % Number of GPS observations

% INS Bias/Scale Factor Paremters
del_biasacc = 1e-2*(2*rand(3,1)-1); % Acc Bias Rate
init_biasacc = 1e-3*(2*rand(3,1)-1); % Acc Bias
scaleacc = 1e-7*(2*rand(3,1)-1);

% INS Bias/Scale Factor Time Constants
beta_a = 1e-2;
beta_ba = 1e-2;
beta_sa = 1e-2;

% INS Noise Standard Deviations
sigma_a = 1e-1;
sigma_ba = 1e-1;
sigma_sa = 1e-1;

% Initial INS state
y0 = zeros(ny,1);
y0(1:3,1) = 6400e3*[1; 0; 0]; % Position at lat / lon = 0,0
y0(6) = 0.25; % Northward velocity
q = [1; 1; 1; 1]; % Set heading, pitch, roll to 0
q = q/norm(q); % Normalize quaternion
% y0(7:10,1) = q; % Set Attitude

% IMU Mechanization function = [ xkp1, latlonalt, attitude, acc, gyr ]
deltaflag = 0;
R_b_m = eye(3);
g1 = @(y, z) [...
    y(1:3) + y(4:6).*dt_IMU + 0.5*(z - y(7:9))./(1 + y(10:12)).*dt_IMU.^2;...
    y(4:6) + (z - y(7:9))./(1 + y(10:12)).*dt_IMU;...
    y(7:9);...
    y(10:12)];
g2 = @(y, z) NaN;
g3 = @(y, z) [0;0;0];
g4 = @(y, z) (z - y(7:9))./(1 + y(10:12));
g5 = @(y, z) zeros(3,1);
g = @(y,z) deal(g1(y,z),g2(y,z),g3(y,z),g4(y,z),g5(y,z));

% F - Dynamics Propagation function (Error State)
f = @(k, x, v, y, zIMU) insErrorDynamicsModel_15state( ...
        dt_IMU, x, v, y, zIMU, beta_ba, 0, beta_sa, 0, nInt );

% H - Measurement Model
h = @(k, x) insErrorMeasurementModel_GNSS_rv( x, zeros(nz_GPS,1) );

% IMU dynamics covariance (noise for manuevering)
Qy = zeros(nvy,nvy);
Qy(1,1) = 1e-7; % Accel
Qy(2,2) = 1e-4; % Accel
Qy(3,3) = 1e-7; % Accel

% R - Measurement Noise Covariance
R_GPS = 1e-1*diag([5;5;7;.1;.1;.2]);
R_IMU = diag([1e-7;1e-7;1e-7]);

% Q - Process Noise Covariance
Q = zeros(nv,nv);
Q(1:3,1:3) = 2*sigma_a^2*beta_a*eye(3);
Q(4:6,4:6) = 2*sigma_ba^2*beta_ba*eye(3);
Q(7:9,7:9) = 2*sigma_sa^2*beta_sa*eye(3);


%% Generate Trajectory time history

% Generate IMU True Manuevering Noise
Sq_vy = chol(Qy)';
vyhist = Sq_vy*randn(nvy,Nx);

% Preallocate
yhist = zeros(ny,Nx);
yhist(:,1) = y0;
zyhist = zeros(nvy,Nx);
hprhist = zeros(3,Nx);
% [ ~, attitude ] = rq2attitude( y0(1:3), y0(7:10) );
% hprhist(:,1) = attitude;

% Generate data
zk = zeros(3,1);
progressbar = waitbar(0, 'Generating Trajectory');
for kp1 = 2:Nx
    k = kp1 - 1;
    
    % Partition Dynamics vector
    yk = yhist(:,k);
    zk = zk + vyhist(:,k);
%     q = yk(7:10);
    zk(1:3) = zk(1:3);% - quat2dircos(q)'*gravitymodel(yk(1:3),q); % Not in free fall
    
    % State Propagation
    [ ykp1, latlonalt, attitude, acc, gyr ] = g(yk,zk);
    
    % Partition k+1 state
    yhist(:,kp1) = ykp1;
    hprhist(:,kp1) = attitude;
    
    % Iterate
    zyhist(:,k) = zk;
    if ~mod(kp1/Nx*100-1,1), waitbar(kp1/Nx, progressbar); end
end
close(progressbar);


%% Generate Measurement Time Historys

% Generate INS True Measurment Noise
Sr_IMU = chol(R_IMU)';
whist_IMU = Sr_IMU*randn(nz_IMU,Nz_IMU);

% Generate INS True Measurment Noise
Sr_GPS = chol(R_GPS)';
whist_GPS = Sr_GPS*randn(nz_GPS,Nz_GPS);

% Generate True Process Noise (Error State)
Sq = chol(Q)';
v0 = Sq*randn(nv,1);
vhist = Sq*randn(nv,Nz_GPS);

% Generate IMU Measurement Noise
Sr_IMU = chol(R_IMU);
wyhist = Sr_IMU *randn(nz_IMU,Nz_IMU);

% Preallocate
zhist_IMU = zeros(nz_IMU,Nz_IMU);
zhist_GPS = zeros(nz_GPS,Nz_GPS);

% Initial values
xk = zeros(nx,1);
wkp1_GPS = whist_GPS(:,1);
zhist_GPS(:,1) = yhist(1:6,1) + wkp1_GPS;

% Generate data
progressbar = waitbar(0, 'Generating Measurements');
kp1gps = 1;
for kp1 = 2:Nx+1
    gpsflag = ~mod(kp1,round(dt_GPS/dt_IMU));
    k = kp1 - 1;
    
    % IMU state components
    yk = yhist(:,k);
    
    % IMU Measurement Noise
    wk = wyhist(:,k);
    
    % IMU measurements, apply parameters and noise
    acc = zyhist(1:3,k) + wk(1:3);
%     gyr = zyhist(4:6,k) + wk(4:6);
    ba = init_biasacc + k*dt_IMU*del_biasacc;
%     bg = init_biasgyr + k*dt_IMU*del_biasgyr;
    sfa = scaleacc;
%     sfg = scalegyr;
    za = (1 + sfa) .* (acc + ba);
%     zg = (1 + sfg) .* (gyr + bg);
    zhist_IMU(:,k) = za;
%     zhist_IMU(:,k) = [za; zg];
    
    % Update true bias / scale factors
%     yk(11:13) = ba;
%     yk(14:16) = sfa;
%     yk(17:19) = bg;
%     yk(20:22) = sfg;
    yk(7:9) = ba;
    yk(10:12) = sfa;
    yhist(:,k) = yk;
    
    % GPS Observation
    if gpsflag
        wkp1_GPS = whist_GPS(:,kp1gps);
        zhist_GPS(:,kp1gps) = yhist(1:6,kp1) + wkp1_GPS;
    end
    
    % Iterate
    if gpsflag, kp1gps = kp1gps + 1; end
    if ~mod(kp1/Nx*100-1,1), waitbar(kp1/Nx, progressbar); end
end
close(progressbar);
y0 = yhist(:,1);


%% Run Filter on Time History

% Dynamic parameters
bai = init_biasacc;
bad = del_biasacc;
% bgi = init_biasgyr;
% bgd = del_biasgyr;
sfa = scaleacc;
% sfg = scalegyr;

% Preallocate data arrays
xhathist = zeros(nx,Nz_GPS);
Phist = zeros(nx,nx,Nx);
epsilonhist = zeros(1,Nz_GPS);
hprhathist = zeros(size(hprhist));

% Initial state estimate - assume we are starting from the first GPS
% measurement
xhathist(:,1) = zeros(nx,1);
yhathist = zeros(size(yhist));
yhathist(:,1) = y0;
zyhathist = zyhist;

% Initial Covariance estimate - This can be set pretty big and it will
% converge fast for a linear system. Don't want it to be too small if you
% think you're not accurate
Phist(:,:,1) = 1e-3*eye(nx);

% Initialize epsilonhist
epsilonhist(1) = 0;

% Initialize loop variables
Qk = Q;

% Initialize INS Kalman Filter
kf_INS = struct(...
    'k', 0, ...
    'xhatk', xhathist(:,1), ...
    'Pk', Phist(:,:,1), ...
    'f', @(k,xk,uk,vk,yk,zkIMU) f(k, xk, vk, yk, zkIMU), ...
    'Qk', Qk, ...
    'zkp1', 0, ...
    'h', h, ...
    'Rkp1', R_GPS, ...
    'uk', 0);

% Run Filter
kp1gps = 1;
progressbar = waitbar(0, 'Running Filter');
for kp1 = 2:Nx % Index by k+1
    k = kp1 - 1;
    gpsflag = ~mod(kp1,round(dt_GPS/dt_IMU));
    
    % Obtain IMU measurements
    zkp1_IMU = zhist_IMU(:,kp1);
    zk_IMU = zhist_IMU(:,k);
    
    % Obtain GPS measurements
    if gpsflag
        zkp1_GPS = zhist_GPS(:,kp1gps);
    end
    
    % Access INS State and Dynamic parameters
    yk = yhathist(:,k);
%     yk(11:13) = bai + k*dt_IMU*bad;
%     yk(14:16) = sfa;
    yk(7:9) = bai + k*dt_IMU*bad;
    yk(10:12) = sfa;
%     yk(17:19) = bgi + k*dt_IMU*bgd;
%     yk(20:22) = sfg;
    
    % IMU State propagation
    [ ykp1, ~, ~, acck, ~ ] = g(yk, zk_IMU);
    
    % Estimate Error State
    if gpsflag
        kf_INS.zkp1 = ykp1(1:6) - zkp1_GPS;
        kf_INS.f = @(k,xk,uk,vk) f(k, xk, vk, yk, acck);
        kf_INS.xhatk = zeros(nx,1);
        [ xhatkp1, Pkp1, Wkp1, epsilonkp1, xbarkp1, Pbarkp1, Fk ] ...
                = kalmaniter_extended( kf_INS );
        
        % Partition INS Error State
        del_r_e = xhatkp1(1:3);         % - Position Error
        del_v_e = xhatkp1(4:6);         % - Velocity Error
%         err_e = xhatkp1(7:9);           % - Misalignment Error
        del_biasacc = xhatkp1(7:9);   % - Acc Bias Drift Error
%         del_biasgyr = xhatkp1(13:15);   % - Byr Bias Drift Error
        biasinitacc_e = xhatkp1(10:12); % - Initial Acc Bias Error
%         biasinitgyr_e = xhatkp1(19:21); % - Initial Gyr Bias Error
        scaleacc_e = xhatkp1(13:15);    % - Acc Scale Factor Error
%         scalegyr_e = xhatkp1(25:27);    % - Gyr Scale Factor Error
        
        % Subtract Errors from parameters
        bai = bai - biasinitacc_e;
        bad = bad - del_biasacc;
%         bgi = bgi - biasinitgyr_e;
%         bgd = bgd - del_biasgyr;
        sfa = sfa - scaleacc_e;
%         sfg = sfg - scalegyr_e;
        
        % Subtract Errors from state
        r = yk(1:3) - del_r_e;
        v = yk(4:6) - del_v_e;
%         q = increQuatWithAngles( yk(7:10), -err_e );
%         ykp1 = [ r; v; q; ykp1(11:end) ];
        ykp1 = [ r; v; ykp1(7:end) ];
    end
    
    % Update true INS state
    yhathist(:,kp1) = ykp1;
%     [ latlonalt, hpr ] = rq2attitude( ykp1(1:3), ykp1(7:10));
%     hprhathist(:,kp1) = hpr;
    
    % Save Error State data
    if gpsflag
        xhathist(:,kp1gps) = xhatkp1;
        Phist(:,:,kp1gps) = Pkp1;
        epsilonhist(kp1gps) = epsilonkp1;
    end
    
    % Iterate
    if gpsflag
        xhatk = xhatkp1;
        Pk = Pkp1;
        kp1gps = kp1gps + 1;
        kf_INS.xhatk = xhatkp1; 
        kf_INS.Pk = Pkp1;
    end
    if ~mod(kp1/Nx*100-1,1), waitbar(kp1/Nx, progressbar); end
end
close(progressbar);


%% Chi-Squared Distribution Test for Filter Consistency

% "1% of points may lie outside these bounds"
alpha = 0.01;

% Weighted average order of chi-squared distribution
nz = (nz_IMU/dt_IMU+(nz_GPS+nz_IMU)/dt_GPS)/(1/dt_IMU+1/dt_GPS);
N = length(epsilonhist);

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
passcount = zeros(Nx,1);
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
plot(yhist(2,:)-yhist(2,1), yhist(3,:)-yhist(3,1), ...
    'k.-', 'linewidth', 1.25, 'markersize', 10);
hold on
plot(yhathist(2,:)-yhist(2,1), yhathist(3,:)-yhist(3,1), ...
    'b.-', 'linewidth', 1.25, 'markersize', 10);
plot(zhist_GPS(2,:)-yhist(2,1), zhist_GPS(3,:)-yhist(3,1), ...
    'ro', 'markersize', 5);
title('Object Ground Track');
xlabel('ECEF Y');
ylabel('ECEF Z');
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

% Attitude Error
figure(3);
hold off
attdiff = hprhathist - hprhist;
for i = 1:size(attdiff,2)
    if attdiff(1,i) > 180
        attdiff(1,i) = attdiff(1,i) - 360;
    elseif attdiff(1,i) < -180
        attdiff(1,i) = attdiff(1,i) + 360;
    end
end
plot(attdiff');
hold on
title('Attitude Error Time History');
legend({'H', 'P', 'R'});
grid on, grid minor

% Acceleration Bias
figure(4);
hold off
plot(yhist(7:9,:)', 'r');
hold on
plot(yhathist(7:9,:)', 'b');
title('Acceleration Bias Time History');
legend({'True', '', '', 'Kalman', '', ''});
grid on, grid minor

% % Gyro Bias
% figure(5);
% hold off
% plot(yhist(17:19,:)', 'r');
% hold on
% plot(yhathist(17:19,:)', 'b');
% title('Gyro Bias Time History');
% legend({'True', '', '', 'Kalman', '', ''});
% grid on, grid minor

% Position Error
figure(6);
hold off
plot((yhathist(1:3,:) - yhist(1:3,:))');
title('Position Error Time History');
grid on, grid minor

% Velocity Error
figure(7);
hold off
plot((yhathist(4:6,:) - yhist(4:6,:))');
title('Velocity Error Time History');
grid on, grid minor

% % Acceleration Error
% figure(8);
% hold off
% plot((yhathist(7:9,:) - yhist(7:9,:))');
% title('Acceleration Error Time History');
% grid on, grid minor

% Acceleration Bias Error
figure(9);
hold off
plot(yhathist(7:9,:)' - yhist(7:9,:)');
title('Acceleration Bias Error Time History');
grid on, grid minor

% % Gyro Bias Error
% figure(10);
% hold off
% plot(yhathist(17:19,:)' - yhist(17:19,:)');
% title('Gyro Bias Error Time History');
% grid on, grid minor

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