function [ xkp1, latlonalt, attitude ] ...
    = imuMechanization( dt, xk, zk, R_b_m, deltaflag )
% Updates IMU State from k time to k+1 time
% Implements the mechanization equations for the INS state transition.
% Based on the INS Mechanization Equations in Section 2.6 of the thesis of
% Eduardo Infante, "Development and Assessment of Loosely-Coupled INS Using
% Smartphone Sensors".
% 
% Measurments are assumed to be acceleration and angular rate, not
% increments.
% 
% The position and velocity updates are different than in Infante. Instead
% of averaging the k+1 state and k state, a zero-order hold model is used.
% Essentially, the velocity and acceleration are constant on the interval
% [ k, k+1 ). This is like a "staircase" time history of velocity and 
% acceleration. 
% 
% @arg
% dt        - double
% xk        - 25 x 1 double matrix
%             INS State vector
% zk        - 6 x 1 double matrix
%             Measurement vector
%             [ ax; ay; az; wx; wy; wz ]
%                 a - Acceleration
%                 t - Angular rate
% R_b_m     - 3 x 3 double matrix
%             Direction Cosines matrix from Sensor Measurement frame to 
%             vehicle body frame. Extrinsics matrix
% deltaflag - bool
%             Flag indicates measurements z are delta angles and delta 
%             velocity if True. If False, then measurements are angular
%             rate and acceleration. Default to false
% 
% @return
% xkp1      - 25 x 1 double matrix
%             State vector at k+1 state
% latlonalt - 3 x 1 double matrix
%             IMU [ Latitude; Longitude; Altitude ]
% attitude  - 3 x 1 double matrix
%             IMU Attitude [ Heading; Pitch; Roll ]
% 
% @author: Matt Marti
% @date: 2019-03-04

% Constants
global OMEGA_EARTH
if isempty(OMEGA_EARTH)
    constants
end

% Check input
if nargin < 5
    deltaflag = 0;
end

% State variables
rk_e   = xk([1,4,7],1);    % Position
vk_e   = xk([2,5,8],1);    % Velocity
ak_e   = xk([3,6,9],1);    % Acceleration
qk_e_b = xk(10:13,1);      % Orientation
betaa  = xk([14,15,16],1); % Accel bias estimate
scalea = xk([17,18,19],1); % Accel Scale Factor error estimate
betag  = xk([20,21,22],1); % Accel bias estimate
scaleg = xk([23,24,25],1); % Accel Scale Factor error estimate

% Measurment variables
dvf_m_raw = zk(1:3); % Raw Specific Force / Acceleration
theta_ib_m_raw = zk(4:6); % Raw Angular Rate
if deltaflag
    dvf_m_raw = dvf_m_raw/dt;
    theta_ib_m_raw = theta_ib_m_raw/dt;
end
dvf_b_raw = R_b_m*dvf_m_raw;
theta_ib_b_raw = R_b_m*theta_ib_m_raw;

% 1. Remove bias from raw measurements
dvf_b = (dvf_b_raw*dt - betaa*dt) ./ (1 + scalea); % (Velocity increment)
dtheta_b_ib = (theta_ib_b_raw*dt - betag*dt) ./ (1 + scaleg); % (Angle increment)

% 2. Remove Earth rotation from gyro measurements
R_e_b = quat2dircos(qk_e_b); % Body to ECEF rotation matrix
R_b_e = R_e_b';
omega_e_ie = [ 0; 0; OMEGA_EARTH ];
theta_b_ie = (R_b_e*omega_e_ie)*dt; % (Angle Increment)
theta_b_eb = dtheta_b_ib - theta_b_ie;

% 3. Quaternion update
qkp1_e_b = increQuatWithAngles( qk_e_b, theta_b_eb );

% 5. Rotate acceleration from Body from to ECEF frame
S_b = skewsym(theta_b_eb);
dvf_e = R_b_e * (eye(3) + 0.5*S_b) * dvf_b;

% 6. Coriolis Effect
Omega_ie_e = skewsym(omega_e_ie);
cf_e = 2 * Omega_ie_e * vk_e;

% 7. Gravity Model and Centrifugal acceleration
gf_e = gravitymodel(rk_e);
wf_e = centrifugalmodel(rk_e);

% 8. Apply coriolis and gravity models to acceleration
dvkp1_e = dvf_e - dt * (gf_e + wf_e - cf_e);
akp1_e = dvkp1_e / dt;

% 9. & 10. Velocity and Position update
vkp1_e = vk_e + ak_e*dt;
rkp1_e = rk_e + vk_e*dt + 0.5*ak_e*dt^2;

% 11. Convert position to Latitude and Longitude
% 12. Local Frame to ECEF Frame Rotation
% 13. Attitude
[ latlonalt, attitude ] = rq2attitude( rkp1_e, qkp1_e_b );

% Update State variables
xkp1 = zeros(25,1);
xkp1([1,4,7],1) = rkp1_e;      % Position
xkp1([2,5,8],1) = vkp1_e;      % Velocity
xkp1([3,6,9],1) = akp1_e;      % Delta Velocity
xkp1(10:13,1)   = qkp1_e_b;    % Orientation
xkp1(14:25,1)   = xk(14:25,1); % Bias and scale factors

end

