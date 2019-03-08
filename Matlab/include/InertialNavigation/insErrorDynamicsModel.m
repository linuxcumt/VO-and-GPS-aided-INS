function [ xkp1, Fk, Gammak ] = insErrorDynamicsModel( dt, ins, xk, vk )
% INS Dead-Reckoning Error State Dynamics Model
% Computes the rate of change of the error state and does numerical
% integration to compute the state and evolution of the partial
% derivatives.
% 
% This models the 21 element error state of the INS, includes modeling
% the accel and gyro scale factor errors.
% 
% INPUTS
% dt     - double
%          Delta time between state updates. Recommended to be INS
%          measurement rate.
% ins    - 21 x 1 double vector
%          Current estimated state by INS
%              [ rk_e;    - Position
%                vk_e;    - Velocity
%                qk_e_b;  - Orientation quaternion
%                betaa;   - Accel bias estimate
%                scalea;  - Accel Scale Factor
%                betag;   - Accel bias estimate
%                scaleg ] - Accel Scale Factor
% xk     - 21 x 1 double vector
%          State vector at time k
%              [ del_r_e;       - Position Error
%                del_v_e;       - Velocity Error
%                err_e;         - Misalignment Error
%                del_biasacc_e; - Acc Bias Drift Error
%                del_biasgyr_e; - Byr Bias Drift Error
%                biasinitacc_e; - Initial Acc Bias Error
%                biasinitgyr_e; - Initial Gyr Bias Error
%                scaleacc_e;    - Acc Scale Factor Error
%                scalegyr_e ]   - Gyr Scale Factor Error
% vk     - 21 x 1 double vector
%          Process Noise vector at time k
%              [ eta_a;   - Acc sensor noise (PSD q_a)
%                eta_g;   - Gyr sensor noise (PSD q_g)
%                eta_ba;  - Acc bias drift noise (PSD q_ba)
%                eta_bg;  - Gyr bias drift noise (PSD q_bg)
%                eta_sa;  - Acc Scale Factor noise (PSD q_sa)
%                eta_sg ] - Gyr Scale Factor noise (PSD q_sg)
% 
% OUTPUTS
% xkp1   - 21 x 1 double vector
%          State vector at time k
% Fk     - 21 x 21 double matrix
%          State Transition partial derivative
% Gammak - 21 x 21 double matrix
%          State Noiise partial derivative
% 
% @author: Matt Marti
% @date: 2019-03-08

% Constants
global MU_EARTH OMEGA_EARTH
nx = 21;
nv = 18;
nIntegrate = 10;

% Partition INS State
rx_e = ins(1); % [m] INS ECEF X position
fx_e = ins(3); % [m/s/s] INS ECEF X acceleration
ry_e = ins(4); % [m] INS ECEF Y position
fy_e = ins(5); % [m/s/s] INS ECEF Y acceleration
rz_e = ins(7); % [m] INS ECEF Z position
fz_e = ins(8); % [m/s/s] INS ECEF Z acceleration
beta_a = ins(11:13); % Acc bias
beta_g = ins(20:22); % Gyr bias

% Gravity model
re = sqrt( rx_e^2 + ry_e^2 + rz_e^2 );  % [m] INS Earth Center Distance
goverr = MU_EARTH / re^3;
omegaEsq = OMEGA_EARTH^2;
N_e = zeros(3,3);
N_e(1,1) = goverr*( 3*(rx_e/re)^2 - 1 ) + omegaEsq;
N_e(1,2) = goverr*( 3*rx_e*ry_e/re^2 );
N_e(1,3) = goverr*( 3*rx_e*rz_e/re^2 );
N_e(2,1) = N_e(1,2);
N_e(2,2) = goverr*( 3*(ry_e/re)^2 - 1 ) + omegaEsq;
N_e(2,3) = goverr*( 3*ry_e*rz_e/re^2 );
N_e(3,1) = N_e(1,3);
N_e(3,2) = N_e(2,3);
N_e(3,3) = goverr*( 3*(rz_e/re)^2 - 1 );

% ECEF Force frame skew semetric matrix
F_e = [ 0, -fz_e, fy_e; fz_e, 0, -fx_e; -fy_e, fx_e, 0 ];

% Earth rotation skew symmetric matrix
Omega_e_ie = skewsym([ 0; 0; OMEGA_EARTH ]);

% Body to ECEF Rotation matrix
qk_e_b = ins(7:10,1);
R_e_b = quat2dircos(qk_e_b); % Body to ECEF rotation matrix
R_b_e = R_e_b';

% State Transition Partial Derivative
Ak = zeros(nx,nx);
Ak(1:3,4:6) = eye(3);
Ak(4:6,1:3) = N_e;
Ak(4:6,4:6) = -2*Omega_e_ie;
Ak(4:6,7:9) = - F_e;
Ak(4:6,10:12) = R_b_e;
Ak(7:9,7:9) = -Omega_e_ie;
Ak(7:9,13:21) = R_b_e;
Ak(10:12,10:12) = -beta_a;
Ak(13:15,13:15) = -beta_g;
Ak(22:24,22:24) = -beta_sf_acc;
Ak(25:27,25:27) = -beta_sf_gyr;

% Noise rate derivative
Dk = zeros(nx,nv);
Dk(4:6,1:3) = R_b_e;
Dk(7:9,4:6) = R_b_e;
Dk(10:12,7:9) = eye(3);
Dk(13:15,10:12) = eye(3);
Dk(16:18,13:15) = eye(3);
Dk(19:21,16:18) = eye(3);

% Integrate state
func = @(delt, x, v, dflag) Ak.*x + Dk.*v;
[xkp1, Fk, Gammak] = trapazoidIntegration(xk, vk, dt, nIntegrate, func, 1);


end

