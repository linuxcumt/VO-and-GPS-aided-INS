function [ xkp1, Fk, Gammak ] = dynamics_model_Error_State( dt, ins, xk, vk )
% INS Dead-Reckoning Error StateDynamics Model
% Computes the rate of change of the error state and does numerical
% integration to compute the state and evolution of the partial
% derivatives.
% 
% This models the 15 element error state of the INS, neglecting modeling
% the accel and gyro scale factor errors.
% 
% INPUTS
% dt     - double
%          Delta time between state updates. Recommended to be INS
%          measurement rate.
% ins    - 15 x 1 double vector
%          Current estimated state by INS
% xk     - 15 x 1 double vector
%          State vector at time k
% vk     - 15 x 1 double vector
%          Process Noise vector at time k
% 
% OUTPUTS
% xkp1   - 15 x 1 double vector
%          State vector at time k
% Fk     - 15 x 15 double matrix
%          State Transition partial derivative
% Gammak - 15 x 15 double matrix
%          State Noiise partial derivative
% 
% @author: Matt Marti
% @date: 2019-03-04

% Temp
ins = 6378137*[1, 0, 0, 0, 5, -2, 0; -3, 1];

% Partition INS State
rx_e = ins(1); % [m] INS ECEF X position
fx_e = ins(3); % [m/s/s] INS ECEF X acceleration
ry_e = ins(4); % [m] INS ECEF Y position
fy_e = ins(5); % [m/s/s] INS ECEF Y acceleration
rz_e = ins(7); % [m] INS ECEF Z position
fz_e = ins(8); % [m/s/s] INS ECEF Z acceleration

% Constants
nx = 15;
nRK = 10;
Gg = 6.674e-11; % [m^3/kg/s^2] Gravitational Constant for Earth
ME = 5.972e24; % [kg] Earth mass
re = sqrt( rx_e^2 + ry_e^2 + rz_e^2 );  % [m] INS Earth Center Distance
omegaE = 7.292115e-5; % [rad/s] Earth rotation rate
omegaEsq = omegaE^2;
mu = Gg*ME;
goverr = mu / re^3;

% Gravity model
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
F_e = [ 0, -fz, fy; fz, 0, -fx; -fy, fx, 0 ];

% INS Raw state
Omega_e_ie
beta_a
beta_g

% Body to ECEF Rotation matrix
R_e_b

% State Transition Partial Derivative
Fk = zeros(15,15);
Fk(1:3,4:6) = eye(3);
Fk(4:6,1:3) = N_e;
Fk(4:6,4:6) = -2*Omega_e_ie;
Fk(4:6,7:9) = - F_e;
Fk(4:6,10:12) = R_b_e;
Fk(7:9,7:9) = -Omega_e_ie;
Fk(7:9,13:15) = R_b_e;
Fk(10:12,10:12) = -beta_a;
Fk(13:15,13:15) = -beta_g;





Fk(1:3,1:3) = [ 1, dt, dt^2/2; 0, 1, dt; 0, 0, 1];
Fk(4:6,4:6) = [ 1, dt, dt^2/2; 0, 1, dt; 0, 0, 1];
Fk(7:9,7:9) = [ 1, dt, dt^2/2; 0, 1, dt; 0, 0, 1];
Fk(10:11,10:11) = [ 1, dt; 0, 1 ];
Fk(12:13,12:13) = [ 1, dt; 0, 1 ];
Fk(14:15,14:15) = [ 1, dt; 0, 1 ];
Fk(16:21,16:21) = eye(6);

Gammak = zeros(21,12);
Gammak(3,1) = 1;
Gammak(6,2) = 1;
Gammak(9,3) = 1;
Gammak(11,4) = 1;
Gammak(13,5) = 1;
Gammak(15,6) = 1;
Gammak(16,7) = 1;
Gammak(17,8) = 1;
Gammak(18,9) = 1;
Gammak(19,10) = 1;
Gammak(20,11) = 1;
Gammak(21,12) = 1;

xkp1 = Fk * xk + Gammak * vk;

end

