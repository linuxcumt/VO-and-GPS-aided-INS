function [ xkp1, Fk, Gammak ] ...
    = insErrorDynamicsModel_15state( dt, xk, vk, yk, f_e, beta_acc, beta_gyr, ...
                             beta_sf_acc, beta_sf_gyr, nInt )
% INS Dead-Reckoning Error State Dynamics Model
% Computes the rate of change of the error state and does numerical
% integration to compute the state and evolution of the partial
% derivatives.
% 
% This models only acceleration errors and does not consider rotating
% Earth.
% 
% INPUTS
% dt          - double
%               Delta time between state updates. Recommended to be INS
%               measurement rate.
% xk          - 15 x 1 double vector
%               State vector at time k
%                   [ del_r_e;       - Position Error
%                     del_v_e;       - Velocity Error
%                     err_e;         - Misalignment Error
%                     del_biasacc;   - Acc Bias Drift Error
%                     del_biasgyr;   - Byr Bias Drift Error
%                     biasinitacc_e; - Initial Acc Bias Error
%                     biasinitgyr_e; - Initial Gyr Bias Error
%                     scaleacc_e;    - Acc Scale Factor Error
%                     scalegyr_e ]   - Gyr Scale Factor Error
% vk          - 18 x 1 double vector
%               Process Noise vector at time k
%                   [ eta_a;   - Acc sensor noise (PSD q_a)
%                     eta_g;   - Gyr sensor noise (PSD q_g)
%                     eta_ba;  - Acc bias drift noise (PSD q_ba)
%                     eta_bg;  - Gyr bias drift noise (PSD q_bg)
%                     eta_sa;  - Acc Scale Factor noise (PSD q_sa)
%                     eta_sg ] - Gyr Scale Factor noise (PSD q_sg)
% ins         - 22 x 1 double vector
%               Current estimated state by INS
%                   [ rk_e;    - Position
%                     vk_e;    - Velocity
%                     qk_e_b;  - Orientation quaternion
%                     betaa;   - Accel bias estimate
%                     scalea;  - Accel Scale Factor
%                     betag;   - Accel bias estimate
%                     scaleg ] - Accel Scale Factor
% beta_acc    - 3 x 3 double matrix
%               Accelerometer bias Power Spectral Density
% beta_gyr    - 3 x 3 double matrix
%               Gyro bias Power Spectral Density
% beta_sf_acc - 3 x 3 double matrix
%               Accelerometer scale factor Power Spectral Density
% beta_sf_gyr - 3 x 3 double matrix
%               Gyro scale factor Power Spectral Density
% nInt        - int
%               Number of steps for integration
% 
% OUTPUTS
% xkp1        - 15 x 1 double vector
%               State vector at time k
% Fk          - 15 x 15 double matrix
%               State Transition partial derivative
% Gammak      - 15 x 18 double matrix
%               State Noiise partial derivative
% 
% @author: Matt Marti
% @date: 2019-03-12

% Constants
% global MUEARTH OMEGA_EARTH
nx = length(xk);
nv = length(vk);

% % Partition INS State
% rx_e = yk(1); % [m] INS ECEF X position
% ry_e = yk(2); % [m] INS ECEF Y position
% rz_e = yk(3); % [m] INS ECEF Z position
% fx_e = f_e(1); % [m/s/s] INS ECEF X acceleration
% fy_e = f_e(2); % [m/s/s] INS ECEF Y acceleration
% fz_e = f_e(3); % [m/s/s] INS ECEF Z acceleration
% 
% % Gravity model
% re = sqrt( rx_e^2 + ry_e^2 + rz_e^2 );  % [m] INS Earth Center Distance
% goverr = MUEARTH / re^3;
% omegaEsq = OMEGA_EARTH^2;
% N_e = zeros(3,3);
% N_e(1,1) = goverr*( 3*(rx_e/re)^2 - 1 ) + omegaEsq;
% N_e(1,2) = goverr*( 3*rx_e*ry_e/re^2 );
% N_e(1,3) = goverr*( 3*rx_e*rz_e/re^2 );
% N_e(2,1) = N_e(1,2);
% N_e(2,2) = goverr*( 3*(ry_e/re)^2 - 1 ) + omegaEsq;
% N_e(2,3) = goverr*( 3*ry_e*rz_e/re^2 );
% N_e(3,1) = N_e(1,3);
% N_e(3,2) = N_e(2,3);
% N_e(3,3) = goverr*( 3*(rz_e/re)^2 - 1 );
% 
% % ECEF Force frame skew semetric matrix
% F_e = [ 0, -fz_e, fy_e; fz_e, 0, -fx_e; -fy_e, fx_e, 0 ];
% 
% % Earth rotation skew symmetric matrix
% Omega_e_ie = skewsym([ 0; 0; OMEGA_EARTH ]);
% 
% % Body to ECEF Rotation matrix
% qk_e_b = yk(7:10,1);
% R_e_b = quat2dircos(qk_e_b); % Body to ECEF rotation matrix
% R_b_e = R_e_b';

% State Transition Partial Derivative
Ak = zeros(nx,nx);
Ak(1:3,4:6) = eye(3);
Ak(4:6,7:9) = -beta_acc;
Ak(13:15,13:15) = -beta_sf_acc;

% Noise rate derivative
Dk = zeros(nx,nv);
Dk(4:12,1:9) = eye(nv);

% Integrate state
func = @(delt, x, v, dflag) deal(Ak*x + Dk*v, Ak, Dk);
[xkp1, Fk, Gammak] = trapazoidIntegration(xk, vk, dt, nInt, func, 1);


end

