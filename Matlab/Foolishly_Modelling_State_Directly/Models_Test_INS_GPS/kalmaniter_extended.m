function [ xhatkp1, Pkp1, Wkp1, epsilonkp1, xbarkp1, Pbarkp1, Fk ] ...
    = kalmaniter_extended( kfdata )
% Computes one iteration of the Extended Kalman Filter
% Given the state estimate xhat(k), covariance matrix, nonlinear 
% measurement and dynamics models, and noise covariance matrices, computes
% the estimate xhat(k+1). This is the implementation of the Extended Kalman
% Filter.
% 
% This function performs the Extended Kalman Filter state estimate for the
% dynamic system model
% 
%   x(k+1)   = f( k, x(k), u(k), v(k) )
%   F(k)     = df / dx(k)
%   Gamma(k) = df / dv(k)
%   z(k+1)   = h( k+1, x(k+1) ) + w(k+1)
%   H(k+1)   = dh / dx(k+1)
% 
% Where v(k) and w(k+1) are Gaussian distributed random noise terms. The
% process noise v(k) has a mean of zero, and a Covariance Q(k). The
% measurement noise w(k+1) has a mean of zero, and a Covariance of R(k+1).
% 
% Variables ending in "kp1" represent the variable at the (k+1) state.
% Variables ending in "k" represent variable at the (k) state.
% 
% INPUTS
% kfdata     - struct
%              Kalman Filter data, including state matrices, xhat(k), etc.
%              A complete list of the required 9 data fields follows:
%              - k     - int
%                        The index of current time in time history
%              - xhatk - nx x 1 double matrix
%                        Previous iteration state estimate vector
%              - Pk    - nx x nx double matrix
%                        Previous state covariance matrix
%              - uk    - nu x 1 double vector
%                        Control input vector
%              - f     - Anonymous function handle
%                        State transition function
%                        Must accept 4 inputs in this order:
%                        - k      - int
%                                   The index of current time in time 
%                                   history
%                        - xk     - nx x 1 double vector
%                                   The state vector at k time
%                        - uk     - nu x 1 double vector
%                                   The control input vector
%                        - vk     - nv x 1 double vector
%                                   Process noise
%                                   Note that for the EKF, this will be
%                                   input as a zero vector.
%                        Must return 3 outputs in this order:
%                        - xkp1   - nx x 1 double vector
%                                   The state vector at k+1 time
%                        - Fk     - nx x nx double matrix
%                                   The partial derivative of state
%                                   transition function with respect to the
%                                   state vector, evaluated at k time
%                        - Gammak - nx x nv double matrix
%                                   The partial derivative of state
%                                   transition function with respect to the
%                                   noise vector, evaluated at k time
%              - Qk    - nv x nv double matrix
%                        Process noise covariance matrix
%              - zkp1  - nz x 1 double vector
%                        Next state measurement vector
%              - h     - Anonymous function handle
%                        Measurement model
%                        Must accept 2 inputs in this order:
%                        - kp1    - int
%                                   The index of current time in time 
%                                   history
%                        - xkp1   - nx x 1 double vector
%                                   The state vector at k time
%                        Must return 3 outputs in this order:
%                        - zkp1   - nz x 1 double vector
%                                   The measurement vector at k+1 time
%                        - Hkp1   - nz x nx double matrix
%                                   The partial derivative of measurement
%                                   function with respect to the state 
%                                   vector, evaluated at k time
%              - Rkp1  - nz x nz double matrix
%                        Measurement noise covariance matrix
% 
% OUTPUTS
% xhatkp1    - nx x 1 double matrix
%              xhat(k+1) next state estimate.
% Pkp1       - nx x nx double matrix
%              P(k+1) State covariance matrix for next state.
% Wkp1       - nx x nx double matrix
%              Kalman Gain matrix
% epsilonkp1 - double scalar
%              Innovation Statistic for study of filter consistency
% xbarkp1    - nx x 1 double matrix
%              A Priori state estimate
% Pbarkp1    - nx x nx double matrix
%              A Priori state covariance
% Fk         - nx x nx double matrix
%              State Transition Function partial derivative with respect to
%              xhatk.
% 
% DEPENDENCIES
%
% @author: Matt Marti
% @date: 2018-12-13


%% Input Assignment

assert(size(fieldnames(kfdata), 1) == 9, ...
    ['Incorrect number of fields in kfdata struct. Check that the ',...
    'field names are correct']);

% Assign variables from struct
k      = kfdata.k;
xhatk  = kfdata.xhatk;
Pk     = kfdata.Pk;
f      = kfdata.f;
uk     = kfdata.uk;
Qk     = kfdata.Qk;
h      = kfdata.h;
zkp1   = kfdata.zkp1;
Rkp1   = kfdata.Rkp1;

% Sizes
nx = size(xhatk, 1);
nz = size(zkp1, 1);
nv = size(Qk,1);

% State propagation
[ xbarkp1, Fk, Gammak ] = f(k, xhatk, uk, zeros(nv,1));
[ zbarkp1, Hkp1 ] = h(k+1, xbarkp1);

% Assertion checking for sizes
assert(size(xbarkp1,1) == nx && size(xbarkp1,2) == 1, ...
    'Incorrect size of f output: xkp1');
assert(size(Pk,1) == nx && size(Pk,2) == nx, ...
    'Incorrect size of Pk');
assert(size(Fk,1) == nx && size(Fk,2) == nx, ...
    'Incorrect size of f output: Fk');
assert(size(Gammak,1) == nx && size(Gammak,2) == nv, ...
    'Incorrect size of f output: Gammak');
assert(size(Qk,1) == nv && size(Qk,2) == nv,...
    'Incorrect size of Qk');
assert(size(Rkp1,1) == nz && size(Rkp1,2) == nz, ...
    'Incorrect size of Rkp1');
assert(size(zbarkp1,1) == nz && size(zbarkp1,2) == 1, ...
    'Incorrect size of h output: zkp1')
assert(size(Hkp1,1) == nz && size(Hkp1,2) == nx, ...
    'Incorrect size of h output: Hkp1');


%% Iteration

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


%% Alternative formulas
% Pkp1 = inv(inv(Pbarkp1) + (Hkp1')*(Rkp1')*Hkp1);
% Wkp1 = Pkp1 * Hkp1' * inv(Rkp1);
% 
% eyemWH = eye(nx) - Wkp1*Hkp1;
% Pkp1 = eyemWH * Pbarkp1 * (eyemWH') + Wkp1*Rkp1*(Wkp1');

end