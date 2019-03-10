function [ xhatkp1, Pkp1, zxkp1, Rxxkp1, epsilonkp1, xbarkp1, Pbarkp1 ] ...
    = srifiter_unscented( kfdata )
% Computes one iteration of the Sigma Points SRIF
% Given the state estimate xhat(k), covariance matrix, nonlinear 
% measurement and dynamics models, and noise covariance matrices, computes
% the estimate xhat(k+1). This is the implementation of the Efficient
% Square Root Information Sigma Points Kalman Filter.
% 
% This function performs the Extended Kalman Filter state estimate for the
% dynamic system model
% 
%   x(k+1)   = f( k, x(k), u(k), v(k) )
%   z(k+1)   = h( k+1, x(k+1) ) + w(k+1)
% 
% Where v(k) and w(k+1) are Gaussian distributed random noise terms. The
% process noise v(k) has a mean of zero, and a Covariance Q(k). The
% measurement noise w(k+1) has a mean of zero, and a Covariance of R(k+1).
% 
% This function is initialized by setting the "Rxxk" and "zxk" values in
% the "kfdata" struct to "NaN". On the first call of this function, the
% cholesky factor of the covariance is computed to initialize the filter.
% For all calls to this function after that, the "zxkp1" and "Rxxkp1"
% should be included in the "kfdata" struct as "zxk" and "Rxxk" fields.
% 
% Variables ending in "kp1" represent the variable at the (k+1) state.
% Variables ending in "k" represent variable at the (k) state.
% 
% INPUTS
% kfdata     - struct
%              Kalman Filter data, including state matrices, xhat(k), etc.
%              A complete list of the required 14 data fields follows:
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
%                                   The process noise vector
%                        Must return 1 output:
%                        - xkp1   - nx x 1 double vector
%                                   The state vector at k+1 time
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
%                        Must return 1 output:
%                        - zkp1   - nz x 1 double vector
%                                   The measurement vector at k+1 time
%              - Rkp1  - nz x nz double matrix
%                        Measurement noise covariance matrix
%              - alpha - double scalar
%                        Central particle weight factor
%              - beta  - double scalar
%                        Central covariance particle weight factor
%              - kappa - double scalar
%                        Central particle weight factor
%              - zxk   - nx x 1 double vector
%                        Factorized State from previous iteration
%                        Initialize the filter by setting this to NaN on
%                        the first iteration of the filter.
%              - Rxxk  - nx x nx double matrix
%                        Factorized Covariance from previous iteration
%                        Initialize the filter by setting this to NaN on
%                        the first iteration of the filter.
% 
% OUTPUTS
% xhatkp1    - nx x 1 double matrix
%              A Posteriori state estimate
% Pkp1       - nx x nx double matrix
%              A Posteriori Covariance matrix
% zxkp1vv    - nx x 1 double vector
%              Factorized A Posteriori state estimate
% Rxxkp1     - nx x nx double matrix
%              Factorized A Posteriori Covariance
% epsilonkp1 - double scalar
%              Innovation Statistic for study of filter consistency
% xbarkp1    - nx x 1 double matrix
%              A Priori state estimate
% Pbarkp1    - nx x nx double matrix
%              A Priori covariance
% 
% DEPENDENCIES
%
% @author: Matt Marti
% @date: 2018-12-05

error('Not implemented');
%% Input Assignment

assert(size(fieldnames(kfdata), 1) == 14, ...
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
alpha  = kfdata.alpha;
kappa  = kfdata.kappa;
beta   = kfdata.beta;
zxk    = kfdata.zxk;
Rxxk   = kfdata.Rxxk;

% Sizes
nx = size(xhatk, 1);
nz = size(zkp1, 1);
nv = size(Qk,1);

% State propagation
xscrbarkp10 = f(k, xhatk, uk, zeros(nv,1));
zscrbarkp10 = h(k+1, xscrbarkp10);

% Assertion checking for sizes
if sum(isnan(zxk))
    assert(isnan(Rxxk), 'Rxxk is not NaN but zxk is');
else
    nx = size(zxk,1);
    assert((size(Rxxk,1) == nx && size(Rxxk,2) == nx), ...
        'Incorrect size of Rxxk');
    assert((size(zxk,1) == nx && size(zxk,2) == 1), ...
        'Incorrect size of zxk');
end
assert(size(xscrbarkp10,1) == nx && size(xscrbarkp10,2) == 1, ...
    'Incorrect size of f output: xkp1');
assert(size(Pk,1) == nx && size(Pk,2) == nx, ...
    'Incorrect size of Pk');
assert(size(Qk,1) == nv && size(Qk,2) == nv,...
    'Incorrect size of Qk');
assert(size(Rkp1,1) == nz && size(Rkp1,2) == nz, ...
    'Incorrect size of Rkp1');
assert(size(zscrbarkp10,1) == nz && size(zscrbarkp10,2) == 1, ...
    'Incorrect size of h output: zkp1')
assert(numel(alpha) == 1, 'Incorrect size of alpha');
assert(numel(beta) == 1, 'Incorrect size of beta');
assert(numel(kappa) == 1, 'Incorrect size of kappa');

% Iteration for Rxxk
if isnan(Rxxk)
    Rxxk = inv(chol(Pk)');
    zxk = Rxxk * xhatk; %#ok
end


%% Iteration

twonxpnv = 2*(nx+nv);

% Compute parameters
lambda = alpha * (nx + nv + kappa) - (nx + nv);
sqrtnxpnvplambda = sqrt(nx + nv + lambda);

% Compute covariance square roots
Sxk = chol(Pk)';
Svk = chol(Qk)';

% Generate the Sigma Points
xscrkmat = zeros(nx, twonxpnv);
vscrkmat = zeros(nv, twonxpnv);
for i = 1:nx
    xscrkmat(:,i) = xhatk + sqrtnxpnvplambda * Sxk(:,i);
end
for i = nx+1:2*nx
    imnx = i - nx;
    xscrkmat(:,i) = xhatk - sqrtnxpnvplambda * Sxk(:,imnx);
end
for i = (2*nx+1):(2*nx+nv)
    im2nx = i - 2*nx;
    xscrkmat(:,i) = xhatk;
    vscrkmat(:,i) = sqrtnxpnvplambda * Svk(:,im2nx);
end
for i = (2*nx+nv+1):twonxpnv
    im2nxmnv = i - 2*nx - nv;
    xscrkmat(:,i) = xhatk;
    vscrkmat(:,i) = - sqrtnxpnvplambda * Svk(:,im2nxmnv);
end

% Propagate sigma points with dynamics model
xscrbarkp1mat = zeros(nx, twonxpnv);
zscrbarkp1mat = zeros(nz, twonxpnv);
for i = 1:twonxpnv
    xscrbarkp1mat(:,i) = f(k, xscrkmat(:,i), uk, vscrkmat(:,i));
    zscrbarkp1mat(:,i) = h(k+1, xscrbarkp1mat(:,i));
end

% Compute the A Priori mean
oneovernxpnvplambda = 1 / (nx + nv + lambda);
oneovertwonxpnvplambda = 0.5 * oneovernxpnvplambda;
wm0 = lambda * oneovernxpnvplambda;
wmvec = oneovertwonxpnvplambda*ones(1,twonxpnv);
xbarkp1 = sum([wm0 wmvec].*[xscrbarkp10 xscrbarkp1mat], 2);
zbarkp1 = sum([wm0 wmvec].*[zscrbarkp10 zscrbarkp1mat], 2);

% Compute the A Priori covariance and cross correlation weights
xscrmxbar0 = xscrbarkp10 - xbarkp1;
zscrmzbar0 = zscrbarkp10 - zbarkp1;
xscrmxbarmat = xscrbarkp1mat - xbarkp1;
zscrmzbarmat = zscrbarkp1mat - zbarkp1;
wc0 = wm0 + 1 - alpha^2 + beta;
wcvec = wmvec;

% Compute the A Priori covariances and cross correlation
Pbarkp1 = wc0*xscrmxbar0*(xscrmxbar0');
Pxzbarkp1 = wc0*xscrmxbar0*(zscrmzbar0');
Pzzbarkp1 = wc0*zscrmzbar0*(zscrmzbar0') + Rkp1;
for i = 1:twonxpnv
    Pbarkp1 = Pbarkp1 + wcvec(i)*xscrmxbarmat(:,i)*(xscrmxbarmat(:,i)');
    Pxzbarkp1 = Pxzbarkp1 + wcvec(i)*xscrmxbarmat(:,i)*(zscrmzbarmat(:,i)');
    Pzzbarkp1 = Pzzbarkp1 + wcvec(i)*zscrmzbarmat(:,i)*(zscrmzbarmat(:,i)');
end

% Kalman Gain calculation
nukp1 = zkp1 - zbarkp1;
invPzzbarkp1 = inv(Pzzbarkp1);
Wkp1 = Pxzbarkp1*invPzzbarkp1; %#ok

% Innovation Statistic
epsilonkp1 = (nukp1')*invPzzbarkp1*nukp1; %#ok

% Measurement update of state covariance
xhatkp1 = xbarkp1 + Wkp1*nukp1;
Pkp1 = Pbarkp1 - Pxzbarkp1*invPzzbarkp1*(Pxzbarkp1'); %#ok


%% Alternative formulas
% Pkp1 = inv(inv(Pbarkp1) + (Hkp1')*(Rkp1')*Hkp1);
% Wkp1 = Pkp1 * Hkp1' * inv(Rkp1);
% 
% eyemWH = eye(nx) - Wkp1*Hkp1;
% Pkp1 = eyemWH * Pbarkp1 * (eyemWH') + Wkp1*Rkp1*(Wkp1');

end