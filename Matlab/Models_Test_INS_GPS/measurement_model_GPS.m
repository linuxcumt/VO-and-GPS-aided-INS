function [ zk, Hk ] = measurement_model_GPS( xk, wk )
% INS Dead-Reckoning Measurement Model
% 
% INPUTS
% kp1     - int
%           time index
% xbarkp1 - 21x1 double vector
%           State vector at time k+1
% wkp1    - 6x1 double vector
%           Measurement Noise vector at time k
% 
% OUTPUTS
% zbarkp1 - 6x1 double vector
%           Measurement prediction
% Hkp1    - 6x21 double matrix
%           Measurement partial derivative
% 
% @author: Matt Marti
% @date: 2019-03-02

% Constants
nx = 21;
nz = 3;

% Pertinent elements of state vector
x1     = xk(1);
x1d    = xk(2); %#ok
x1dd   = xk(3); %#ok
x2     = xk(4);
x2d    = xk(5); %#ok
x2dd   = xk(6); %#ok
x3     = xk(7);
x3d    = xk(8); %#ok
x3dd   = xk(9); %#ok
psi    = xk(10); %#ok
psid   = xk(11); %#ok
theta  = xk(12); %#ok
thetad = xk(13); %#ok
phi    = xk(14); %#ok
phid   = xk(15); %#ok
beta1  = xk(16); %#ok
beta2  = xk(17); %#ok
beta3  = xk(18); %#ok
gamma1 = xk(19); %#ok
gamma2 = xk(20); %#ok
gamma3 = xk(21); %#ok

% Measurement prediction
zk = zeros(nz,1);
zk(1:3) = [x1; x2; x3];
zk = zk + wk;

% Measurement Partial Derivative
Hk = zeros(nz,nx);
Hk(1,1) = 1;
Hk(2,4) = 1;
Hk(3,7) = 1;

end

