function [ zk, Hk ] = measurement_model_INS( xk, wk )
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
% zk      - 6x1 double vector
%           Measurement prediction
% Hk      - 6x21 double matrix
%           Measurement partial derivative
% 
% @author: Matt Marti
% @date: 2019-03-01

% Constants
nx = 21;
nz = 6;

% Pertinent elements of state vector
x1     = xk(1);
x1d    = xk(2);
x1dd   = xk(3);
x2     = xk(4);
x2d    = xk(5);
x2dd   = xk(6);
x3     = xk(7);
x3d    = xk(8);
x3dd   = xk(9);
psi    = xk(10);
psid   = xk(11);
theta  = xk(12);
thetad = xk(13);
phi    = xk(14);
phid   = xk(15);
beta1  = xk(16);
beta2  = xk(17);
beta3  = xk(18);
gamma1 = xk(19);
gamma2 = xk(20);
gamma3 = xk(21);

zk = [ x1; x1d; x2; x2d; x3; x3d ];
zk = zk + wk;

Hk = zeros(nz,nx);
Hk(1,1) = 1;
Hk(2,2) = 1;
Hk(3,4) = 1;
Hk(4,5) = 1;
Hk(5,7) = 1;
Hk(6,8) = 1;

end

