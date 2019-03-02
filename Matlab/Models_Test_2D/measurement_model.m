function [ z, Hk ] = measurement_model( k, xk, wk )
% INS Dead-Reckoning Measurement Model for 2D system
% 
% @arg
% kp1       - int
%             time index
% xkp1      - 21x1 double vector
%             State vector at time k+1
% wkp1      - 6x1 double vector
%             Measurement Noise vector at time k
% 
% @return
% zkp1      - 6x1 double vector
%             Measurement prediction
% Hkp1      - 6x21 double matrix
%             Measurement partial derivative
% 
% @author: Matt Marti
% @date: 2019-03-01

nx = 11;
nz = 3;

% Pertinent elements of state vector
x1dd  = xk(3);
x2dd  = xk(6);
psi   = xk(7);
psid  = xk(8);
beta1 = xk(9);
beta2 = xk(10);
gamma = xk(11);

% Direction Cosine Matrix
invR0 = [0 -1; 1 0];
spsi = sin(psi);
cpsi = cos(psi);
Rpsi = [cpsi, spsi; -spsi, cpsi];

% Rotate acceleration to body frame
acc = invR0 * Rpsi * [x1dd; x2dd];

% Measurement prediction
z = zeros(nz,1);
z(1:2) = acc + [beta1; beta2];
z(3) = psid + gamma;
z = z + wk;


%% Partial Derivatives for Measurements

% Decision making matrix
Hk = zeros(nz, nx);
switchmatrix = eye(nx);

% Derivative loop
for i = 1:nx

    % Activate partial derivatives
    dx1_ds    = switchmatrix(i,1); %#ok
    dx1d_ds   = switchmatrix(i,2); %#ok
    dx1dd_ds  = switchmatrix(i,3);
    dx2_ds    = switchmatrix(i,4); %#ok
    dx2d_ds   = switchmatrix(i,5); %#ok
    dx2dd_ds  = switchmatrix(i,6);
    dpsi_ds   = switchmatrix(i,7);
    dpsid_ds  = switchmatrix(i,8);
    dbeta1_ds = switchmatrix(i,9);
    dbeta2_ds = switchmatrix(i,10);
    dgamma_ds = switchmatrix(i,11);

    % Direction Cosine Matrix
    dspsi_ds = cpsi * dpsi_ds;
    dcpsi_ds = - spsi * dpsi_ds;
    dRpsi_ds = [dcpsi_ds, dspsi_ds; -dspsi_ds, dcpsi_ds];

    % Rotate acceleration to body frame
    dacc_ds = invR0 * ( dRpsi_ds * [x1dd; x2dd] ...
                      + Rpsi * [dx1dd_ds; dx2dd_ds] );

    % Assign Derivative
    Hk(1:2,i) = dacc_ds + [dbeta1_ds; dbeta2_ds];
    Hk(3,i) = dpsid_ds + dgamma_ds;
end

end

