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
g0 = -9.81;

% Pertinent elements of state vector
x1     = xk(1);
x1d    = xk(2); %#ok
x1dd   = xk(3);
x2     = xk(4);
x2d    = xk(5); %#ok
x2dd   = xk(6);
x3     = xk(7);
x3d    = xk(8); %#ok
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

% ECEF to ENV approximate Rotation Matrix
lonapprox = atan2(x2, x1);
slon = sin(lonapprox);
clon = cos(lonapprox);
Rlon = [ clon, slon, 0; -slon, clon, 0; 0, 0, 1 ];
ecef = [x1; x2; x3];
sumecefsq = x1^2 + x2^2 + x3^2;
normecef = sqrt(sumecefsq);
oneovernormecef = 1 / normecef;
ecefunit = ecef * oneovernormecef;
latapprox = asin(ecefunit(3));
slat = sin(latapprox);
clat = cos(latapprox);
Rlat = [ clat, 0, -slat; 0, 1, 0; slat, 0, clat ];
Renvframe = [ 0, -1, 0; 1, 0, 0; 0, 0, 1 ];
R_ecef_env = Renvframe * Rlat * Rlon;

% ENV to Body Rotation Matrix
spsi = sin(psi);
cpsi = cos(psi);
Rpsi = [cpsi, spsi, 0; -spsi, cpsi, 0; 0, 0, 1];
stheta = sin(theta);
ctheta = cos(theta);
Rtheta = [ctheta, 0, -stheta; 0, 1, 0; stheta, 0, ctheta];
sphi = sin(phi);
cphi = cos(phi);
Rphi = [1, 0, 0; 0, cphi, sphi; 0, -sphi, cphi ];
R_env_body = Rphi * Rtheta * Rpsi;

% Full Direction Cosine Matrix to Body Frame
Rbodyframe = [0 -1 0; 1 0 0; 0 0 1];
R_ecef_body = Rbodyframe * R_env_body * R_ecef_env;

% Accels
g = Rbodyframe*R_env_body*[0;0;g0];
acc = R_ecef_body * [x1dd; x2dd; x3dd];
beta = [beta1; beta2; beta3];

% Gyros
gyr = [psid; thetad; phid];
gamma = [gamma1; gamma2; gamma3];

% Measurement prediction
zk = zeros(nz,1);
zk(1:3) = acc + beta + g;
zk(4:6) = gyr + gamma;
zk = zk + wk;


%% Partial Derivatives for Measurements

% Decision making matrix
Hk = zeros(nz, nx);
switchmatrix = eye(nx);

% Derivative loop
for i = 1:nx

    % Activate partial derivatives
    dx1_ds     = switchmatrix(1,i);
    dx1d_ds    = switchmatrix(2,i); %#ok
    dx1dd_ds   = switchmatrix(3,i);
    dx2_ds     = switchmatrix(4,i);
    dx2d_ds    = switchmatrix(5,i); %#ok
    dx2dd_ds   = switchmatrix(6,i);
    dx3_ds     = switchmatrix(7,i);
    dx3d_ds    = switchmatrix(8,i); %#ok
    dx3dd_ds   = switchmatrix(9,i);
    dpsi_ds    = switchmatrix(10,i);
    dpsid_ds   = switchmatrix(11,i);
    dtheta_ds  = switchmatrix(12,i);
    dthetad_ds = switchmatrix(13,i);
    dphi_ds    = switchmatrix(14,i);
    dphid_ds   = switchmatrix(15,i);
    dbeta1_ds  = switchmatrix(16,i);
    dbeta2_ds  = switchmatrix(17,i);
    dbeta3_ds  = switchmatrix(18,i);
    dgamma1_ds = switchmatrix(19,i);
    dgamma2_ds = switchmatrix(20,i);
    dgamma3_ds = switchmatrix(21,i);
    
    % ECEF to ENV approximate Rotation Matrix
    xsqpysq = x1^2 + x2^2;
    dlonapprox_ds = (-x2*dx1_ds+x1*dx2_ds)/xsqpysq;
    dslon_ds = clon * dlonapprox_ds;
    dclon_ds = -slon * dlonapprox_ds;
    dRlon_ds = [ dclon_ds, dslon_ds, 0; -dslon_ds, dclon_ds, 0; 0, 0, 0 ];
    decef_ds = [dx1_ds; dx2_ds; dx3_ds];
    dsumecefsq_ds = 2*(x1*dx1_ds + x2*dx2_ds + x3*dx3_ds);
    dnormecef_ds = 0.5/sqrt(sumecefsq)*dsumecefsq_ds;
    doneovernormecef_ds = - oneovernormecef^2 * dnormecef_ds;
    decefunit_ds = decef_ds*oneovernormecef + ecef*doneovernormecef_ds;
    dlatapprox_ds = decefunit_ds(3)/sqrt(1-ecefunit(3)^2);
    dslat_ds = clat * dlatapprox_ds;
    dclat_ds = -slat * dlatapprox_ds;
    dRlat_ds = [ dclat_ds, 0, -dslat_ds; 0, 0, 0; dslat_ds, 0, dclat_ds ];
    dR_ecef_env_ds = Renvframe*(Rlat*dRlon_ds+dRlat_ds*Rlon);

    % ENV to Body Rotation Matrix
    dspsi_ds = cpsi * dpsi_ds;
    dcpsi_ds = -spsi * dpsi_ds;
    dRpsi_ds = [dcpsi_ds, dspsi_ds, 0; -dspsi_ds, dcpsi_ds, 0; 0, 0, 0];
    dstheta_ds = ctheta * dtheta_ds;
    dctheta_ds = -stheta * dtheta_ds;
    dRtheta_ds = [dctheta_ds, 0, -dstheta_ds; 0, 0, 0; dstheta_ds, 0, dctheta_ds];
    dsphi_ds = cphi * dphi_ds;
    dcphi_ds = -sphi * dphi_ds;
    dRphi_ds = [0, 0, 0; 0, dcphi_ds, dsphi_ds; 0, -dsphi_ds, dcphi_ds ];
    dR_env_body_ds = dRphi_ds*Rtheta*Rpsi + Rphi*dRtheta_ds*Rpsi + Rphi*Rtheta*dRpsi_ds;

    % Full Direction Cosine Matrix to Body Frame
    Rbodyframe = [0 -1 0; 1 0 0; 0 0 1];
    dR_ecef_body_ds = Rbodyframe*(dR_env_body_ds*R_ecef_env + R_env_body*dR_ecef_env_ds);

    % Accels
    dg_ds = Rbodyframe*dR_env_body_ds*[0;0;g0];
    dacc_ds = dR_ecef_body_ds*[x1dd; x2dd; x3dd] + R_ecef_body*[dx1dd_ds; dx2dd_ds; dx3dd_ds];
    dbeta_ds = [dbeta1_ds; dbeta2_ds; dbeta3_ds];

    % Gyros
    dgyr_ds = [dpsid_ds; dthetad_ds; dphid_ds];
    dgamma_ds = [dgamma1_ds; dgamma2_ds; dgamma3_ds];

    % Assign Derivative
    Hk(1:3,i) = dacc_ds + dbeta_ds + dg_ds;
    Hk(4:6,i) = dgyr_ds + dgamma_ds;
end

end

