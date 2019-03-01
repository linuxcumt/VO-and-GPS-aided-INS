function [ zbarkp1, Hkp1 ] = measurement_model( kp1, xbarkp1, wkp1 )
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
% @date: 2019-02-25

nx = 21;
nz = 6;

% Direction Cosine Matrix
invR0 = [0 -1 0; 1 0 0; 0 0 1];
smphi = sin(-xbarkp1(14));
cmphi = cos(-xbarkp1(14));
invRphi = [1, 0, 0; 0, cmphi, smphi; 0, -smphi, cmphi ];
smtheta = sin(-xbarkp1(12));
cmtheta = cos(-xbarkp1(12));
invRtheta = [cmtheta, 0, -smtheta; 0, 1, 0; smtheta, 0, cmtheta];
smpsi = sin(-xbarkp1(10));
cmpsi = cos(-xbarkp1(10));
invRpsi = [cmpsi, smpsi, 0; -smpsi, cmpsi, 0; 0, 0, 1];

% Rotate acceleration to body frame
acc = invR0 * invRphi * invRtheta * invRpsi * xbarkp1([3,6,9],1);

% Measurement prediction
zbarkp1 = zeros(nz,1);
zbarkp1(1:3) = acc + xbarkp1([16,17,18]);
zbarkp1(4:6) = xbarkp1([11,13,15]) + xbarkp1([19,20,21]);
zbarkp1 = zbarkp1 + wkp1;

Hkp1 = zeros(nz,nx);
Hkp1(1,3) = 1;
Hkp1(2,6) = 1;
Hkp1(3,9) = 1;
Hkp1(4,11) = 1;
Hkp1(5,13) = 1;
Hkp1(6,15) = 1;
Hkp1(1,16) = 1;
Hkp1(2,17) = 1;
Hkp1(3,18) = 1;
Hkp1(4,19) = 1;
Hkp1(5,20) = 1;
Hkp1(6,21) = 1;

end

