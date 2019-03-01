function [ zbarkp1, Hkp1 ] = measurement_model( kp1, xbarkp1, wkp1 )
% INS Dead-Reckoning Measurement Model for 2D system
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

nx = 11;
nz = 3;

% Direction Cosine Matrix
invR0 = [0 -1; 1 0];
smpsi = sin(-xbarkp1(10));
cmpsi = cos(-xbarkp1(10));
invRpsi = [cmpsi, smpsi; -smpsi, cmpsi];

% Rotate acceleration to body frame
acc = invR0 * invRpsi * xbarkp1([3,6],1);

% Measurement prediction
zbarkp1 = zeros(nz,1);
zbarkp1(1:2) = acc + xbarkp1([9,10]);
zbarkp1(3:3) = xbarkp1([8]) + xbarkp1([11]);
zbarkp1 = zbarkp1 + wkp1;

Hkp1 = zeros(nz,nx);
Hkp1(1,3) = 1;
Hkp1(2,6) = 1;
Hkp1(3,8) = 1;
Hkp1(1,9) = 1;
Hkp1(2,10) = 1;
Hkp1(3,11) = 1;

end

