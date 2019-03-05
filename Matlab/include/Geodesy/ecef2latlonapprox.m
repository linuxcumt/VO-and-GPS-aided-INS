function [ latlon, dlatlon_dr ] = ecef2latlonapprox( r, dflag )
% Gives approximate Latitude and Longitude from ECEF coordinates
% Used to get Lat and Lon angles to make a rotation matrix. Also can
% compute derivatives of the approximate Lat Lon values with respect to
% ECEF coordinate vector r.
% 
% @arg
% r          - 3 x 1 double matrix
%              ECEF coordinates
% dflag      - bool
%              Compute derivatives flag. True to compute derivatives
% 
% @return
% latlon     - 2 x 1 double matrix
%              [ Lat; Lon ]
% dlatlon_dr - 2 x 3 double matrix
%              Latitude and Longitude partial derivatives
% 
% @author: Matt Marti
% @date: 2019-03-04

% Input check
if nargin == 1
    dflag = 0;
    dlatlon_dr = NaN;
end

% Input assignment
r1 = r(1);
r2 = r(2);
r3 = r(3);

% ECEF to ENV approximate Rotation Matrix
lonapprox = atan2(r2, r1);
sumecefsq = r1^2 + r2^2 + r3^2;
normecef = sqrt(sumecefsq);
oneovernormecef = 1 / normecef;
runit = r * oneovernormecef;
latapprox = asin(runit(3));
latlon = [latapprox; lonapprox];


%% Partial Derivatives for Measurements

% Decision making matrix
dlatlon_dr = zeros(2, 3);
switchmatrix = eye(3);

% Derivative loop
for i = 1:3

    % Activate partial derivatives
    dr1_ds = switchmatrix(1,i);
    dr2_ds = switchmatrix(2,i);
    dr3_ds = switchmatrix(3,i);
    
    % ECEF to ENV approximate Rotation Matrix
    rxsqprysq = r1^2 + r2^2;
    dlonapprox_ds = (-r2*dr1_ds+r1*dr2_ds)/rxsqprysq;
    
    dr_ds = [dr1_ds; dr2_ds; dr3_ds];
    dsumecefsq_ds = 2*(r1*dr1_ds + r2*dr2_ds + r3*dr3_ds);
    dnormecef_ds = 0.5/sqrt(sumecefsq)*dsumecefsq_ds;
    doneovernormecef_ds = - oneovernormecef^2 * dnormecef_ds;
    decefunit_ds = dr_ds*oneovernormecef + r*doneovernormecef_ds;
    dlatapprox_ds = decefunit_ds(3)/sqrt(1-runit(3)^2);

    % Assign Derivative
    dlatlon_dr(1,i) = dlonapprox_ds;
    dlatlon_dr(2,i) = dlatapprox_ds;
end

end

