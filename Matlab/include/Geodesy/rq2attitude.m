function [ latlonalt, attitude ] = rq2attitude( r_e, q_e_b )
% Converts ECEF position vector and body frame quaternion to attitude
% Computes the latitude, longitude, altitude, heading, pitch, and roll
% using the ECEF position and quaternion vectors.
% 
% @arg
% r_e       - 3 x 1 double matrix
%             ECEF position vector
% q_e_b     - 4 x 1 double matrix
%             Body frame to ECEF rotation quaternion vector
% 
% @return
% latlonalt - 3 x 1 double matrix
%             [ Latitude; Longitude; Altitude ]
% attitude  - 3 x 1 double matrix
%             [ Heading; Pitch; Roll ]
% 
% @author: Matt Marti
% @date: 2019-03-06

% ECEF to body frame rotation vector
R_e_b = quat2dircos(q_e_b); % Body to ECEF rotation matrix
R_b_e = R_e_b';

% 11. Convert position to Latitude and Longitude
latlonalt = ecef2latlon(r_e);

% 12. Local Frame to ECEF Frame Rotation
lat = latlonalt(1);
slat = sind(-lat);
clat = cosd(-lat);
Rlat = [ clat, 0, -slat; 0, 1, 0; slat, 0, clat ];
lon = latlonalt(2);
slon = sind(lon);
clon = cosd(lon);
Rlon = [ clon, slon, 0; -slon, clon, 0; 0, 0, 1 ];
% Renvframe = [ 0, 1, 0; -1, 0, 0; 0, 0, 1 ];
Renvframe = [ 0, 1, 0; 0, 0, 1; 1, 0, 0 ];
R_l_e = Renvframe * Rlat * Rlon;
R_e_l = R_l_e';
R_b_l = R_b_e * R_e_l;

% head = mod( - real(atan2d( - R_l_b(3,1), R_l_b(3,3) )), 360);
% roll = - real(atand( - R_l_b(1,2) / R_l_b(2,2) ));
% pitc = real(asind( R_l_b(3,2)));

% 13. Attitude
head = mod( - real(atan2d( - R_b_l(2,1), R_b_l(2,2) )), 360);
pitc = real(asind( R_b_l(2,3) ));
roll = real(atan2d( - R_b_l(1,3), R_b_l(3,3) ));
attitude = [ head; pitc; roll ];

end

