function ECEFxyz = latlon2ecef(latlonalt)
% Converts WGS-84 latitude-longitude-altitude to ECEF
%
% @arg 
% latlonalt - 3x1 double array
%             vector which contains a latitude-longitude-altitude 
%             position specified by lat(deg),long(deg, alt(m)
%             FORMAT: [ latitude; longitude; altitude ]
%
% @return
% ecefxyz - 3x1 double array
%           vector which contains the position specified in the ECEF 
%           coordinate frame (meters)
%           FORMAT: [ ECEFx; ECEFy; ECEFz ]
%
% @author: Matt Marti
% @date: 2018-08-05

% Constants
global A_EARTH B_EARTH
if isempty(A_EARTH)
    constants;
end

% Input checking
if size(latlonalt,1) == 3
    transposeflag = 0;
elseif size(latlonalt,1) == 1 && size(latlonalt,2) == 3
    transposeflag = 1;
else
    error('Incorrect size of input');
end
if transposeflag
    latlonalt = latlonalt';
end

% get lat-long-alt location to be converted to ECEF coordinates
lat = latlonalt(1,:);
lon = latlonalt(2,:);
alt = latlonalt(3,:);
   
% computes the ECEF coordinates 
esq = 1 - (B_EARTH/A_EARTH)^2;
NN = A_EARTH^2./sqrt(A_EARTH^2.*cosd(lat).^2 + B_EARTH^2.*sind(lat).^2);
x = (NN+alt) .* cosd(lat).*cosd(lon); % meters
y = (NN+alt) .* cosd(lat).*sind(lon); % meters
z = (NN + alt - NN*esq) .* sind(lat); % meters

% return location in ECEF coordinates
if ~transposeflag
    ECEFxyz = [ x; y; z ];
else
    ECEFxyz = [ x, y, z ];
end

return