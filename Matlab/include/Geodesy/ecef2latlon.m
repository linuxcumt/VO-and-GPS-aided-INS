function latlonalt = ecef2latlon(ecef)
% Converts WGS-84 latitude-longitude-altitude to ECEF
%
% @arg 
% ecefxyz - 3x1 double array
%           vector which contains the position specified in the ECEF 
%           coordinate frame (meters)
%           FORMAT: [ ECEFx; ECEFy; ECEFz ]
%
% @return
% latlonalt - 3x1 double array
%             vector which contains a latitude-longitude-altitude 
%             position specified by lat(deg),long(deg, alt(m)
%             FORMAT: [ latitude; longitude; altitude ]
%
% @author: Matt Marti
% @date: 2018-08-05

% define physical constants
global A_EARTH B_EARTH ECEF2LATLON_MAXITER ECEF2LATLON_PRECISION DEGRAD
if isempty(A_EARTH)
    constants;
end

% Input checking
if size(ecef,1) == 3
    transposeflag = 0;
elseif size(ecef,1) == 1 && size(ecef,2) == 3
    transposeflag = 1;
else
    error('Incorrect size of input');
end
if transposeflag
    ecef = ecef';
end

% get ECEF location to be converted to latitude-longitude-altitude coords
x = ecef(1,:);
y = ecef(2,:);
z = ecef(3,:);

% compute the longitude which is an exact calculation
lon = atan2( y, x ); % radians

% compute the latitude using iteration
p = sqrt(x.^2 + y.^2);

% compute approximate latitude
lat = atan2(z, p);

% Loop
iter = 0;
while (iter < ECEF2LATLON_MAXITER)
    N0 = A_EARTH^2./sqrt(A_EARTH^2.*cos(lat).^2 + B_EARTH^2.*sin(lat).^2);
    alt = p./cos(lat) - N0; % meters
	
	% calculate improved latitude
    e_sqr = 1 - (B_EARTH/A_EARTH)^2 ;
    lat_new  = atan2(z./p, 1 - (e_sqr*N0./(N0+alt))); % radians
	
	% check if result is close enough, 
    if (abs(lat_new - lat) < ECEF2LATLON_PRECISION)
        lat = lat_new;
        break;
    end
    lat = lat_new;
    iter = iter + 1;
end

% convert the latitude and longitude to degrees
lat = lat / DEGRAD; % degrees
lon = lon / DEGRAD; % degrees

% return location in latitude-longitude-altitude coordinates
if ~transposeflag
    latlonalt = [ lat; lon; alt ];
else
    latlonalt = [ lat, lon, alt ];
end


return;