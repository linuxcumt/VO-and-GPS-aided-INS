function [ a, da_dr ] = centrifugalmodel( r, dflag )
% Computes acceleration due to gravity
% 
% @arg
% r     - 3 x 1 double matrix
%         Position vector in ECEF coordinates
% dflag - bool
%         Compute derivative flag
% 
% @return
% a     - 3 x 1 double matrix
%         Gravity acceleration vector
% da_dr - 3 x 3 double matrix
%         Partial derivative matrix with respect to input r
% 
% @author: Matt Marti
% @date: 2019-03-04

% Constants
global OMEGA_EARTH
if isempty(OMEGA_EARTH)
    constants;
end

% Check input
if nargin == 1
    dflag = 0;
    da_dr = NaN;
end

% Centrifugal force
omegasq = OMEGA_EARTH^2;
a = omegasq * [r(1); r(2); 0];


%% Partial derivative calculation

if dflag
    
    da_dr = omegasq * [ 1, 0, 0; 0, 1, 0; 0, 0, 0 ];
    
end

end

