function [ xkp1, Fk, Gammak ] = dynamics_model( k, xk, vk )
% INS Dead-Reckoning Dynamics Model
% 
% @arg
% k         - int
%             time index
% xk        - 21x1 double vector
%             State vector at time k
% vk        - 21x1 double vector
%             Process Noise vector at time k
% derivflag - bool
%             True to enable calculation of derivative
% 
% @return
% xkp1      - 11x1 double vector
%             State vector at time k
% Fk        - 11x11 double matrix
%             State Transition partial derivative
% Gammak    - 11x11 double matrix
%             State Noiise partial derivative
% 
% @author: Matt Marti
% @date: 2019-03-01

global dt

% State Transition Partial Derivative
Fk = zeros(11,11);
Fk(1:3,1:3) = [ 1, dt, dt^2/2; 0, 1, dt; 0, 0, 1];
Fk(4:6,4:6) = [ 1, dt, dt^2/2; 0, 1, dt; 0, 0, 1];
Fk(7:8,7:8) = [ 1, dt; 0, 1 ];
Fk(9:11,9:11) = eye(3);

Gammak = zeros(11,6);
Gammak(3,1) = 1;
Gammak(6,2) = 1;
Gammak(8,3) = 1;
Gammak(9,4) = 1;
Gammak(10,5) = 1;
Gammak(11,6) = 1;

xkp1 = Fk * xk + Gammak * vk;

end

