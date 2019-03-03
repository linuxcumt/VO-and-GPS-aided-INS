function [ xkp1, Fk, Gammak ] = dynamics_model( dt, xk, vk )
% INS Dead-Reckoning Dynamics Model
% 
% INPUTS
% k      - double
%          Delta time between state updates. Recommended to be INS
%          measurement rate.
% xk     - 21x1 double vector
%          State vector at time k
% vk     - 21x1 double vector
%          Process Noise vector at time k
% 
% OUTPUTS
% xkp1   - 21x1 double vector
%          State vector at time k
% Fk     - 21x21 double matrix
%          State Transition partial derivative
% Gammak - 21x21 double matrix
%          State Noiise partial derivative
% 
% @author: Matt Marti
% @date: 2019-02-25

% State Transition Partial Derivative
Fk = zeros(21,21);
Fk(1:3,1:3) = [ 1, dt, dt^2/2; 0, 1, dt; 0, 0, 1];
Fk(4:6,4:6) = [ 1, dt, dt^2/2; 0, 1, dt; 0, 0, 1];
Fk(7:9,7:9) = [ 1, dt, dt^2/2; 0, 1, dt; 0, 0, 1];
Fk(10:11,10:11) = [ 1, dt; 0, 1 ];
Fk(12:13,12:13) = [ 1, dt; 0, 1 ];
Fk(14:15,14:15) = [ 1, dt; 0, 1 ];
Fk(16:21,16:21) = eye(6);

Gammak = zeros(21,12);
Gammak(3,1) = 1;
Gammak(6,2) = 1;
Gammak(9,3) = 1;
Gammak(11,4) = 1;
Gammak(13,5) = 1;
Gammak(15,6) = 1;
Gammak(16,7) = 1;
Gammak(17,8) = 1;
Gammak(18,9) = 1;
Gammak(19,10) = 1;
Gammak(20,11) = 1;
Gammak(21,12) = 1;

xkp1 = Fk * xk + Gammak * vk;

end

