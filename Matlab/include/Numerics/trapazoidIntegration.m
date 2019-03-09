function [ xkp1, Fk, Gammak ] ...
    = trapazoidIntegration( xk, vk, dt, n, func, dflag )
% Integrates the time derivative of a state vector according to a function
% Solves a first order, time dependent, ordinary differential equation
% 
% @arg
% xk     - nx x 1 double matrix
%          State vector at time k
% vk     - nv x 1 double matrix
%          Process noise vector
% dt     - double
%          Delta time for integration duration
% n      - int
%          Number of integration steps
% func   - Anonymous function handle
%          Function handle for state rate of change
%          Arguments: func(delt, xk, vk, dflag)
% dflag  - bool
%          Flag to compute derivatives
% 
% @return
% xkp1   - nx x 1 double matrix
%          State vector at time k+1
% Fk     - nx x nx double matrix
%          State transition matrix
% Gammak - nx x nv double matrix
%          Noise partial derivative matrix
% 
% @author: Matt Marti
% @date: 2019-03-04

% Parameters
nx = length(xk);
nv = length(vk);
delt = dt / n;

% Initialize partial derivatives
Fk = eye(nx);
Gammak = zeros(nx,nv);

% Integration Loop
xkp1 = xk;
t = 0;
[ xdota, ~, ~ ] = func(delt, xk, vk, 1);
for i = 1:n
    
    % Next state approximation
    [ xdotb, dFb, dGammab ] = func(delt, xkp1, vk, dflag);
    xkp1 = xkp1 + 0.5*delt*(xdotb + xdota);
    
    % Update partials
    if dflag
        Fk = Fk + dFb*Fk*delt;
        Gammak = Gammak + (dFb*Gammak + dGammab)*delt;
    end
    
    % Iterate
    xdota = xdotb;
    t = t + delt;
end



end

