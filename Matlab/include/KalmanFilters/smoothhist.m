function [xstarhist,Pstarhist,epsstarhist] = ...
    smoothhist(Fhist,xhathist,Phist,xbarhist,Pbarhist)
% First Order Smoothing Filter implementation
% Given the time history from either a Linear or Extended Kalman Filter, 
% computes the smoothed estimates of the state x and associtated 
% covariances for each measurement time. This smoother uses first order
% approximations of the covariances for non-linear systems.
% 
% This function performs the first order smoothed Kalman Filter state 
% estimate for the dynamic system model
% 
%   x(k+1)   = f( k, x(k), u(k), v(k) )
%   F(k)     = df / dx(k)
%   Gamma(k) = df / dv(k)
%   z(k+1)   = h( k+1, x(k+1) ) + w(k+1)
%   H(k+1)   = dh / dx(k+1)
% 
% Where v(k) and w(k+1) are Gaussian distributed random noise terms. The
% process noise v(k) has a mean of zero, and a Covariance Q(k). The
% measurement noise w(k+1) has a mean of zero, and a Covariance of R(k+1).
% 
% Variables ending in "_kp1" represent the variable at the (k+1) state.
% Variables ending in "_k" represent variable at the (k) state.
% 
% Note that for dimensions that have size K, if a tensor with K = 1 is
% passed as an argument, this smoother funnction will assume that the
% system is time invariant.
% 
% INPUTS
% Fhist       - nx x nx x K double tensor
%               State Transition Matrix at each k state
%               For a time invariant system, set to a nx x nx double matrix
% xhathist    - nx x K double matrix
%               A Posteriori state estimate time history
% Phist       - nx x nx x K double tensor
%               A Posteriori state covariance time history
% xbarhist    - nx x K double matrix
%               A Prioir state estimate time history
% Pbarhist    - nx x nx x K double tensor
%               A Priori state covariance time history
% 
% OUTPUTS     
% xstarhist   - nx x K double matrix
%               Smoothed State estimates for each state
% Pstarhist   - nx x nx x K double tensor
%               Smoothed State Covariance matrix for each state.
% epsstarhist - K x 1 double vector
%               Smoothed Innovation Statistic time history
% xhathist    - nx x K double matrix
%               Forward pass State estimates
% Phist       - nx x nx x K double tensor
%               Forward pass State Covariance matrix for each state.
% epshist     - K x 1 double vector
%               Forward pass Innovation Statistic
% 
% DEPENDENCIES
% kalmaniter.m v 2018-10-29
% 
% @author: Matt Marti
% @date: 2018-11-28

% Measure sizes
nx = size(xhathist, 1);
K  = size(xhathist, 2);
tiflag = size(Fhist, 3) == 1;


%% Input checking

% Fhist - nx x nx x K double tensor
assert(size(Fhist,1) == nx && ...
    size(Fhist,2) == nx && ...
    ((size(Fhist,3)==K && ~tiflag) || ((size(Fhist,3)==1) && tiflag)), ...
    'Incorrect size of Fhist');

% xhathist - nx x K double matrix
assert(size(xhathist,1) == nx && ...
    size(xhathist,2) == K && ...
    size(xhathist,3) == 1  ,...
    'Incorrect size of xhathist');

% Phist - nx x nx x K double tensor
assert(size(Phist,1) == nx && ...
    size(Phist,2) == nx && ...
    size(Phist,3) == K  ,...
    'Incorrect size of Phist');

% xbarhist - nx x K double matrix
assert(size(xbarhist,1) == nx && ...
    size(xbarhist,2) == K && ...
    size(xbarhist,3) == 1  ,...
    'Incorrect size of xbarhist');

% Pbarhist - nx x nx x K double tensor
assert(size(Pbarhist,1) == nx && ...
    size(Pbarhist,2) == nx && ...
    size(Pbarhist,3) == K  ,...
    'Incorrect size of Pbarhist');


%% Computation

% Preallocate data arrays
xstarhist   = zeros(nx,K);
Pstarhist   = zeros(nx,nx,K);
epsstarhist = zeros(1,K);

% Run backwards smoothing filter
xstarhist(:,K) = xhathist(:,K);
Pstarhist(:,:,K) = Phist(:,:,K);
for k = K-1:-1:1
    
    % Iterate
    kp1 = k + 1;
    xbarkp1 = xbarhist(:,kp1);
    Pbarkp1 = Pbarhist(:,:,kp1);
    invPbarkp1 = inv(Pbarkp1);
    xhatk = xhathist(:,k);
    Pk = Phist(:,:,k);
    xstarkp1 = xstarhist(:,kp1);
    Pstarkp1 = Pstarhist(:,:,kp1);
    
    % Assign Fk based on time-invarient system
    if ~tiflag
        Fk = Fhist(:,:,kp1);
    else
        Fk = Fhist;
    end
    
    % Compute backwards step
    Ck = Pk*(Fk')*invPbarkp1; %#ok
    innovkp1 = xstarkp1 - xbarkp1;
    xstark = xhatk + Ck*innovkp1;
    Pstark = Pk - Ck*(Pbarkp1 - Pstarkp1)*(Ck');
    
    % Innovation statistic
    xdiffk = xstark - xhatk;
    epsilonk = (xdiffk')*inv(Pk)*xdiffk; %#ok
    
    % Save data
    xstarhist(:,k) = xstark;
    Pstarhist(:,:,k) = Pstark;
    epsstarhist(k) = epsilonk;
end

