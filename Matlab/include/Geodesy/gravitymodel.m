function [ g, dg_dr ] = gravitymodel( r, dflag )
% Computes acceleration due to gravity
% 
% @arg
% r     - 3 x 1 double matrix
%         Position vector in ECEF coordinates
% dflag - bool
%         Compute derivative flag
% 
% @return
% g     - 3 x 1 double matrix
%         Gravity acceleration vector
% dg_dr - 3 x 3 double matrix
%         Partial derivative matrix with respect to input r
% 
% @author: Matt Marti
% @date: 2019-03-04

% Constants
global MUEARTH
if isempty(MUEARTH)
    constants;
end

% Check input
if nargin == 1
    dflag = 0;
    dg_dr = NaN;
end

% Unit vector
r1 = r(1);
r2 = r(2);
r3 = r(3);
sumrsq = r1^2 + r2^2 + r3^2;
r = sqrt(sumrsq);
oneoverr = 1 / r;
runit = [ r1; r2; r3 ] * oneoverr;

% Gravity scalar force
g0 = - MUEARTH / sumrsq;
g = g0 * runit;


%% Partial derivative calculation

if dflag
    
    % Decision matrix and initialize partial
    smat = eye(3,3);
    dg_dr = zeros(3,3);
    
    % Loop
    for i = 1:3
        
        % Activate partials
        dr1_ds = smat(1,i);
        dr2_ds = smat(2,i);
        dr3_ds = smat(3,i);
        
        % Unit vector
        dsumrsq_ds = 2*r1*dr1_ds + 2*r2*dr2_ds + 2*r3*dr3_ds;
        dr_ds = -0.5/r * dsumrsq_ds;
        doneoverr_ds = - oneoverr^2 * dr_ds;
        drunit_ds = [dr1_ds;dr2_ds;dr3_ds]*oneoverr + [r1;r2;r3]*doneoverr_ds;

        % Gravity scalar force
        dg0_ds = MUEARTH / sumrsq^2 * dsumrsq_ds;
        dg_ds = dg0_ds*runit + g0*drunit_ds;
        
        % Assign output
        dg_dr(:,i) = dg_ds;
    end

end

