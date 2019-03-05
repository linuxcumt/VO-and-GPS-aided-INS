function [ R, dR_dq ] = quat2dircos( q, dflag )
% Quaternion to direction cosines matrix function
% Input quaternion vector, return direction cosines matrix. Based on the
% formula presented in Savage, P., "Strapdown Analytics", 2000.
% 
% @arg
% q     - 4 x 1 double matrix
%         Quaternion vector
% dflag - bool
%         Compute derivatives flag. True to compute derivatives
% 
% @return
% R     - 3 x 3 double matrix
%         Directions cosine matrix
% dR_dq - 3 x 3 x 4 double tensor
%         Partial derivative of direction cosine matrix with respect to
%         quanternion vector elements.
% 
% @author: Matt Marti
% @date: 2019-03-04

% Input check
if nargin == 1
    dflag = 0;
    dR_dq = NaN;
end
assert(abs(sum(q.^2) - 1) < 1e-9, 'Quaternion does not have unit norm');

% Variables
a = q(1);
b = q(2);
c = q(3);
d = q(4);

% Direction cosines
R = zeros(3,3);
R(1,1) = a^2 + b^2 - c^2 - d^2;
R(1,2) = 2*(b*c - a*d);
R(1,3) = 2*(b*d + a*c);
R(2,1) = 2*(b*c + a*d);
R(2,2) = a^2 - b^2 + c^2 - d^2;
R(2,3) = 2*(c*d - a*b);
R(3,1) = 2*(b*d - a*c);
R(3,2) = 2*(c*d + a*b);
R(3,3) = a^2 - b^2 - c^2 + d^2;


%% Partial derivative

if dflag
    n = 4;
    
    % Decision making matrix
    dR_dq = zeros(3, 3, n);
    switchmatrix = eye(n);

    % Derivative loop
    for i = 1:n

        % Activate partial derivatives
        da_ds = switchmatrix(1,i);
        db_ds = switchmatrix(2,i);
        dc_ds = switchmatrix(3,i);
        dd_ds = switchmatrix(4,i);

        % Direction cosines matrix
        r = zeros(3,3);
        
%         R(1,1) = a^2 + b^2 - c^2 - d^2;
        r(1,1) = 2*a*da_ds + 2*b*db_ds - 2*c*dc_ds - 2*d*dd_ds;
        
%         R(1,2) = 2*(b*c - a*d);
        r(1,2) = 2*((db_ds*c + b*dc_ds) - (da_ds*d + a*dd_ds));
        
%         R(1,3) = 2*(b*d + a*c);
        r(1,3) = 2*((db_ds*d + b*dd_ds) + (da_ds*c + a*dc_ds));
        
%         R(2,1) = 2*(b*c + a*d);
        r(2,1) = 2*((db_ds*c + b*dc_ds) + (da_ds*d + a*dd_ds));
        
%         R(2,2) = a^2 - b^2 + c^2 - d^2;
        r(2,2) = 2*a*da_ds - 2*b*db_ds + 2*c*dc_ds - 2*d*dd_ds;
        
%         R(2,3) = 2*(c*d - a*b);
        r(2,3) = 2*((dc_ds*d + c*dd_ds) - (da_ds*b + a*db_ds));
        
%         R(3,1) = 2*(b*d - a*c);
        r(3,1) = 2*((db_ds*d + b*dd_ds) - (da_ds*c + a*dc_ds));
        
%         R(3,2) = 2*(c*d + a*b);
        r(3,2) = 2*((dc_ds*d + c*dd_ds) + (da_ds*b + a*db_ds));
        
%         R(3,3) = a^2 - b^2 - c^2 + d^2;
        r(3,3) = 2*a*da_ds - 2*b*db_ds - 2*c*dc_ds + 2*d*dd_ds;

        % Assign Derivative
        dR_dq(1:3,1:3,i) = r;
    end
end
    
end

