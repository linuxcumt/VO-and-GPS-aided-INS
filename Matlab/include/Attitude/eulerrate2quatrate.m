function [ qdot, dqdot_ds ] = eulerrate2quatrate( omega, q, dflag )
% Quaternion to direction cosines matrix function
% Input quaternion vector, return direction cosines matrix
% 
% @arg
% omega    - 3 x 1 double matrix
%            Rotation rate for Euler angles
%            Angles measured in radians!
% q        - 4 x 1 double matrix
%            Quaternion vector
% 
% @return
% dq_dt    - 4 x 1 double matrix
%            Quaternion rate of change vector
% dqdot_ds - 4 x 7 double matrix
%            Partial derivative with respect to inputs. 
%            Takes the form: [ dqdot_domega(4x3); dqdot_dq(4x4) ]
% 
% @author: Matt Marti
% @date: 2019-03-04

% Input check
if nargin == 1
    dflag = 0;
    dqdot_ds = NaN;
end
assert(abs(sum(q.^2) - 1) < 1e-9, 'Quaternion does not have unit norm');

omega


end

