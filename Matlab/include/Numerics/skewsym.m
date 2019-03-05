function [ X, dX_dx ] = skewsym( x, dflag )
% Creates the Skew Symetric form of the vector
% 
% @arg
% x     - 3 x 1 double matrix
%         Input vector
% dflag - bool
%         Flag to compute derivative
% 
% @return
% X     - 3 x 1 double matrix
%         Skew symmetric matrix
% dX_dx - 3 x 3 x 3 double tensor
%         Derivative of skew symmetric matrix with respect to input x
% 
% @author: Matt Marti
% @2019-03-04

if nargin == 1
    dflag = 0;
    dX_dx = NaN;
end

X = [ 0, -x(3), x(2); x(3), 0, -x(1); -x(2), x(1), 0 ];

if dflag
    dX_dx = zeros(3,3,3);
    dX_dx(2,3,1) = -1;
    dX_dx(3,2,1) = 1;
    dX_dx(1,3,2) = 1;
    dX_dx(3,1,2) = -1;
    dX_dx(1,2,3) = -1;
    dX_dx(2,1,3) = 1;
end

end

