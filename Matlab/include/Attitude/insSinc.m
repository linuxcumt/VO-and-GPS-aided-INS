function [ sa, dsa_dx ] = insSinc(x, dflag)
% Unscaled sinc function
% sinc(x) = sin(x)/x for x ~= 0, = 1 otherwise
% 
% @arg
% x      - nx x 1 double matrix
%          Input vector
% dflag  - bool
%          Flag to calculate derivative. True to calculate derivative
% @return
% sa     - nx x 1 double matrix
%          Output of sinc function ( y in y = sin(x)/x )
% dsa_dx - nx x nx double matrix
%          Derivative of output sa with respect to the input vector x
% 
% @author: Matt Marti
% @date: 2019-03-05

% Input check
if nargin == 1
    dflag = 0;
    dsa_dx = NaN;
end
assert(size(x,2) == 1, 'argument x must be column vector and is not');
nx = size(x,1);

% Function
sa = zeros(nx,1);
if dflag
    dsa_dx = zeros(nx,nx);
end
for i = 1:nx
    if x(i) ~= 0
        sinx = sin(x(i));
        oneoverx = 1 / x(i);
        sa(i) = sinx * oneoverx;
        if dflag
            dsinx_dx = cos(x(i));
            doneoverx_dx = - 1 / x(i)^2;
            dsa_dx(i,i) = dsinx_dx*oneoverx + sinx*doneoverx_dx;
        end
    else
        sa(i) = 1;
        if dflag
            dsa_dx(i,i) = 0;
        end
    end
end

end

