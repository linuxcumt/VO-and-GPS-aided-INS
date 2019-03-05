%% Test_insSinc
%
% Tests the custom sinc function
%
% @author: Matt Marti
% @date: 2019-03-04

clear;


%% Zero
% x0 = 0;
% [ sa, dsa_dx ] = insSinc( x0, 1 );
% assert(sa == 1, 'Bad value at x = 0');
% assert(dsa_dx == 0, 'Bad derivative at x = 0');


%% Test

x = (-1:.01:1)';
[ sa, dsa_dx ] = insSinc( x, 1 );

% True
satrue = sin(x) ./ x;
dsa_dxtrue = cos(x) ./ x - sin(x) ./ x.^2;
i = 1;
while x(i) && i < length(x)
    i = i + 1;
end
assert(i ~= length(x), 'Couldn''t find zero x');
satrue(i) = 1;

% Check value
for i = 1:length(x)
    assert(abs(satrue(i) - sa(i)) < 1e-12, 'failed value');
end

% Numerical derivative
dx = 1e-4;
D_num = zeros(length(x));
for i = 1:length(x)
    xp = x;
    xm = x;
    xp(i) = xp(i) + dx;
    xm(i) = xm(i) - dx;
    sap = insSinc( xp );
    sam = insSinc( xm );
    D_num(:,i) = (sap - sam) / (2*dx);
end

% Assertions
p = 2e-6;
D_diff = dsa_dx - D_num;
Dcheck = abs(D_diff) <= p;
assert(prod(prod(Dcheck)) == 1, 'Bad Derivative');


%% Output
fprintf('PASSED: Test_insSinc\n');