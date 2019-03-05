%% Test_centrifugalmodel
%
% Tests the centrifugal force model
%
% @author: Matt Marti
% @date: 2019-03-04

clear;


%% Test

ecef = [1061985.2311866253;-4837817.3795761978;4005245.7246151697];
[ a, da_dr ] = centrifugalmodel( ecef, 1 );

% Not sure what value it should be, but it's low
assert(norm(a) < 0.1 && norm(a) > 0);

% Numerical derivative
dr = 1e1;
D_num = zeros(3,3);
for i = 1:3
    rp = ecef;
    rm = ecef;
    rp(i) = rp(i) + dr;
    rm(i) = rm(i) - dr;
    ap = centrifugalmodel( rp );
    am = centrifugalmodel( rm );
    D_num(:,i) = (ap - am) / (2*dr);
end

% Assertions
p = 2e-6;
D_diff = da_dr - D_num;
Dcheck = abs(D_diff) <= p;
assert(prod(prod(Dcheck)) == 1, 'Bad Derivative');


%% Output
fprintf('PASSED: Test_centrifugalmodel\n');