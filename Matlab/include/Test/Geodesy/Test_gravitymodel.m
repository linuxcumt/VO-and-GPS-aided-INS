%% Test_gravitymodel
%
% Tests the gravity model
%
% @author: Matt Marti
% @date: 2019-03-04

clear;


%% Test

ecef = [1061985.2311866253;-4837817.3795761978;4005245.7246151697];
[ g, dg_dr ] = gravitymodel( ecef, 1 );

% Check gravity value
assert(abs(norm(g) - 9.81) < 2e-2, 'Bad gravity norm')
assert(prod( abs(g/norm(g) + ecef/norm(ecef)) < 1e-12 ) == 1, ...
    'Wrong direction unit vector');

% Numerical derivative
dr = 1e1;
D_num = zeros(3,3);
for i = 1:3
    rp = ecef;
    rm = ecef;
    rp(i) = rp(i) + dr;
    rm(i) = rm(i) - dr;
    gp = gravitymodel( rp );
    gm = gravitymodel( rm );
    D_num(:,i) = (gp - gm) / (2*dr);
end

% Assertions
p = 2e-6;
D_diff = dg_dr - D_num;
Dcheck = abs(D_diff) <= p;
assert(prod(prod(Dcheck)) == 1, 'Bad Derivative');


%% Output
fprintf('PASSED: Test_gravitymodel\n');