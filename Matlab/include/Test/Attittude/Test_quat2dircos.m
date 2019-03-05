%% Test_quat2dircos.m
% 
% Script to test partial derivatives of the Quaternion to Direction Cosines
% matrix
% 
% @author: Matt Marti
% @date: 2019-03-04

clear


%% Test
    
% State
q = [-5; 4; 3; 1];
q = q/norm(q);

% Run function for exact values
[ R ] = quat2dircos( q, 0 );

% Gradient descent test for derivatives b/c numerical derivatives get
% distorted due to quaternion normalization
q_perturb = q + [.2; -.6; -.3; .8].*q;
q_perturb = q_perturb / norm(q_perturb);
err = @(qq) sqrt(sum((qq - q).^2));
iter = 0;
while err(q_perturb) > 1e-12 && iter < 1000
    
    [ r, dr ] = quat2dircos( q_perturb, 1 );
    de_dq = zeros(4,1);
    for k = 1:4
        for i = 1:3
            for j = 1:3
                de_dq(k) = de_dq(k) + 2 * (r(i,j)-R(i,j))*dr(i,j,k);
            end
        end
    end
    q_perturb = q_perturb - (1/(iter+1)^2)*de_dq;
    q_perturb = q_perturb ./ norm(q_perturb);
    err(q_perturb);
    
    iter = iter + 1;
end
assert(iter < 1000, 'Correct q not found using derivative');
assert(logical(prod(abs(q - q_perturb) < 1e-12)), 'Incorrect q');

% Assertions for rotation matrix
p = 1e-9;
Rcheck = abs(inv(R)' - R) < p;
assert(prod(prod(Rcheck)) == 1, 'Failed rotation check');


%% Output
fprintf('PASSED: Test_quat2dircos\n');