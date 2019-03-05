%% Test_skewsym
% 
% Brief test to make sure that the skew symmetric matrix function works as
% intended.

clear


%% Test

% Data
x = [1; 3; 5];
y = [6; 3; 6];

% Function
[X, dX] = skewsym(x, 1);

% Correctness test
ztrue = cross(x,y);
z = X * y;
for i = 1:length(z)
    assert(abs(z(i)-ztrue(i)) < 1e-9, 'Bad result');
end

% Derivative check
p = 1e-12;
d1 = [ 0, 0, 0; 0, 0, -1; 0, 1, 0 ];
d2 = [ 0, 0, 1; 0, 0, 0; -1, 0, 0 ];
d3 = [ 0, -1, 0; 1, 0, 0; 0, 0, 0 ];
check1 = prod(prod(abs(d1 - squeeze(dX(:,:,1))) < p)) == 1;
check2 = prod(prod(abs(d2 - squeeze(dX(:,:,2))) < p)) == 1;
check3 = prod(prod(abs(d3 - squeeze(dX(:,:,3))) < p)) == 1;
assert(check1, 'Failed derivative check')
assert(check3, 'Failed derivative check')
assert(check2, 'Failed derivative check')



%% Out
fprintf('PASSED: Test_skewsym\n');