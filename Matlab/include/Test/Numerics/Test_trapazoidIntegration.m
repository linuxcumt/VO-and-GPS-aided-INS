function Test_trapazoidIntegration()
%% Test_trapazoidIntegration.m
% 
% Test script for trapazoidal method based numerical integration script for
% first order ODEs.
% 
% @author: Matt Marti
% @date: 2019-03-04
clc
% Run the integration function
xk = 5*ones(3,1);
vk = 4*ones(3,1);
dt = 2;
n = 20000;
func = @(dt,x,v,dflag) fscript(dt,x,v,dflag);
[ xkp1, Fk, Gammak ] = ...
             trapazoidIntegration( xk, vk, dt, n, func, 1 );

% Compute the state at tkp1 using analytical integration
Fk_true = [ 1, 2, 8; 0, 1, 2; 0, 0, 1 ];
Gammak_true = ...
    [ -0.4333333333333334, 0.3333333333333334, 0.1800000000000000; ...
      -0.0800000000000000, 0.0400000000000000, 0.0800000000000000; ...
      -0.1400000000000000, 0.0800000000000000, 0.0600000000000000];
xkp1true = Fk_true*xk + Gammak_true*vk;

% State check
diff = xkp1 - xkp1true;
for i = 1:length(xkp1)
    assert(abs(diff(i)) < 1e-3*xkp1true(i), 'Bad x');
end

% F check
diff = Fk - Fk_true;
comparenum = 1e-4*max(max(Fk));
for i = 1:size(Gammak,1)
    for j = 1:size(Gammak,2)
        assert(abs(diff(i,j)) < comparenum, 'Bad F');
    end
end

% Gamma check
diff = Gammak - Gammak_true;
comparenum = 1e-4*max(max(Gammak));
for i = 1:size(Gammak,1)
    for j = 1:size(Gammak,2)
        assert(abs(diff(i,j)) < comparenum, 'Bad Gamma');
    end
end


%% Out
fprintf('PASSED: Test_trapazoidIntegration\n');


end


function [xkp1, A, D] = fscript(dt, xk, vk, dflag) %#ok

A = [0, 1, 3; 0, 0, 1; 0, 0, 0];
D = 1e-2*[1, 4, -3; 3, -2, 1; -7, 4, 3];

xkp1 = A*xk + D*vk;


end